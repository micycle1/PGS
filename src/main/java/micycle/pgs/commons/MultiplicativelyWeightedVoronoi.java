package micycle.pgs.commons;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.tinspin.index.Index.PointEntryKnn;
import org.tinspin.index.PointMap;
import org.tinspin.index.kdtree.KDTree;

import net.jafama.FastMath;
import processing.core.PVector;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * <b>Multiplicatively Weighted Voronoi Diagrams</b>
 * <p>
 * MWVDs are a generalisation of Voronoi diagrams where each site has a weight
 * associated with it. These weights influence the boundaries between cells in
 * the diagram. Instead of being equidistant from generator points, the
 * boundaries are defined by the <b>ratio</b> of distances to the weighted
 * generator points. This results in characteristically curved cell boundaries,
 * unlike the straight line boundaries seen in standard Voronoi diagrams.
 * <p>
 * For two sites a,b the region of a's "dominance" over b is defined as the
 * Apollonius circle made from a,b. This circle's perimeter represents the set
 * of all points where the ratio of distances to a and b is equal to the ratio
 * of their weights. The MWV cell for a site is then formed by the geometric
 * intersection of all the Apollonius circles it forms with every other site in
 * the diagram. In other words, it is the region where that site "dominates" (in
 * terms of distance-weighting) all other sites.
 * <p>
 * The intersection of Apollonius circles provides a very good approximation to
 * the geometrically "true" MWVD; this approach is considerably simpler both
 * conceptually and in implementation than alternative approaches such as
 * wavefront propagation.
 * 
 * @author Michael Carleton
 */
public class MultiplicativelyWeightedVoronoi {

	private final static GeometryFactory geometryFactory = new GeometryFactory();

	private MultiplicativelyWeightedVoronoi() {
	}

	/**
	 * Computes Multiplicatively Weighted Voronoi (MWV) diagram from a list of
	 * weighted points.
	 *
	 * @param vectors List of PVectors where: x, y are coordinates of the site, and
	 *                z is the weight.
	 * @param bounds  [minX, minY, maxX, maxY] array defining the rectangular
	 *                extent. for the voronoi diagram. This should enclose all
	 *                points.
	 * @return List of Polygons representing MWV regions.
	 */
	public static List<Geometry> getMWVFromPVectors(List<PVector> vectors, double[] bounds) {
		List<Coordinate> sites = new ArrayList<>();
		List<Double> weights = new ArrayList<>();

		for (PVector vector : vectors) {
			sites.add(new Coordinate(vector.x, vector.y));
			weights.add((double) vector.z);
		}

		Envelope extent = new Envelope(bounds[0], bounds[2], bounds[1], bounds[3]); // NOTE x,x,y,y

		boolean allWeightsDifferent = weights.stream().distinct().count() == weights.size();

		return allWeightsDifferent ? getMWVDFast(sites, weights, extent) : getMWVD(sites, weights, extent);
	}

	/**
	 * Optimised implementation. Reduces redundant geometric boolean operations.
	 * Requires all weights to be different (doesn't handle bisectors to remain
	 * simpler).
	 */
	private static List<Geometry> getMWVDFast(List<Coordinate> sites, List<Double> weights, Envelope extent) {
		// The MWV cell for each site is equivalent to:
		// intersection(all ap_circles containing s)-union(all ap_circles not containing
		// s)

		final Geometry extentGeometry = geometryFactory.toGeometry(extent);

		return sites.parallelStream().map(site -> {
			int i = sites.indexOf(site); // Get index of current site

			// s2 dominates s1 (hence s1 contained in apollo circle)
			List<double[]> inCircleData = new ArrayList<>();
			List<double[]> exCircleData = new ArrayList<>();

			// prepare apollonian circles
			for (int j = 0; j < sites.size(); j++) {
				if (i == j) {
					continue;
				}
				Coordinate site2 = sites.get(j);
				var circle = calculateWeightedApollonius(site, site2, weights.get(i), weights.get(j));
				Coordinate center = new Coordinate(circle[0], circle[1]);
				double radius = circle[2];
				if (radius == 0) { // ideally not!
					continue;
				}
				double distance = site.distance(center);

				if (distance < radius) {
					inCircleData.add(circle);
				} else {
					exCircleData.add(circle);
				}
			}

			List<Geometry> inCircles = new ArrayList<>();
			List<Geometry> exCircles = new ArrayList<>();

			if (inCircleData.isEmpty()) {
				// if no incircles, region is defined by difference between plane and exCircles
				exCircleData.forEach(c -> exCircles.add(createCircle(c[0], c[1], c[2])));
				var outerDom = UnaryUnionOp.union(exCircles);
				return extentGeometry.difference(outerDom);
			}
			/*
			 * NOTE optimisation: The sites's dominance area can be no bigger than its
			 * smallest circle. Therefore any exCircles that do not intersect it do not need
			 * to be carved away. Will usually reduce the number of outcircles to union
			 * together.
			 */
			inCircleData.sort((a, b) -> Double.compare(a[2], b[2])); // find smallest circle
			var maxDominanceRegion = inCircleData.get(0);
			exCircleData = exCircleData.stream().filter(c -> {
				double dx = maxDominanceRegion[0] - c[0];
				double dy = maxDominanceRegion[1] - c[1];
				double distanceSquared = (dx * dx) + (dy * dy);
				double radiusSum = maxDominanceRegion[2] + c[2];
				return distanceSquared <= radiusSum * radiusSum;
			}).toList();

			exCircleData.forEach(c -> exCircles.add(createCircle(c[0], c[1], c[2])));

			/*
			 * NOTE optimisation: filter out redundant circles that completely cover other
			 * circles to reduce number of required intersection operations. List is already
			 * sorted above.
			 */
			List<double[]> essentialCircles = new ArrayList<>(inCircleData);
			for (int k = 0; k < inCircleData.size(); k++) {
				double[] smaller = inCircleData.get(k);

				// Remove any larger circles that completely contain this one
				essentialCircles.removeIf(larger -> {
					if (larger[2] <= smaller[2])
						return false; // Skip if not larger

					double dx = smaller[0] - larger[0];
					double dy = smaller[1] - larger[1];
					double distanceSquared = (dx * dx) + (dy * dy);
					return Math.sqrt(distanceSquared) + smaller[2] <= larger[2];
				});
			}
			essentialCircles.forEach(c -> inCircles.add(createCircle(c[0], c[1], c[2])));

			// intersect all inCircles to find the dominant region for this site
			var localDominance = inCircles.stream().reduce((geom1, geom2) -> geom1.intersection(geom2));

			if (exCircles.isEmpty()) {
				return localDominance.get().intersection(extentGeometry);
			} else {
				var outerDom = UnaryUnionOp.union(exCircles);
				return localDominance.get().difference(outerDom).intersection(extentGeometry);
			}
		}).toList();
	}

	/**
	 * Vanilla implementation. Handles equal-weighted pairs.
	 */
	private static List<Geometry> getMWVD(List<Coordinate> sites, List<Double> weights, Envelope extent) {
		if (sites.size() != weights.size()) {
			throw new IllegalArgumentException("The number of sites and weights must be equal.");
		}

		final PointMap<Point> tree = KDTree.create(2);
		for (int i = 0; i < sites.size(); i++) {
			Point p = geometryFactory.createPoint(sites.get(i));
			p.setUserData(weights.get(i));
			tree.insert(new double[] { p.getX(), p.getY() }, p);
		}

		Geometry extentGeometry = geometryFactory.toGeometry(extent);

		/*
		 * NOTE optimisation. Look at first 30 (or N/3, if larger, but up to 60)
		 * neighbors only, since subsequent neighbors usually have negligible effect.
		 * This makes producing MWVD for hundreds of points a lot more feasible.
		 */
		final int n = Math.min(Math.max(sites.size() / 3, Math.min(sites.size(), 30)), 60);

		List<Geometry> polygons = sites.stream().map(site -> {
			var query = tree.queryKnn(new double[] { site.x, site.y }, n);
			var me = query.next().value(); // first query is always the site itself
			Double w = (Double) me.getUserData();
			Geometry dominance = extentGeometry; // dominance area begins with whole plane

			List<PointEntryKnn<Point>> neighbors = new ArrayList<>();
			query.forEachRemaining(item -> {
				neighbors.add(item);
			});

			/*
			 * The dominance region of a site is formed by taking the boolean AND
			 * (intersection) of all the Apollonius circles it forms with every other site.
			 */
			for (var otherSite : neighbors) {
				Coordinate other = otherSite.value().getCoordinate();
				Double wOther = (Double) otherSite.value().getUserData();
				Geometry localDominanceCircle = apolloniusCircle(site, other, w, wOther, extentGeometry);
				dominance = dominance.intersection(localDominanceCircle);
			}

			return dominance != null && !dominance.isEmpty() ? dominance : null;
		}).filter(dominance -> dominance != null) // Filter out empty geoms
				.toList();

		return polygons;
	}

	private static Geometry apolloniusCircle(Coordinate s1, Coordinate s2, double w1, double w2, Geometry extentG) {
		Geometry localDominanceCircle;
		if (w1 == w2 || s1.distance(s2) < 1e-7) { // Regular Voronoi (Perpendicular Bisector)
			localDominanceCircle = calculatePerpendicularBisector(s1, s2, extentG);
		} else { // Weighted Voronoi (Circle)
			double[] circle = calculateWeightedApollonius(s1, s2, w1, w2);
			Coordinate center = new Coordinate(circle[0], circle[1]); // often lies outside bounds
			double radius = circle[2];
			double distance = s1.distance(center);
			localDominanceCircle = createCircle(center.x, center.y, radius);
			/*
			 * The circle will either enclose site 1 (i.e. it's dominated by site 2), or
			 * will bend away from site 1 (enclosing and dominating site 2).
			 */
			if (distance < radius) { // circle encloses site 1
				return localDominanceCircle.intersection(extentG);
			} else { // circle encloses site 2
				// site 1's dominance is the rest of the plane
				var out = extentG.difference(localDominanceCircle);
				return out;
			}
		}

		Point s1Geom = geometryFactory.createPoint(s1);

		if (s1Geom.intersects(localDominanceCircle)) {
			return localDominanceCircle.intersection(extentG);
		} else {
			return extentG.difference(localDominanceCircle);
		}
	}

	private static Geometry calculatePerpendicularBisector(Coordinate s1, Coordinate s2, Geometry extent) {
		Coordinate midPoint = new Coordinate((s1.x + s2.x) / 2, (s1.y + s2.y) / 2);
		double dx = s2.x - s1.x;
		double dy = s2.y - s1.y;
		double length = extent.getEnvelopeInternal().getWidth() + extent.getEnvelopeInternal().getHeight();
		Coordinate endPoint1, endPoint2;

		if (dy == 0) {
			// Handle the case where dy is zero (horizontal line)
			// The bisector will be a vertical line passing through the midpoint.
			endPoint1 = new Coordinate(midPoint.x, midPoint.y + length / 2);
			endPoint2 = new Coordinate(midPoint.x, midPoint.y - length / 2);
		} else {
			// regular case (non-horizontal line)
			double angle = Math.atan2(dy, dx) + Math.PI / 2;
			endPoint1 = new Coordinate(midPoint.x + length * Math.cos(angle), midPoint.y + length * Math.sin(angle));
			endPoint2 = new Coordinate(midPoint.x - length * Math.cos(angle), midPoint.y - length * Math.sin(angle));
		}

		// create the bisector line
		LineString bisectorLine = geometryFactory.createLineString(new Coordinate[] { endPoint1, endPoint2 });

		// split the extent with the bisector
		Geometry nodedLine = extent.getBoundary().union(bisectorLine);
		Polygonizer polygonizer = new Polygonizer();
		polygonizer.add(nodedLine);

		Collection<Geometry> polygons = polygonizer.getPolygons();
		Geometry[] polygonArray = polygons.toArray(new Geometry[0]);

		// both [0] and [1] seem to work.
		return polygonArray[0];
	}

	private static double[] calculateWeightedApollonius(Coordinate s1, Coordinate s2, double w1, double w2) {
		double s1x = s1.x;
		double s1y = s1.y;
		double s2x = s2.x;
		double s2y = s2.y;

		double den = 1.0 / (w1 * w1 - w2 * w2);
		double cx = (w1 * w1 * s2x - w2 * w2 * s1x) * den;
		double cy = (w1 * w1 * s2y - w2 * w2 * s1y) * den;
		double d = Math.sqrt(((s1x - s2x) * (s1x - s2x) + (s1y - s2y) * (s1y - s2y)));
		double r = w1 * w2 * d * den;
		if (r < 0) {
			r = r * -1;
		}
		// NOTE r can be huge (as circle may tend towards straight line)
		return new double[] { cx, cy, r };
	}

	private static Polygon createCircle(double x, double y, double r) {
		final double maxDeviation = 0.5;
		// Calculate the number of points based on the radius and maximum deviation.
		int nPts = (int) Math.ceil(2 * Math.PI / Math.acos(1 - maxDeviation / r));
		nPts = Math.max(nPts, 21); // min of 21 points for tiny circles
		final int circumference = (int) (Math.PI * r * 2);
		if (nPts > circumference * 2) {
			// AT MOST 1 point every half pixel, hard limit=100000
			nPts = Math.min(100000, Math.abs(circumference * 2));
		}

		Coordinate[] pts = new Coordinate[nPts + 1];
		for (int i = 0; i < nPts; i++) {
			double ang = i * (2 * Math.PI / nPts);
			double px = r * FastMath.cos(ang) + x;
			double py = r * FastMath.sin(ang) + y;
			pts[i] = new Coordinate(px, py);
		}
		pts[nPts] = new Coordinate(pts[0]); // Close the circle

		return geometryFactory.createPolygon(pts);
	}
}