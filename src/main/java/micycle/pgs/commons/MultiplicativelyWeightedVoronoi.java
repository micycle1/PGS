package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.overlayng.RingClipper;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.tinspin.index.Index.PointEntryKnn;
import org.tinspin.index.PointMap;
import org.tinspin.index.kdtree.KDTree;

import net.jafama.FastMath;
import processing.core.PVector;

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

		for (PVector vector : vectors) {
			sites.add(new Coordinate(vector.x, vector.y, vector.z));
		}

		Envelope extent = new Envelope(bounds[0], bounds[2], bounds[1], bounds[3]); // NOTE x,x,y,y

		boolean allWeightsDifferent = sites.stream().map(s -> s.z).distinct().count() == sites.size();

		return allWeightsDifferent ? getMWVDFast(sites, extent) : getMWVD(sites, extent);
	}

	/**
	 * Optimised implementation. Reduces redundant geometric boolean operations.
	 * Requires all weights to be different (doesn't handle bisectors arising from
	 * equal pairwise weights to remain simpler).
	 */
	private static List<Geometry> getMWVDFast(List<Coordinate> sites, Envelope extent) {
		// The MWV cell for each site is equivalent to:
		// intersection(all ap_circles containing s)-union(all ap_circles not containing
		// s)

		sites.sort((s1, s2) -> Double.compare(s1.z, s2.z));
		final Geometry extentGeometry = geometryFactory.toGeometry(extent);
		final RingClipper rc = new RingClipper(extent);

		return sites.parallelStream().map(site -> { // NOTE parallel
			// s2 dominates s1 (hence s1 contained in apollo circle)
			List<double[]> inCircleData = new ArrayList<>();
			List<double[]> exCircleData = new ArrayList<>();

			// prepare apollonian circles
			for (int j = 0; j < sites.size(); j++) {
				Coordinate site2 = sites.get(j);
				if (site == site2) {
					continue;
				}
				var circle = calculateWeightedApollonius(site, site2);
				double radius = circle[2];
				if (radius == 0) { // ideally not!
					continue;
				}

				if (site.z < site2.z) { // equivalent to distance(s1,s2) < radius
					inCircleData.add(circle);
				} else {
					exCircleData.add(circle);
				}
			}

			List<Geometry> inCircles = new ArrayList<>();
			List<Geometry> exCircles = new ArrayList<>();

			if (inCircleData.isEmpty()) {
				// if no incircles, cell is simply defined by difference between plane and
				// exCircles
				exCircleData.forEach(c -> exCircles.add(createClippedCircle(c[0], c[1], c[2], rc)));
				var outerDom = UnaryUnionOp.union(exCircles);
				return extentGeometry.difference(outerDom);
			}

			/*
			 * NOTE optimisation: filter out redundant circles that completely cover other
			 * circles to reduce number of required intersection operations.
			 */
			inCircleData.sort((a, b) -> Double.compare(a[2], b[2])); // sort by radius, smallest first
			List<double[]> essentialCircles = new ArrayList<>(inCircleData);
			for (double[] smaller : inCircleData) {
				// Remove any larger circles that completely contain this one
				essentialCircles.removeIf(larger -> {
					if (larger[2] <= smaller[2]) {
						return false; // Skip if not larger
					}

					double dx = smaller[0] - larger[0];
					double dy = smaller[1] - larger[1];
					double distanceSquared = (dx * dx) + (dy * dy);
					return Math.sqrt(distanceSquared) + smaller[2] <= larger[2];
				});
			}

			essentialCircles.forEach(c -> inCircles.add(createClippedCircle(c[0], c[1], c[2], rc)));
			// intersect all inCircles to find the dominant region for this site
			var localDominance = inCircles.stream().reduce((geom1, geom2) -> OverlayNG.overlay(geom1, geom2, OverlayNG.INTERSECTION)).get();

			/*
			 * NOTE optimisation: only circles that intersect the MBC of the localDominance
			 * area need subtracting from it. Will usually reduce the number of outcircles
			 * to union together.
			 */
			MinimumBoundingCircle mbc = new MinimumBoundingCircle(localDominance);
			var mbcP = mbc.getCentre();
			var maxDominanceRegion = new double[] { mbcP.x, mbcP.y, mbc.getRadius() };
//			var maxDominanceRegion = inCircleData.get(0); // smallest circle (but mbc even smaller...)
			exCircleData = exCircleData.stream().filter(c -> {
				double dx = maxDominanceRegion[0] - c[0];
				double dy = maxDominanceRegion[1] - c[1];
				double distanceSquared = (dx * dx) + (dy * dy);
				double radiusSum = maxDominanceRegion[2] + c[2];
				// keep if intersects
				return distanceSquared < radiusSum * radiusSum;
			}).collect(Collectors.toList());

			/*
			 * NOTE optimisation, similar to inCircleData optimisation, but this time,
			 * remove any SMALLER circles that are contained by another.
			 */
			exCircleData = filterContainedCircles(exCircleData);
			exCircleData.forEach(c -> exCircles.add(createClippedCircle(c[0], c[1], c[2], rc)));

			if (exCircles.isEmpty()) {
				return localDominance;
			} else {
				var outerDom = UnaryUnionOp.union(exCircles);
				return localDominance.difference(outerDom);
			}
		}).toList();
	}

	/**
	 * Filters the circles, removing any circles that are fully contained by
	 * another.
	 */
	private static List<double[]> filterContainedCircles(List<double[]> exCircleData) {
		if (exCircleData.size() <= 1) {
			return exCircleData;
		}

		// Sort by radius (descending) to process largest first
		exCircleData.sort((a, b) -> Double.compare(b[2], a[2]));

		// Use boolean array to track removed circles - faster than repeated list
		// operations
		boolean[] isRemoved = new boolean[exCircleData.size()];

		// Check each circle against smaller ones
		for (int i = 0; i < exCircleData.size(); i++) {
			if (isRemoved[i]) {
				continue;
			}

			double[] larger = exCircleData.get(i);
			double largerX = larger[0];
			double largerY = larger[1];
			double largerRadius = larger[2];

			// Only need to check against smaller circles (after i in sorted list)
			for (int j = i + 1; j < exCircleData.size(); j++) {
				if (isRemoved[j]) {
					continue;
				}

				double[] smaller = exCircleData.get(j);
				// Quick radius check
				if (smaller[2] >= largerRadius) {
					continue;
				}

				double dx = smaller[0] - largerX;
				double dy = smaller[1] - largerY;
				double distanceSquared = dx * dx + dy * dy;

				// Check if smaller circle is completely contained
				// distance + smaller_radius <= larger_radius
				// Optimize by comparing squares first
				double radiusDiff = largerRadius - smaller[2];
				if (distanceSquared <= radiusDiff * radiusDiff) {
					isRemoved[j] = true;
				}
			}
		}

		// Build final list in one pass
		List<double[]> filteredCircles = new ArrayList<>();
		for (int i = 0; i < exCircleData.size(); i++) {
			if (!isRemoved[i]) {
				filteredCircles.add(exCircleData.get(i));
			}
		}

		return filteredCircles;
	}

	/**
	 * Vanilla implementation. Handles equal-weighted pairs.
	 */
	private static List<Geometry> getMWVD(List<Coordinate> sites, Envelope extent) {
		final PointMap<Point> tree = KDTree.create(2);
		sites.forEach(s -> {
			Point p = geometryFactory.createPoint(s);
			p.setUserData(s.z);
			tree.insert(new double[] { p.getX(), p.getY() }, p);

		});

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
			double[] circle = calculateWeightedApollonius(s1, s2);
			Coordinate center = new Coordinate(circle[0], circle[1]); // often lies outside bounds
			double radius = circle[2];
			double distance = s1.distance(center);
			localDominanceCircle = createClippedCircle(center.x, center.y, radius, null);
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

	/**
	 * Find the apollonius circle representing the bisector between two sites.
	 */
	private static double[] calculateWeightedApollonius(Coordinate s1, Coordinate s2) {
		final double w1 = s1.z;
		final double w2 = s2.z;
		final double s1x = s1.x;
		final double s1y = s1.y;
		final double s2x = s2.x;
		final double s2y = s2.y;

		final double den = 1.0 / (w1 * w1 - w2 * w2);
		final double cx = (w1 * w1 * s2x - w2 * w2 * s1x) * den;
		final double cy = (w1 * w1 * s2y - w2 * w2 * s1y) * den;
		final double d = Math.sqrt(((s1x - s2x) * (s1x - s2x) + (s1y - s2y) * (s1y - s2y)));
		// NOTE r can be huge (as circle may tend towards straight line)
		final double r = Math.abs(w1 * w2 * d * den);

		return new double[] { cx, cy, r };
	}

	private static Polygon createClippedCircle(double x, double y, double r, RingClipper rc) {
		final double maxDeviation = 0.49;
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

		// NOTE clip the circle to bounds now. slightly speeds up 2d boolean ops later
		// on.
		return geometryFactory.createPolygon(rc == null ? pts : rc.clip(pts));
	}
}