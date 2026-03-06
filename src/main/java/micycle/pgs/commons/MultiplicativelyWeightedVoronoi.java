package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.overlayng.RingClipper;
import org.locationtech.jts.operation.union.UnaryUnionOp;
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

	private static final double HALF_PLANE_EPS = 1e-12;
	/**
	 * Approximates a near-linear Apollonius arc by a perpendicular bisector when
	 * its maximum deviation from straightness across the bounds is below this
	 * fraction of the bounds diagonal.
	 */
	private static final double BISECTOR_APPROX_MAX_DEVIATION_RATIO = 0.01; // 1%
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
		return getMWVDFast(sites, extent);
	}

	/**
	 * Optimised implementation.
	 *
	 * Handles: - unequal-weight pairs via Apollonius circles - equal-weight pairs
	 * via clipped bisector half-planes
	 * <p>
	 * Cell(site) = extent ∩ all equal-weight bisector half-planes containing site ∩
	 * all Apollonius circles that contain site \ union(all Apollonius circles that
	 * exclude site)
	 */
	private static List<Geometry> getMWVDFast(List<Coordinate> sites, Envelope extent) {
		sites.sort((s1, s2) -> Double.compare(s1.z, s2.z));

		final Geometry extentGeometry = geometryFactory.toGeometry(extent);
		final RingClipper rc = new RingClipper(extent);

		return sites.parallelStream().map(site -> {
			List<double[]> inCircleData = new ArrayList<>();
			List<double[]> exCircleData = new ArrayList<>();

			// Equal-weight constraints are linear half-planes, clipped against extent.
			// We keep them as a convex polygon in coordinate form to avoid expensive
			// polygonisation/splitting.
			List<Coordinate> bisectorCell = null;

			for (Coordinate other : sites) {
				if (site == other) {
					continue;
				}

				int cmp = Double.compare(site.z, other.z);

				// Exact equal-weight case: true bisector
				if (cmp == 0) {
					if (bisectorCell == null) {
						bisectorCell = extentPolygon(extent);
					}
					bisectorCell = clipPolygonToBisectorHalfPlane(bisectorCell, site, other);

					if (bisectorCell.size() < 3) {
						return null;
					}
					continue;
				}

				double[] circle = calculateWeightedApollonius(site, other);
				double radius = circle[2];
				if (radius == 0 || !Double.isFinite(radius)) {
					continue;
				}

				// Near-equal / huge-radius case: approximate the circle by the bisector
				if (shouldApproximateAsBisector(radius, extent)) {
					if (bisectorCell == null) {
						bisectorCell = extentPolygon(extent);
					}
					bisectorCell = clipPolygonToBisectorHalfPlane(bisectorCell, site, other);

					if (bisectorCell.size() < 3) {
						return null;
					}
					continue;
				}

				if (cmp < 0) {
					// site is weaker, so its dominance lies INSIDE the circle
					inCircleData.add(circle);
				} else {
					// site is stronger, so its dominance lies OUTSIDE the circle
					exCircleData.add(circle);
				}
			}

			// Fast trivial case: no bisectors and no in-circles
			if (bisectorCell == null && inCircleData.isEmpty()) {
				if (exCircleData.isEmpty()) {
					return extentGeometry;
				}

				exCircleData = filterContainedCircles(exCircleData);
				List<Geometry> exCircles = new ArrayList<>(exCircleData.size());
				exCircleData.forEach(c -> exCircles.add(createClippedCircle(c[0], c[1], c[2], rc)));

				Geometry outerDom = UnaryUnionOp.union(exCircles);
				return extentGeometry.difference(outerDom);
			}

			Geometry localDominance = (bisectorCell == null) ? extentGeometry : toPolygon(bisectorCell);
			if (localDominance.isEmpty()) {
				return null;
			}

			/*
			 * Optimisation: remove any larger in-circle that completely contains a smaller
			 * one, because intersecting with the larger adds no information.
			 */
			if (!inCircleData.isEmpty()) {
				List<double[]> essentialInCircles = filterRedundantInCircles(inCircleData);

				for (double[] c : essentialInCircles) {
					Geometry inCircle = createClippedCircle(c[0], c[1], c[2], rc);

					if (!localDominance.getEnvelopeInternal().intersects(inCircle.getEnvelopeInternal())) {
						return null;
					}

					localDominance = OverlayNG.overlay(localDominance, inCircle, OverlayNG.INTERSECTION);
					if (localDominance.isEmpty()) {
						return null;
					}
				}
			}

			if (exCircleData.isEmpty()) {
				return localDominance;
			}

			/*
			 * Only ex-circles that can intersect the current local dominance region matter.
			 * Use the MBC as a cheap coarse filter.
			 */
			MinimumBoundingCircle mbc = new MinimumBoundingCircle(localDominance);
			Coordinate mbcC = mbc.getCentre();
			double mbcR = mbc.getRadius();

			exCircleData = exCircleData.stream().filter(c -> {
				double dx = mbcC.x - c[0];
				double dy = mbcC.y - c[1];
				double radiusSum = mbcR + c[2];
				return dx * dx + dy * dy < radiusSum * radiusSum;
			}).collect(Collectors.toList());

			exCircleData = filterContainedCircles(exCircleData);

			if (exCircleData.isEmpty()) {
				return localDominance;
			}

			List<Geometry> exCircles = new ArrayList<>(exCircleData.size());
			exCircleData.forEach(c -> exCircles.add(createClippedCircle(c[0], c[1], c[2], rc)));

			Geometry outerDom = UnaryUnionOp.union(exCircles);
			return localDominance.difference(outerDom);

		}).filter(g -> g != null && !g.isEmpty()).toList();
	}

	private static List<double[]> filterRedundantInCircles(List<double[]> inCircleData) {
		if (inCircleData.size() <= 1) {
			return inCircleData;
		}

		inCircleData.sort((a, b) -> Double.compare(a[2], b[2])); // smallest first
		List<double[]> essentialCircles = new ArrayList<>(inCircleData);

		for (double[] smaller : inCircleData) {
			essentialCircles.removeIf(larger -> {
				if (larger == smaller || larger[2] <= smaller[2]) {
					return false;
				}
				return circleContains(larger, smaller);
			});
		}

		return essentialCircles;
	}

	private static boolean circleContains(double[] outer, double[] inner) {
		double dx = inner[0] - outer[0];
		double dy = inner[1] - outer[1];
		double radiusDiff = outer[2] - inner[2];
		return radiusDiff >= 0 && dx * dx + dy * dy <= radiusDiff * radiusDiff;
	}

	private static List<Coordinate> extentPolygon(Envelope extent) {
		List<Coordinate> poly = new ArrayList<>(4);
		poly.add(new Coordinate(extent.getMinX(), extent.getMinY()));
		poly.add(new Coordinate(extent.getMaxX(), extent.getMinY()));
		poly.add(new Coordinate(extent.getMaxX(), extent.getMaxY()));
		poly.add(new Coordinate(extent.getMinX(), extent.getMaxY()));
		return poly;
	}

	/**
	 * Returns true if the Apollonius circle is so large that, across the current
	 * bounds, it is effectively indistinguishable from a straight line.
	 *
	 * We measure this using the sagitta (max arc deviation from its chord) across
	 * the extent diagonal.
	 */
	private static boolean shouldApproximateAsBisector(double radius, Envelope extent) {
		if (!Double.isFinite(radius) || radius <= 0) {
			return false;
		}

		double span = extent.getDiameter();
		if (span <= HALF_PLANE_EPS) {
			return true;
		}

		double halfChord = span * 0.5;

		// If the radius isn't even big enough to support such a chord, then the arc is
		// definitely not near-linear over the full bounds.
		if (radius <= halfChord) {
			return false;
		}

		// Sagitta = r - sqrt(r^2 - (L/2)^2)
		// Use the approximation for very large radii to avoid cancellation error.
		double deviation;
		double t = halfChord / radius;
		if (t < 1e-3) {
			deviation = (halfChord * halfChord) / (2.0 * radius); // ~ L^2 / (8r)
		} else {
			deviation = radius - Math.sqrt(radius * radius - halfChord * halfChord);
		}

		return deviation <= BISECTOR_APPROX_MAX_DEVIATION_RATIO * span;
	}

	/**
	 * Clips a convex polygon against the half-plane of the perpendicular bisector
	 * that contains s1, i.e. points p such that dist(p, s1) <= dist(p, s2).
	 */
	private static List<Coordinate> clipPolygonToBisectorHalfPlane(List<Coordinate> polygon, Coordinate s1, Coordinate s2) {
		if (polygon.isEmpty()) {
			return polygon;
		}

		double nx = s2.x - s1.x;
		double ny = s2.y - s1.y;

		// Degenerate coincident case: undefined bisector, leave polygon unchanged.
		if (Math.abs(nx) < HALF_PLANE_EPS && Math.abs(ny) < HALF_PLANE_EPS) {
			return polygon;
		}

		double mx = (s1.x + s2.x) * 0.5;
		double my = (s1.y + s2.y) * 0.5;

		List<Coordinate> out = new ArrayList<>(polygon.size() + 2);

		Coordinate prev = polygon.get(polygon.size() - 1);
		double prevSide = bisectorSide(prev, mx, my, nx, ny);
		boolean prevInside = prevSide <= HALF_PLANE_EPS;

		for (Coordinate curr : polygon) {
			double currSide = bisectorSide(curr, mx, my, nx, ny);
			boolean currInside = currSide <= HALF_PLANE_EPS;

			if (prevInside && currInside) {
				addDistinct(out, curr);
			} else if (prevInside) {
				addDistinct(out, bisectorIntersection(prev, curr, prevSide, currSide));
			} else if (currInside) {
				addDistinct(out, bisectorIntersection(prev, curr, prevSide, currSide));
				addDistinct(out, curr);
			}

			prev = curr;
			prevSide = currSide;
			prevInside = currInside;
		}

		if (out.size() > 1 && samePoint(out.get(0), out.get(out.size() - 1))) {
			out.remove(out.size() - 1);
		}

		return out;
	}

	private static double bisectorSide(Coordinate p, double mx, double my, double nx, double ny) {
		return (p.x - mx) * nx + (p.y - my) * ny;
	}

	private static Coordinate bisectorIntersection(Coordinate a, Coordinate b, double fa, double fb) {
		double denom = fa - fb;
		if (Math.abs(denom) < HALF_PLANE_EPS) {
			return new Coordinate(b);
		}

		double t = fa / denom;
		t = Math.max(0.0, Math.min(1.0, t));

		return new Coordinate(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
	}

	private static void addDistinct(List<Coordinate> coords, Coordinate c) {
		if (c == null) {
			return;
		}
		if (coords.isEmpty() || !samePoint(coords.get(coords.size() - 1), c)) {
			coords.add(c);
		}
	}

	private static boolean samePoint(Coordinate a, Coordinate b) {
		return a.distanceSq(b) <= HALF_PLANE_EPS * HALF_PLANE_EPS;
	}

	private static Polygon toPolygon(List<Coordinate> coords) {
		if (coords == null || coords.size() < 3) {
			return geometryFactory.createPolygon();
		}

		List<Coordinate> clean = new ArrayList<>(coords.size());
		for (Coordinate c : coords) {
			addDistinct(clean, c);
		}

		if (clean.size() < 3) {
			return geometryFactory.createPolygon();
		}

		if (samePoint(clean.get(0), clean.get(clean.size() - 1))) {
			clean.remove(clean.size() - 1);
		}

		if (clean.size() < 3) {
			return geometryFactory.createPolygon();
		}

		Coordinate[] ring = new Coordinate[clean.size() + 1];
		for (int i = 0; i < clean.size(); i++) {
			ring[i] = clean.get(i);
		}
		ring[clean.size()] = new Coordinate(clean.get(0));

		return geometryFactory.createPolygon(ring);
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