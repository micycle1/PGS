package micycle.pgs.utility;

import static java.lang.Math.floorMod;

import java.util.Arrays;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;

/**
 * Computes the Minimum Bounding Triangle (MBT) for the points in a Geometry.
 * The MBT is the smallest triangle which covers all the input points (this is
 * also known as the Smallest Enclosing Triangle).
 * <p>
 * The implementation of the algorithm is based on O'Rourke's 'An optimal
 * algorithm for finding minimal enclosing triangles' and Klee & Laskowski's
 * 'Finding the smallest triangles containing a given convex polygon'. O'Rourke
 * provides a θ(n) algorithm for finding the minimal enclosing triangle of a 2D
 * convex polygon with n vertices. However, the overall complexity for the
 * computation is O(nlog(n)) because a convex hull must first be computed for
 * the input geometry.
 * 
 * @author Python implementation by Charlie Marsh
 * @author Java port by Michael Carleton
 *
 */
public class MinimumBoundingTriangle {

	// TODO copy docs from openCV implementation
	// port of https://github.com/crm416/point-location/blob/master/min_triangle.py

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));
	private static final double EPSILON = 0.01; // Account for floating-point errors

	private final int n;
	private final Coordinate[] points;

	/**
	 * Creates a new instance of a Maximum Inscribed Triangle computation.
	 * 
	 * @param shape an areal geometry
	 */
	public MinimumBoundingTriangle(Geometry shape) {
		shape = shape.convexHull();
		points = Arrays.copyOfRange(shape.getCoordinates(), 0, shape.getCoordinates().length - 1); // treat coordinates as unclosed
		n = points.length;
	}

	/**
	 * Gets a geometry which represents the Minimium Bounding Triangle.
	 * 
	 * @return a triangle Geometry representing the Minimum Bounding Triangle
	 */
	public Geometry getTriangle() {
		int a = 1;
		int b = 2;

		double minArea = Double.MAX_VALUE;
		Polygon minAreaTriangle = null;

		for (int i = 0; i < n; i++) {
			triangleForIndex tForIndex = new triangleForIndex(i, a, b);
			Polygon triangle = tForIndex.triangle;
			a = tForIndex.aOut;
			b = tForIndex.bOut;
			if (triangle != null) {
				double area = triangle.getArea();
				if (minAreaTriangle == null || area < minArea) {
					minArea = area;
					minAreaTriangle = triangle;
				}
			}
		}

		return minAreaTriangle;
	}

	private Side side(final int i) {
		return new Side(points[floorMod(i - 1, n)], points[i % n]);
	}

	/**
	 * Checks that a midpoint touches the polygon on the appropriate side.
	 */
	private boolean validateMidpoint(Coordinate midpoint, int index) {
		Side s = side(index);

		if (s.vertical) {
			if (midpoint.x != s.p1.x) {
				return false;
			}
			double maxY = Math.max(s.p1.y, s.p2.y) + EPSILON;
			double minY = Math.min(s.p1.y, s.p2.y) - EPSILON;
			if (!(midpoint.y <= maxY && midpoint.y >= minY)) {
				return false;
			}
			return true;
		} else {
			double maxX = Math.max(s.p1.x, s.p2.x) + EPSILON;
			double minX = Math.min(s.p1.x, s.p2.x) - EPSILON;
			// Must touch polygon
			if (!(midpoint.x <= maxX && midpoint.x >= minX)) {
				return false;
			}

			if ((s.atX(midpoint.x).distance(midpoint) > 0.01)) {
				return false;
			}

			return true;
		}
	}

	private static Coordinate midpoint(Coordinate a, Coordinate b) {
		return new Coordinate((a.x + b.x) / 2, (a.y + b.y) / 2);
	}

	/**
	 * Computes the minimal triangle with edge C flush to vertex c.
	 * 
	 * Abstracted into class (during Java port) to better structure the many
	 * methods.
	 */
	private class triangleForIndex {

		// return values
		final int aOut, bOut;
		final Polygon triangle;

		private final Side sideC;
		private final int c;
		private Side sideA, sideB;

		triangleForIndex(int c, int a, int b) {
			a = Math.max(a, c + 1) % n;
			b = Math.max(b, c + 2) % n;
			sideC = side(c);
			this.c = c;

			while (onLeftChain(b)) { // Increment b while low
				b = (b + 1) % n;
			}

			while (dist(b, sideC) > dist(a, sideC)) { // Increment a if low, b if high
				int[] ab = incrementLowHigh(a, b, c);
				a = ab[0];
				b = ab[1];
			}

			while (tangency(a, b)) { // Search for b tangency
				b = (b + 1) % n;
			}

			Coordinate gammaB = gamma(points[b], side(a), sideC);
			// Adjust if necessary
			if (low(a, b, c, gammaB) || dist(b, sideC) < dist((a - 1) % n, sideC)) {
				sideB = side(b);
				sideA = side(a);
				sideB = new Side(sideC.intersection(sideB), sideA.intersection(sideB));

				if (dist(sideB.midpoint(), sideC) < dist(floorMod(a - 1, n), sideC)) {
					Coordinate gammaA = gamma(points[floorMod(a - 1, n)], sideB, sideC);
					sideA = new Side(gammaA, points[floorMod(a - 1, n)]);
				}
			} else {
				gammaB = gamma(points[b], side(a), sideC);
				sideB = new Side(gammaB, points[b]);
				sideA = new Side(gammaB, points[floorMod(a - 1, n)]);
			}

			// Calculate final intersections
			final Coordinate vertexA = sideC.intersection(sideB);
			final Coordinate vertexB = sideC.intersection(sideA);
			final Coordinate vertexC = sideA.intersection(sideB);

			// Check if triangle is valid local minimum
			if (!isValidTriangle(vertexA, vertexB, vertexC, a, b, c)) {
				triangle = null;
			} else {
				triangle = GEOM_FACTORY.createPolygon(new Coordinate[] { vertexA, vertexB, vertexC, vertexA });
			}

			aOut = a;
			bOut = b;
		}

		/**
		 * Computes the distance from the point (specified by its index) to the side.
		 */
		private double dist(int point, Side side) {
			return side.distance(points[floorMod(point, points.length)]);
		}

		/**
		 * Computes the distance from the point to the side.
		 */
		private double dist(Coordinate point, Side side) {
			return side.distance(point);
		}

		/**
		 * Calculate the point on 'on' that is twice as far from 'base' as 'point'.
		 */
		private Coordinate gamma(Coordinate point, Side on, Side base) {
			Coordinate intersection = on.intersection(base);
			if (intersection != null) {
				double dist = 2 * dist(point, base);
				// Calculate differential change in distance
				if (on.vertical) {
					double ddist = dist(new Coordinate(intersection.x, intersection.y + 1), base);
					Coordinate guess = new Coordinate(intersection.x, intersection.y + dist / ddist);
					if (ccw(base.p1, base.p2, guess) != ccw(base.p1, base.p2, point)) {
						guess = new Coordinate(intersection.x, intersection.y - dist / ddist);
					}
					return guess;
				} else {
					double ddist = dist(on.atX(intersection.x + 1), base);
					Coordinate guess = on.atX(intersection.x + dist / ddist);
					if (ccw(base.p1, base.p2, guess) != ccw(base.p1, base.p2, point)) {
						guess = on.atX(intersection.x - dist / ddist);
					}
					return guess;
				}
			}
			return intersection;
		}

		private boolean high(int a, int b, int c, Coordinate gammaB) {
			// Test if two adjacent vertices are on same side of line (implies tangency)
			if (ccw(gammaB, points[b], points[floorMod(b - 1, n)]) == ccw(gammaB, points[b], points[(b + 1) % n])) {
				return false;
			}

			// Test if Gamma and B are on same side of line from adjacent vertices
			if (ccw(points[floorMod(b - 1, n)], points[(b + 1) % n], gammaB) == ccw(points[floorMod(b - 1, n)], points[(b + 1) % n],
					points[b])) {
				return dist(gammaB, sideC) > dist(b, sideC);
			} else {
				return false;
			}
		}

		private boolean low(int a, int b, int c, Coordinate gammaB) {
			// Test if two adjacent vertices are on same side of line (implies tangency)
			if (ccw(gammaB, points[b], points[floorMod(b - 1, n)]) == ccw(gammaB, points[b], points[(b + 1) % n])) {
				return false;
			}

			// Test if Gamma and B are on same side of line from adjacent vertices
			if (ccw(points[floorMod(b - 1, n)], points[(b + 1) % n], gammaB) == ccw(points[floorMod(b - 1, n)], points[(b + 1) % n],
					points[b])) {
				return false;
			} else {
				return dist(gammaB, sideC) > dist(b, sideC);
			}
		}

		private boolean onLeftChain(int b) {
			return dist((b + 1) % n, sideC) >= dist(b, sideC);
		}

		/**
		 * @return [a, b] tuple
		 */
		private int[] incrementLowHigh(int a, int b, int c) {
			Coordinate gammaA = gamma(points[a], side(a), sideC);

			if (high(a, b, c, gammaA)) {
				b = (b + 1) % n;
			} else {
				a = (a + 1) % n;
			}
			return new int[] { a, b };
		}

		private boolean tangency(int a, int b) {
			Coordinate gammaB = gamma(points[b], side(a), sideC);
			return dist(b, sideC) >= dist((a - 1) % n, sideC) && high(a, b, c, gammaB);
		}

		/**
		 * Tests whether the line formed by A, B, and C is ccw.
		 */
		private boolean ccw(Coordinate a, Coordinate b, Coordinate c) {
			return (b.x - a.x) * (c.y - a.y) > (b.y - a.y) * (c.x - a.x);
		}

		/**
		 * Checks that a triangle composed of the given vertices is a valid local
		 * minimum (entails that all midpoints of the triangle should touch the
		 * polygon).
		 */
		private boolean isValidTriangle(Coordinate vertexA, Coordinate vertexB, Coordinate vertexC, int a, int b, int c) {
			if (vertexA == null || vertexB == null || vertexC == null) {
				return false;
			}
			Coordinate midpointA = midpoint(vertexC, vertexB);
			Coordinate midpointB = midpoint(vertexA, vertexC);
			Coordinate midpointC = midpoint(vertexA, vertexB);
			return (validateMidpoint(midpointA, a) && validateMidpoint(midpointB, b) && validateMidpoint(midpointC, c));
		}
	}

	private class Side {

		final Coordinate p1, p2;
		final double slope, intercept;
		final boolean vertical;

		Side(Coordinate p1, Coordinate p2) {
			this.p1 = p1;
			this.p2 = p2;
			slope = (p2.y - p1.y) / (p2.x - p1.x);
			intercept = p1.y - slope * p1.x;
			vertical = p1.x == p2.x;
		}

		private double sqrDistance(Coordinate p) {
			double numerator = (p2.x - p1.x) * (p1.y - p.y) - (p1.x - p.x) * (p2.y - p1.y);
			numerator *= numerator;
			double denominator = (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y);
			return numerator / denominator;
		}

		/**
		 * @return Returns the distance of p from the line
		 */
		private double distance(Coordinate p) {
			return Math.sqrt(sqrDistance(p));
		}

		private Coordinate atX(double x) {
			if (vertical) {
				return p1; // NOTE rather return null
			}
			return new Coordinate(x, slope * x + intercept);
		}

		private Coordinate intersection(Side that) {
			if (that.slope == slope) {
				return null;
			}

			if (vertical) {
				return that.atX(p1.x);
			} else if (that.vertical) {
				return atX(that.p1.x);
			}

			double x = (intercept - that.intercept) / (that.slope - slope);
			return atX(x);
		}

		private Coordinate midpoint() {
			return new Coordinate((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
		}
	}

}
