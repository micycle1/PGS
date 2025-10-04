package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;

/**
 * Implementation of the Bayazit convex decomposition algorithm for simple
 * polygons.
 * <p>
 * This algorithm is a O(nr) complexity algorithm where n is the number of input
 * vertices and r is the number of output convex polygons. This algorithm can
 * achieve optimal decompositions, however this is not guaranteed.
 *
 * @author William Bittle
 * @author Refactored for JTS by Michael Carleton
 * @see <a href= "https://mpen.ca/406/bayazit">Mark Bayazits Algorithm </a>
 */
public class PolygonDecomposition {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	private PolygonDecomposition() {
	}

	public static List<Polygon> decompose(Polygon polygon) {
		// Use only the exterior ring; algorithm expects a simple polygon boundary.
		Coordinate[] closed = polygon.getExteriorRing().getCoordinates(); // closed ring
		boolean ccw = Orientation.isCCW(closed);

		Coordinate[] open = openRing(closed); // remove duplicate closing coord
		if (!ccw) {
			reverse(open);
		}
		return decompose(open);
	}

	private static List<Polygon> decompose(Coordinate... points) {
		if (points == null) {
			throw new NullPointerException("Points are null.");
		}
		int size = points.length;
		if (size < 3) {
			throw new IllegalArgumentException("Points have invalid size (<3 open vertices).");
		}

		List<Coordinate> polygon = new ArrayList<>();
		Collections.addAll(polygon, points);

		List<Polygon> polygons = new ArrayList<>();
		decomposePolygon(polygon, polygons);
		return polygons;
	}

	private static void decomposePolygon(final List<Coordinate> polygon, final List<Polygon> polygons) {
		final int size = polygon.size();

		Coordinate upperIntersection = null;
		Coordinate lowerIntersection = null;
		double upperDistance = Double.MAX_VALUE;
		double lowerDistance = Double.MAX_VALUE;
		double closestDistance = Double.MAX_VALUE;
		int upperIndex = 0;
		int lowerIndex = 0;
		int closestIndex = 0;

		final List<Coordinate> lower = new ArrayList<>();
		final List<Coordinate> upper = new ArrayList<>();

		for (int i = 0; i < size; i++) {
			final Coordinate p = polygon.get(i);
			final Coordinate p0 = polygon.get(i - 1 < 0 ? size - 1 : i - 1);
			final Coordinate p1 = polygon.get(i + 1 == size ? 0 : i + 1);

			if (isReflex(p0, p, p1)) {
				for (int j = 0; j < size; j++) {
					final Coordinate q = polygon.get(j);
					final Coordinate q0 = polygon.get(j - 1 < 0 ? size - 1 : j - 1);
					final Coordinate q1 = polygon.get(j + 1 == size ? 0 : j + 1);

					// extend the previous edge: infinite lines p0-p with q-q0
					if (left(p0, p, q) && rightOn(p0, p, q0)) {
						final Coordinate s = lineLineIntersection(p0, p, q, q0);
						if (s != null && right(p1, p, s)) {
							final double dist = p.distanceSq(s);
							if (dist < lowerDistance) {
								lowerDistance = dist;
								lowerIntersection = s;
								lowerIndex = j;
							}
						}
					}

					// extend the next edge: infinite lines p1-p with q-q1
					if (left(p1, p, q1) && rightOn(p1, p, q)) {
						final Coordinate s = lineLineIntersection(p1, p, q, q1);
						if (s != null && left(p0, p, s)) {
							final double dist = p.distanceSq(s);
							if (dist < upperDistance) {
								upperDistance = dist;
								upperIntersection = s;
								upperIndex = j;
							}
						}
					}
				}

				if (lowerIndex == (upperIndex + 1) % size) {
					// create a Steiner point
					final Coordinate s = midpoint(upperIntersection, lowerIntersection);
					// guard in case intersections weren’t found due to degeneracy
					if (s == null) {
						// fall back to skipping this reflex (should be rare)
						return;
					}

					if (i < upperIndex) {
						lower.addAll(polygon.subList(i, upperIndex + 1));
						lower.add(s);
						upper.add(s);
						if (lowerIndex != 0) {
							upper.addAll(polygon.subList(lowerIndex, size));
						}
						upper.addAll(polygon.subList(0, i + 1));
					} else {
						if (i != 0) {
							lower.addAll(polygon.subList(i, size));
						}
						lower.addAll(polygon.subList(0, upperIndex + 1));
						lower.add(s);
						upper.add(s);
						upper.addAll(polygon.subList(lowerIndex, i + 1));
					}
				} else {
					if (lowerIndex > upperIndex) {
						upperIndex += size;
					}

					closestIndex = lowerIndex;
					for (int j = lowerIndex; j <= upperIndex; j++) {
						final int jmod = j % size;
						final Coordinate q = polygon.get(jmod);
						if (coordsEqual(q, p) || coordsEqual(q, p0) || coordsEqual(q, p1)) {
							continue;
						}
						final double dist = p.distanceSq(q);
						if (dist < closestDistance) {
							if (isVisible(polygon, i, jmod)) {
								closestDistance = dist;
								closestIndex = jmod;
							}
						}
					}

					if (i < closestIndex) {
						lower.addAll(polygon.subList(i, closestIndex + 1));
						if (closestIndex != 0) {
							upper.addAll(polygon.subList(closestIndex, size));
						}
						upper.addAll(polygon.subList(0, i + 1));
					} else {
						if (i != 0) {
							lower.addAll(polygon.subList(i, size));
						}
						lower.addAll(polygon.subList(0, closestIndex + 1));
						upper.addAll(polygon.subList(closestIndex, i + 1));
					}
				}

				if (lower.size() < upper.size()) {
					decomposePolygon(lower, polygons);
					decomposePolygon(upper, polygons);
				} else {
					decomposePolygon(upper, polygons);
					decomposePolygon(lower, polygons);
				}
				return;
			}
		}

		// No reflex vertices => polygon is convex
		if (polygon.size() < 3) {
			return;
		}
		final Coordinate[] vertices = polygon.toArray(new Coordinate[0]);
		final Coordinate[] jtsCoords = new Coordinate[vertices.length + 1];
		for (int j = 0; j < vertices.length; j++) {
			jtsCoords[j] = new Coordinate(vertices[j].x, vertices[j].y);
		}
		jtsCoords[vertices.length] = jtsCoords[0];
		polygons.add(GEOM_FACTORY.createPolygon(jtsCoords));
	}

	private static boolean isReflex(final Coordinate p0, final Coordinate p, final Coordinate p1) {
		return right(p1, p0, p);
	}

	private static boolean left(final Coordinate a, final Coordinate b, final Coordinate p) {
		return Orientation.index(a, b, p) > 0;
	}

	private static boolean leftOn(final Coordinate a, final Coordinate b, final Coordinate p) {
		return Orientation.index(a, b, p) >= 0;
	}

	private static boolean right(final Coordinate a, final Coordinate b, final Coordinate p) {
		return Orientation.index(a, b, p) < 0;
	}

	private static boolean rightOn(final Coordinate a, final Coordinate b, final Coordinate p) {
		return Orientation.index(a, b, p) <= 0;
	}

	private static Coordinate lineLineIntersection(final Coordinate a1, final Coordinate a2, final Coordinate b1, final Coordinate b2) {
		final LineSegment la = new LineSegment(a1, a2);
		final LineSegment lb = new LineSegment(b1, b2);
		// lineIntersection returns the intersection of the supporting lines (not
		// segments), or null if parallel/collinear
		return la.lineIntersection(lb);
	}

	private static boolean isVisible(final List<Coordinate> polygon, final int i, final int j) {
		final int s = polygon.size();
		final Coordinate iv0 = polygon.get(i == 0 ? s - 1 : i - 1);
		final Coordinate iv = polygon.get(i);
		final Coordinate iv1 = polygon.get(i + 1 == s ? 0 : i + 1);

		final Coordinate jv0 = polygon.get(j == 0 ? s - 1 : j - 1);
		final Coordinate jv = polygon.get(j);
		final Coordinate jv1 = polygon.get(j + 1 == s ? 0 : j + 1);

		if (isReflex(iv0, iv, iv1)) {
			if (leftOn(iv, iv0, jv) && rightOn(iv, iv1, jv)) {
				return false;
			}
		} else {
			if (rightOn(iv, iv1, jv) || leftOn(iv, iv0, jv)) {
				return false;
			}
		}
		if (isReflex(jv0, jv, jv1)) {
			if (leftOn(jv, jv0, iv) && rightOn(jv, jv1, iv)) {
				return false;
			}
		} else {
			if (rightOn(jv, jv1, iv) || leftOn(jv, jv0, iv)) {
				return false;
			}
		}

		final LineSegment segA = new LineSegment(iv, jv);
		for (int k = 0; k < s; k++) {
			final int ki1 = k + 1 == s ? 0 : k + 1;
			if (k == i || k == j || ki1 == i || ki1 == j) {
				continue;
			}
			final Coordinate k1 = polygon.get(k);
			final Coordinate k2 = polygon.get(ki1);

			final LineSegment segB = new LineSegment(k1, k2);
			final Coordinate in = segA.intersection(segB);
			if (in != null) {
				return false;
			}
		}

		return true;
	}

	private static void reverse(final Coordinate[] pts) {
		for (int i = 0, j = pts.length - 1; i < j; i++, j--) {
			final Coordinate tmp = pts[i];
			pts[i] = pts[j];
			pts[j] = tmp;
		}
	}

	private static Coordinate[] openRing(final Coordinate[] coords) {
		if (coords.length >= 2 && coords[0].equals2D(coords[coords.length - 1])) {
			final Coordinate[] open = new Coordinate[coords.length - 1];
			System.arraycopy(coords, 0, open, 0, coords.length - 1);
			return open;
		}
		return coords;
	}

	private static Coordinate midpoint(final Coordinate a, final Coordinate b) {
		if (a == null || b == null) {
			return null;
		}
		return new Coordinate((a.x + b.x) * 0.5, (a.y + b.y) * 0.5);
	}

	private static boolean coordsEqual(final Coordinate a, final Coordinate b) {
		return a.equals2D(b);
	}
}
