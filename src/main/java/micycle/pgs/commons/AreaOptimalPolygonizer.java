package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;

/**
 * <p>
 * Computes approximate area-optimal (minimum- and maximum-area) polygonizations
 * of a given set of 2D points. The approach combines a triangle-based
 * constructive heuristic with a triangle-swap local-search procedure.
 * </p>
 *
 * <p>
 * Based on Natanael Ramos, Raí C. de Jesus, Pedro J. de Rezende, Cid C. de
 * Souza, and Fábio L. Usberti. "Triangle-Based Heuristics for Area Optimal
 * Polygonizations".
 * </p>
 *
 * @author Michael Carleton
 */
public final class AreaOptimalPolygonizer {

	public enum AreaObjective {
		MAXIMIZE, MINIMIZE
	}

	static LineIntersector INTERSECTOR = new RobustLineIntersector();

	/**
	 * Area-optimal-heuristic polygonization using default settings.
	 *
	 * <p>
	 * Constructs an approximate area-optimal polygon visiting all distinct input
	 * points using the proximity heuristic. Duplicates and nulls in
	 * {@code inputPoints} are ignored. The pivot is chosen automatically (the first
	 * unique point in iteration order). Triangle-swap local-search is applied by
	 * default and a new {@link GeometryFactory} is used to construct the resulting
	 * polygon.
	 * </p>
	 *
	 * @param inputPoints list of points (duplicates allowed; duplicates are
	 *                    removed, preserving iteration order)
	 * @param objective   area objective (MAXIMIZE or MINIMIZE)
	 * @return a simple, non-degenerate {@link Polygon} visiting the given points
	 */
	public static Polygon polygonize(List<Coordinate> inputPoints, AreaObjective objective) {
		return polygonize(inputPoints, null, objective, true, new GeometryFactory());
	}

	/**
	 * <p>
	 * Constructs an approximate area-optimal polygon visiting all distinct input
	 * points using a heuristic. Optionally, the polygon may be refined by applying
	 * the triangle-swap local-search (SLS) procedure.
	 * </p>
	 *
	 * <p>
	 * Duplicates and nulls in inputPoints are ignored. If pivot is null, the first
	 * unique point (in input iteration order) is used as the pivot. The returned
	 * polygon is oriented CCW and is guaranteed to be simple and non-degenerate
	 * (otherwise an exception is thrown).
	 * </p>
	 *
	 * @param inputPoints                list of input points (duplicates allowed;
	 *                                   duplicates are removed preserving order)
	 * @param pivot                      optional pivot; if null, the first unique
	 *                                   point is used
	 * @param objective                  area objective (MAXIMIZE or MINIMIZE)
	 * @param useTriangleSwapLocalSearch if true, apply triangle-swap local-search
	 *                                   as a post-optimization step
	 * @param gf                         geometry factory used to construct the
	 *                                   resulting polygon
	 * @return a simple, non-degenerate Polygon visiting the given points
	 */
	public static Polygon polygonize(Collection<Coordinate> inputPoints, Coordinate pivot, AreaObjective objective, boolean postOptimize, GeometryFactory gf) {
		Objects.requireNonNull(inputPoints, "inputPoints");
		Objects.requireNonNull(objective, "objective");
		Objects.requireNonNull(gf, "gf");

		List<Coordinate> pts = dedupePreserveOrder(inputPoints);
		if (pts.size() < 3) {
			throw new IllegalArgumentException("Need at least 3 distinct points.");
		}

		Coordinate piv = (pivot != null) ? pivot : pts.get(0);

		// Sort by non-decreasing distance to pivot
		List<Coordinate> sorted = new ArrayList<>(pts);
		sorted.sort(Comparator.comparingDouble((Coordinate c) -> c.distanceSq(piv)).thenComparingDouble(c -> c.x).thenComparingDouble(c -> c.y));

		// Build initial non-degenerate polygon (handles initial collinearity)
		List<Coordinate> ring = buildInitialRingHandlingCollinear(sorted);

		// Track used points
		Set<Coordinate> used = new HashSet<>(ring.size() * 2);
		for (Coordinate c : ring) {
			used.add(c);
		}

		// Reusable candidate buffers (avoid allocations per point)
		final int maxN = pts.size();
		final int[] candEdge = new int[maxN];
		final double[] candScore = new double[maxN];

		// Insert remaining points
		for (Coordinate p : sorted) {
			if (used.contains(p)) {
				continue;
			}

			final int m = ring.size();

			// 1) Score all candidate edges (cheap)
			int c = 0;
			for (int i = 0; i < m; i++) {
				Coordinate q = ring.get(i);
				Coordinate r = ring.get((i + 1) % m);

				double area2 = triangleArea2(p, q, r);
				if (area2 == 0.0) {
					continue; // degenerate triangle
				}

				candEdge[c] = i;
				candScore[c] = area2;
				c++;
			}

			// 2) Sort candidates by objective score (best-first)
			sortCandidatesByScore(candEdge, candScore, c, objective);

			// 3) Test candidates in that order until first valid (expensive part)
			int bestEdge = -1;
			for (int k = 0; k < c; k++) {
				int i = candEdge[k];
				if (isInsertionNonCrossingFast(ring, i, p)) {
					bestEdge = i;
					break;
				}
			}

			if (bestEdge < 0) {
				throw new IllegalStateException("No valid insertion edge found for point " + p);
			}

			ring.add(bestEdge + 1, new Coordinate(p));
			used.add(p);
		}

		// Ensure CCW (JTS doesn't require, but usually preferred)
		if (!isCCW(ring)) {
			Collections.reverse(ring);
		}

		Polygon poly = toPolygon(ring, gf);
		if (!poly.isValid() || poly.getArea() == 0.0) {
			throw new IllegalStateException("Constructed polygon is invalid/degenerate. Try a different pivot.");
		}

		if (postOptimize) {
			return TriangleSwapLocalSearch.improveWithSLS(poly, ring, objective, gf);
		} else {
			return poly;
		}
	}

	/**
	 * Returns true if (s1-s2) intersects (e1-e2) in a way that is NOT allowed. The
	 * only allowed intersection point is 'allowedEndpoint' (exactly).
	 */
	private static boolean segmentsIntersectDisallowed(Coordinate s1, Coordinate s2, Coordinate e1, Coordinate e2, Coordinate allowedEndpoint) {
		INTERSECTOR.computeIntersection(s1, s2, e1, e2);
		if (!INTERSECTOR.hasIntersection()) {
			return false;
		}

		// Proper intersection (interior-interior crossing) is not allowed
		if (INTERSECTOR.isProper()) {
			return true;
		}

		// Non-proper intersection: must be exactly at allowedEndpoint (and only that)
		int n = INTERSECTOR.getIntersectionNum();
		if (n == LineIntersector.COLLINEAR) {
			return true;
		}
		for (int i = 0; i < n; i++) {
			Coordinate ip = INTERSECTOR.getIntersection(i);
			if (!ip.equals2D(allowedEndpoint)) {
				return true;
			}
		}
		return false;
	}

	private static boolean isInsertionNonCrossingFast(List<Coordinate> ring, int edgeIndex, Coordinate p) {
		int m = ring.size();
		Coordinate q = ring.get(edgeIndex);
		Coordinate r = ring.get((edgeIndex + 1) % m);

		// Adjacent edges can only meet at endpoints; skip them to avoid needless tests
		int prevEdge = (edgeIndex - 1 + m) % m; // (... -> q)
		int nextEdge = (edgeIndex + 1) % m; // (r -> ...)

		// Check (q,p) against all edges except replaced edge (edgeIndex) and prevEdge
		// Check (p,r) against all edges except replaced edge (edgeIndex) and nextEdge
		if (segmentHitsAnyEdge(ring, q, p, edgeIndex, prevEdge, q) || segmentHitsAnyEdge(ring, p, r, edgeIndex, nextEdge, r)) {
			return false;
		}

		return true;
	}

	private static boolean segmentHitsAnyEdge(List<Coordinate> ring, Coordinate s1, Coordinate s2, int skipEdge1, int skipEdge2, Coordinate allowedEndpoint) {
		int m = ring.size();

		// Segment AABB
		double sMinX = Math.min(s1.x, s2.x), sMaxX = Math.max(s1.x, s2.x);
		double sMinY = Math.min(s1.y, s2.y), sMaxY = Math.max(s1.y, s2.y);

		for (int j = 0; j < m; j++) {
			if (j == skipEdge1 || j == skipEdge2) {
				continue;
			}

			Coordinate a = ring.get(j);
			Coordinate b = ring.get((j + 1) % m);

			// Edge AABB and quick reject
			double eMinX = Math.min(a.x, b.x), eMaxX = Math.max(a.x, b.x);
			double eMinY = Math.min(a.y, b.y), eMaxY = Math.max(a.y, b.y);

			if (sMaxX < eMinX || eMaxX < sMinX || sMaxY < eMinY || eMaxY < sMinY) {
				continue;
			}

			if (segmentsIntersectDisallowed(s1, s2, a, b, allowedEndpoint)) {
				return true;
			}
		}
		return false;
	}

	private static void sortCandidatesByScore(int[] edge, double[] score, int n, AreaObjective obj) {
		boolean desc = (obj == AreaObjective.MAXIMIZE);
		if (n > 1) {
			quicksort(edge, score, 0, n - 1, desc);
		}
	}

	/**
	 * Paper suggests handling initial collinear prefix; this builds a
	 * non-degenerate starting polygon.
	 */
	private static List<Coordinate> buildInitialRingHandlingCollinear(List<Coordinate> sorted) {
		if (sorted.size() < 3) {
			throw new IllegalArgumentException();
		}

		Coordinate a = sorted.get(0);
		Coordinate b = sorted.get(1);

		int idx = 2;
		while (idx < sorted.size() && index(a, b, sorted.get(idx)) == Orientation.COLLINEAR) {
			idx++;
		}
		if (idx >= sorted.size()) {
			throw new IllegalArgumentException("All points appear collinear; no simple polygonization exists.");
		}

		// Collinear prefix: sorted[0..idx-1], next non-collinear point is r =
		// sorted[idx]
		List<Coordinate> col = new ArrayList<>(sorted.subList(0, idx));
		Coordinate r = sorted.get(idx);

		// Pick extremes along the line to be endpoints p,q
		CollinearOrderComparator order = lineOrder(a, b);
		col.sort(order);

		Coordinate p = col.get(0);
		Coordinate q = col.get(col.size() - 1);

		// Start with triangle p-q-r (ensure CCW later)
		List<Coordinate> ring = new ArrayList<>();
		ring.add(new Coordinate(p));

		// Insert intermediate collinear points between p and q
		for (int i = 1; i < col.size() - 1; i++) {
			ring.add(new Coordinate(col.get(i)));
		}

		ring.add(new Coordinate(q));
		ring.add(new Coordinate(r));
		return ring;
	}

	private static Polygon toPolygon(List<Coordinate> ring, GeometryFactory gf) {
		Coordinate[] coords = new Coordinate[ring.size() + 1];
		for (int i = 0; i < ring.size(); i++) {
			coords[i] = ring.get(i);
		}
		coords[coords.length - 1] = ring.get(0); // close
		LinearRing shell = gf.createLinearRing(coords);
		return gf.createPolygon(shell);
	}

	/** Twice the triangle area (absolute cross product). */
	private static double triangleArea2(Coordinate p, Coordinate q, Coordinate r) {
		return Math.abs((q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x));
	}

	private static boolean isCCW(List<Coordinate> ring) {
		// Signed area (shoelace). Positive => CCW.
		double a = 0.0;
		int n = ring.size();
		for (int i = 0; i < n; i++) {
			Coordinate c1 = ring.get(i);
			Coordinate c2 = ring.get((i + 1) % n);
			a += (c1.x * c2.y) - (c2.x * c1.y);
		}
		return a > 0.0;
	}

	private static boolean pointOnSegment(Coordinate p, Coordinate a, Coordinate b) {
		if (index(a, b, p) != Orientation.COLLINEAR) {
			return false;
		}
		double minX = Math.min(a.x, b.x), maxX = Math.max(a.x, b.x);
		double minY = Math.min(a.y, b.y), maxY = Math.max(a.y, b.y);
		return p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY;
	}

	private static List<Coordinate> dedupePreserveOrder(Collection<Coordinate> pts) {
		Map<Coordinate, Coordinate> map = new LinkedHashMap<>();
		for (Coordinate c : pts) {
			if (c == null) {
				continue;
			}
			map.putIfAbsent(c, new Coordinate(c));
		}
		return new ArrayList<>(map.values());
	}

	/**
	 * Orders points along a line (choose x or y depending on direction), to pick
	 * extremes of a collinear set.
	 */
	private static CollinearOrderComparator lineOrder(Coordinate a, Coordinate b) {
		double dx = Math.abs(b.x - a.x);
		double dy = Math.abs(b.y - a.y);
		if (dx >= dy) {
			return new CollinearOrderComparator(true);
		} else {
			return new CollinearOrderComparator(false);
		}
	}

	private static int index(Coordinate p1, Coordinate p2, Coordinate q) {
		return Orientation.index(p1, p2, q);
//		double dx1 = p2.x - p1.x;
//		double dy1 = p2.y - p1.y;
//		double dx2 = q.x - p1.x;
//		double dy2 = q.y - p1.y;
//
//		double cross = dx1 * dy2 - dy1 * dx2;
//
//		if (cross > 0.0) {
//			return 1; // COUNTERCLOCKWISE / LEFT
//		}
//		if (cross < 0.0) {
//			return -1; // CLOCKWISE / RIGHT
//		}
//		return 0; // COLLINEAR / STRAIGHT
	}

	private static void quicksort(int[] edge, double[] score, int lo, int hi, boolean desc) {
		int i = lo, j = hi;
		int mid = (lo + hi) >>> 1;
		double pivotScore = score[mid];
		int pivotEdge = edge[mid];

		while (i <= j) {
			if (desc) {
				while (score[i] > pivotScore || (score[i] == pivotScore && edge[i] < pivotEdge)) {
					i++;
				}
				while (score[j] < pivotScore || (score[j] == pivotScore && edge[j] > pivotEdge)) {
					j--;
				}
			} else {
				while (score[i] < pivotScore || (score[i] == pivotScore && edge[i] < pivotEdge)) {
					i++;
				}
				while (score[j] > pivotScore || (score[j] == pivotScore && edge[j] > pivotEdge)) {
					j--;
				}
			}

			if (i <= j) {
				double ts = score[i];
				score[i] = score[j];
				score[j] = ts;
				int te = edge[i];
				edge[i] = edge[j];
				edge[j] = te;
				i++;
				j--;
			}
		}

		if (lo < j) {
			quicksort(edge, score, lo, j, desc);
		}
		if (i < hi) {
			quicksort(edge, score, i, hi, desc);
		}
	}

	private static final class CollinearOrderComparator implements Comparator<Coordinate> {
		private final boolean byX;

		CollinearOrderComparator(boolean byX) {
			this.byX = byX;
		}

		@Override
		public int compare(Coordinate o1, Coordinate o2) {
			if (byX) {
				int c = Double.compare(o1.x, o2.x);
				if (c != 0) {
					return c;
				}
				return Double.compare(o1.y, o2.y);
			} else {
				int c = Double.compare(o1.y, o2.y);
				if (c != 0) {
					return c;
				}
				return Double.compare(o1.x, o2.x);
			}
		}
	}

	final class TriangleSwapLocalSearch {

		static Polygon improveWithSLS(Polygon start, List<Coordinate> allPoints, AreaOptimalPolygonizer.AreaObjective objective, GeometryFactory gf) {
			List<Coordinate> ring = exteriorRingToList(start);
			// Ensure CCW for reflex/convex tests
			if (!isCCW(ring)) {
				Collections.reverse(ring);
			}

			List<Coordinate> pts = dedupePreserveOrder(allPoints);

			boolean improved;
			do {
				improved = false;
				int n = ring.size();
				for (int j = 0; j < n; j++) {
					int jm2 = mod(j - 2, n);
					int jm1 = mod(j - 1, n);
					int jp1 = mod(j + 1, n);
					int jp2 = mod(j + 2, n);

					if (!isReflex(ring, j) || !isConvex(ring, jm1) || !isConvex(ring, jp1)) {
						continue;
					}

					Coordinate a = ring.get(jm1);
					Coordinate b = ring.get(j);
					Coordinate c = ring.get(jp1);

					if (triangleArea2(a, b, c) == 0.0) {
						continue;
					}
					if (!triangleEmptyWrtPoints(pts, a, b, c)) {
						continue;
					}

					double ins = triangleArea2(a, b, c);

					// left remove candidate: (jm2, jm1, j), swap jm1 <-> j
					Coordinate l0 = ring.get(jm2), l1 = ring.get(jm1), l2 = ring.get(j);
					boolean leftEmpty = triangleArea2(l0, l1, l2) != 0.0 && triangleEmptyWrtPoints(pts, l0, l1, l2);
					double leftArea = triangleArea2(l0, l1, l2);

					// right remove candidate: (j, jp1, jp2), swap j <-> jp1
					Coordinate r0 = ring.get(j), r1 = ring.get(jp1), r2 = ring.get(jp2);
					boolean rightEmpty = triangleArea2(r0, r1, r2) != 0.0 && triangleEmptyWrtPoints(pts, r0, r1, r2);
					double rightArea = triangleArea2(r0, r1, r2);

					int swapI = -1, swapK = -1; // swap indices
					if (objective == AreaOptimalPolygonizer.AreaObjective.MAXIMIZE) {
						boolean leftOk = leftEmpty && leftArea < ins;
						boolean rightOk = rightEmpty && rightArea < ins;

						if (leftOk && rightOk) {
							// remove smaller triangle
							if (leftArea <= rightArea) {
								swapI = jm1;
								swapK = j;
							} else {
								swapI = j;
								swapK = jp1;
							}
						} else if (leftOk) {
							swapI = jm1;
							swapK = j;
						} else if (rightOk) {
							swapI = j;
							swapK = jp1;
						}
					} else { // MIN_AREA
						boolean leftOk = leftEmpty && leftArea > ins;
						boolean rightOk = rightEmpty && rightArea > ins;

						if (leftOk && rightOk) {
							// remove larger triangle
							if (leftArea >= rightArea) {
								swapI = jm1;
								swapK = j;
							} else {
								swapI = j;
								swapK = jp1;
							}
						} else if (leftOk) {
							swapI = jm1;
							swapK = j;
						} else if (rightOk) {
							swapI = j;
							swapK = jp1;
						}
					}

					if (swapI >= 0) {
						if (wouldBeSimpleAfterAdjacentSwap(ring, swapI, swapK)) {
							Collections.swap(ring, swapI, swapK);
							improved = true;
							// ring changed; restart scan (optional but tends to behave better)
							break;
						}
					}
				}
			} while (improved);

			if (!isCCW(ring)) {
				Collections.reverse(ring);
			}
			return toPolygon(ring, gf);
		}

		private static boolean wouldBeSimpleAfterAdjacentSwap(List<Coordinate> ring, int i, int k) {
			int n = ring.size();
			// ensure i and k are adjacent cyclically
			if ((n < 4) || !(mod(i + 1, n) == k || mod(k + 1, n) == i)) {
				return false;
			}

			// normalize so k = i+1 mod n
			if (mod(i + 1, n) != k) {
				int tmp = i;
				i = k;
				k = tmp;
			}
			int prev = mod(i - 1, n);
			int next = mod(k + 1, n);

			Coordinate vPrev = ring.get(prev);
			Coordinate vi = ring.get(i);
			Coordinate vk = ring.get(k);
			Coordinate vNext = ring.get(next);

			// New edges after swap:
			// (vPrev - vk), (vk - vi), (vi - vNext)
			// Edge (vk-vi) is same as old (vi-vk) reversed, but check anyway for collinear
			// overlaps.
			return edgeDoesNotCrossBoundary(ring, vPrev, vk, Set.of(prev, i)) && // replaces (vPrev-vi)
					edgeDoesNotCrossBoundary(ring, vk, vi, Set.of(i, k)) && // replaces (vi-vk)
					edgeDoesNotCrossBoundary(ring, vi, vNext, Set.of(k, next)); // replaces (vk-vNext)
		}

		/**
		 * Checks segment (s1-s2) does not intersect polygon boundary edges, except at
		 * shared endpoints. skipEdgeIndices is a set of edge start indices to ignore.
		 * (Edge index e means ring[e] -> ring[e+1].)
		 */
		private static boolean edgeDoesNotCrossBoundary(List<Coordinate> ring, Coordinate s1, Coordinate s2, Set<Integer> skipEdgeIndices) {
			int n = ring.size();

			for (int e = 0; e < n; e++) {
				if (skipEdgeIndices.contains(e)) {
					continue;
				}
				Coordinate a = ring.get(e);
				Coordinate b = ring.get((e + 1) % n);

				INTERSECTOR.computeIntersection(s1, s2, a, b);
				if (!INTERSECTOR.hasIntersection()) {
					continue;
				}

				if (INTERSECTOR.isProper()) {
					return false;
				}

				// only allow intersections at shared endpoints
				int m = INTERSECTOR.getIntersectionNum();
				for (int t = 0; t < m; t++) {
					Coordinate ip = INTERSECTOR.getIntersection(t);
					boolean shared = ip.equals2D(s1) || ip.equals2D(s2) ? (ip.equals2D(a) || ip.equals2D(b)) : false;
					if (!shared) {
						return false;
					}
				}
			}
			return true;
		}

		private static boolean triangleEmptyWrtPoints(List<Coordinate> pts, Coordinate a, Coordinate b, Coordinate c) {
			// Empty w.r.t. S: no point (other than a,b,c) inside OR on boundary.
			for (Coordinate p : pts) {
				if (p.equals2D(a) || p.equals2D(b) || p.equals2D(c)) {
					continue;
				}
				if (pointInTriangleOrOnEdge(p, a, b, c)) {
					return false;
				}
			}
			return true;
		}

		private static boolean pointInTriangleOrOnEdge(Coordinate p, Coordinate a, Coordinate b, Coordinate c) {
			int o1 = index(a, b, p);
			int o2 = index(b, c, p);
			int o3 = index(c, a, p);

			boolean hasPos = (o1 > 0) || (o2 > 0) || (o3 > 0);
			boolean hasNeg = (o1 < 0) || (o2 < 0) || (o3 < 0);
			// if not both signs, point is inside or on boundary
			if (!(hasPos && hasNeg)) {
				// also ensure within bounding box (handles collinear outside segment range)
				double minX = Math.min(a.x, Math.min(b.x, c.x));
				double maxX = Math.max(a.x, Math.max(b.x, c.x));
				double minY = Math.min(a.y, Math.min(b.y, c.y));
				double maxY = Math.max(a.y, Math.max(b.y, c.y));
				return p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY;
			}
			return false;
		}

		private static boolean isReflex(List<Coordinate> ring, int i) {
			int n = ring.size();
			Coordinate prev = ring.get(mod(i - 1, n));
			Coordinate cur = ring.get(i);
			Coordinate next = ring.get(mod(i + 1, n));
			return index(prev, cur, next) < 0; // CW turn in CCW polygon
		}

		private static boolean isConvex(List<Coordinate> ring, int i) {
			int n = ring.size();
			Coordinate prev = ring.get(mod(i - 1, n));
			Coordinate cur = ring.get(i);
			Coordinate next = ring.get(mod(i + 1, n));
			return index(prev, cur, next) > 0;
		}

		private static int mod(int x, int n) {
			int r = x % n;
			return (r < 0) ? r + n : r;
		}

		private static List<Coordinate> exteriorRingToList(Polygon p) {
			Coordinate[] coords = p.getExteriorRing().getCoordinates();
			// coords is closed; drop last
			List<Coordinate> ring = new ArrayList<>(coords.length - 1);
			for (int i = 0; i < coords.length - 1; i++) {
				ring.add(new Coordinate(coords[i]));
			}
			return ring;
		}
	}
}