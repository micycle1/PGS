package micycle.pgs.commons;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Objects;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;

/**
 * Keil & Snoeyink optimal convex partitioner.
 *
 * <p>
 * Computes an optimal convex partition of a polygon (possibly with holes) using
 * the dynamic-programming algorithm originally exposed as ConvexPartition_OPT
 * in the polypartition library. The algorithm minimizes the number of added
 * diagonals (hence the number of convex pieces) and reconstructs an optimal
 * partition of the input.
 * <p>
 * Complexity: worst-case O(n^3) time.
 * </p>
 * 
 * @author Michael Carleton
 */
public final class KeilSnoeyinkConvexPartitioner {

	private KeilSnoeyinkConvexPartitioner() {
	}

	/**
	 * Computes an optimal convex partition of the given polygon.
	 *
	 * <p>
	 * This is the public entry point for the Keil &amp; Snoeyink
	 * dynamic-programming partition (the implementation corresponds to the
	 * ConvexPartition_OPT strategy). Holes in the input polygon are bridged to
	 * produce simple polygons prior to the DP on each simple polygon.
	 * </p>
	 *
	 * @param input non-null Polygon (may contain holes)
	 * @return a non-null List of Polygons. Each returned polygon is simple (no
	 *         holes), oriented CCW, and the union of the returned polygons equals
	 *         the input polygon.
	 */
	public static List<Polygon> convexPartition(Polygon input) {
		Objects.requireNonNull(input, "input");
		if (input.isEmpty()) {
			return List.of();
		}

		GeometryFactory gf = input.getFactory();
		PrecisionModel pm = gf.getPrecisionModel();

		// Remove holes by bridging (same as original C++ library).
		List<List<Coordinate>> simplePolys = removeHoles(input, pm);

		List<Polygon> out = new ArrayList<>();
		for (List<Coordinate> polyPts : simplePolys) {
			List<Coordinate> pts = cleanupOpenRing(polyPts);
			pts = removeConsecutiveDuplicates(pts);
			if (pts.size() < 3) {
				continue;
			}

			pts = ensureCCW(pts);

			if (isConvexWeakly(pts)) {
				out.add(toJTSPolygon(pts, gf, pm));
				continue;
			}

			List<List<Coordinate>> parts = convexPartitionOptSimple(pts);
			for (List<Coordinate> part : parts) {
				List<Coordinate> p = cleanupOpenRing(part);
				p = removeConsecutiveDuplicates(p);
				if (p.size() < 3) {
					continue;
				}
				p = ensureCCW(p);
				out.add(toJTSPolygon(p, gf, pm));
			}
		}
		return out;
	}

	private static List<List<Coordinate>> convexPartitionOptSimple(List<Coordinate> ptsCCW) {
		final int n = ptsCCW.size();
		if (n < 3) {
			throw new IllegalArgumentException("Polygon has < 3 vertices");
		}
		if (n == 3) {
			return List.of(List.copyOf(ptsCCW));
		}

		// vertices[]
		final V[] vertices = new V[n];
		for (int i = 0; i < n; i++) {
			vertices[i] = new V(ptsCCW.get(i));
		}
		for (int i = 0; i < n; i++) {
			vertices[i].prev = (i == 0) ? (n - 1) : (i - 1);
			vertices[i].next = (i == n - 1) ? 0 : (i + 1);
		}
		for (int i = 1; i < n; i++) {
			vertices[i].isConvex = !isReflex(vertices[vertices[i].prev].p, vertices[i].p, vertices[vertices[i].next].p);
		}
		vertices[0].isConvex = false; // by convention (as in C++)

		final boolean[][] visible = new boolean[n][n];
		final int[][] weight = new int[n][n];
		final DiagonalDeque[][] pairs = new DiagonalDeque[n][n];

		// Initialize visibility and base weights.
		for (int i = 0; i < n - 1; i++) {
			final Coordinate p1 = vertices[i].p;
			for (int j = i + 1; j < n; j++) {
				visible[i][j] = true;
				weight[i][j] = (j == i + 1) ? 0 : Integer.MAX_VALUE;

				if (j != i + 1) {
					final Coordinate p2 = vertices[j].p;

					if (!inCone(vertices, i, p2)) {
						visible[i][j] = false;
						continue;
					}
					if (!inCone(vertices, j, p1)) {
						visible[i][j] = false;
						continue;
					}

					for (int k = 0; k < n; k++) {
						final Coordinate p3 = vertices[k].p;
						final Coordinate p4 = vertices[(k == n - 1) ? 0 : (k + 1)].p;
						if (intersects(p1, p2, p3, p4)) {
							visible[i][j] = false;
							break;
						}
					}
				}
			}
		}

		// Triangles (gap=2) that are visible have weight 0 and one "pair".
		for (int i = 0; i < n - 2; i++) {
			int j = i + 2;
			if (visible[i][j]) {
				weight[i][j] = 0;
				pairs[i][j] = new DiagonalDeque();
				pairs[i][j].addLast(i + 1, i + 1);
			}
		}
		visible[0][n - 1] = true;

		// DP
		for (int gap = 3; gap < n; gap++) {
			for (int i = 0; i < n - gap; i++) {
				if (vertices[i].isConvex) {
					continue;
				}
				int k = i + gap;
				if (!visible[i][k]) {
					continue;
				}

				if (!vertices[k].isConvex) {
					for (int j = i + 1; j < k; j++) {
						typeA(i, j, k, vertices, visible, weight, pairs);
					}
				} else {
					for (int j = i + 1; j < k - 1; j++) {
						if (vertices[j].isConvex) {
							continue;
						}
						typeA(i, j, k, vertices, visible, weight, pairs);
					}
					typeA(i, k - 1, k, vertices, visible, weight, pairs);
				}
			}

			for (int k = gap; k < n; k++) {
				if (vertices[k].isConvex) {
					continue;
				}
				int i = k - gap;
				if (!vertices[i].isConvex) {
					continue;
				}
				if (!visible[i][k]) {
					continue;
				}

				typeB(i, i + 1, k, vertices, visible, weight, pairs);
				for (int j = i + 2; j < k; j++) {
					if (vertices[j].isConvex) {
						continue;
					}
					typeB(i, j, k, vertices, visible, weight, pairs);
				}
			}
		}

		// Recover solution (first pass prunes pair-lists to chosen solution).
		boolean ok = true;
		ArrayDeque<Diagonal> diagonals = new ArrayDeque<>();
		diagonals.addFirst(new Diagonal(0, n - 1));

		while (!diagonals.isEmpty()) {
			Diagonal d = diagonals.removeFirst();
			if (d.b - d.a <= 1) {
				continue;
			}

			DiagonalDeque pd = pairs[d.a][d.b];
			if (pd == null || pd.isEmpty()) {
				ok = false;
				break;
			}

			if (!vertices[d.a].isConvex) {
				Diagonal chosen = pd.peekLast();
				int j = chosen.b;

				diagonals.addFirst(new Diagonal(j, d.b));
				if (j - d.a > 1) {
					if (chosen.a != chosen.b) {
						DiagonalDeque pd2 = pairs[d.a][j];
						while (true) {
							if (pd2 == null || pd2.isEmpty()) {
								ok = false;
								break;
							}
							Diagonal t = pd2.peekLast();
							if (chosen.a != t.a) {
								pd2.removeLast();
							} else {
								break;
							}
						}
						if (!ok) {
							break;
						}
					}
					diagonals.addFirst(new Diagonal(d.a, j));
				}
			} else {
				Diagonal chosen = pd.peekFirst();
				int j = chosen.a;

				diagonals.addFirst(new Diagonal(d.a, j));
				if (d.b - j > 1) {
					if (chosen.a != chosen.b) {
						DiagonalDeque pd2 = pairs[j][d.b];
						while (true) {
							if (pd2 == null || pd2.isEmpty()) {
								ok = false;
								break;
							}
							Diagonal t = pd2.peekFirst();
							if (chosen.b != t.b) {
								pd2.removeFirst();
							} else {
								break;
							}
						}
						if (!ok) {
							break;
						}
					}
					diagonals.addFirst(new Diagonal(j, d.b));
				}
			}
		}

		if (!ok) {
			throw new IllegalStateException("ConvexPartition_OPT failed to recover a solution.");
		}

		// Recover actual polygons (second pass).
		List<List<Coordinate>> parts = new ArrayList<>();
		diagonals.clear();
		diagonals.addFirst(new Diagonal(0, n - 1));

		while (!diagonals.isEmpty()) {
			Diagonal root = diagonals.removeFirst();
			if (root.b - root.a <= 1) {
				continue;
			}

			ArrayDeque<Diagonal> diagonals2 = new ArrayDeque<>();
			diagonals2.addFirst(root);

			IntList indices = new IntList();
			indices.add(root.a);
			indices.add(root.b);

			while (!diagonals2.isEmpty()) {
				Diagonal d = diagonals2.removeFirst();
				if (d.b - d.a <= 1) {
					continue;
				}

				DiagonalDeque pd = pairs[d.a][d.b];

				boolean ijReal = true, jkReal = true;
				int j;

				if (!vertices[d.a].isConvex) {
					Diagonal chosen = pd.peekLast();
					j = chosen.b;
					if (chosen.a != chosen.b) {
						ijReal = false;
					}
				} else {
					Diagonal chosen = pd.peekFirst();
					j = chosen.a;
					if (chosen.a != chosen.b) {
						jkReal = false;
					}
				}

				Diagonal ij = new Diagonal(d.a, j);
				if (ijReal) {
					diagonals.addLast(ij);
				} else {
					diagonals2.addLast(ij);
				}

				Diagonal jk = new Diagonal(j, d.b);
				if (jkReal) {
					diagonals.addLast(jk);
				} else {
					diagonals2.addLast(jk);
				}

				indices.add(j);
			}

			int[] idx = indices.toSortedArray();
			List<Coordinate> poly = new ArrayList<>(idx.length);
			for (int v : idx) {
				poly.add(vertices[v].p);
			}
			parts.add(poly);
		}

		return parts;
	}

	private static void updateState(int a, int b, int w, int i, int j, boolean[][] visible, int[][] weight, DiagonalDeque[][] pairs) {
		if (!visible[a][b]) {
			return;
		}

		int w2 = weight[a][b];
		if (w > w2) {
			return;
		}

		DiagonalDeque q = pairs[a][b];
		if (q == null) {
			q = new DiagonalDeque();
			pairs[a][b] = q;
		}

		if (w < w2) {
			q.clear();
			q.addFirst(i, j);
			weight[a][b] = w;
		} else {
			if (!q.isEmpty() && i <= q.peekFirst().a) {
				return;
			}
			while (!q.isEmpty() && q.peekFirst().b >= j) {
				q.removeFirst();
			}
			q.addFirst(i, j);
		}
	}

	private static void typeA(int i, int j, int k, V[] vertices, boolean[][] visible, int[][] weight, DiagonalDeque[][] pairs) {
		if (!visible[i][j]) {
			return;
		}

		int top = j;
		int w = weight[i][j];

		if (k - j > 1) {
			if (!visible[j][k]) {
				return;
			}
			w += weight[j][k] + 1;
		}

		if (j - i > 1) {
			DiagonalDeque q = pairs[i][j];

			int lastPos = -1;
			if (q != null && !q.isEmpty()) {
				for (int pos = q.size() - 1; pos >= 0; pos--) {
					int idx2 = q.get(pos).b;
					if (!isReflex(vertices[idx2].p, vertices[j].p, vertices[k].p)) {
						lastPos = pos;
					} else {
						break;
					}
				}
			}

			if (lastPos < 0) {
				w++;
			} else {
				int idx1 = q.get(lastPos).a;
				if (isReflex(vertices[k].p, vertices[i].p, vertices[idx1].p)) {
					w++;
				} else {
					top = idx1;
				}
			}
		}

		updateState(i, k, w, top, j, visible, weight, pairs);
	}

	private static void typeB(int i, int j, int k, V[] vertices, boolean[][] visible, int[][] weight, DiagonalDeque[][] pairs) {
		if (!visible[j][k]) {
			return;
		}

		int top = j;
		int w = weight[j][k];

		if (j - i > 1) {
			if (!visible[i][j]) {
				return;
			}
			w += weight[i][j] + 1;
		}

		if (k - j > 1) {
			DiagonalDeque q = pairs[j][k];

			if (q != null && !q.isEmpty() && !isReflex(vertices[i].p, vertices[j].p, vertices[q.peekFirst().a].p)) {

				int lastPos = 0;
				for (int pos = 0; pos < q.size(); pos++) {
					int idx1 = q.get(pos).a;
					if (!isReflex(vertices[i].p, vertices[j].p, vertices[idx1].p)) {
						lastPos = pos;
					} else {
						break;
					}
				}
				int idx2 = q.get(lastPos).b;
				if (isReflex(vertices[idx2].p, vertices[k].p, vertices[i].p)) {
					w++;
				} else {
					top = idx2;
				}
			} else {
				w++;
			}
		}

		updateState(i, k, w, j, top, visible, weight, pairs);
	}

	// -------------------- Geometry primitives (ported) --------------------

	private static final class V {
		final Coordinate p; // treat as immutable (we only create our own)
		int prev, next;
		boolean isConvex;

		V(Coordinate p) {
			this.p = p;
		}
	}

	private record Diagonal(int a, int b) {
	}

	private static boolean isConvex(Coordinate p1, Coordinate p2, Coordinate p3) {
		return Orientation.index(p1, p2, p3) == Orientation.COUNTERCLOCKWISE;
	}

	private static boolean isReflex(Coordinate p1, Coordinate p2, Coordinate p3) {
		return Orientation.index(p1, p2, p3) == Orientation.CLOCKWISE;
	}

	private static boolean inCone(V[] vertices, int vi, Coordinate p) {
		V v = vertices[vi];
		Coordinate p1 = vertices[v.prev].p;
		Coordinate p2 = v.p;
		Coordinate p3 = vertices[v.next].p;
		return inCone(p1, p2, p3, p);
	}

	private static boolean inCone(Coordinate p1, Coordinate p2, Coordinate p3, Coordinate p) {
		boolean convex = isConvex(p1, p2, p3);
		if (convex) {
			if (!isConvex(p1, p2, p)) {
				return false;
			}
			if (!isConvex(p2, p3, p)) {
				return false;
			}
			return true;
		} else {
			if (isConvex(p1, p2, p)) {
				return true;
			}
			if (isConvex(p2, p3, p)) {
				return true;
			}
			return false;
		}
	}

	// C++ Intersects() exact port (true if segments intersect, excluding shared
	// endpoints).
	private static boolean intersects(Coordinate p11, Coordinate p12, Coordinate p21, Coordinate p22) {
		if (p11.equals2D(p21) || p11.equals2D(p22) || p12.equals2D(p21) || p12.equals2D(p22)) {
			return false;
		}

		double v1ortx = p12.y - p11.y;
		double v1orty = p11.x - p12.x;

		double v2ortx = p22.y - p21.y;
		double v2orty = p21.x - p22.x;

		double vx, vy;

		vx = p21.x - p11.x;
		vy = p21.y - p11.y;
		double dot21 = vx * v1ortx + vy * v1orty;

		vx = p22.x - p11.x;
		vy = p22.y - p11.y;
		double dot22 = vx * v1ortx + vy * v1orty;

		vx = p11.x - p21.x;
		vy = p11.y - p21.y;
		double dot11 = vx * v2ortx + vy * v2orty;

		vx = p12.x - p21.x;
		vy = p12.y - p21.y;
		double dot12 = vx * v2ortx + vy * v2orty;

		if (dot11 * dot12 > 0) {
			return false;
		}
		if (dot21 * dot22 > 0) {
			return false;
		}
		return true;
	}

	private static final class PolyHole {
		final boolean hole;
		final List<Coordinate> pts; // open ring

		PolyHole(boolean hole, List<Coordinate> pts) {
			this.hole = hole;
			this.pts = pts;
		}
	}

	private static List<List<Coordinate>> removeHoles(Polygon input, PrecisionModel pm) {
		List<PolyHole> in = new ArrayList<>();

		List<Coordinate> shell = ringToCoords(input.getExteriorRing().getCoordinates(), pm);
		shell = cleanupOpenRing(shell);
		shell = ensureCCW(shell);
		in.add(new PolyHole(false, shell));

		for (int i = 0; i < input.getNumInteriorRing(); i++) {
			List<Coordinate> hole = ringToCoords(input.getInteriorRingN(i).getCoordinates(), pm);
			hole = cleanupOpenRing(hole);
			hole = ensureCW(hole);
			in.add(new PolyHole(true, hole));
		}

		boolean hasHoles = false;
		for (PolyHole ph : in) {
			if (ph.hole) {
				hasHoles = true;
				break;
			}
		}
		if (!hasHoles) {
			return List.of(List.copyOf(shell));
		}

		List<PolyHole> polys = new ArrayList<>(in);

		while (true) {
			// Find the hole point with the largest x.
			int holeIdx = -1;
			int holePointIndex = 0;

			for (int pi = 0; pi < polys.size(); pi++) {
				PolyHole ph = polys.get(pi);
				if (!ph.hole) {
					continue;
				}
				if (holeIdx < 0) {
					holeIdx = pi;
					holePointIndex = 0;
				}
				for (int i = 0; i < ph.pts.size(); i++) {
					if (ph.pts.get(i).x > polys.get(holeIdx).pts.get(holePointIndex).x) {
						holeIdx = pi;
						holePointIndex = i;
					}
				}
			}
			if (holeIdx < 0) {
				break;
			}

			Coordinate holePoint = polys.get(holeIdx).pts.get(holePointIndex);

			boolean pointFound = false;
			int polyIdx = -1;
			int polyPointIndex = -1;
			Coordinate bestPolyPoint = null;

			for (int pi = 0; pi < polys.size(); pi++) {
				PolyHole ph = polys.get(pi);
				if (ph.hole) {
					continue;
				}

				int m = ph.pts.size();
				for (int i = 0; i < m; i++) {
					Coordinate candidate = ph.pts.get(i);
					if (candidate.x <= holePoint.x) {
						continue;
					}

					Coordinate prev = ph.pts.get((i + m - 1) % m);
					Coordinate next = ph.pts.get((i + 1) % m);
					if (!inCone(prev, candidate, next, holePoint)) {
						continue;
					}

					if (pointFound) {
						double d1 = dist2(holePoint, candidate);
						double d2 = dist2(holePoint, bestPolyPoint);
						if (d2 < d1) {
							continue;
						}
					}

					boolean visible = true;
					for (int pj = 0; pj < polys.size() && visible; pj++) {
						PolyHole ph2 = polys.get(pj);
						if (ph2.hole) {
							continue;
						}

						int mm = ph2.pts.size();
						for (int e = 0; e < mm; e++) {
							Coordinate a = ph2.pts.get(e);
							Coordinate b = ph2.pts.get((e + 1) % mm);
							if (intersects(holePoint, candidate, a, b)) {
								visible = false;
								break;
							}
						}
					}

					if (visible) {
						pointFound = true;
						bestPolyPoint = candidate;
						polyIdx = pi;
						polyPointIndex = i;
					}
				}
			}

			if (!pointFound) {
				throw new IllegalStateException("RemoveHoles failed: no visible bridge found.");
			}

			PolyHole hole = polys.get(holeIdx);
			PolyHole poly = polys.get(polyIdx);

			List<Coordinate> newPts = new ArrayList<>(hole.pts.size() + poly.pts.size() + 2);
			for (int i = 0; i <= polyPointIndex; i++) {
				newPts.add(poly.pts.get(i));
			}
			for (int i = 0; i <= hole.pts.size(); i++) {
				newPts.add(hole.pts.get((i + holePointIndex) % hole.pts.size()));
			}
			for (int i = polyPointIndex; i < poly.pts.size(); i++) {
				newPts.add(poly.pts.get(i));
			}

			int a = Math.max(holeIdx, polyIdx);
			int b = Math.min(holeIdx, polyIdx);
			polys.remove(a);
			polys.remove(b);
			polys.add(new PolyHole(false, newPts));
		}

		List<List<Coordinate>> out = new ArrayList<>();
		for (PolyHole ph : polys) {
			out.add(List.copyOf(ph.pts));
		}
		return out;
	}

	private static double dist2(Coordinate a, Coordinate b) {
		double dx = b.x - a.x, dy = b.y - a.y;
		return dx * dx + dy * dy;
	}

	/**
	 * Small deque for Diagonal pairs per DP state
	 */
	private static final class DiagonalDeque {
		private int[] a = new int[4];
		private int[] b = new int[4];
		private int head = 0;
		private int size = 0;

		boolean isEmpty() {
			return size == 0;
		}

		int size() {
			return size;
		}

		void clear() {
			head = 0;
			size = 0;
		}

		Diagonal peekFirst() {
			int i = head;
			return new Diagonal(a[i], b[i]);
		}

		Diagonal peekLast() {
			int i = idx(size - 1);
			return new Diagonal(a[i], b[i]);
		}

		Diagonal get(int pos) {
			int i = idx(pos);
			return new Diagonal(a[i], b[i]);
		}

		void addFirst(int ia, int ib) {
			ensureCap(size + 1);
			head = (head - 1 + a.length) % a.length;
			a[head] = ia;
			b[head] = ib;
			size++;
		}

		void addLast(int ia, int ib) {
			ensureCap(size + 1);
			int i = idx(size);
			a[i] = ia;
			b[i] = ib;
			size++;
		}

		void removeFirst() {
			if (size == 0) {
				throw new NoSuchElementException();
			}
			head = (head + 1) % a.length;
			size--;
		}

		void removeLast() {
			if (size == 0) {
				throw new NoSuchElementException();
			}
			size--;
		}

		private int idx(int pos) {
			return (head + pos) % a.length;
		}

		private void ensureCap(int cap) {
			if (cap <= a.length) {
				return;
			}
			int newCap = Math.max(cap, a.length * 2);
			int[] na = new int[newCap];
			int[] nb = new int[newCap];
			for (int i = 0; i < size; i++) {
				int j = idx(i);
				na[i] = a[j];
				nb[i] = b[j];
			}
			a = na;
			b = nb;
			head = 0;
		}
	}

	// -------------------- Utilities / JTS interop --------------------

	private static List<Coordinate> ringToCoords(Coordinate[] coords, PrecisionModel pm) {
		int len = coords.length;
		if (len >= 2 && coords[0].equals2D(coords[len - 1])) {
			len--;
		}

		List<Coordinate> out = new ArrayList<>(len);
		for (int i = 0; i < len; i++) {
			Coordinate c = new Coordinate(coords[i].x, coords[i].y);
			pm.makePrecise(c);
			if (c.x == 0.0) {
				c.x = 0.0; // normalize -0.0
			}
			if (c.y == 0.0) {
				c.y = 0.0;
			}
			out.add(c);
		}
		return out;
	}

	private static Polygon toJTSPolygon(List<Coordinate> pts, GeometryFactory gf, PrecisionModel pm) {
		List<Coordinate> p = removeConsecutiveDuplicates(cleanupOpenRing(pts));
		if (p.size() < 3) {
			throw new IllegalArgumentException("Part has < 3 vertices after cleanup");
		}

		p = ensureCCW(p);

		Coordinate[] cs = new Coordinate[p.size() + 1];
		for (int i = 0; i < p.size(); i++) {
			Coordinate c = new Coordinate(p.get(i).x, p.get(i).y);
			pm.makePrecise(c);
			cs[i] = c;
		}
		cs[p.size()] = new Coordinate(cs[0]);
		return gf.createPolygon(cs);
	}

	private static List<Coordinate> cleanupOpenRing(List<Coordinate> pts) {
		if (pts.isEmpty()) {
			return List.of();
		}
		int n = pts.size();
		if (n >= 2 && pts.get(0).equals2D(pts.get(n - 1))) {
			return List.copyOf(pts.subList(0, n - 1));
		}
		return List.copyOf(pts);
	}

	private static List<Coordinate> removeConsecutiveDuplicates(List<Coordinate> pts) {
		if (pts.size() < 2) {
			return pts;
		}

		List<Coordinate> out = new ArrayList<>(pts.size());
		Coordinate prev = null;
		for (Coordinate c : pts) {
			if (prev == null || !c.equals2D(prev)) {
				out.add(c);
			}
			prev = c;
		}
		if (out.size() >= 2 && out.get(0).equals2D(out.get(out.size() - 1))) {
			out.remove(out.size() - 1);
		}
		return List.copyOf(out);
	}

	private static List<Coordinate> ensureCCW(List<Coordinate> pts) {
		if (signedArea2(pts) < 0.0) {
			List<Coordinate> rev = new ArrayList<>(pts);
			Collections.reverse(rev);
			return List.copyOf(rev);
		}
		return pts;
	}

	private static List<Coordinate> ensureCW(List<Coordinate> pts) {
		if (signedArea2(pts) > 0.0) {
			List<Coordinate> rev = new ArrayList<>(pts);
			Collections.reverse(rev);
			return List.copyOf(rev);
		}
		return pts;
	}

	private static double signedArea2(List<Coordinate> pts) {
		double a2 = 0.0;
		for (int i = 0, n = pts.size(); i < n; i++) {
			Coordinate p = pts.get(i);
			Coordinate q = pts.get((i + 1) % n);
			a2 += p.x * q.y - p.y * q.x;
		}
		return a2;
	}

	private static boolean isConvexWeakly(List<Coordinate> pts) {
		int n = pts.size();
		if (n < 3) {
			return false;
		}

		double a2 = signedArea2(pts);
		if (a2 == 0.0) {
			return false;
		}
		boolean ccw = a2 > 0.0;

		for (int i = 0; i < n; i++) {
			Coordinate pPrev = pts.get((i + n - 1) % n);
			Coordinate p = pts.get(i);
			Coordinate pNext = pts.get((i + 1) % n);

			int o = Orientation.index(pPrev, p, pNext);
			if (ccw) {
				if (o == Orientation.CLOCKWISE) {
					return false;
				}
			} else {
				if (o == Orientation.COUNTERCLOCKWISE) {
					return false;
				}
			}
		}
		return true;
	}

	private static final class IntList {
		private int[] a = new int[8];
		private int size = 0;

		void add(int v) {
			if (size == a.length) {
				a = Arrays.copyOf(a, a.length * 2);
			}
			a[size++] = v;
		}

		int[] toSortedArray() {
			int[] r = Arrays.copyOf(a, size);
			Arrays.sort(r);
			return r;
		}
	}
}