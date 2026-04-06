package micycle.pgs.commons;

import java.util.List;

import processing.core.PVector;

/**
 * Hyper-optimized 2-opt uncrossing (untangling) for self-crossing polygon
 * tours.
 *
 * <h2>What It Does</h2>
 * <p>
 * Given a closed polygon tour (cyclic sequence of vertices), this routine
 * removes self-intersections to produce a simple (non-self-intersecting)
 * polygon by repeatedly applying a 2-opt move: when two non-adjacent edges
 * cross, it reverses the vertex subsequence between them to eliminate that
 * crossing.
 * </p>
 *
 * <h2>How It Works (Key Idea)</h2>
 * <p>
 * A 2-opt reversal removes exactly two geometric edges and creates exactly two
 * new geometric edges. All other segments in the polygon remain the same (just
 * traversed in reverse order), so their crossing status does not change.
 * Therefore, after a swap only the <em>two new boundary edges</em> can
 * introduce new crossings.
 * </p>
 *
 * <h2>Algorithm</h2>
 * <ol>
 * <li><strong>Dirty-edge stack initialization:</strong> Treat each edge (by its
 * start index {@code i}, representing {@code (i -> i+1)} with wraparound) as
 * "dirty" and push all indices onto a stack.</li>
 * <li><strong>Edge validation:</strong> Pop one dirty edge {@code i} and scan
 * it against all non-adjacent edges to find any crossing partner {@code j}. If
 * none is found, edge {@code i} is "clean" and is not revisited unless it
 * becomes a new boundary edge of a later 2-opt move.</li>
 * <li><strong>2-opt reversal:</strong> If a crossing is found between edges
 * {@code (i,i+1)} and {@code (j,j+1)}, reverse the appropriate contiguous
 * vertex range to remove the crossing.</li>
 * <li><strong>Incremental rechecking:</strong> After the reversal, push only
 * the two new boundary-edge indices back onto the dirty stack (because only
 * those two edges have changed geometrically).</li>
 * </ol>
 *
 * <h2>Performance Optimizations</h2>
 * <ul>
 * <li><strong>No global rescans / no crossing graph:</strong> Avoids
 * {@code O(n^2)} rescans per swap and avoids hash-map/set maintenance
 * overhead.</li>
 * <li><strong>Array-based coordinates:</strong> Uses primitive {@code float[]}
 * for cache-friendly access in tight loops.</li>
 * <li><strong>Low-branch intersection test:</strong> Inlined proper
 * intersection using XOR sign checks, matching the common "strict" crossing
 * definition.</li>
 * <li><strong>Modulo avoidance:</strong> Wraparound edges are handled
 * explicitly to minimize {@code % n} in hot loops.</li>
 * <li><strong>Optional AABB reject:</strong> A fast bounding-box overlap test
 * can be enabled to prune most non-crossing candidates before orientation
 * math.</li>
 * </ul>
 *
 * <h2>Complexity</h2>
 * <ul>
 * <li><strong>Time:</strong> Typically {@code O((n + s) * n)} where {@code n}
 * is the number of vertices and {@code s} is the number of 2-opt swaps
 * performed. Each swap triggers rechecks of only two edges, each checked
 * against {@code O(n)} candidates. In the worst case this remains
 * {@code O(s n)} with {@code s} potentially large, but it avoids the practical
 * {@code n^3}-like behavior of "restart from scratch" implementations.</li>
 * <li><strong>Space:</strong> {@code O(n)} for coordinate arrays and the
 * dirty-edge stack. No {@code O(c)} storage of crossing pairs is used.</li>
 * </ul>
 *
 * <h2>Assumptions / Semantics</h2>
 * <ul>
 * <li>Input is a closed tour (last vertex connects back to first).</li>
 * <li>Polygon has at least 4 vertices.</li>
 * <li>Uses a "proper" intersection test (touching at endpoints or collinear
 * overlaps are treated as non-crossings). If you need to treat those as
 * crossings, the predicate must be adjusted.</li>
 * <li>This implementation reorders the {@code PVector} references in the list
 * (it reverses ranges of the vertex sequence), rather than shuffling
 * coordinates among fixed objects.</li>
 * </ul>
 *
 * @author Michael Carleton
 * @see <a href="https://en.wikipedia.org/wiki/2-opt">2-opt algorithm
 *      (Wikipedia)</a>
 */
public final class Uncrossing2Opt {

	// Toggle: cheap AABB overlap check before orientation math
	private static final boolean USE_AABB = !true;

	public static void uncross(List<PVector> seq) {
		final int n = seq.size();
		if (n < 4) {
			return;
		}

		// Work on arrays (fast) and write back at end
		final PVector[] arr = seq.toArray(new PVector[n]);
		final float[] x = new float[n];
		final float[] y = new float[n];
		for (int i = 0; i < n; i++) {
			x[i] = arr[i].x;
			y[i] = arr[i].y;
		}

		// Dirty-edge stack: indices of edge-starts i (edge is i -> i+1, and n-1 -> 0)
		final IntStack stack = new IntStack(n * 2);
		for (int i = 0; i < n; i++) {
			stack.push(i);
		}

		while (!stack.isEmpty()) {
			final int i = stack.pop();

			final int j = findAnyCrossingPartner(i, x, y, n);
			if (j < 0) {
				continue;
			}

			// Apply the 2-opt reversal in the appropriate contiguous range
			// and mark ONLY the two new boundary edges dirty.
			if (i == n - 1) {
				// crossing between closing edge (n-1 -> 0) and (j -> j+1): reverse [0..j]
				reverseRange(arr, x, y, 0, j);

				// Only changed geometric edges are now at indices (n-1) and (j)
				stack.push(n - 1);
				stack.push(j);

			} else if (j == n - 1) {
				// crossing between (i -> i+1) and closing edge: reverse [i+1 .. n-1]
				reverseRange(arr, x, y, i + 1, n - 1);

				// Only changed geometric edges are now at indices (i) and (n-1)
				stack.push(i);
				stack.push(n - 1);

			} else {
				// general case: normalize so p < q, reverse [p+1..q]
				int p = i, q = j;
				if (p > q) {
					int t = p;
					p = q;
					q = t;
				}

				reverseRange(arr, x, y, p + 1, q);

				// Only changed geometric edges are now at indices p and q
				stack.push(p);
				stack.push(q);
			}
		}

		// write back reordered points (fast)
		for (int i = 0; i < n; i++) {
			seq.set(i, arr[i]);
		}
	}

	/**
	 * Returns any non-adjacent edge index j such that edge(i) crosses edge(j), else
	 * -1.
	 *
	 * Heavily optimized: - avoids abs/modulo in the main candidate loops - handles
	 * wrap edges explicitly - optional AABB reject - inlined proper intersection
	 * (strict) using XOR sign tests
	 */
	private static int findAnyCrossingPartner(int i, float[] x, float[] y, int n) {
		final float ax, ay, bx, by;

		if (i == n - 1) {
			ax = x[n - 1];
			ay = y[n - 1];
			bx = x[0];
			by = y[0];
		} else {
			ax = x[i];
			ay = y[i];
			bx = x[i + 1];
			by = y[i + 1];
		}

		final float abx = bx - ax;
		final float aby = by - ay;

		final float minAx = (ax < bx) ? ax : bx;
		final float maxAx = (ax > bx) ? ax : bx;
		final float minAy = (ay < by) ? ay : by;
		final float maxAy = (ay > by) ? ay : by;

		// Candidate edges depend on i to avoid adjacency checks.
		if (i == 0) {
			// Edge (0->1): skip adjacent edges (n-1->0) and (1->2).
			for (int k = 2; k <= n - 2; k++) {
				final int k1 = k + 1;
				if (segmentsCrossFast(ax, ay, bx, by, abx, aby, minAx, maxAx, minAy, maxAy, x[k], y[k], x[k1], y[k1])) {
					return k;
				}
			}
			return -1;
		}

		if (i == n - 1) {
			// Closing edge (n-1->0): skip adjacent edges (n-2->n-1) and (0->1).
			for (int k = 1; k <= n - 3; k++) {
				final int k1 = k + 1;
				if (segmentsCrossFast(ax, ay, bx, by, abx, aby, minAx, maxAx, minAy, maxAy, x[k], y[k], x[k1], y[k1])) {
					return k;
				}
			}
			return -1;
		}

		// General i in [1..n-2]:
		// Check k in [0..i-2] and [i+2..n-2], plus possibly k=n-1 (closing edge) if not
		// adjacent.
		for (int k = 0; k <= i - 2; k++) {
			final int k1 = k + 1;
			if (segmentsCrossFast(ax, ay, bx, by, abx, aby, minAx, maxAx, minAy, maxAy, x[k], y[k], x[k1], y[k1])) {
				return k;
			}
		}

		for (int k = i + 2; k <= n - 2; k++) {
			final int k1 = k + 1;
			if (segmentsCrossFast(ax, ay, bx, by, abx, aby, minAx, maxAx, minAy, maxAy, x[k], y[k], x[k1], y[k1])) {
				return k;
			}
		}

		// Check closing edge k = n-1 (n-1 -> 0) unless adjacent (only adjacent when i
		// == n-2)
		if (i != n - 2) {
			if (segmentsCrossFast(ax, ay, bx, by, abx, aby, minAx, maxAx, minAy, maxAy, x[n - 1], y[n - 1], x[0], y[0])) {
				return n - 1;
			}
		}

		return -1;
	}

	// Proper intersection test (strict): excludes collinear/touching cases.
	private static boolean segmentsCrossFast(float ax, float ay, float bx, float by, float abx, float aby, float minAx, float maxAx, float minAy, float maxAy,
			float cx, float cy, float dx, float dy) {
		if (USE_AABB) {
			final float minCx = (cx < dx) ? cx : dx;
			final float maxCx = (cx > dx) ? cx : dx;
			if (maxAx < minCx || maxCx < minAx) {
				return false;
			}

			final float minCy = (cy < dy) ? cy : dy;
			final float maxCy = (cy > dy) ? cy : dy;
			if (maxAy < minCy || maxCy < minAy) {
				return false;
			}
		}

		// o1 and o2: C and D on opposite sides of AB?
		final float acx = cx - ax, acy = cy - ay;
		final float adx = dx - ax, ady = dy - ay;

		final float o1 = abx * acy - aby * acx;
		final float o2 = abx * ady - aby * adx;

		// strict opposite sign (o1==0 or o2==0 => reject)
		if (!((o1 > 0) ^ (o2 > 0))) {
			return false;
		}

		// o3 and o4: A and B on opposite sides of CD?
		final float cdx = dx - cx, cdy = dy - cy;
		final float cax = ax - cx, cay = ay - cy;
		final float cbx = bx - cx, cby = by - cy;

		final float o3 = cdx * cay - cdy * cax;
		final float o4 = cdx * cby - cdy * cbx;

		return ((o3 > 0) ^ (o4 > 0));
	}

	private static void reverseRange(PVector[] arr, float[] x, float[] y, int l, int r) {
		while (l < r) {
			PVector tp = arr[l];
			arr[l] = arr[r];
			arr[r] = tp;

			float tx = x[l];
			x[l] = x[r];
			x[r] = tx;
			float ty = y[l];
			y[l] = y[r];
			y[r] = ty;

			l++;
			r--;
		}
	}

	// Tiny int stack (no boxing)
	private static final class IntStack {
		private int[] a;
		private int sz;

		IntStack(int cap) {
			a = new int[Math.max(8, cap)];
		}

		boolean isEmpty() {
			return sz == 0;
		}

		void push(int v) {
			if (sz == a.length) {
				int[] b = new int[a.length << 1];
				System.arraycopy(a, 0, b, 0, a.length);
				a = b;
			}
			a[sz++] = v;
		}

		int pop() {
			return a[--sz];
		}
	}
}