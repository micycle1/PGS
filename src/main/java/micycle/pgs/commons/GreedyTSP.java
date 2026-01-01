package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.ToDoubleBiFunction;

/**
 * Heuristic Traveling Salesman Problem (TSP) tour builder for a complete,
 * weighted graph.
 * <p>
 * Given a list of vertices and a distance function, this class computes a fast,
 * high-quality <em>approximate</em> Hamiltonian cycle (tour) that visits every
 * vertex exactly once and returns to the start. The solution is not guaranteed
 * to be optimal.
 * </p>
 *
 * <h2>Input / Graph Model</h2>
 * <ul>
 * <li>The input is treated as a complete graph over the provided vertices.</li>
 * <li>Distances are assumed to be symmetric. This implementation precomputes a
 * symmetric distance table by evaluating {@code distFunc} only for
 * {@code i < j} and mirroring the result. If {@code distFunc} is asymmetric,
 * the effective distance used will be {@code d(i,j)=d(j,i)=distFunc(i,j)} for
 * {@code i<j} (i.e., it will be implicitly symmetrized by the evaluation
 * order).</li>
 * <li>{@code distFunc} should be deterministic, side-effect free, and should
 * not return NaN. Returning NaN or extreme values may lead to poor tours or
 * undefined behavior.</li>
 * </ul>
 *
 * <h2>Output</h2>
 * <ul>
 * <li>{@link #getTour()} returns a <strong>closed</strong> tour: a
 * {@code List<V>} of length {@code n+1} where the first vertex is repeated at
 * the end.</li>
 * <li>The returned tour is anchored to vertex index {@code 0} (the first input
 * vertex) as the start/end point.</li>
 * </ul>
 *
 * <h2>Determinism and Thread Safety</h2>
 * <ul>
 * <li>For a fixed vertex order and deterministic {@code distFunc}, the produced
 * tour is deterministic.</li>
 * <li>Instances are effectively immutable after construction. Concurrent calls
 * to {@link #getTour()} are safe provided {@code distFunc} itself is
 * thread-safe and has no side effects.</li>
 * </ul>
 *
 * @param <V> vertex type
 * @author Michael Carleton
 */
public final class GreedyTSP<V> {

	// Tuning knobs (good defaults)
	private static final int CANDIDATES_K = 24; // 16..40 common; higher => better, slower
	private static final int RESTARTS_SMALL_N = 8; // more restarts => better, slower
	private static final int RESTARTS_LARGE_N = 4;
	private static final double EPS = 1e-12;

	private final Object[] verts;
	private final ToDoubleBiFunction<V, V> distFunc;

	private final int n;
	private final double[] dist; // flat n*n
	private final int[] rowBase;

	// Candidate lists: cand[i*k + t] = t-th nearest neighbor of i (sorted by
	// distance)
	private final int k;
	private final int[] cand;

	public GreedyTSP(Collection<V> vertices, ToDoubleBiFunction<V, V> distFunc) {
		if (vertices == null || vertices.isEmpty()) {
			throw new IllegalArgumentException("Vertex list must not be null or empty");
		}
		this.verts = vertices.toArray();
		this.distFunc = distFunc;
		this.n = verts.length;

		this.rowBase = new int[n];
		for (int i = 0; i < n; i++) {
			rowBase[i] = i * n;
		}

		this.dist = new double[n * n];
		initDistanceTable();

		this.k = Math.min(CANDIDATES_K, Math.max(0, n - 1));
		this.cand = (k == 0) ? new int[0] : buildCandidateLists(k);
	}

	@SuppressWarnings("unchecked")
	private void initDistanceTable() {
		for (int i = 0; i < n; i++) {
			final int ri = rowBase[i];
			dist[ri + i] = 0.0;
			final V vi = (V) verts[i];
			for (int j = i + 1; j < n; j++) {
				final double dij = distFunc.applyAsDouble(vi, (V) verts[j]);
				dist[ri + j] = dij;
				dist[rowBase[j] + i] = dij;
			}
		}
	}

	/**
	 * Returns a CLOSED tour (first vertex repeated at end).
	 */
	@SuppressWarnings("unchecked")
	public List<V> getTour() {
		if (n == 1) {
			final V a = (V) verts[0];
			return List.of(a, a);
		}
		if (n == 2) {
			final V a = (V) verts[0];
			final V b = (V) verts[1];
			return List.of(a, b, a);
		}

		final int restarts = (n <= 2000) ? RESTARTS_SMALL_N : RESTARTS_LARGE_N;

		final int[] bestNext = new int[n];
		double bestLen = Double.POSITIVE_INFINITY;

		// Working buffers reused across restarts (minimize allocations)
		final int[] next = new int[n];
		final int[] prev = new int[n];
		final int[] order = new int[n];

		final int[] visitedStamp = new int[n];
		int stamp = 1;

		final byte[] dlb2 = new byte[n];
		final byte[] dlbR = new byte[n];

		// Deterministic seed set (keeps runs reproducible)
		final int[] seeds = makeSeeds(restarts);

		for (int seed : seeds) {
			if (++stamp == 0) { // extremely unlikely, but keep safe
				Arrays.fill(visitedStamp, 0);
				stamp = 1;
			}

			buildNearestNeighborTour(seed, order, next, prev, visitedStamp, stamp);
			localSearch(next, prev, dlb2, dlbR);

			final double len = tourLength(next);
			if (len < bestLen) {
				bestLen = len;
				System.arraycopy(next, 0, bestNext, 0, n);
			}
		}

		// Materialize as CLOSED tour starting from vertex 0
		final ArrayList<V> out = new ArrayList<>(n + 1);
		int cur = 0;
		for (int i = 0; i < n; i++) {
			out.add((V) verts[cur]);
			cur = bestNext[cur];
		}
		out.add((V) verts[0]);
		return out;
	}

	private int[] makeSeeds(int restarts) {
		final int[] seeds = new int[Math.min(restarts, n)];

		// Always include 0
		seeds[0] = 0;
		int count = 1;
		if (count == seeds.length) {
			return seeds;
		}

		// Add farthest-from-0 (often a good diversification)
		int far = 1;
		double farD = dist[rowBase[0] + 1];
		for (int i = 2; i < n; i++) {
			final double d = dist[rowBase[0] + i];
			if (d > farD) {
				farD = d;
				far = i;
			}
		}
		seeds[count++] = far;
		if (count == seeds.length) {
			return seeds;
		}

		// Add farthest-from-far
		int far2 = 0;
		double far2D = dist[rowBase[far] + 0];
		for (int i = 1; i < n; i++) {
			final double d = dist[rowBase[far] + i];
			if (d > far2D) {
				far2D = d;
				far2 = i;
			}
		}
		seeds[count++] = far2;
		if (count == seeds.length) {
			return seeds;
		}

		// Fill remaining deterministically (hash-like spread)
		for (int i = count; i < seeds.length; i++) {
			final long x = (i * 0x9E3779B97F4A7C15L);
			seeds[i] = (int) Long.remainderUnsigned(x, n);
		}
		return seeds;
	}

	/**
	 * Builds an NN tour; uses candidate list first, falls back to full scan if
	 * needed.
	 */
	private void buildNearestNeighborTour(int seed, int[] order, int[] next, int[] prev, int[] visitedStamp, int stamp) {

		order[0] = seed;
		visitedStamp[seed] = stamp;

		int cur = seed;
		for (int pos = 1; pos < n; pos++) {
			int best = -1;
			double bestD = Double.POSITIVE_INFINITY;

			// Fast attempt: search among k nearest candidates
			if (k != 0) {
				final int base = cur * k;
				final int rc = rowBase[cur];
				for (int t = 0; t < k; t++) {
					final int candNode = cand[base + t];
					if (visitedStamp[candNode] == stamp) {
						continue;
					}
					final double d = dist[rc + candNode];
					best = candNode;
					bestD = d;
					break; // candidates are sorted by distance
				}
			}

			// Fallback: exact NN (full scan)
			if (best < 0) {
				final int rc = rowBase[cur];
				for (int j = 0; j < n; j++) {
					if (visitedStamp[j] == stamp) {
						continue;
					}
					final double d = dist[rc + j];
					if (d < bestD) {
						bestD = d;
						best = j;
					}
				}
			}

			order[pos] = best;
			visitedStamp[best] = stamp;
			cur = best;
		}

		// Convert order[] into next/prev cycle
		for (int i = 0; i < n - 1; i++) {
			final int a = order[i];
			final int b = order[i + 1];
			next[a] = b;
			prev[b] = a;
		}
		final int first = order[0];
		final int last = order[n - 1];
		next[last] = first;
		prev[first] = last;
	}

	private void localSearch(int[] next, int[] prev, byte[] dlb2, byte[] dlbR) {
		Arrays.fill(dlb2, (byte) 0);
		Arrays.fill(dlbR, (byte) 0);

		boolean improved;
		do {
			improved = false;

			// 2-opt phase
			boolean changed2;
			do {
				changed2 = false;
				for (int a = 0; a < n; a++) {
					if (dlb2[a] != 0) {
						continue;
					}
					if (tryTwoOptAt(a, next, prev, dlb2)) {
						changed2 = true;
						improved = true;
					} else {
						dlb2[a] = 1;
					}
				}
			} while (changed2);

			// Relocation phase (Or-opt-1)
			boolean changedR;
			do {
				changedR = false;
				for (int x = 0; x < n; x++) {
					if (dlbR[x] != 0) {
						continue;
					}
					if (tryRelocateAt(x, next, prev, dlbR)) {
						changedR = true;
						improved = true;
					} else {
						dlbR[x] = 1;
					}
				}
			} while (changedR);

			// After relocation, allow 2-opt bits to re-activate a bit
			if (improved) {
				Arrays.fill(dlb2, (byte) 0);
			}
		} while (improved);
	}

	private boolean tryTwoOptAt(int a, int[] next, int[] prev, byte[] dlb) {
		final int b = next[a];
		final int ra = rowBase[a];
		final int rb = rowBase[b];

		final double dab = dist[ra + b];

		// Search c among candidate neighbors of a (sorted nearest-first)
		if (k == 0) {
			return false;
		}
		final int base = a * k;

		for (int t = 0; t < k; t++) {
			final int c = cand[base + t];
			if (c == a || c == b) {
				continue;
			}

			final int d = next[c];
			if (d == a || d == b)
			 {
				continue; // edges share an endpoint -> invalid 2-opt
			}

			final double delta = (dist[ra + c] + dist[rb + d]) - (dab + dist[rowBase[c] + d]);
			if (delta < -EPS) {
				twoOptSwap(a, b, c, d, next, prev);
				clearDlbAround(dlb, a, b, c, d, next, prev);
				return true;
			}
		}
		return false;
	}

	/**
	 * 2-opt: remove (a,b) and (c,d), add (a,c) and (b,d) reversing segment [b..c].
	 */
	private static void twoOptSwap(int a, int b, int c, int d, int[] next, int[] prev) {
		// Reverse pointers along the path from b to c following next[]
		int x = b;
		while (true) {
			final int nx = next[x];
			final int px = prev[x];
			next[x] = px;
			prev[x] = nx;
			if (x == c) {
				break;
			}
			x = nx;
		}

		// Reconnect endpoints
		next[a] = c;
		prev[c] = a;

		next[b] = d;
		prev[d] = b;
	}

	/**
	 * Or-opt-1: remove node x and insert it after a (between a and b=next[a]).
	 */
	private boolean tryRelocateAt(int x, int[] next, int[] prev, byte[] dlb) {
		final int p = prev[x];
		final int q = next[x];

		// If x is the only node? not possible here, but keep structure safe
		if (p == x || q == x) {
			return false;
		}

		final int rx = rowBase[x];
		final int rp = rowBase[p];

		final double dpx = dist[rp + x];
		final double dxq = dist[rx + q];
		final double dpq = dist[rp + q];

		if (k == 0) {
			return false;
		}
		final int base = x * k;

		for (int t = 0; t < k; t++) {
			final int a = cand[base + t];
			if (a == x || a == p) {
				continue;
			}

			final int b = next[a];
			if (b == x || b == q)
			 {
				continue; // would reinsert into same place / adjacent issues
			}

			// delta = (p,q) + (a,x) + (x,b) - (p,x) - (x,q) - (a,b)
			final double delta = dpq + dist[rowBase[a] + x] + dist[rx + b] - dpx - dxq - dist[rowBase[a] + b];

			if (delta < -EPS) {
				// remove x
				next[p] = q;
				prev[q] = p;

				// insert x after a
				next[a] = x;
				prev[x] = a;
				next[x] = b;
				prev[b] = x;

				clearDlbAround(dlb, x, p, q, a, next, prev);
				return true;
			}
		}

		return false;
	}

	private static void clearDlbAround(byte[] dlb, int a, int b, int c, int d, int[] next, int[] prev) {
		// Clear a few affected nodes + their immediate neighbors (cheap and effective)
		clear(dlb, a);
		clear(dlb, b);
		clear(dlb, c);
		clear(dlb, d);
		clear(dlb, next[a]);
		clear(dlb, prev[a]);
		clear(dlb, next[b]);
		clear(dlb, prev[b]);
		clear(dlb, next[c]);
		clear(dlb, prev[c]);
		clear(dlb, next[d]);
		clear(dlb, prev[d]);
	}

	private static void clear(byte[] dlb, int i) {
		if (i >= 0 && i < dlb.length) {
			dlb[i] = 0;
		}
	}

	private double tourLength(int[] next) {
		double sum = 0.0;
		int cur = 0;
		for (int i = 0; i < n; i++) {
			final int nx = next[cur];
			sum += dist[rowBase[cur] + nx];
			cur = nx;
		}
		return sum;
	}

	/**
	 * Build k-nearest candidate lists for each node (sorted nearest-first).
	 */
	private int[] buildCandidateLists(int k) {
		final int[] out = new int[n * k];
		final int[] bestIdx = new int[k];
		final double[] bestD = new double[k];

		for (int i = 0; i < n; i++) {
			Arrays.fill(bestIdx, -1);
			Arrays.fill(bestD, Double.POSITIVE_INFINITY);

			int maxPos = 0;
			double maxVal = Double.POSITIVE_INFINITY;

			final int ri = rowBase[i];
			for (int j = 0; j < n; j++) {
				if (j == i) {
					continue;
				}
				final double d = dist[ri + j];
				if (d < maxVal) {
					bestD[maxPos] = d;
					bestIdx[maxPos] = j;

					// recompute current worst
					maxPos = 0;
					maxVal = bestD[0];
					for (int t = 1; t < k; t++) {
						final double v = bestD[t];
						if (v > maxVal) {
							maxVal = v;
							maxPos = t;
						}
					}
				}
			}

			// sort bestIdx by bestD (small k => insertion sort is fine)
			for (int a = 1; a < k; a++) {
				final double kd = bestD[a];
				final int ki = bestIdx[a];
				int b = a - 1;
				while (b >= 0 && bestD[b] > kd) {
					bestD[b + 1] = bestD[b];
					bestIdx[b + 1] = bestIdx[b];
					b--;
				}
				bestD[b + 1] = kd;
				bestIdx[b + 1] = ki;
			}

			final int base = i * k;
			for (int t = 0; t < k; t++) {
				out[base + t] = bestIdx[t];
			}
		}

		return out;
	}
}