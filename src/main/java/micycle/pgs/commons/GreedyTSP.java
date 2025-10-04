package micycle.pgs.commons;

import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.function.ToDoubleBiFunction;
import java.util.stream.Collectors;

/**
 * A high-performance implementation of the Traveling Salesman Problem (TSP)
 * using a greedy construction heuristic followed by 2-opt local search
 * improvement.
 *
 * <h2>Algorithm Overview</h2>
 * <p>
 * This implementation uses a two-phase approach:
 * </p>
 * <ol>
 * <li><strong>Greedy Construction:</strong> Builds an initial tour by
 * repeatedly selecting the shortest available edge that doesn't violate TSP
 * constraints (no cycles except the final one, maximum degree 2 per
 * vertex)</li>
 * <li><strong>2-opt Improvement:</strong> Iteratively improves the tour by
 * swapping edges until no further improvement is possible</li>
 * </ol>
 *
 * @param <V> the type of vertices in the graph. Can be any type for which
 *            distances can be computed.
 *
 * @author Michael Carleton
 */
public class GreedyTSP<V> {

	private final List<V> vertices;
	private final ToDoubleBiFunction<V, V> distFunc;
	private final double[][] allDist;

	public GreedyTSP(List<V> vertices, ToDoubleBiFunction<V, V> distFunc) {
		if (vertices == null || vertices.isEmpty()) {
			throw new IllegalArgumentException("Vertex list must not be null or empty");
		}
		this.vertices = List.copyOf(vertices);
		this.distFunc = distFunc;
		this.allDist = initDistanceTable();
	}

	/**
	 * Build the full symmetric distance matrix.
	 */
	private double[][] initDistanceTable() {
		int n = vertices.size();
		double[][] d = new double[n][n];
		for (int i = 0; i < n; i++) {
			d[i][i] = 0;
			for (int j = i + 1; j < n; j++) {
				double dij = distFunc.applyAsDouble(vertices.get(i), vertices.get(j));
				d[i][j] = dij;
				d[j][i] = dij;
			}
		}
		return d;
	}

	/**
	 * Runs greedy construction heuristic, then improves with 2-opt, and returns a
	 * CLOSED tour (first vertex repeated at end).
	 */
	public List<V> getTour() {
		int n = vertices.size();
		if (n == 1) {
			return List.of(vertices.get(0), vertices.get(0));
		}
		if (n == 2) {
			V a = vertices.get(0), b = vertices.get(1);
			return List.of(a, b, a);
		}

		// 1) Build tour using greedy edge selection
		int[] tour = buildGreedyTour();

		// 2) Improve with 2-opt
		improve(tour);

		// 3) Map back to V
		return Arrays.stream(tour).mapToObj(vertices::get).collect(Collectors.toList());
	}

	/**
	 * Edge record for efficient immutable edge representation.
	 */
	private record Edge(int u, int v, double weight) implements Comparable<Edge> {
		@Override
		public int compareTo(Edge other) {
			return Double.compare(this.weight, other.weight);
		}
	}

	/**
	 * Build tour using greedy edge selection with optimizations.
	 */
	private int[] buildGreedyTour() {
		int n = vertices.size();

		// Pre-allocate exact capacity
		int edgeCount = n * (n - 1) / 2;
		Edge[] edges = new Edge[edgeCount];
		int idx = 0;

		// Create edges array directly (avoid List overhead)
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				edges[idx++] = new Edge(i, j, allDist[i][j]);
			}
		}
		Arrays.sort(edges);

		// Use byte array for degrees (max degree is 2)
		byte[] degree = new byte[n];
		UnionFind uf = new UnionFind(n);

		// Pre-allocate adjacency lists with exact capacity (2)
		int[][] adj = new int[n][2];
		for (int i = 0; i < n; i++) {
			adj[i][0] = adj[i][1] = -1;
		}

		int edgesAdded = 0;
		for (Edge e : edges) {
			// Fast degree check
			// Fast cycle check (skip for last edge)
			if (degree[e.u] == 2 || degree[e.v] == 2 || (edgesAdded < n - 1 && uf.connected(e.u, e.v))) {
				continue;
			}

			// Add edge to adjacency (no list needed, max 2 neighbors)
			adj[e.u][degree[e.u]] = e.v;
			adj[e.v][degree[e.v]] = e.u;
			degree[e.u]++;
			degree[e.v]++;
			uf.union(e.u, e.v);

			if (++edgesAdded == n) {
				break;
			}
		}

		// Convert adjacency representation to tour array
		return buildTourFromAdjacency(adj);
	}

	/**
	 * Convert adjacency array representation to tour array. Optimized to avoid list
	 * operations.
	 */
	private int[] buildTourFromAdjacency(int[][] adj) {
		int n = vertices.size();
		int[] tour = new int[n + 1];

		// Use bitset for visited tracking (more cache-friendly)
		BitSet visited = new BitSet(n);

		tour[0] = 0;
		visited.set(0);
		int current = 0;
		int prev = -1;

		// Follow the path (each vertex has exactly 2 neighbors)
		for (int i = 1; i < n; i++) {
			int next = adj[current][0];
			if (next == prev || visited.get(next)) {
				next = adj[current][1];
			}
			tour[i] = next;
			visited.set(next);
			prev = current;
			current = next;
		}

		tour[n] = 0; // close the tour
		return tour;
	}

	/**
	 * Optimized Union-Find with path compression and union by rank.
	 */
	private static class UnionFind {
		private final int[] parent;
		private final byte[] rank; // rank never exceeds log(n)

		UnionFind(int n) {
			parent = new int[n];
			rank = new byte[n];
			for (int i = 0; i < n; i++) {
				parent[i] = i;
			}
		}

		int find(int x) {
			int root = x;
			// Find root
			while (parent[root] != root) {
				root = parent[root];
			}
			// Path compression
			while (x != root) {
				int next = parent[x];
				parent[x] = root;
				x = next;
			}
			return root;
		}

		boolean connected(int x, int y) {
			return find(x) == find(y);
		}

		void union(int x, int y) {
			int px = find(x);
			int py = find(y);
			if (px == py) {
				return;
			}

			// Union by rank
			if (rank[px] < rank[py]) {
				parent[px] = py;
			} else if (rank[px] > rank[py]) {
				parent[py] = px;
			} else {
				parent[py] = px;
				rank[px]++;
			}
		}
	}

	/**
	 * Improve tour with 2-opt. Optimized with early termination and better cache
	 * patterns.
	 */
	private void improve(int[] tour) {
		int N = tour.length - 1;
		double minImprovement = 1e-9;
		int stallCount = 0;
		int maxStalls = 3; // stop after 3 rounds with tiny improvements

		while (true) {
			double bestDelta = 0;
			int bestI = -1, bestJ = -1;

			// Cache-friendly iteration pattern
			for (int i = 0; i < N - 2; i++) {
				int ci = tour[i], ci1 = tour[i + 1];
				double currentEdge = allDist[ci][ci1];

				// Start j from i+2 to avoid adjacent edges
				for (int j = i + 2; j < N; j++) {
					int cj = tour[j], cj1 = tour[j + 1];

					// Quick calculation with early exit
					double newEdges = allDist[ci][cj] + allDist[ci1][cj1];
					double oldEdges = currentEdge + allDist[cj][cj1];
					double delta = newEdges - oldEdges;

					if (delta < bestDelta) {
						bestDelta = delta;
						bestI = i;
						bestJ = j;
					}
				}
			}

			if (bestDelta < -minImprovement) {
				// Apply the improvement
				reverse(tour, bestI + 1, bestJ);
				stallCount = 0;
			} else if (bestDelta < 0) {
				// Very small improvement
				reverse(tour, bestI + 1, bestJ);
				if (++stallCount >= maxStalls) {
					break;
				}
			} else {
				// No improvement found
				break;
			}
		}
	}

	/**
	 * Optimized in-place reverse using XOR swap for primitives.
	 */
	private void reverse(int[] tour, int from, int to) {
		while (from < to) {
			// XOR swap (avoids temp variable)
			tour[from] ^= tour[to];
			tour[to] ^= tour[from];
			tour[from] ^= tour[to];
			from++;
			to--;
		}
	}
}