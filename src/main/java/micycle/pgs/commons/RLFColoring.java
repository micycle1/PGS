package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import org.jgrapht.Graph;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm;
import org.jgrapht.alg.util.NeighborCache;

/**
 * The Recursive Largest First (RLF) algorithm for graph coloring.
 * <p>
 * The RLF algorithm was originally designed by F. Leighton (1979) in <a href=
 * "https://nvlpubs.nist.gov/nistpubs/jres/84/jresv84n6p489_a1b.pdf"><i>A Graph
 * Coloring Algorithm for Large Scheduling Problems</i></a>, in part for use in
 * constructing solutions to large timetabling problems. It sequentially builds
 * color classes on the basis of greedy choices. In particular, the first vertex
 * placed in a color class C is one with a maximum number of uncolored
 * neighbors, and the next vertices placed in C are chosen so that they have as
 * many uncolored neighbors which cannot be placed in C.
 * <p>
 * This implementation is based on the original algorithm pseudocode provided in
 * 'A new efficient RLF-like Algorithm for the Vertex Coloring Problem' : "for
 * practical purposes, the RLF algorithm, if programmed properly, exhibits an
 * O(n<sup>2</sup>) time dependence for many applications".
 * <p>
 * RLF exhibits similar chromatic performance compared to DSATUR. In <i>'A
 * Performance Comparison of Graph Coloring Algorithms'</i> RLF tended to
 * produce the best colorings (as measured by color number), marginally ahead of
 * DSATUR.
 * <p>
 * Improved drop-in replacement heuristics for RLF are explored in <i>'A new
 * efficient RLF-like Algorithm for the Vertex Coloring Problem'</i> (though
 * most increase runtime complexity).
 * 
 * @author Michael Carleton
 * 
 * @param <V> the graph vertex type
 * @param <E> the graph edge type
 */
public class RLFColoring<V, E> implements VertexColoringAlgorithm<V> {

	private final Graph<V, E> graph;
	private final List<V> vertexList;
	private final Map<V, Integer> vertexIndex;
	private final int n;

	private final BitSet U; // uncolored vertices
	private final BitSet W; // uncolored vertices with neighbor in C
	private final Map<V, Integer> C; // colored vertices

	/**
	 * An array storing the number of uncolored neighbors for each vertex (AU(x) in
	 * the RLF algorithm).
	 */
	private final int[] AU;
	/**
	 * An array storing the number of neighbors in the set W for each vertex (AW(x)
	 * in the RLF algorithm).
	 */
	private final int[] AW;

	// Pre-computed adjacency lists as indices
	private final int[][] adjacency;
	private final NeighborCache<V, E> neighborCache;

	private int activeColor;

	public RLFColoring(Graph<V, E> graph, long seed) {
		this.graph = Objects.requireNonNull(graph);
		this.n = graph.vertexSet().size();
		this.neighborCache = new NeighborCache<>(graph);

		// Create vertex indexing
		vertexList = new ArrayList<>(graph.vertexSet());
		Collections.shuffle(vertexList, new Random(seed));
		vertexIndex = new HashMap<>(n);
		for (int i = 0; i < n; i++) {
			vertexIndex.put(vertexList.get(i), i);
		}

		U = new BitSet(n);
		U.set(0, n); // all initially uncolored
		W = new BitSet(n);

		AU = new int[n];
		AW = new int[n];
		C = new HashMap<>(n);

		// Pre-compute adjacency as indices
		adjacency = new int[n][];
		for (int i = 0; i < n; i++) {
			V v = vertexList.get(i);
			Set<V> neighbors = neighborCache.neighborsOf(v);
			adjacency[i] = new int[neighbors.size()];
			int j = 0;
			for (V neighbor : neighbors) {
				adjacency[i][j++] = vertexIndex.get(neighbor);
			}
			AU[i] = adjacency[i].length;
		}
	}

	@Override
	public Coloring<V> getColoring() {
		activeColor = -1;

		while (C.size() < n) {
			// Recompute U
			U.set(0, n);
			for (V v : C.keySet()) {
				U.clear(vertexIndex.get(v));
			}

			activeColor++;

			// Find vertex with largest AU value
			int maxAU = Integer.MIN_VALUE;
			V nextVertex = null;
			for (int i = U.nextSetBit(0); i >= 0; i = U.nextSetBit(i + 1)) {
				if (AU[i] > maxAU) {
					maxAU = AU[i];
					nextVertex = vertexList.get(i);
				}
			}

			createColorClass(nextVertex);
		}

		return new ColoringImpl<>(C, activeColor + 1);
	}

	/**
	 * Constructs the color class Cv (using current color value) and assigns color k
	 * to all vertices in Cv.
	 * 
	 * @param vertex inital vertex of color class
	 */
	private void createColorClass(V initialVertex) {
		W.clear();
		Arrays.fill(AW, 0); // Reset AW for this color class

		color(initialVertex);

		while (U.cardinality() > 0) {
			V candidate = findNextCandidate();
			if (candidate == null)
				break;
			color(candidate);
		}
	}

	/**
	 * Find next vertex that should belong to the current color class. The next
	 * vertex to be moved from U to C is one having the largest number of neighbors
	 * in W (the set of uncolored vertices with at least one neighbor in C).
	 */
	private V findNextCandidate() {
		V candidate = null;
		int maxAW = -1;

		// Find vertex in U with maximum AW value
		for (int i = U.nextSetBit(0); i >= 0; i = U.nextSetBit(i + 1)) {
			if (AW[i] > maxAW) {
				maxAW = AW[i];
				candidate = vertexList.get(i);
			}
		}

		return candidate;
	}

	/**
	 * Colors the given vertex with the current color class.
	 */
	private void color(V vertex) {
		int vIdx = vertexIndex.get(vertex);

		// Get uncolored neighbors efficiently
		BitSet uncoloredNeighbors = new BitSet(n);
		for (int nIdx : adjacency[vIdx]) {
			if (U.get(nIdx)) {
				uncoloredNeighbors.set(nIdx);
			}
		}

		// Move uncolored neighbors to W
		for (int nIdx = uncoloredNeighbors.nextSetBit(0); nIdx >= 0; nIdx = uncoloredNeighbors.nextSetBit(nIdx + 1)) {

			W.set(nIdx);
			U.clear(nIdx);

			// Update AW and AU for neighbors of this moved vertex
			for (int nnIdx : adjacency[nIdx]) {
				if (U.get(nnIdx)) {
					AW[nnIdx]++;
					AU[nnIdx]--;
				}
			}

			// Update AU for the moved vertex itself
			AU[nIdx]--;
		}

		// Move vertex from U to C
		U.clear(vIdx);
		C.put(vertex, activeColor);

		// Update AU for remaining neighbors
		for (int nIdx : adjacency[vIdx]) {
			if (U.get(nIdx)) {
				AU[nIdx]--;
			}
		}
	}
}