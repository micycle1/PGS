package micycle.pgs.utility;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import org.jgrapht.Graph;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.util.CollectionUtil;

/**
 * The Recursive Largest First (RLF) Algorithm, from A Graph Coloring Algorithm
 * for Large Scheduling Problems, Frank Thomson Leighton: paper @
 * https://nvlpubs.nist.gov/nistpubs/jres/84/jresv84n6p489_a1b.pdf
 * <p>
 * The Recursive Largest First (RLF) algorithm was proposed in 1979 by F.
 * Leighton
 * <p>
 * It sequentially builds color classes on the basis of greedy choices. In
 * particular the first vertex placed in a color class C is one with a maximum
 * number of uncolored neighbors, and the next vertices placed in C are chosen
 * so that they have as many uncolored neighbors which cannot be placed in C. T
 * <p>
 * As was true with the SLI (smallest last with interchange) algorithm, the RLF
 * algorithm, in general, requires 0(n3) time and 0(n2) space to color an n node
 * graph. Unlike the SLI algorithm, however, the RLF algorithm requires only 0
 * (n2) time to color graphs for which k·e = n 2 where k is the number of colors
 * used to color the graph, e is the number of edges in the graph, and n is the
 * number of nodes in the graph.
 * <p>
 * Also see https://www.gerad.ca/~alainh/RLFPaper.pdf
 * <p>
 * Improved drop-in replacement heuristics for RLF are explored in 'A new
 * efficient RLF-like Algorithm for the Vertex Coloring Problem' (though most
 * increase runtime).
 * <p>
 * This implementation of this class is based on the algorithm desctiption in 'A
 * new efficient RLF-like Algorithm for the Vertex Coloring Problem'.
 * <p>
 * for practical purposes, the RLF algorithm, if programmed properly, exhibits
 * an 0(n2 ) time depend ence for many applications.
 * <p>
 * RLF exhibits similar (if not slightly better) chromatic performance compared
 * to DSATUR. In 'A Performance Comparison of Graph Coloring Algorithms Murat
 * Aslan* 1 , Nurdan Akhan Baykan' RLF performs tends to produce best colorings
 * (as measured by color number) (just marginally ahead of DSATUR).
 * 
 * @author Michael Carleton
 * 
 *         file:///C:/Users/micyc/Downloads/A_Performance_Comparison_of_Graph_Coloring_Algorit.pdf
 *
 * @param <V> the graph vertex type
 * @param <E> the graph edge type
 */
public class RLFColoring<V, E> implements VertexColoringAlgorithm<V> {

	/*
	 * TODO see 'Efficiency issues in the RLF heuristic for graph coloring' for data
	 * structure and better time complexity.
	 */

	private final Graph<V, E> graph;
	/** let U denote the set of uncolored vertices */
	private final Set<V> U;
	/**
	 * let W be the set (initially empty) of uncolored vertices with at least one
	 * neighbor in C
	 */
	private final Set<V> W;
	final Map<V, Integer> C;
	private final NeighborCache<V, E> adjacency;

	/**
	 * Use these maps to store the numbers AU (x) and AW (x) each time a vertex is
	 * removed from U (rather than calling aU lots).
	 */
	final Map<V, Integer> AU;
	final Map<V, Integer> AW;

	private int activeColor;

	public RLFColoring(Graph<V, E> graph) {
		this.graph = Objects.requireNonNull(graph, "Graph cannot be null");
		U = new HashSet<>(graph.vertexSet());
		W = new HashSet<>(graph.vertexSet().size());
		C = CollectionUtil.newHashMapWithExpectedSize(graph.vertexSet().size());
		adjacency = new NeighborCache<>(graph);
		AU = new HashMap<>(graph.vertexSet().size());
		AW = new HashMap<>(graph.vertexSet().size());
	}

	/**
	 * {@inheritDoc} NOTE notation follows https://www.gerad.ca/~alainh/RLFPaper.pdf
	 */
	@Override
	public Coloring<V> getColoring() {
		final int n = graph.vertexSet().size();

		V activeVertex = null;
		int maxDegree = 0;
		activeColor = -1; // current color class (RLF builds color classes sequentially)

		for (V v : U) {
			AU.put(v, adjacency.neighborsOf(v).size());
			AW.put(v, 0);
		}

		while (C.size() < n) { // while G contains uncolored vertices
			/*
			 * Recompute U (set of uncolored vertices).
			 */
			U.addAll(graph.vertexSet());
			U.removeAll(C.keySet());

			activeColor++;

			// Choose a vertex v ∈ U with largest value AU(v).
			// Select the uncolored vertex which has the largest degree for coloring.
			maxDegree = -1;
			for (V v : U) {
				int d = AU.get(v);
				if (d > maxDegree) {
					maxDegree = d;
					activeVertex = v;
				}
			}
			cV(activeVertex);
		}

		return new ColoringImpl<>(C, activeColor + 1);
	}

	/**
	 * Finds the neighbours of u in U (the set of uncolored vertices).
	 * 
	 * @return
	 */
	private Set<V> aU(V u) {
		// neighborsOf(v) does not include v (which is desirable)
		final Set<V> aU = new HashSet<>(adjacency.neighborsOf(u));
		aU.retainAll(U);
		return aU;
	}

	/**
	 * Constructs the color class Cv (using current color value) and assigns color k
	 * to all vertices in Cv.
	 * 
	 * @param vertex
	 */
	void cV(V vertex) {
		// Initialize W as the set of vertices in U adjacent to v
		W.clear();

		moveFromUnseen(vertex);

		while (!U.isEmpty()) {
			// Select a vertex u ∈ U with largest value AW(u)
			V candidate = findNextCandidate();

			/*
			 * Every time a vertex in U is chosen to be moved to C, all its neighbors in U
			 * are moved from U to W
			 */

			moveFromUnseen(candidate);
		}
	}

	private void moveFromUnseen(V vertex) {
		// move all neighbors w ∈ U of u to W
		final Set<V> uncoloredNeighbours = aU(vertex);
		W.addAll(uncoloredNeighbours);
		U.removeAll(uncoloredNeighbours);
		U.remove(vertex);

		C.put(vertex, activeColor); // Move u from U to C

		uncoloredNeighbours.forEach(n -> {
			final Set<V> neighbours = adjacency.neighborsOf(n); // NOTE use graph adjacency (not aU())
			/*
			 * Each time a vertex w is moved from U to W, AW(x) is incremented by one unit
			 * and AU(x) is decreased by one unit for all neighbors x ∈ U of w.
			 */
			neighbours.forEach(n2 -> {
				AW.merge(n2, 1, Integer::sum);
				AU.merge(n2, -1, Integer::sum);
			});
			/*
			 * When a vertex u ∈ U is moved from U to Cv, AU(x) is decreased by one unit for
			 * all neighbors x ∈ U of u.
			 */
			AU.merge(n, -1, Integer::sum);
		});
	}

	/**
	 * Find next vertex that should belong to the current color class. The next
	 * vertex to be moved from U to C is one having the largest number of neighbors
	 * in W (the set of uncolored vertices with at least one neighbor in C).
	 * 
	 * @return
	 */
	private V findNextCandidate() {
		V candidate = null;

		int maxDegree = -1;
		// TODO optimise this O(n) loop (treeset?).
		for (V v : U) {
			int d = AW.get(v);
			if (d > maxDegree) {
				maxDegree = d;
				candidate = v;
			}
		}
		/*
		 * NOTE Leighton specifies tie-breaking at this stage: "Ties are, if possible,
		 * broken by choosing a vertex with the smallest number of neighbors in U", but
		 * I've found that such tie-breaking consistently leads to worse colorings
		 * (higher chromatic number), so tie-breaking has not been included.
		 */

		return candidate;
	}
}
