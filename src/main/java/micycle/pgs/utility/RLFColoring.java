package micycle.pgs.utility;

import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.util.CollectionUtil;

/**
 * The Recursive Largest First (RLF) algorithm for graph coloring.
 * 
 * <p>
 * 
 * The RLF algorithm was originally designed by F. Leighton (1979) in <a href=
 * "https://nvlpubs.nist.gov/nistpubs/jres/84/jresv84n6p489_a1b.pdf"><i>A Graph
 * Coloring Algorithm for Large Scheduling Problems</i></a>, in part for use in
 * constructing solutions to large timetabling problems. It sequentially builds
 * color classes on the basis of greedy choices. In particular, the first vertex
 * placed in a color class C is one with a maximum number of uncolored
 * neighbors, and the next vertices placed in C are chosen so that they have as
 * many uncolored neighbors which cannot be placed in C.
 * <p>
 * This implementation is based on the algorithm description in 'A new efficient
 * RLF-like Algorithm for the Vertex Coloring Problem' : "for practical
 * purposes, the RLF algorithm, if programmed properly, exhibits an
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
	private final NeighborCache<V, E> neighborCache;

	/**
	 * These maps to store the numbers AU(x) and AW(x). Each time a vertex is
	 * removed from U (and its uncolored neighbours added to W) they are updated
	 * (provided an efficient computing the value each call).
	 */
	final Map<V, Integer> AU;
	final Map<V, Integer> AW;

	private int activeColor;

	public RLFColoring(Graph<V, E> graph) {
		this.graph = Objects.requireNonNull(graph, "Graph cannot be null");
		U = new HashSet<>(graph.vertexSet());
		W = new HashSet<>(graph.vertexSet().size());
		C = CollectionUtil.newHashMapWithExpectedSize(graph.vertexSet().size());
		neighborCache = new NeighborCache<>(graph);
		AU = CollectionUtil.newHashMapWithExpectedSize(graph.vertexSet().size());
		AW = CollectionUtil.newHashMapWithExpectedSize(graph.vertexSet().size());
	}

	@Override
	public Coloring<V> getColoring() {
		final int n = graph.vertexSet().size();

		V nextClassVertex = null;
		int maxDegree = 0;
		activeColor = -1; // current color class (RLF builds color classes sequentially)

		for (V v : U) {
			AU.put(v, neighborCache.neighborsOf(v).size());
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
			maxDegree = Integer.MIN_VALUE;
			for (V v : U) {
				int d = AU.get(v);
				if (d > maxDegree) {
					maxDegree = d; // may be negative (which is fine)
					nextClassVertex = v;
				}
			}

			createColorClass(nextClassVertex);
		}

		return new ColoringImpl<>(C, activeColor + 1);
	}

	/**
	 * Constructs the color class Cv (using current color value) and assigns color k
	 * to all vertices in Cv.
	 * 
	 * @param vertex inital vertex of color class
	 */
	void createColorClass(V vertex) {
		// Initialize W as the set of vertices in U adjacent to v
		W.clear();

		color(vertex);

		while (!U.isEmpty()) {
			V candidate = findNextCandidate(); // select a vertex u ∈ U with largest value AW(u)
			color(candidate);
		}
	}

	/**
	 * Colors the given vertex with the current color class.
	 */
	private void color(V vertex) {
		// Move all neighbors w ∈ U of u to W
		final Set<V> uncoloredNeighbours = findUncoloredNeighbours(vertex);
		W.addAll(uncoloredNeighbours);
		U.removeAll(uncoloredNeighbours);

		// Move u from U to C
		U.remove(vertex);
		C.put(vertex, activeColor);

		uncoloredNeighbours.forEach(n -> {
			final Set<V> neighbours = neighborCache.neighborsOf(n); // NOTE use graph adjacency (not uncolored neighbours)
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
		// TODO optimise this O(n) loop
		// look at https://imada.sdu.dk/~marco/Publications/Files/MIC2011-ChiGalGua.pdf
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

	/**
	 * Finds the neighbours of u in U (the set of uncolored vertices).
	 * 
	 * @return set of uncolored neighbouring vertices of u
	 */
	private Set<V> findUncoloredNeighbours(V u) {
		// neighborsOf(v) does not include v (which is desirable)
		final Set<V> aU = new HashSet<>(neighborCache.neighborsOf(u));
		aU.retainAll(U);
		return aU;
	}
}
