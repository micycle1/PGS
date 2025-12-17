package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.jgrapht.Graph;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm;

/**
 * DBLAC (Degree-Based Largest Adjacency Count) graph coloring.
 * <p>
 * <strong>DBLAC selection rule</strong>:
 * </p>
 * <p>
 * Repeatedly select an uncolored vertex <code>v</code> that maximises
 * </p>
 * 
 * <pre>
 * LAC(v) = number of already-colored neighbors of v
 * </pre>
 * <p>
 * Tie-break by larger static degree (original degree), then by the shuffled
 * index.
 * </p>
 * <p>
 * Color the selected vertex using first-fit (smallest feasible color).
 * </p>
 *
 * @author Michael Carleton
 */
public class DBLACColoring<V, E> implements VertexColoringAlgorithm<V> {

	/*-
	 * Implementation written for maximum performance (not readability):
	 *  - CSR/flat adjacency (offsets + neighbors array)
	 *  - adjacency built from edgeSet() (no NeighborCache/Set allocations)
	 *  - max-heap with increase-key (no O(n) scans)
	 *  - first-fit using 64-bit mask for <=64 colors, else BitSet + touched-list
	 */

	private final Random rnd;

	private final List<V> vertexList; // shuffled
	private final Map<V, Integer> vertexIndex;
	private final int n;

	// CSR adjacency: neighbors in nbrs[off[v]..off[v+1])
	private final int[] off;
	private final int[] nbrs;

	// static degree (unique neighbors)
	private final int[] degTotal;

	// dynamic: #colored neighbors
	private final int[] lac;

	// color[v] = -1 if uncolored
	private final int[] color;

	// Max-heap over uncolored vertices by (lac desc, degTotal desc, index asc)
	private final int[] heap;
	private final int[] posInHeap; // -1 if removed
	private int heapSize;

	// first-fit for many colors
	private final BitSet forbidden = new BitSet();
	private int[] touched; // colors set in forbidden for current vertex
	private int touchedSize;

	private Coloring<V> cached;

	public DBLACColoring(Graph<V, E> graph, long seed) {
		this.rnd = new Random(seed);

		this.n = graph.vertexSet().size();

		// Vertex indexing (shuffle for random tie-breaking)
		this.vertexList = new ArrayList<>(graph.vertexSet());
		Collections.shuffle(vertexList, rnd);

		this.vertexIndex = new HashMap<>(Math.max(16, n * 2));
		for (int i = 0; i < n; i++) {
			vertexIndex.put(vertexList.get(i), i);
		}

		// Build CSR adjacency from edgeSet (2-pass), then deduplicate per vertex
		// (stamp-based).
		// Pass 1: count degrees (with duplicates)
		int[] degDup = new int[n];
		for (E e : graph.edgeSet()) {
			int s = vertexIndex.get(graph.getEdgeSource(e));
			int t = vertexIndex.get(graph.getEdgeTarget(e));
			degDup[s]++;
			degDup[t]++;
		}

		int[] offDup = new int[n + 1];
		for (int i = 0; i < n; i++) {
			offDup[i + 1] = offDup[i] + degDup[i];
		}

		int[] nbrsDup = new int[offDup[n]];
		int[] cur = offDup.clone();

		// Pass 2: fill adjacency (with duplicates)
		for (E e : graph.edgeSet()) {
			int s = vertexIndex.get(graph.getEdgeSource(e));
			int t = vertexIndex.get(graph.getEdgeTarget(e));
			nbrsDup[cur[s]++] = t;
			nbrsDup[cur[t]++] = s;
		}

		// Deduplicate per vertex without sorting (O(n+m)) using stamps
		int[] seen = new int[n];
		int stamp = 1;

		int[] degUniq = new int[n];
		for (int v = 0; v < n; v++) {
			stamp++;
			if (stamp == 0) {
				Arrays.fill(seen, 0);
				stamp = 1;
			} // ultra-defensive
			int a = offDup[v], b = offDup[v + 1];
			int cnt = 0;
			for (int p = a; p < b; p++) {
				int nb = nbrsDup[p];
				if (seen[nb] != stamp) {
					seen[nb] = stamp;
					cnt++;
				}
			}
			degUniq[v] = cnt;
		}

		int[] offUniq = new int[n + 1];
		for (int i = 0; i < n; i++) {
			offUniq[i + 1] = offUniq[i] + degUniq[i];
		}

		int[] nbrsUniq = new int[offUniq[n]];
		stamp = 1;
		for (int v = 0; v < n; v++) {
			stamp++;
			if (stamp == 0) {
				Arrays.fill(seen, 0);
				stamp = 1;
			}
			int write = offUniq[v];
			int a = offDup[v], b = offDup[v + 1];
			for (int p = a; p < b; p++) {
				int nb = nbrsDup[p];
				if (seen[nb] != stamp) {
					seen[nb] = stamp;
					nbrsUniq[write++] = nb;
				}
			}
		}

		this.off = offUniq;
		this.nbrs = nbrsUniq;
		this.degTotal = degUniq;

		int maxDeg = 0;
		for (int d : degTotal) {
			maxDeg = Math.max(maxDeg, d);
		}

		this.lac = new int[n];
		this.color = new int[n];
		Arrays.fill(color, -1);

		// Heap init: all vertices uncolored
		this.heap = new int[n];
		this.posInHeap = new int[n];
		for (int i = 0; i < n; i++) {
			heap[i] = i;
			posInHeap[i] = i;
		}
		this.heapSize = n;
		for (int i = (heapSize >>> 1) - 1; i >= 0; i--) {
			siftDown(i);
		}

		this.touched = new int[Math.max(16, maxDeg)]; // typical touched size ~ degree
	}

	public DBLACColoring(Graph<V, E> graph) {
		this(graph, System.nanoTime());
	}

	@Override
	public Coloring<V> getColoring() {
		if (cached != null) {
			return cached;
		}

		if (n == 0) {
			cached = new ColoringImpl<>(Collections.emptyMap(), 0);
			return cached;
		}

		int numColors = 0;

		while (heapSize > 0) {
			int v = extractMax();

			int c = chooseSmallestAvailableColor(v, numColors);
			if (c == numColors) {
				numColors++;
			}
			color[v] = c;

			// update LAC for uncolored neighbors and fix heap keys (increase-key)
			for (int p = off[v], end = off[v + 1]; p < end; p++) {
				int nb = nbrs[p];
				if (color[nb] == -1) {
					lac[nb]++;
					increaseKey(nb);
				}
			}
		}

		Map<V, Integer> colorMap = new HashMap<>(Math.max(16, n * 2));
		for (int i = 0; i < n; i++) {
			colorMap.put(vertexList.get(i), color[i]);
		}

		cached = new ColoringImpl<>(colorMap, numColors);
		return cached;
	}

	// ---- first-fit color choice ----

	private int chooseSmallestAvailableColor(int v, int numColors) {
		if (numColors <= 64) {
			long forb = 0L;
			for (int p = off[v], end = off[v + 1]; p < end; p++) {
				int c = color[nbrs[p]];
				if (c >= 0) {
					forb |= (1L << c);
				}
			}
			long avail = ~forb;
			int c = Long.numberOfTrailingZeros(avail); // 0..64 (64 means none in [0..63])
			return (c < numColors) ? c : numColors;
		}

		// BitSet + touched list: O(deg(v)) to mark and clear; nextClearBit finds
		// smallest available.
		touchedSize = 0;
		for (int p = off[v], end = off[v + 1]; p < end; p++) {
			int c = color[nbrs[p]];
			if (c >= 0 && !forbidden.get(c)) {
				forbidden.set(c);
				if (touchedSize == touched.length) {
					touched = Arrays.copyOf(touched, touched.length << 1);
				}
				touched[touchedSize++] = c;
			}
		}

		int chosen = forbidden.nextClearBit(0);
		// clear only what we set
		for (int i = 0; i < touchedSize; i++) {
			forbidden.clear(touched[i]);
		}

		return (chosen < numColors) ? chosen : numColors;
	}

	// ---- heap (max by lac, then degree, then shuffled index) ----

	private boolean better(int a, int b) {
		int la = lac[a], lb = lac[b];
		if (la != lb) {
			return la > lb;
		}

		int da = degTotal[a], db = degTotal[b];
		if (da != db) {
			return da > db;
		}

		return a < b; // shuffled index tie-break
	}

	private void swapHeap(int i, int j) {
		int vi = heap[i], vj = heap[j];
		heap[i] = vj;
		heap[j] = vi;
		posInHeap[vi] = j;
		posInHeap[vj] = i;
	}

	private void siftUp(int i) {
		while (i > 0) {
			int p = (i - 1) >>> 1;
			if (better(heap[p], heap[i])) {
				break;
			}
			swapHeap(p, i);
			i = p;
		}
	}

	private void siftDown(int i) {
		for (;;) {
			int l = (i << 1) + 1;
			if (l >= heapSize) {
				return;
			}
			int r = l + 1;

			int bestChild = (r < heapSize && better(heap[r], heap[l])) ? r : l;
			if (better(heap[i], heap[bestChild])) {
				return;
			}

			swapHeap(i, bestChild);
			i = bestChild;
		}
	}

	private int extractMax() {
		int v = heap[0];
		posInHeap[v] = -1;

		int last = heap[--heapSize];
		if (heapSize > 0) {
			heap[0] = last;
			posInHeap[last] = 0;
			siftDown(0);
		}
		return v;
	}

	private void increaseKey(int v) {
		int p = posInHeap[v];
		if (p >= 0) {
			siftUp(p);
		}
	}
}