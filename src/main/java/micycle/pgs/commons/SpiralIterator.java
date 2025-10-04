package micycle.pgs.commons;

import java.util.*;

import net.jafama.FastMath;
import processing.core.PShape;
import processing.core.PVector;

/**
 * An iterator that walks the faces of a polygonal mesh in a concentric “spiral”
 * order from a given face.
 * 
 * @author Michael Carleton
 */

public class SpiralIterator implements Iterator<PShape> {

	/*-
	 * Build the next BFS‐ring in an “edge‐aware CW‐spiral” by:
	 *
	 * 1) Gathering all unvisited faces that touch any vertex of the current ring.
	 * 2) Sorting them globally by descending polar‐angle around the start‐face centroid
	 *    (i.e. a single CW sort).
	 * 3) Splitting that sorted list into connected components under edge‐adjacency
	 *    (so faces sharing an edge stay together).
	 * 4) Locating the component which contains a face sharing an edge with the
	 *    last‐emitted face, and rotating the component‐list so it comes first.
	 * 5) Rotating that first component internally so its very first face
	 *    actually shares an edge with the last‐emitted face.
	 * 6) Finally, flattening all components in order—each component still in CW order.
	 *
	 * This guarantees:
	 *   - maximal edge‐continuity (no artificial “gaps” inside a component),
	 *   - strictly CW progression (no backwards jumps),
	 *   - and a single global angle sort per ring for efficiency.
	 */

	// ─── immutable data ──────────────────────────────────────────────────────────
	private final List<PShape> faces; // idx → PShape
	private final int F; // number of faces
	private final int[][] faceVerts; // face → its vertex‐IDs
	private final List<Integer>[] vertToFaces; // vertex → touching faces
	private final int[][] edgeNbrs; // face → edge‐adjacent face‐IDs
	private final double[] angle; // face → global CW‐polar angle

	// ─── mutable iteration state ────────────────────────────────────────────────
	private final boolean[] visited; // face visited?
	private List<Integer> currentRing; // face‐IDs in current ring
	private Iterator<Integer> ringIter;
	private int lastEmitted; // last face‐ID we handed out

	// ───────────────────────────────────────────────────────────────────────────────
	/**
	 * Compute the spiral ordering of {@code allFaces}, beginning at
	 * {@code startFace}. The returned list contains every face exactly once: first
	 * the start face, then all faces 1 hop away (in a continuous edge‐first walk,
	 * with angle‐based tie‐breaking), then all faces 2 hops away, etc.
	 *
	 * @param startFace the face from which to begin the spiral; must be an element
	 *                  of {@code allFaces}
	 * @param allFaces  the complete list of mesh faces (triangles or polygons)
	 * @return a new List of faces in spiral order
	 */
	public static List<PShape> spiral(PShape startFace, List<PShape> allFaces) {
		SpiralIterator it = new SpiralIterator(startFace, allFaces);
		List<PShape> out = new ArrayList<>(allFaces.size());
		while (it.hasNext()) {
			out.add(it.next());
		}
		return out;
	}

	/**
	 * Construct a new spiral iterator over the mesh faces. Performs all
	 * preprocessing (indexing vertices, building adjacency, computing angles) so
	 * that subsequent {@link #hasNext()} and {@link #next()} calls run in amortized
	 * O(1) per face plus one sort/reorder per BFS ring.
	 *
	 * @param startFace one of the faces in {@code allFaces}; this will be returned
	 *                  first by the iterator
	 * @param allFaces  the list of all mesh faces to visit
	 * @throws IllegalArgumentException if {@code startFace} is not in
	 *                                  {@code allFaces}
	 */
	@SuppressWarnings("unchecked")
	public SpiralIterator(PShape startFace, List<PShape> allFaces) {
		this.faces = allFaces;
		this.F = allFaces.size();
		Map<PShape, Integer> fidx = new HashMap<>(F);
		for (int i = 0; i < F; i++)
			fidx.put(allFaces.get(i), i);
		Integer startIdx = fidx.get(startFace);
		if (startIdx == null)
			throw new IllegalArgumentException("startFace not in allFaces");

		// 1) build per‐face vertex‐lists + global vertex‐ID map
		Map<VertexKey, Integer> vmap = new HashMap<>();
		int nextVid = 0;
		faceVerts = new int[F][];
		for (int f = 0; f < F; f++) {
			PShape shp = allFaces.get(f);
			int vc = shp.getVertexCount();
			int[] vids = new int[vc];
			for (int j = 0; j < vc; j++) {
				PVector v = shp.getVertex(j);
				VertexKey k = new VertexKey(v);
				Integer vid = vmap.get(k);
				if (vid == null) {
					vid = nextVid++;
					vmap.put(k, vid);
				}
				vids[j] = vid;
			}
			faceVerts[f] = vids;
		}

		// 2) build vertex→faces incidence
		vertToFaces = new List[nextVid];
		for (int v = 0; v < nextVid; v++) {
			vertToFaces[v] = new ArrayList<>();
		}
		for (int f = 0; f < F; f++) {
			for (int v : faceVerts[f]) {
				vertToFaces[v].add(f);
			}
		}

		// 3) build face→edge‐neighbors via an undirected‐edge map
		Map<EdgeKey, List<Integer>> e2f = new HashMap<>();
		for (int f = 0; f < F; f++) {
			int[] vids = faceVerts[f];
			int vc = vids.length;
			for (int i = 0; i < vc; i++) {
				EdgeKey ek = new EdgeKey(vids[i], vids[(i + 1) % vc]);
				e2f.computeIfAbsent(ek, __ -> new ArrayList<>()).add(f);
			}
		}
		edgeNbrs = new int[F][];
		for (int f = 0; f < F; f++) {
			Set<Integer> nb = new HashSet<>();
			int[] vids = faceVerts[f];
			for (int i = 0; i < vids.length; i++) {
				EdgeKey ek = new EdgeKey(vids[i], vids[(i + 1) % vids.length]);
				for (int g : e2f.get(ek)) {
					if (g != f)
						nb.add(g);
				}
			}
			edgeNbrs[f] = nb.stream().mapToInt(x -> x).toArray();
		}

		// 4) precompute global CW‐angle around startFace
		angle = new double[F];
		PVector c0 = computeCentroid(startFace);
		for (int f = 0; f < F; f++) {
			PVector cf = computeCentroid(allFaces.get(f));
			angle[f] = FastMath.atan2(cf.y - c0.y, cf.x - c0.x);
		}

		// 5) init BFS rings
		visited = new boolean[F];
		visited[startIdx] = true;
		lastEmitted = startIdx;
		currentRing = List.of(startIdx);
		ringIter = currentRing.iterator();
	}

	@Override
	public boolean hasNext() {
		if (ringIter.hasNext()) {
			return true;
		}
		// build the next vertex‐adjacent ring
		Set<Integer> Rset = new LinkedHashSet<>();
		for (int f : currentRing) {
			for (int v : faceVerts[f]) {
				for (int g : vertToFaces[v]) {
					if (!visited[g]) {
						visited[g] = true;
						Rset.add(g);
					}
				}
			}
		}
		if (Rset.isEmpty()) {
			return false;
		}

		// reorder so that edge‐connected components remain together,
		// the component touching lastEmitted comes first, and each
		// component is globally CW‐sorted & rotated to hug the edge.
		currentRing = reorderRing(Rset, lastEmitted);
		ringIter = currentRing.iterator();
		return ringIter.hasNext();
	}

	@Override
	public PShape next() {
		if (!hasNext())
			throw new NoSuchElementException();
		int f = ringIter.next();
		lastEmitted = f;
		return faces.get(f);
	}

	// ─────────────────────────────────────────────────────────────────────────────
	/**
	 * Given the new BFS‐ring Rset and the face‐ID lastEmitted, partition Rset into
	 * edge‐connected components, sort each component CW by global angle, rotate the
	 * “seed” component to start at its member sharing an edge with lastEmitted,
	 * then concatenate.
	 */
	private List<Integer> reorderRing(Set<Integer> Rset, int seedFace) {
		// 1) globally CW‐sort all ring‐members by angle[] descending
		List<Integer> ringSorted = new ArrayList<>(Rset);
		ringSorted.sort((a, b) -> Double.compare(angle[b], angle[a]));

		// 2) extract edge‐connected components in ringSorted order
		Set<Integer> seen = new HashSet<>();
		List<List<Integer>> comps = new ArrayList<>();
		for (int f : ringSorted)
			if (!seen.contains(f)) {
				// flood‐fill
				List<Integer> comp = new ArrayList<>();
				Deque<Integer> stack = new ArrayDeque<>();
				stack.push(f);
				seen.add(f);
				while (!stack.isEmpty()) {
					int u = stack.pop();
					comp.add(u);
					for (int g : edgeNbrs[u]) {
						if (Rset.contains(g) && !seen.contains(g)) {
							seen.add(g);
							stack.push(g);
						}
					}
				}
				comps.add(comp);
			}

		// 3) find which component touches seedFace by an edge
		Set<Integer> seedNbrs = new HashSet<>();
		for (int g : edgeNbrs[seedFace]) {
			if (Rset.contains(g))
				seedNbrs.add(g);
		}
		int firstIdx = 0;
		for (int i = 0; i < comps.size(); i++) {
			for (int f : comps.get(i)) {
				if (seedNbrs.contains(f)) {
					firstIdx = i;
					break;
				}
			}
		}
		// rotate so that comps[firstIdx] is first
		Collections.rotate(comps, -firstIdx);

		// 4) for each component, CW‐sort by global angle
		for (List<Integer> comp : comps) {
			comp.sort((a, b) -> Double.compare(angle[b], angle[a]));
		}

		// 5) build the output list
		List<Integer> out = new ArrayList<>(Rset.size());
		for (int ci = 0; ci < comps.size(); ci++) {
			List<Integer> comp = comps.get(ci);

			// for the first (seeded) comp, rotate to start at the best seed‐nbr
			if (ci == 0 && !seedNbrs.isEmpty()) {
				// pick the comp‐member connected to seedFace
				int best = comp.get(0);
				double bestD = Double.POSITIVE_INFINITY;
				for (int f : comp) {
					if (seedNbrs.contains(f)) {
						double d = cwDist(angle[seedFace], angle[f]);
						if (d < bestD) {
							bestD = d;
							best = f;
						}
					}
				}
				// rotate comp so best is at index 0
				int k = comp.indexOf(best);
				Collections.rotate(comp, -k);
			}

			out.addAll(comp);
		}
		return out;
	}

	/** CW‐distance from angle a1 to a2, in [0,2π). */
	private double cwDist(double a1, double a2) {
		double d = a1 - a2;
		if (d < 0)
			d += Math.PI * 2;
		return d;
	}

	/** centroid of a face */
	private PVector computeCentroid(PShape f) {
		float cx = 0, cy = 0;
		int vc = f.getVertexCount();
		for (int i = 0; i < vc; i++) {
			PVector v = f.getVertex(i);
			cx += v.x;
			cy += v.y;
		}
		return new PVector(cx / vc, cy / vc);
	}

	// ── small helpers for hashing PVector/Edges ─────────────────────────────────
	private static class VertexKey {
		final int xh, yh, zh;

		VertexKey(PVector v) {
			xh = Float.floatToIntBits(v.x);
			yh = Float.floatToIntBits(v.y);
			zh = Float.floatToIntBits(v.z);
		}

		@Override
		public int hashCode() {
			return Objects.hash(xh, yh, zh);
		}

		@Override
		public boolean equals(Object o) {
			if (!(o instanceof VertexKey))
				return false;
			VertexKey k = (VertexKey) o;
			return xh == k.xh && yh == k.yh && zh == k.zh;
		}
	}

	private static class EdgeKey {
		final int a, b;

		EdgeKey(int x, int y) {
			if (x < y) {
				a = x;
				b = y;
			} else {
				a = y;
				b = x;
			}
		}

		@Override
		public int hashCode() {
			return Objects.hash(a, b);
		}

		@Override
		public boolean equals(Object o) {
			if (!(o instanceof EdgeKey))
				return false;
			EdgeKey e = (EdgeKey) o;
			return a == e.a && b == e.b;
		}
	}
}