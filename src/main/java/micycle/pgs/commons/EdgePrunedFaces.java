package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import org.jgrapht.alg.interfaces.VertexColoringAlgorithm.Coloring;
import org.jgrapht.alg.spanning.GreedyMultiplicativeSpanner;
import org.jgrapht.graph.AbstractBaseGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.utils.TriangleCollector;
import org.tinspin.index.PointMap;
import org.tinspin.index.kdtree.KDTree;

import micycle.pgs.PGS_Conversion;
import micycle.pgs.PGS_Triangulation;
import processing.core.PShape;

/**
 * Utilities for extracting polygonal faces from a TIN by pruning TIN edges
 * according to a rule, grouping triangles across the pruned edges, and tracing
 * group boundaries.
 *
 * <p>
 * Core idea
 * </p>
 * <ul>
 * <li>Start with a Delaunay TIN ({@link IIncrementalTin}).</li>
 * <li>Drop a subset of base edges according to a rule (Urquhart, Gabriel, RNG,
 * k-spanner, etc.).</li>
 * <li>Merge triangles across any dropped edge (DSU / union-find).</li>
 * <li>For each merged component, collect only its boundary half-edges (those
 * not dropped and not shared by two triangles of the component), then sequence
 * them to form a ring and emit a polygon.</li>
 * </ul>
 *
 * <p>
 * Benefits
 * </p>
 * <ul>
 * <li>Linear-time over triangles/edges for the TIN portion; avoids
 * polygonization / noding / unions.</li>
 * <li>Robust topology, since edges are taken directly from the TIN; boundaries
 * are sequenced only.</li>
 * <li>Pluggable rules via a simple {@code DropRule} interface.</li>
 * <li>Consistent perimeter handling: when {@code preservePerimeter} is true,
 * constraint/hull borders are never dropped.</li>
 * </ul>
 *
 * <p>
 * Implemented rules
 * </p>
 * <ul>
 * <li>Urquhart: drop the longest base edge of each triangle.</li>
 * <li>Gabriel: drop uv if the midpoint of uv has a nearest vertex that is
 * neither u nor v.</li>
 * <li>Relative Neighborhood (RNG): drop uv if there exists w with max(d(u,w),
 * d(v,w)) &lt; d(u,v).</li>
 * <li>k-Spanner: keep only edges chosen by a greedy multiplicative spanner on
 * the TIN graph (SimpleGraph&lt;Vertex, IQuadEdge&gt;), drop the rest.</li>
 * <li>Edge-collapse quadrangulation: 3-color the three edges of each triangle
 * (via graph coloring) and drop all edges with color &ge; 2 (or equivalently,
 * keep colors 0–1) to form mostly quadrilateral faces; honors
 * perimeter/constraint borders when requested.</li>
 * </ul>
 * 
 * @author Michael Carleton
 */
public class EdgePrunedFaces {

	// thin wrappers around the common pipelines

	public static PShape urquhartFaces(final IIncrementalTin tin, final boolean preservePerimeter) {
		return facesSequencedCommon(tin, preservePerimeter, URQUHART_RULE);
	}

	public static PShape gabrielFaces(final IIncrementalTin tin, final boolean preservePerimeter) {
		return facesSequencedCommon(tin, preservePerimeter, GABRIEL_RULE);
	}

	public static PShape relativeNeighborFaces(final IIncrementalTin tin, final boolean preservePerimeter) {
		return facesSequencedCommon(tin, preservePerimeter, RELATIVE_NEIGHBOR_RULE);
	}

	public static PShape spannerFaces(final IIncrementalTin tin, int k, final boolean preservePerimeter) {
		return facesSequencedCommon(tin, preservePerimeter, spannerDropRule(tin, k));
	}

	public static PShape edgeCollapseQuadrangulation(final IIncrementalTin tin, final boolean preservePerimeter) {
		return facesSequencedCommon(tin, preservePerimeter, edgeCollapseQuadrangulationDropRule(tin));
	}

	/**
	 * Common pipeline: collect mesh, mark dropped edges via rule, group via DSU,
	 * sequence boundaries to polygons.
	 * 
	 * @param tin
	 * @param preservePerimeter
	 * @param rule
	 * @return
	 */
	private static PShape facesSequencedCommon(final IIncrementalTin tin, final boolean preservePerimeter, final DropRule rule) {
		final GeometryFactory gf = new GeometryFactory();
		final boolean notConstrained = tin.getConstraints().isEmpty();

		// Collect triangles, adjacency, and unique vertices
		final List<TriRec> tris = new ArrayList<>();
		final Map<IQuadEdge, int[]> edgeAdj = new IdentityHashMap<>();
		final Set<Vertex> vset = new HashSet<>();

		TriangleCollector.visitSimpleTriangles(tin, t -> {
			final IConstraint c = t.getContainingRegion();
			if (!(notConstrained || (c != null && c.definesConstrainedRegion()))) {
				return;
			}

			final int tid = tris.size();
			final IQuadEdge ea = t.getEdgeA();
			final IQuadEdge eb = t.getEdgeB();
			final IQuadEdge ec = t.getEdgeC();
			tris.add(new TriRec(ea, eb, ec));

			addAdj(edgeAdj, ea.getBaseReference(), tid);
			addAdj(edgeAdj, eb.getBaseReference(), tid);
			addAdj(edgeAdj, ec.getBaseReference(), tid);

			vset.add(t.getVertexA());
			vset.add(t.getVertexB());
			vset.add(t.getVertexC());
		});

		final int nTri = tris.size();
		if (nTri == 0) {
			return new PShape();
		}

		final MeshCtx ctx = new MeshCtx(tris, edgeAdj, new ArrayList<>(vset));

		// Mark dropped edges according to the rule
		final Set<IQuadEdge> dropped = Collections.newSetFromMap(new IdentityHashMap<>());
		rule.markDropped(ctx, preservePerimeter, dropped);

		// DSU across any dropped base edge
		final DSU dsu = new DSU(nTri);
		for (IQuadEdge e : dropped) {
			final int[] inc = edgeAdj.get(e);
			if (inc != null && inc[0] >= 0 && inc[1] >= 0) {
				dsu.union(inc[0], inc[1]);
			}
		}

		// Group triangles
		final Map<Integer, List<Integer>> comp = new HashMap<>();
		for (int t = 0; t < nTri; t++) {
			comp.computeIfAbsent(dsu.find(t), k -> new ArrayList<>()).add(t);
		}

		// Build faces: collect boundary half-edges, sequence into one ring, emit
		// polygon
		var faces = comp.values().parallelStream().map(triIds -> {
			final List<IQuadEdge> boundary = new ArrayList<>();

			for (int tid : triIds) {
				final TriRec tr = tris.get(tid);
				for (int i = 0; i < 3; i++) {
					final IQuadEdge he = tr.e[i];
					final IQuadEdge base = he.getBaseReference();

					// Skip dropped inner seams
					if (dropped.contains(base)) {
						continue;
					}

					// If neighbor is in same group, edge is interior, skip
					final int[] inc = edgeAdj.get(base);
					final int nb = (inc == null) ? -1 : (inc[0] == tid ? inc[1] : inc[0]);
					if (nb >= 0 && dsu.find(nb) == dsu.find(tid)) {
						continue;
					}

					// Boundary half-edge, interior on the left
					boundary.add(he);
				}
			}

			if (boundary.isEmpty()) {
				return null;
			}

			final Coordinate[] ring = sequenceSingleLoop(boundary);
			if (ring == null) {
				return null;
			}

			return PGS_Conversion.toPShape(gf.createPolygon(ring));
		});

		return PGS_Conversion.flatten(faces.filter(Objects::nonNull).toList());
	}

	private static Coordinate[] sequenceSingleLoop(List<IQuadEdge> edges) {
		if (edges.isEmpty()) {
			return null;
		}

		final Map<Vertex, IQuadEdge> out = new IdentityHashMap<>(edges.size() * 2);
		for (IQuadEdge e : edges) {
			out.put(e.getA(), e); // assumes at most one outgoing per vertex on the boundary
		}

		final IQuadEdge start = edges.get(0);
		final Vertex startV = start.getA();

		final List<Coordinate> coords = new ArrayList<>(edges.size() + 1);
		coords.add(new Coordinate(start.getA().x, start.getA().y));

		IQuadEdge cur = start;
		for (int i = 0; i < edges.size(); i++) {
			coords.add(new Coordinate(cur.getB().x, cur.getB().y));
			if (cur.getB() == startV) {
				break; // closed the loop
			}
			final IQuadEdge next = out.get(cur.getB());
			if (next == null) {
				return null; // unexpected: dangling
			}
			cur = next;
		}

		// Ensure closed
		if (!coords.get(0).equals2D(coords.get(coords.size() - 1))) {
			coords.add(new Coordinate(coords.get(0)));
		}
		if (coords.size() < 4) {
			return null;
		}

		return coords.toArray(new Coordinate[0]);
	}

	private static final DropRule URQUHART_RULE = (ctx, preservePerimeter, dropped) -> {
		// Urquhart rule: drop the longest base edge of each triangle
		for (TriRec tr : ctx.tris) {
			IQuadEdge e0 = tr.e[0].getBaseReference();
			IQuadEdge e1 = tr.e[1].getBaseReference();
			IQuadEdge e2 = tr.e[2].getBaseReference();

			IQuadEdge max = e0;
			double m2 = len2(e0);
			double l1 = len2(e1);
			if (l1 > m2) {
				m2 = l1;
				max = e1;
			}
			double l2 = len2(e2);
			if (l2 > m2) {
				max = e2;
			}

			if (!preservePerimeter || !max.isConstraintRegionBorder()) {
				dropped.add(max);
			}
		}
	};

	private static final DropRule GABRIEL_RULE = (ctx, preservePerimeter, dropped) -> {
		// Gabriel rule: drop edges whose midpoint’s nearest vertex is neither endpoint
		// Build KD-tree of all vertices
		final PointMap<Vertex> tree = KDTree.create(2);
		for (Vertex v : ctx.vertices) {
			tree.insert(new double[] { v.x, v.y }, v);
		}

		for (IQuadEdge base : ctx.edgeAdj.keySet()) {
			if (preservePerimeter && base.isConstraintRegionBorder()) {
				continue;
			}

			final double mx = 0.5 * (base.getA().x + base.getB().x);
			final double my = 0.5 * (base.getA().y + base.getB().y);
			final Vertex nn = tree.query1nn(new double[] { mx, my }).value();

			if (nn != base.getA() && nn != base.getB()) {
				dropped.add(base);
			}
		}
	};

	private static final DropRule RELATIVE_NEIGHBOR_RULE = (ctx, preservePerimeter, dropped) -> {
		// Relative Neighborhood Graph: drop uv if exists w in N(u) or N(v) with
		// max(d(u,w), d(v,w)) < d(u,v)
		// Build 1-ring vertex adjacency from base edges
		final IdentityHashMap<Vertex, Set<Vertex>> nbrs = new IdentityHashMap<>();
		for (IQuadEdge base : ctx.edgeAdj.keySet()) {
			final Vertex a = base.getA();
			final Vertex b = base.getB();
			nbrs.computeIfAbsent(a, k -> Collections.newSetFromMap(new IdentityHashMap<>())).add(b);
			nbrs.computeIfAbsent(b, k -> Collections.newSetFromMap(new IdentityHashMap<>())).add(a);
		}

		// Test each base edge against RNG condition, using squared distances
		for (IQuadEdge base : ctx.edgeAdj.keySet()) {
			if (preservePerimeter && base.isConstraintRegionBorder()) {
				continue;
			}

			final Vertex a = base.getA();
			final Vertex b = base.getB();
			final double l2 = dist2(a, b);

			boolean drop = false;

			// Check neighbors of a
			final Set<Vertex> Na = nbrs.getOrDefault(a, Collections.emptySet());
			for (Vertex w : Na) {
				if (w == b) {
					continue;
				}
				if (Math.max(dist2(w, a), dist2(w, b)) < l2) {
					drop = true;
					break;
				}
			}
			// If not dropped, check neighbors of b
			if (!drop) {
				final Set<Vertex> Nb = nbrs.getOrDefault(b, Collections.emptySet());
				for (Vertex w : Nb) {
					if (w == a) {
						continue;
					}
					if (Math.max(dist2(w, a), dist2(w, b)) < l2) {
						drop = true;
						break;
					}
				}
			}

			if (drop) {
				dropped.add(base);
			}
		}
	};

	private static DropRule spannerDropRule(final IIncrementalTin tin, final int kParam) {
		final int k = Math.max(2, kParam);
		return (ctx, preservePerimeter, dropped) -> {
			final SimpleGraph<Vertex, IQuadEdge> g = PGS_Triangulation.toTinfourGraph(tin);
			if (g.edgeSet().isEmpty()) {
				return;
			}

			final GreedyMultiplicativeSpanner<Vertex, IQuadEdge> sp = new GreedyMultiplicativeSpanner<>(g, k);

			// Build identity set of kept base refs from the spanner result
			final Set<IQuadEdge> kept = Collections.newSetFromMap(new IdentityHashMap<>());
			for (IQuadEdge e : sp.getSpanner()) {
				kept.add(e.getBaseReference());
			}

			// Drop every base edge not in the spanner (unless perimeter preserved)
			for (IQuadEdge base : ctx.edgeAdj.keySet()) {
				if (preservePerimeter && base.isConstraintRegionBorder()) {
					continue;
				}
				if (!kept.contains(base)) {
					dropped.add(base);
				}
			}
		};
	}

	private static DropRule edgeCollapseQuadrangulationDropRule(final IIncrementalTin tin) {
		/*-
		 * From 'Fast unstructured quadrilateral mesh generation'.
		 * A better coloring approach is given in 'Face coloring in unstructured CFD codes'.
		 * 
		 * First partition the edges of the triangular mesh into three groups such that
		 * no triangle has two edges of the same color (find groups by reducing to a
		 * graph-coloring).
		 * Then obtain an all-quadrilateral mesh by removing all edges of *one* 
		 * particular color.
		 */
		final boolean unconstrained = tin.getConstraints().isEmpty();

		// Collect unconstrained perimeter edges (base refs) if applicable
		final Set<IQuadEdge> perimeterBaseRefs = Collections.newSetFromMap(new IdentityHashMap<>());
		if (unconstrained) {
			for (IQuadEdge e : tin.getPerimeter()) {
				perimeterBaseRefs.add(e.getBaseReference());
			}
		}

		return (ctx, preservePerimeter, dropped) -> {
			// Build a graph where each vertex is a base edge; connect the 3 edges of each
			// triangle
			final AbstractBaseGraph<IQuadEdge, DefaultEdge> g = new SimpleGraph<>(DefaultEdge.class);
			for (TriRec tr : ctx.tris) {
				final IQuadEdge a = tr.e[0].getBaseReference();
				final IQuadEdge b = tr.e[1].getBaseReference();
				final IQuadEdge c = tr.e[2].getBaseReference();

				g.addVertex(a);
				g.addVertex(b);
				g.addVertex(c);

				g.addEdge(a, b);
				g.addEdge(a, c);
				g.addEdge(b, c);
			}
			if (g.vertexSet().isEmpty()) {
				return;
			}

			// 3-color the "edge graph" so no triangle has two edges of the same color
			final Coloring<IQuadEdge> coloring = new RLFColoring<>(g, 1337).getColoring();

			// Mark all edges of the chosen color as dropped, honoring perimeter
			// preservation
			for (Map.Entry<IQuadEdge, Integer> e : coloring.getColors().entrySet()) {
				final IQuadEdge base = e.getKey();
				final int color = e.getValue();

				/*
				 * NOTE 4-colorings are possible, so some triangles may have two or three edges
				 * with color >= 2, and yield faces larger than quads once edges are dropped.
				 */
				if (color < 2) {
					continue;
				}

				if (preservePerimeter) {
					// Preserve constraint borders, and unconstrained outer perimeter
					if (base.isConstraintRegionBorder() || perimeterBaseRefs.contains(base)) {
						continue;
					}
				}

				dropped.add(base);
			}
		};
	}

	private static void addAdj(Map<IQuadEdge, int[]> adj, IQuadEdge base, int tid) {
		int[] a = adj.get(base);
		if (a == null) {
			a = new int[] { -1, -1 };
			adj.put(base, a);
		}
		if (a[0] < 0) {
			a[0] = tid;
		} else {
			a[1] = tid;
		}
	}

	private static double dist2(Vertex p, Vertex q) {
		return p.getDistanceSq(q);
	}

	private static double len2(IQuadEdge e) {
		final double dx = e.getA().x - e.getB().x;
		final double dy = e.getA().y - e.getB().y;
		return dx * dx + dy * dy;
	}

	/**
	 * Strategy for selecting TIN base edges to drop prior to face extraction.
	 * <p>
	 * The pipeline will union triangles across every dropped base edge, then trace
	 * the remaining edges to form polygon boundaries. A rule inspects the immutable
	 * mesh context and adds base edges (getBaseReference()) to the provided set.
	 * </p>
	 * <p>
	 * Typical usage idea:
	 * </p>
	 * 
	 * <pre>
	 * DropRule rule = (ctx, preservePerimeter, dropped) -> {
	 * 	for (IQuadEdge base : ctx.edgeAdj.keySet()) {
	 * 		if (preservePerimeter && isPerimeter(ctx, base))
	 * 			continue; // edge has &lt;2 incident triangles
	 * 		if (shouldDropAccordingToRule(base, ctx))
	 * 			dropped.add(base);
	 * 	}
	 * };
	 * </pre>
	 * <p>
	 * Notes:
	 * </p>
	 * <ul>
	 * <li>Add base edges only (not half-edges); identity semantics apply.</li>
	 * <li>If preservePerimeter is true, do not drop perimeter/constraint
	 * borders.</li>
	 * <li>Do not mutate ctx; only populate the dropped set.</li>
	 * <li>Deterministic selection is recommended.</li>
	 * </ul>
	 */
	private interface DropRule {
		/**
		 * Marks base edges to be removed.
		 *
		 * @param ctx               immutable mesh data (triangles, base-edge adjacency,
		 *                          vertices)
		 * @param preservePerimeter when true, perimeter/constraint borders must be kept
		 * @param dropped           identity set to populate with base edges to drop
		 */
		void markDropped(MeshCtx ctx, boolean preservePerimeter, Set<IQuadEdge> dropped);
	}

	// Collected mesh context from the TIN
	private static final class MeshCtx {
		final List<TriRec> tris;
		final Map<IQuadEdge, int[]> edgeAdj; // base edge -> up to 2 incident triangle ids
		final List<Vertex> vertices; // unique vertices seen in accepted triangles

		MeshCtx(List<TriRec> tris, Map<IQuadEdge, int[]> edgeAdj, List<Vertex> vertices) {
			this.tris = tris;
			this.edgeAdj = edgeAdj;
			this.vertices = vertices;
		}
	}

	// Triangle record with oriented half-edges (triangle on left)
	private static record TriRec(IQuadEdge[] e) {
		public TriRec {
			if (e == null)
				throw new NullPointerException("e");
			if (e.length != 3)
				throw new IllegalArgumentException("array must be length 3");
			e = e.clone(); // defensive copy before assignment
		}

		public TriRec(IQuadEdge a, IQuadEdge b, IQuadEdge c) {
			this(new IQuadEdge[] { a, b, c });
		}

		// override accessor to return a copy so callers can't mutate internal array
		@Override
		public IQuadEdge[] e() {
			return e.clone();
		}
	}

	// Simple DSU
	private static final class DSU {
		final int[] p, r;

		DSU(int n) {
			p = new int[n];
			r = new int[n];
			for (int i = 0; i < n; i++) {
				p[i] = i;
			}
		}

		int find(int x) {
			return p[x] == x ? x : (p[x] = find(p[x]));
		}

		void union(int a, int b) {
			a = find(a);
			b = find(b);
			if (a == b) {
				return;
			}
			if (r[a] < r[b]) {
				int t = a;
				a = b;
				b = t;
			}
			p[b] = a;
			if (r[a] == r[b]) {
				r[a]++;
			}
		}
	}

}
