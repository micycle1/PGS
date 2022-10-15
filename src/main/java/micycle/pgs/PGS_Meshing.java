package micycle.pgs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jgrapht.alg.interfaces.VertexColoringAlgorithm.Coloring;
import org.jgrapht.alg.spanning.GreedyMultiplicativeSpanner;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.graph.AbstractBaseGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.strtree.STRtree;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.utils.TriangleCollector;
import org.tinspin.index.PointIndex;
import org.tinspin.index.kdtree.KDTree;

import micycle.pgs.color.RGB;
import micycle.pgs.commons.ClipperPGS;
import micycle.pgs.commons.IncrementalTinDual;
import micycle.pgs.commons.PEdge;
import micycle.pgs.commons.RLFColoring;
import micycle.pgs.commons.SpiralQuadrangulation;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Mesh generation (excluding triangulation).
 * <p>
 * Many of the methods within this class process an existing Delaunay
 * triangulation; you may first generate such a triangulation from a shape using
 * the
 * {@link PGS_Triangulation#delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
 * delaunayTriangulationMesh()} method
 * 
 * @author Michael Carleton
 * @since 1.2.0
 */
public class PGS_Meshing {

	private PGS_Meshing() {
	}

	/**
	 * Generates a shape consisting of polygonal faces of an <i>Urquhart graph</i>.
	 * An Urquhart graph is obtained by removing the longest edge from each triangle
	 * in a triangulation.
	 * <p>
	 * In practice this is a way to tessellate a shape into polygons (with the
	 * resulting tessellation being in between a
	 * {@link PGS_Triangulation#delaunayTriangulation(PShape) triangulation} and a
	 * {@link micycle.pgs.PGS_Processing#convexPartition(PShape) partition}).
	 * <p>
	 * Note that this method processes a Delaunay triangulation. Process a shape
	 * using
	 * {@link PGS_Triangulation#delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * delaunayTriangulationMesh()} first and then feed it to this method.
	 * <p>
	 * The <i>Urquhart graph</i> is a good approximation to the
	 * {@link #relativeNeighborFaces(IIncrementalTin, boolean) <i>relative
	 * neighborhood</i>} graph (having only about 2% additional edges).
	 * 
	 * @param triangulation     a triangulation mesh
	 * @param preservePerimeter whether to retain/preserve edges on the perimeter
	 *                          even if they should be removed according to the
	 *                          urquhart condition
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.1.0
	 * @see #gabrielFaces(IIncrementalTin, boolean)
	 */
	public static PShape urquhartFaces(final IIncrementalTin triangulation, final boolean preservePerimeter) {
		final HashSet<IQuadEdge> edges = PGS.makeHashSet(triangulation.getMaximumEdgeAllocationIndex());
		final HashSet<IQuadEdge> uniqueLongestEdges = PGS.makeHashSet(triangulation.getMaximumEdgeAllocationIndex());

		final boolean notConstrained = triangulation.getConstraints().isEmpty();

		TriangleCollector.visitSimpleTriangles(triangulation, t -> {
			final IConstraint constraint = t.getContainingRegion();
			if (notConstrained || (constraint != null && constraint.definesConstrainedRegion())) {
				edges.add(t.getEdgeA().getBaseReference());
				edges.add(t.getEdgeB().getBaseReference());
				edges.add(t.getEdgeC().getBaseReference());
				final IQuadEdge longestEdge = findLongestEdge(t).getBaseReference();
				if (!preservePerimeter || (preservePerimeter && !longestEdge.isConstrainedRegionBorder())) {
					uniqueLongestEdges.add(longestEdge);
				}
			}
		});

		edges.removeAll(uniqueLongestEdges);

		final Collection<PEdge> meshEdges = new ArrayList<>(edges.size());
		edges.forEach(edge -> meshEdges.add(new PEdge(edge.getA().x, edge.getA().y, edge.getB().x, edge.getB().y)));

		return PGS.polygonizeEdges(meshEdges);
	}

	/**
	 * Generates a shape consisting of polygonal faces of a <i>Gabriel graph</i>. A
	 * Gabriel graph is obtained by removing each edge E from a triangulation if a
	 * vertex lies within a circle of diameter = length(E), centered on the midpoint
	 * of E.
	 * <p>
	 * In practice this is a way to tessellate a shape into polygons (with the
	 * resulting tessellation being reminiscent of shattering the shape as if it
	 * were glass).
	 * <p>
	 * Note that this method processes a Delaunay triangulation. Process a shape
	 * using
	 * {@link PGS_Triangulation#delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * delaunayTriangulationMesh()} first and then feed it to this method.
	 * 
	 * @param triangulation     a triangulation mesh
	 * @param preservePerimeter whether to retain/preserve edges on the perimeter
	 *                          even if they should be removed according to the
	 *                          gabriel condition
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.1.0
	 * @see #urquhartFaces(IIncrementalTin, boolean)
	 */
	public static PShape gabrielFaces(final IIncrementalTin triangulation, final boolean preservePerimeter) {
		final HashSet<IQuadEdge> edges = new HashSet<>();
		final HashSet<Vertex> vertices = new HashSet<>();

		final boolean notConstrained = triangulation.getConstraints().isEmpty();
		TriangleCollector.visitSimpleTriangles(triangulation, t -> {
			final IConstraint constraint = t.getContainingRegion();
			if (notConstrained || (constraint != null && constraint.definesConstrainedRegion())) {
				edges.add(t.getEdgeA().getBaseReference()); // add edge to set
				edges.add(t.getEdgeB().getBaseReference()); // add edge to set
				edges.add(t.getEdgeC().getBaseReference()); // add edge to set
				vertices.add(t.getVertexA());
				vertices.add(t.getVertexB());
				vertices.add(t.getVertexC());
			}
		});

		final PointIndex<Vertex> tree = KDTree.create(2, (p1, p2) -> {
			final double deltaX = p1[0] - p2[0];
			final double deltaY = p1[1] - p2[1];
			return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
		});
		vertices.forEach(v -> tree.insert(new double[] { v.x, v.y }, v));

		final HashSet<IQuadEdge> nonGabrielEdges = new HashSet<>(); // base references to edges that should be removed
		edges.forEach(edge -> {
			final double[] midpoint = midpoint(edge);
			final Vertex near = tree.query1NN(midpoint).value();
			if (near != edge.getA() && near != edge.getB()) {
				if (!preservePerimeter || (preservePerimeter && !edge.isConstrainedRegionBorder())) {
					nonGabrielEdges.add(edge); // base reference
				}
			}
		});
		edges.removeAll(nonGabrielEdges);

		final Collection<PEdge> meshEdges = new ArrayList<>(edges.size());
		edges.forEach(edge -> meshEdges.add(new PEdge(edge.getA().x, edge.getA().y, edge.getB().x, edge.getB().y)));

		return PGS.polygonizeEdges(meshEdges);
	}

	/**
	 * Generates a shape consisting of polygonal faces of a <i>Relative neighborhood
	 * graph</i> (RNG).
	 * <p>
	 * An RNG is obtained by removing each edge E from a triangulation if any vertex
	 * is nearer to both vertices of E than the length of E.
	 * <p>
	 * The RNG is a subgraph of the {@link #urquhartFaces(IIncrementalTin, boolean)
	 * urquhart} graph, having only slightly fewer edges.
	 * 
	 * @param triangulation     a triangulation mesh
	 * @param preservePerimeter whether to retain/preserve edges on the perimeter
	 *                          even if they should be removed according to the
	 *                          relative neighbor condition
	 * @return
	 * @since 1.3.0
	 */
	public static PShape relativeNeighborFaces(final IIncrementalTin triangulation, final boolean preservePerimeter) {
		SimpleGraph<Vertex, IQuadEdge> graph = PGS_Triangulation.toTinfourGraph(triangulation);
		NeighborCache<Vertex, IQuadEdge> cache = new NeighborCache<>(graph);

		Set<IQuadEdge> edges = new HashSet<>(graph.edgeSet());

		/*
		 * If any vertex is nearer to both vertices of an edge, than the length of the
		 * edge, this edge does not belong in the RNG.
		 */
		graph.edgeSet().forEach(e -> {
			double l = e.getLength();
			cache.neighborsOf(e.getA()).forEach(n -> {
				if (Math.max(n.getDistance(e.getA()), n.getDistance(e.getB())) < l) {
					if (!preservePerimeter || (preservePerimeter && !e.isConstrainedRegionBorder())) {
						edges.remove(e);
					}
				}
			});
			cache.neighborsOf(e.getB()).forEach(n -> {
				if (Math.max(n.getDistance(e.getA()), n.getDistance(e.getB())) < l) {
					if (!preservePerimeter || (preservePerimeter && !e.isConstrainedRegionBorder())) {
						edges.remove(e);
					}
				}
			});
		});

		List<PEdge> edgesOut = edges.stream().map(PGS_Triangulation::toPEdge).collect(Collectors.toList());

		if (preservePerimeter) {
			return PGS.polygonizeEdgesRobust(edgesOut);
		} else {
			return PGS.polygonizeEdges(edgesOut);
		}

	}

	/**
	 * Generates a shape consisting of polygonal faces formed by edges returned by a
	 * greedy sparse spanner applied to a triangulation.
	 * 
	 * @param triangulation     a triangulation mesh
	 * @param k                 the order of the spanner. Should be at least 1.
	 *                          Higher numbers collapse more edges resulting in
	 *                          larger faces, until a single face remains
	 * @param preservePerimeter whether to retain/preserve edges on the perimeter
	 *                          even if they should be removed according to the
	 *                          spanner condition
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.3.0
	 */
	public static PShape spannerFaces(final IIncrementalTin triangulation, int k, final boolean preservePerimeter) {
		SimpleGraph<PVector, PEdge> graph = PGS_Triangulation.toGraph(triangulation);
		if (graph.edgeSet().isEmpty()) {
			return new PShape();
		}

		k = Math.max(2, k); // min(2) since k=1 returns triangulation
		GreedyMultiplicativeSpanner<PVector, PEdge> spanner = new GreedyMultiplicativeSpanner<>(graph, k);
		List<PEdge> spannerEdges = spanner.getSpanner().stream().collect(Collectors.toList());
		if (preservePerimeter) {
			if (triangulation.getConstraints().isEmpty()) { // does not have constraints
				spannerEdges.addAll(triangulation.getPerimeter().stream().map(PGS_Triangulation::toPEdge).collect(Collectors.toList()));
			} else { // has constraints
				spannerEdges.addAll(triangulation.getEdges().stream().filter(IQuadEdge::isConstrainedRegionBorder)
						.map(PGS_Triangulation::toPEdge).collect(Collectors.toList()));
			}
		}

		return PGS.polygonizeEdgesRobust(spannerEdges);
	}

	/**
	 * Generates a (mesh-like) shape consisting of polygonal faces of the dual graph
	 * of the given triangulation.
	 * <p>
	 * In practice, the resulting dual mesh has hexagonal-like cells.
	 * <p>
	 * Note that this method processes a Delaunay triangulation. Process a shape
	 * using
	 * {@link PGS_Triangulation#delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * delaunayTriangulationMesh()} first and then feed it to this method.
	 * <p>
	 * If the input has been generated from a PShape, consider generating the
	 * triangulation with <code>refinements > 1</code> for better dual mesh results.
	 * 
	 * @param triangulation a triangulation mesh
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.2.0
	 */
	public static PShape dualFaces(final IIncrementalTin triangulation) {
		// TODO SEE HOT: Hodge-Optimized Triangulations - use voronoi dual / HOT
		final IncrementalTinDual dual = new IncrementalTinDual(triangulation);
		final PShape dualMesh = dual.getMesh();
		PGS_Conversion.setAllFillColor(dualMesh, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(dualMesh, RGB.PINK, 2);
		return dualMesh;
	}

	/**
	 * Produces a quadrangulation from a triangulation, by splitting each triangle
	 * into three quadrangles (using the <i>Catmull and Clark</i> technique). A
	 * quadrangulation is a mesh where every face is a quadrangle.
	 * <p>
	 * Since this method employs a very simple technique to produce a
	 * quadrangulation, the result is poor-quality, containing many helix-like
	 * structures (it's not at all "regular").
	 * <p>
	 * Note that this method processes a Delaunay triangulation. Process a shape
	 * using
	 * {@link PGS_Triangulation#delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * delaunayTriangulationMesh()} first and then feed it to this method.
	 * 
	 * @param triangulation a triangulation mesh
	 * @return a GROUP PShape, where each child shape is one quadrangle
	 * @since 1.2.0
	 */
	public static PShape splitQuadrangulation(final IIncrementalTin triangulation) {
		// https://www.cs.mcgill.ca/~cs507/projects/1998/rachelp/welcome.html
		final PShape quads = new PShape(PConstants.GROUP);

		/*-
		 * 1. Insert a Steiner point along the interior of every edge of each 
		 * triangle.
		 * 2. Insert an extra Steiner point in the interior of each triangle.
		 * 3. Connect the Steiner point inside each triangle to the Steiner points on the
		 * edges of that triangle.
		 * Each triangle is converted into three quadrangles.
		 */
		final boolean unconstrained = triangulation.getConstraints().isEmpty();

		TriangleCollector.visitSimpleTriangles(triangulation, t -> {
			final IConstraint constraint = t.getContainingRegion();
			if (unconstrained || (constraint != null && constraint.definesConstrainedRegion())) {
				final PVector p1 = PGS_Triangulation.toPVector(t.getVertexA());
				final PVector p2 = PGS_Triangulation.toPVector(t.getVertexB());
				final PVector p3 = PGS_Triangulation.toPVector(t.getVertexC());
				final PVector sA = PVector.add(p1, p2).div(2); // steiner point p1-p2
				final PVector sB = PVector.add(p2, p3).div(2); // steiner point p2-p3
				final PVector sC = PVector.add(p3, p1).div(2); // steiner point p3-p1

				// compute ?barycenter? of triangle = interior steiner point
				final PVector cSeg = PVector.add(sA, sB).div(2);
				final PVector sI = PVector.add(cSeg, sC).div(2); // interior steiner point

				// anti-clockwise, starting at original vertex
				quads.addChild(PGS_Conversion.fromPVector(p1, sC, sI, sA));
				quads.addChild(PGS_Conversion.fromPVector(p2, sA, sI, sB));
				quads.addChild(PGS_Conversion.fromPVector(p3, sB, sI, sC));
			}
		});

		/*-
		 * Now ideally "regularize" the mesh using techniques explored here:
		 * https://acdl.mit.edu/ESP/Publications/AIAApaper2019-1988.pdf
		 * https://acdl.mit.edu/ESP/Publications/IMR28.pdf
		 */

		PGS_Conversion.setAllFillColor(quads, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(quads, RGB.PINK, 2);

		return quads;
	}

	/**
	 * Generates a quadrangulation from a triangulation by selectively removing (or
	 * "collapsing") the edges shared by neighboring triangles (via edge coloring).
	 * <p>
	 * This method may be slow on large inputs (as measured by vertex count), owing
	 * to the graph coloring it performs.
	 * 
	 * @param triangulation     a triangulation mesh
	 * @param preservePerimeter whether to preserve the perimeter of the input
	 *                          triangulation; when true, retains edges that lie on
	 *                          the perimeter of the triangulation mesh that would
	 *                          have otherwise been removed (this results in some
	 *                          triangles being included in the output).
	 * @return a GROUP PShape, where each child shape is one quadrangle
	 * @since 1.2.0
	 */
	public static PShape edgeCollapseQuadrangulation(final IIncrementalTin triangulation, final boolean preservePerimeter) {
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
		final boolean unconstrained = triangulation.getConstraints().isEmpty();
		final AbstractBaseGraph<IQuadEdge, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
		TriangleCollector.visitSimpleTriangles(triangulation, t -> {
			final IConstraint constraint = t.getContainingRegion();
			if (unconstrained || (constraint != null && constraint.definesConstrainedRegion())) {
				graph.addVertex(t.getEdgeA().getBaseReference());
				graph.addVertex(t.getEdgeB().getBaseReference());
				graph.addVertex(t.getEdgeC().getBaseReference());

				graph.addEdge(t.getEdgeA().getBaseReference(), t.getEdgeB().getBaseReference());
				graph.addEdge(t.getEdgeA().getBaseReference(), t.getEdgeC().getBaseReference());
				graph.addEdge(t.getEdgeB().getBaseReference(), t.getEdgeC().getBaseReference());
			}
		});

		Coloring<IQuadEdge> coloring = new RLFColoring<>(graph).getColoring();

		final HashSet<IQuadEdge> perimeter = new HashSet<>(triangulation.getPerimeter());
		if (!unconstrained) {
			perimeter.clear(); // clear, the perimeter of constrained tin is unaffected by the constraint
		}

		final Collection<PEdge> meshEdges = new ArrayList<>();
		coloring.getColors().forEach((edge, color) -> {
			/*
			 * "We can remove the edges of any one of the colors, however a convenient
			 * choice is the one that leaves the fewest number of unmerged boundary
			 * triangles". -- ideal, but not implemented here...
			 */
			// NOTE could now apply Topological optimization, as given in paper.
			if ((color < 2) || (preservePerimeter && (edge.isConstrainedRegionBorder() || perimeter.contains(edge)))) {
				meshEdges.add(new PEdge(edge.getA().x, edge.getA().y, edge.getB().x, edge.getB().y));
			}
		});

		PShape quads = PGS.polygonizeEdges(meshEdges);
		if (triangulation.getConstraints().size() < 2) { // assume constraint 1 is the boundary (not a hole)
			return quads;
		} else {
			return removeHoles(quads, triangulation);
		}
	}

	/**
	 * Generates a quadrangulation from a triangulation by "inverting" triangles
	 * (for each triangle, create edges joining its centroid to each of its
	 * vertices).
	 * <p>
	 * This approach tends to create a denser quad mesh than
	 * {@link #edgeCollapseQuadrangulation(IIncrementalTin, boolean)
	 * <code>edgeCollapseQuadrangulation()</code>} on the same input.
	 * 
	 * @param triangulation     a triangulation mesh
	 * @param preservePerimeter whether to preserve the perimeter of the input
	 *                          triangulation; when true, retains edges that lie on
	 *                          the perimeter of the triangulation mesh that would
	 *                          have otherwise been removed (this results in some
	 *                          triangles being included in the output).
	 * @return a GROUP PShape, where each child shape is one quadrangle
	 * @since 1.2.0
	 */
	public static PShape centroidQuadrangulation(final IIncrementalTin triangulation, final boolean preservePerimeter) {
		final boolean unconstrained = triangulation.getConstraints().isEmpty();
		final HashSet<PEdge> edges = new HashSet<>();
		TriangleCollector.visitSimpleTriangles(triangulation, t -> {
			final IConstraint constraint = t.getContainingRegion();
			if (unconstrained || (constraint != null && constraint.definesConstrainedRegion())) {
				Vertex centroid = centroid(t);
				edges.add(new PEdge(centroid.getX(), centroid.getY(), t.getVertexA().x, t.getVertexA().y));
				edges.add(new PEdge(centroid.getX(), centroid.getY(), t.getVertexB().x, t.getVertexB().y));
				edges.add(new PEdge(centroid.getX(), centroid.getY(), t.getVertexC().x, t.getVertexC().y));
			}
		});

		if (preservePerimeter) {
			List<IQuadEdge> perimeter = triangulation.getPerimeter();
			triangulation.edges().forEach(edge -> {
				if (edge.isConstrainedRegionBorder() || (unconstrained && perimeter.contains(edge))) {
					edges.add(new PEdge(edge.getA().x, edge.getA().y, edge.getB().x, edge.getB().y));
				}
			});
		}

		final PShape quads = PGS.polygonizeEdges(edges);
		if (triangulation.getConstraints().size() < 2) { // assume constraint 1 is the boundary (not a hole)
			return quads;
		} else {
			return removeHoles(quads, triangulation);
		}
	}

	/**
	 * Removes (what should be) holes from a polygonized quadrangulation.
	 * <p>
	 * When the polygonizer is applied to the collapsed triangles of a
	 * triangulation, it cannot determine which collapsed regions represent holes in
	 * the quadrangulation and will consequently fill them in. The subroutine below
	 * restores holes/topology, detecting which polygonized face(s) are original
	 * holes. Note the geometry of the original hole/constraint and its associated
	 * polygonized face are different, since quads are polygonized, not triangles
	 * (hence an overlap metric is used to match candidates).
	 * 
	 * @param faces         faces of the quadrangulation
	 * @param triangulation
	 * @return
	 */
	private static PShape removeHoles(PShape faces, IIncrementalTin triangulation) {
		List<IConstraint> holes = new ArrayList<>(triangulation.getConstraints()); // copy list
		holes = holes.subList(1, holes.size()); // slice off perimeter constraint (not a hole)

		STRtree tree = new STRtree();
		holes.stream().map(constraint -> constraint.getVertices()).iterator().forEachRemaining(vertices -> {
			CoordinateList coords = new CoordinateList(); // coords of constraint
			vertices.forEach(v -> coords.add(new Coordinate(v.x, v.y)));
			coords.closeRing();

			if (!Orientation.isCCWArea(coords.toCoordinateArray())) { // triangulation holes are CW
				Polygon polygon = PGS.GEOM_FACTORY.createPolygon(coords.toCoordinateArray());
				tree.insert(polygon.getEnvelopeInternal(), polygon);
			}
		});

		List<PShape> nonHoles = PGS_Conversion.getChildren(faces).parallelStream().filter(quad -> {
			/*
			 * If quad overlaps with a hole detect whether it *is* that hole via Hausdorff
			 * Similarity.
			 */
			final Geometry g = PGS_Conversion.fromPShape(quad);

			@SuppressWarnings("unchecked")
			List<Polygon> matches = tree.query(g.getEnvelopeInternal());

			for (Polygon m : matches) {
				try {
					// PGS_ShapePredicates.overlap() inlined here
					Geometry overlap = ClipperPGS.intersectPolygons(m, g);
					double a1 = g.getArea();
					double a2 = m.getArea();
					double total = a1 + a2;
					double aOverlap = overlap.getArea();
					double w1 = a1 / total;
					double w2 = a2 / total;

					double similarity = w1 * (aOverlap / a1) + w2 * (aOverlap / a2);
					if (similarity > 0.33) { // magic constant, unsure what the best value is
						return false; // is hole; keep=false
					}
				} catch (Exception e) { // catch occasional noded error
					continue;
				}

			}
			return true; // is not hole; keep=true
		}).collect(Collectors.toList());

		return PGS_Conversion.flatten(nonHoles);
	}

	/**
	 * Produces a quadrangulation from a point set. The resulting quadrangulation
	 * has a characteristic spiral pattern.
	 * 
	 * @param points
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.2.0
	 */
	public static PShape spiralQuadrangulation(List<PVector> points) {
		SpiralQuadrangulation sq = new SpiralQuadrangulation(points);
		return PGS.polygonizeEdges(sq.getQuadrangulationEdges());
	}

	/**
	 * Calculate the longest edge of a given triangle.
	 */
	private static IQuadEdge findLongestEdge(final SimpleTriangle t) {
		if (t.getEdgeA().getLength() > t.getEdgeB().getLength()) {
			if (t.getEdgeC().getLength() > t.getEdgeA().getLength()) {
				return t.getEdgeC();
			} else {
				return t.getEdgeA();
			}
		} else {
			if (t.getEdgeC().getLength() > t.getEdgeB().getLength()) {
				return t.getEdgeC();
			} else {
				return t.getEdgeB();
			}
		}
	}

	private static double[] midpoint(final IQuadEdge edge) {
		final Vertex a = edge.getA();
		final Vertex b = edge.getB();
		return new double[] { (a.x + b.x) / 2d, (a.y + b.y) / 2d };
	}

	private static Vertex centroid(final SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new Vertex(x, y, 0);
	}

}
