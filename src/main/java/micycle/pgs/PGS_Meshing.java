package micycle.pgs;

import static micycle.pgs.PGS_Conversion.toPShape;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.polygonize.QuickPolygonizer;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;
import org.tinspin.index.PointIndex;
import org.tinspin.index.kdtree.KDTree;

import micycle.pgs.color.RGB;
import micycle.pgs.utility.IncrementalTinDual;
import micycle.pgs.utility.SpiralQuadrangulation;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Mesh generation (excluding triangulation).
 * <p>
 * So far, each method within this class processes an existing Delaunay
 * triangulation. You must first process a shape using
 * {@link PGS_Triangulation#delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
 * delaunayTriangulationMesh()} and then feed it to methods within PGS_Meshing.
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
	 * {@link PGS_Triangulation#delaunayTriangulation(PShape, List, boolean, int, boolean)
	 * triangulation} and a {@link micycle.pgs.PGS_Processing#partition(PShape)
	 * partition}).
	 * <p>
	 * Note that this method processes a Delaunay triangulation. Process a shape
	 * using
	 * {@link PGS_Triangulation#delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * delaunayTriangulationMesh()} first and then feed it to this method.
	 * 
	 * @param triangulation     a triangulation mesh
	 * @param preservePerimeter whether to preserve the perimeter of the input
	 *                          triangulation; when true, retains edges that lie on
	 *                          the perimeter of the triangulation mesh that would
	 *                          have otherwise been removed according to the
	 *                          urquhart condition.
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.1.0
	 * @see #gabrielFaces(IIncrementalTin, boolean)
	 */
	@SuppressWarnings("unchecked")
	public static PShape urquhartFaces(final IIncrementalTin triangulation, final boolean preservePerimeter) {
		final HashSet<IQuadEdge> edges = new HashSet<>();
		final HashSet<IQuadEdge> uniqueLongestEdges = new HashSet<>();

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

		final Polygonizer polygonizer = new QuickPolygonizer(false);
		polygonizer.setCheckRingsValid(false);
		edges.removeAll(uniqueLongestEdges);
		edges.forEach(edge -> polygonizer
				.add(PGS.GEOM_FACTORY.createLineString(new Coordinate[] { toCoord(edge.getA()), toCoord(edge.getB()) })));

		final PShape out = new PShape(PConstants.GROUP);
		polygonizer.getPolygons().forEach(p -> {
			final PShape face = toPShape((Polygon) p);
			face.setStrokeWeight(3);
			out.addChild(face);
		});

		return out;
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
	 * @param preservePerimeter whether to preserve the perimeter of the input
	 *                          triangulation; when true, retains edges that lie on
	 *                          the perimeter of the triangulation mesh that would
	 *                          have otherwise been removed according to the
	 *                          urquhart condition.
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.1.0
	 * @see #urquhartFaces(IncrementalTin, boolean)
	 */
	@SuppressWarnings("unchecked")
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

		final Polygonizer polygonizer = new QuickPolygonizer(false);
		polygonizer.setCheckRingsValid(false);
		edges.forEach(edge -> polygonizer
				.add(PGS.GEOM_FACTORY.createLineString(new Coordinate[] { toCoord(edge.getA()), toCoord(edge.getB()) })));

		final PShape out = new PShape(PConstants.GROUP);
		polygonizer.getPolygons().forEach(p -> {
			final PShape face = toPShape((Polygon) p);
			face.setStrokeWeight(3);
			out.addChild(face);
		});
		return out;
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
		// PGS_Conversion.disableAllFill(dualMesh);
		// PGS_Conversion.setAllStrokeColor(dualMesh, RGB.composeColor(0, 50, 255), 2);
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
	 * Produces a quadrangulation from a point set. The resulting quadrangulation
	 * has a characteristic spiral pattern.
	 * 
	 * @param points
	 * @return a GROUP PShape where each child shape is a single face
	 * @since 1.2.0
	 */
	public static PShape spiralQuadrangulation(List<PVector> points) {
		SpiralQuadrangulation sq = new SpiralQuadrangulation(points);
		List<SegmentString> segments = new ArrayList<>(sq.getQuadrangulationEdges().size());
		sq.getQuadrangulationEdges().forEach(e -> segments.add(PGS.createSegmentString(e.a, e.b)));
		return PGS_Conversion.toPShape(PGS.polygonizeSegments(segments));
	}

	private static Coordinate toCoord(final Vertex v) {
		return new Coordinate(v.x, v.y);
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

}
