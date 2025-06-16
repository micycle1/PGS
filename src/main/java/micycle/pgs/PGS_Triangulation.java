package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.triangulate.polygon.PolygonTriangulator;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.edge.QuadEdge;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.HilbertSort;
import org.tinfour.utils.TriangleCollector;

import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.color.Colors;
import micycle.pgs.commons.Nullable;
import micycle.pgs.commons.PEdge;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Delaunay and earcut triangulation of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Triangulation {

	private PGS_Triangulation() {
	}

	/**
	 * Generates a constrained Delaunay Triangulation from the given shape.
	 * 
	 * @param shape the shape whose vertices to generate a triangulation from
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 */
	public static PShape delaunayTriangulation(PShape shape) {
		return delaunayTriangulation(shape, null, true, 0, true);
	}

	/**
	 * Generates a Delaunay Triangulation from the given shape. The triangulation
	 * can be both constrained (meaning the triangulation is masked by the original
	 * shape) and refined (meaning additional points are inserted, usually leading
	 * to more uniform triangle shapes and sizes).
	 * 
	 * @param shape         the shape whose vertices to generate a triangulation
	 *                      from
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to maintain the Delaunay nature when
	 *                      constraining the triangulation, and whether to check
	 *                      that centroid locations lie within the shape during
	 *                      refinement. When pretty=true, triangles in the
	 *                      triangulation may be slightly more regular in
	 *                      shape/size. There is a small performance overhead which
	 *                      becomes more considerable at higher refinement levels.
	 *                      When constrain=false and refinements=0, this argument
	 *                      has no effect.
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #delaunayTriangulationPoints(PShape, Collection, boolean, int, boolean)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 */
	public static PShape delaunayTriangulation(PShape shape, @Nullable Collection<PVector> steinerPoints, boolean constrain, int refinements, boolean pretty) {
		final IIncrementalTin tin = delaunayTriangulationMesh(shape, steinerPoints, constrain, refinements, pretty);
		return toPShape(tin);
	}

	/**
	 * Generates a Delaunay Triangulation from a collection of points.
	 * 
	 * @param points the point collection to triangulate
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 * @since 1.1.0
	 */
	public static PShape delaunayTriangulation(Collection<PVector> points) {
		return delaunayTriangulation(null, points, false, 0, false);
	}

	/**
	 * Generates a Delaunay Triangulation having a shape constraint from a
	 * collection of points.
	 * 
	 * @param points          the collection of points to triangulate
	 * @param shapeConstraint a shape that defines the boundary constraint on the
	 *                        triangulation
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 * @since 2.0
	 */
	public static PShape delaunayTriangulation(Collection<PVector> points, PShape shapeConstraint) {
		final IIncrementalTin tin = delaunayTriangulationMesh(shapeConstraint, points, true, 0, true, false);
		return toPShape(tin);
	}

	/**
	 * Generates a constrained Delaunay Triangulation from a collection of points.
	 * <p>
	 * This method returns the triangulation as a list of points, rather than a
	 * PShape.
	 * 
	 * @param shape the shape whose vertices to generate a triangulation from
	 * @return List of PVector coordinates, where each consecutive triplet of
	 *         coordinates are the 3 vertices belonging to one triangle
	 */
	public static List<PVector> delaunayTriangulationPoints(PShape shape) {
		return delaunayTriangulationPoints(shape, null, true, 0, true);
	}

	/**
	 * Generates a Delaunay Triangulation from the given shape. The triangulation
	 * can be both constrained (meaning the triangulation is masked by the original
	 * shape) and refined (meaning additional points are inserted, usually leading
	 * to more uniform triangle shapes and sizes).
	 * <p>
	 * This method returns the triangulation as a list of points, rather than a
	 * PShape.
	 * 
	 * @param shape         the shape whose vertices to generate a triangulation of
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to maintain the Delaunay nature when
	 *                      constraining the triangulation, and whether to check
	 *                      that centroid locations lie within the shape during
	 *                      refinement. When pretty=true, triangles in the
	 *                      triangulation may be slightly more regular in
	 *                      shape/size. There is a small performance overhead which
	 *                      becomes more considerable at higher refinement levels.
	 *                      When constrain=false and refinements=0, this argument
	 *                      has no effect.
	 * @return List of PVector coordinates, where each consecutive triplet of
	 *         coordinates are the 3 vertices belonging to one triangle
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 */
	public static List<PVector> delaunayTriangulationPoints(PShape shape, @Nullable Collection<PVector> steinerPoints, boolean constrain, int refinements,
			boolean pretty) {
		final IIncrementalTin tin = delaunayTriangulationMesh(shape, steinerPoints, constrain, refinements, pretty);

		final ArrayList<PVector> triangles = new ArrayList<>();
		final Consumer<Vertex[]> triangleVertexConsumer = t -> {
			triangles.add(toPVector(t[0]));
			triangles.add(toPVector(t[1]));
			triangles.add(toPVector(t[2]));
		};
		if (constrain) {
			TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		} else {
			TriangleCollector.visitTriangles(tin, triangleVertexConsumer);
		}
		return triangles;
	}

	/**
	 * Generates a Delaunay Triangulation from a collection of points.
	 * <p>
	 * This method returns the triangulation as a list of points, rather than a
	 * PShape.
	 * 
	 * @param points the point collection to triangulate
	 * @return List of PVector coordinates, where each consecutive triplet of
	 *         coordinates are the 3 vertices belonging to one triangle
	 * @see #delaunayTriangulationPoints(PShape, Collection, boolean, int, boolean)
	 * @since 1.1.0
	 */
	public static List<PVector> delaunayTriangulationPoints(Collection<PVector> points) {
		return delaunayTriangulationPoints(null, points, false, 0, false);
	}

	/**
	 * Generates a constrained Delaunay Triangulation from the given shape.
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 * 
	 * @param shape the shape whose vertices to generate a triangulation from
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 */
	public static IIncrementalTin delaunayTriangulationMesh(PShape shape) {
		return delaunayTriangulationMesh(shape, null, true, 0, true);
	}

	/**
	 * Generates a Delaunay Triangulation of the given shape and returns it in raw
	 * form as a Triangulated Irregular Network (mesh).
	 * <p>
	 * The triangulation can be constrained to the shape's boundary and refined by
	 * adding additional points, resulting in more uniform triangle shapes and
	 * sizes.
	 * 
	 * @param shape         the shape to generate a triangulation from. <b>Can be
	 *                      null</b>.
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     whether to constrain the triangulation to the shape's
	 *                      boundary. If using a shape, it is recommended to set
	 *                      this to true.
	 * @param refinements   The number of times to subdivide the triangulation by
	 *                      inserting the centroid of each triangle. Should be 0 or
	 *                      greater, typically no more than 5.
	 * @param pretty        Whether to maintain Delaunay nature when constraining
	 *                      the triangulation and check that centroid locations are
	 *                      within the shape during refinement. This can result in
	 *                      more regular triangle shapes and sizes, but with a
	 *                      performance overhead that increases with higher
	 *                      refinement levels. Has no effect if
	 *                      <code>constrain=false</code> and
	 *                      <code>refinements=0</code>.
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 * @see #delaunayTriangulationPoints(PShape, Collection, boolean, int, boolean)
	 */
	public static IIncrementalTin delaunayTriangulationMesh(@Nullable PShape shape, @Nullable Collection<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		return delaunayTriangulationMesh(shape, steinerPoints, constrain, refinements, pretty, true);

	}

	/**
	 * Generates a Delaunay Triangulation from a collection of points.
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 * 
	 * @param points the point collection to triangulate
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * @since 1.1.0
	 */
	public static IIncrementalTin delaunayTriangulationMesh(Collection<PVector> points) {
		return delaunayTriangulationMesh(null, points, false, 0, false);
	}

	/**
	 * Generates a Delaunay Triangulation having a shape constraint from a
	 * collection of points.
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 * 
	 * @param points          the collection of points to triangulate
	 * @param shapeConstraint a shape that defines the boundary constraint on the
	 *                        triangulation
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * @since 2.0
	 */
	public static IIncrementalTin delaunayTriangulationMesh(Collection<PVector> points, PShape shapeConstraint) {
		return delaunayTriangulationMesh(shapeConstraint, points, true, 0, false, false);
	}

	/**
	 * @param shape               the shape to generate a triangulation from. <b>Can
	 *                            be null</b>.
	 * @param steinerPoints       A list of additional points to insert into the
	 *                            triangulation in addition to the vertices of the
	 *                            input shape. <b>Can be null</b>.
	 * @param constrain           whether to constrain the triangulation to the
	 *                            shape's boundary. If using a shape, it is
	 *                            recommended to set this to true.
	 * @param refinements         The number of times to subdivide the triangulation
	 *                            by inserting the centroid of each triangle. Should
	 *                            be 0 or greater, typically no more than 5.
	 * @param pretty              Whether to maintain Delaunay nature when
	 *                            constraining the triangulation and check that
	 *                            centroid locations are within the shape during
	 *                            refinement. This can result in more regular
	 *                            triangle shapes and sizes, but with a performance
	 *                            overhead that increases with higher refinement
	 *                            levels. Has no effect if
	 *                            <code>constrain=false</code> and
	 *                            <code>refinements=0</code>.
	 * @param insertShapeVertices Determines how input shape vertices are treated:
	 *                            as part of the triangulation (=true), or as a
	 *                            boundary constraint only (=false).
	 */
	private static IIncrementalTin delaunayTriangulationMesh(@Nullable PShape shape, @Nullable Collection<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty, boolean insertShapeVertices) {
		Geometry g = shape == null ? PGS.GEOM_FACTORY.createEmpty(2) : fromPShape(shape);
		final IncrementalTin tin = new IncrementalTin(10);

		final List<Vertex> vertices = new ArrayList<>();
		final Coordinate[] coords = g.getCoordinates();
		int vIndex = 0;
		if (insertShapeVertices) {
			for (vIndex = 0; vIndex < coords.length; vIndex++) {
				vertices.add(new Vertex(coords[vIndex].x, coords[vIndex].y, Double.NaN, vIndex));
			}
		}

		if (steinerPoints != null) {
			for (PVector v : steinerPoints) { // add steiner points
				vertices.add(new Vertex(v.x, v.y, Double.NaN, vIndex++));
			}
		}

		HilbertSort hs = new HilbertSort();
		hs.sort(vertices); // prevent degenerate insertion
		tin.add(vertices, null); // initial triangulation

		if (refinements > 0) {

			final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
			final ArrayList<Vertex> refinementVertices = new ArrayList<>();

			/*
			 * All vertices must be added to Tin before constraints are added (hence the
			 * need for pointLocator, since visitSimpleTriangles() visits triangles that lie
			 * outside the shape at this stage.
			 */
			for (int i = 0; i < refinements; i++) {
				/*
				 * Must be iterated like this so that triangles to refine are delaunay at each
				 * stage (rather than recurse the original triangles).
				 */
				refinementVertices.clear();
				TriangleCollector.visitSimpleTriangles(tin, t -> {
					if (t.getArea() > 99) { // don't refine small triangles
						final Coordinate center = centroid(t); // use centroid rather than circumcircle center
						if (pretty || pointLocator.locate(center) != Location.EXTERIOR) {
							refinementVertices.add(new Vertex(center.x, center.y, Double.NaN));
						}
					}
				});
				tin.add(refinementVertices, null); // add refinement (steiner) points
			}
		}

		if (constrain) {
			// If geom is a point set, constrain tin using its concave hull.
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOINT)) {
				g = fromPShape(PGS_Hull.concaveHullBFS2(PGS_Conversion.toPVector(shape), 0.3));
			}
			List<IConstraint> constraints = new ArrayList<>();
			for (int n = 0; n < g.getNumGeometries(); n++) {
				boolean exterior = true;

				if (g instanceof Polygonal) {
					LinearRingIterator lri = new LinearRingIterator(g.getGeometryN(n));
					for (LinearRing ring : lri) {
						final List<Vertex> points = new ArrayList<>();
						final Coordinate[] c = ring.getCoordinates();
						if (c.length == 0) {
							exterior = false;
							continue;
						}

						for (Coordinate element : c) {
							points.add(new Vertex(element.x, element.y, Double.NaN));
						}
						/*
						 * In Tinfour, the shape exterior must be CCW and the holes must be CW. This is
						 * true for most PShapes, but some shapes (like those created from fonts) may
						 * have the rings orientated the other way, which needs to be corrected.
						 */
						if ((exterior && !Orientation.isCCWArea(c)) || (!exterior && Orientation.isCCWArea(c))) {
							Collections.reverse(points);
						}
						constraints.add(new PolygonConstraint(points));
						exterior = false;
					}
				}
			}
			if (!constraints.isEmpty()) {
				tin.addConstraints(constraints, pretty);
			}
		}

		return tin;
	}

	/**
	 * Creates a Delaunay triangulation of the shape where additional steiner
	 * points, populated by poisson sampling, are included.
	 * 
	 * @param shape
	 * @param spacing (Minimum) spacing between poisson points
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #poissonTriangulationPoints(PShape, double)
	 */
	public static PShape poissonTriangulation(PShape shape, double spacing) {
		final IIncrementalTin tin = poissonTriangulationMesh(shape, spacing);

		final PShape triangulation = new PShape(PConstants.GROUP);

		TriangleCollector.visitTrianglesConstrained(tin, t -> {
			final PShape triangle = new PShape(PShape.PATH);
			triangle.beginShape();
			triangle.vertex((float) t[0].x, (float) t[0].y);
			triangle.vertex((float) t[1].x, (float) t[1].y);
			triangle.vertex((float) t[2].x, (float) t[2].y);
			triangle.endShape(PConstants.CLOSE);
			triangulation.addChild(triangle);
		});

		PGS_Conversion.setAllFillColor(triangulation, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(triangulation, Colors.PINK, 2);

		return triangulation;
	}

	/**
	 * Creates a Delaunay triangulation of the shape where additional steiner
	 * points, populated by poisson sampling, are included.
	 * 
	 * @param shape
	 * @param spacing (Minimum) spacing between poisson points
	 * @return list of PVectors, where each successive triplet of PVectors
	 *         correspond to the 3 vertices of one triangle
	 * @see #poissonTriangulation(PShape, double)
	 */
	public static List<PVector> poissonTriangulationPoints(PShape shape, double spacing) {
		final IIncrementalTin tin = poissonTriangulationMesh(shape, spacing);

		final ArrayList<PVector> triangles = new ArrayList<>();
		TriangleCollector.visitTrianglesConstrained(tin, t -> {
			triangles.add(toPVector(t[0]));
			triangles.add(toPVector(t[1]));
			triangles.add(toPVector(t[2]));
		});
		return triangles;
	}

	/**
	 * Creates a Delaunay triangulation of the shape where additional steiner
	 * points, populated by poisson sampling, are included.
	 * <p>
	 * This method returns the triangulation in its raw form: a
	 * TriangulatedIrregular Network (mesh).
	 * 
	 * @param shape
	 * @param spacing (Minimum) spacing between poisson points
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #poissonTriangulation(PShape, double)
	 */
	public static IIncrementalTin poissonTriangulationMesh(PShape shape, double spacing) {
		final Envelope e = fromPShape(shape).getEnvelopeInternal();

		final List<PVector> poissonPoints = PGS_PointSet.poisson(e.getMinX(), e.getMinY(), e.getMinX() + e.getWidth(), e.getMinY() + e.getHeight(), spacing, 0);

		final IIncrementalTin tin = delaunayTriangulationMesh(shape, poissonPoints, true, 0, false);
		return tin;
	}

	/**
	 * Computes a triangulation of the shape according to the ear clipping
	 * ("earcut") method. The triangulation is constrained to the shape outline.
	 * 
	 * @param shape shape whose vertices to triangulate
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @since 1.1.0
	 * @since Supports holes since 1.3.0
	 */
	public static PShape earCutTriangulation(PShape shape) {
		PolygonTriangulator pt = new PolygonTriangulator(fromPShape(shape));
		return PGS_Conversion.toPShape(pt.getResult());
	}

	/**
	 * Converts a triangulated mesh object to a PShape representing the
	 * triangulation.
	 * 
	 * @param triangulation the IIncrementalTin object to convert
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @since 1.4.0
	 */
	public static PShape toPShape(IIncrementalTin triangulation) {
		final PShape out = new PShape(PConstants.GROUP);

		final Consumer<Vertex[]> triangleVertexConsumer = t -> {
			final PShape triangle = new PShape(PShape.PATH);
			triangle.beginShape();
			triangle.vertex((float) t[0].x, (float) t[0].y);
			triangle.vertex((float) t[1].x, (float) t[1].y);
			triangle.vertex((float) t[2].x, (float) t[2].y);
			triangle.endShape(PConstants.CLOSE);
			out.addChild(triangle);
		};

		if (!triangulation.getConstraints().isEmpty()) {
			TriangleCollector.visitTrianglesConstrained(triangulation, triangleVertexConsumer);
		} else {
			TriangleCollector.visitTriangles(triangulation, triangleVertexConsumer);
		}

		PGS_Conversion.setAllFillColor(out, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(out, Colors.PINK, 2);

		return out;
	}

	/**
	 * Finds the graph equivalent to a triangulation. Graph vertices are
	 * triangulation vertices; graph edges are triangulation edges.
	 * <p>
	 * The output is an undirected weighted graph of Processing primitives; edge
	 * weights are their euclidean length of their triangulation equivalent.
	 * 
	 * @param triangulation triangulation mesh
	 * @return
	 * @since 1.3.0
	 * @see #toTinfourGraph(IIncrementalTin)
	 * @see #toDualGraph(IIncrementalTin)
	 */
	public static SimpleGraph<PVector, PEdge> toGraph(IIncrementalTin triangulation) {
		final SimpleGraph<PVector, PEdge> graph = new SimpleWeightedGraph<>(PEdge.class);
		final boolean notConstrained = triangulation.getConstraints().isEmpty();
		triangulation.edges().forEach(e -> {
//			if (isEdgeOnPerimeter(e)) {
//				return; // skip to next triangle
//			}
			if (notConstrained || e.isConstraintRegionMember()) {
				final IQuadEdge base = e.getBaseReference();
				PVector a = toPVector(base.getA());
				PVector b = toPVector(base.getB());
				PEdge edge = new PEdge(a, b);
				graph.addVertex(a);
				graph.addVertex(b);
				graph.addEdge(a, b, edge);
				graph.setEdgeWeight(edge, edge.length());
			}
		});
		return graph;
	}

	/**
	 * Finds the graph equivalent to a triangulation. Graph vertices are
	 * triangulation vertices; graph edges are triangulation edges.
	 * <p>
	 * The output is an undirected weighted graph of Tinfour primtives; edge weights
	 * are their euclidean length of their triangulation equivalent.
	 * 
	 * @param triangulation triangulation mesh
	 * @return
	 * @since 1.3.0
	 * @see #toGraph(IIncrementalTin)
	 * @see #toDualGraph(IIncrementalTin)
	 */
	public static SimpleGraph<Vertex, IQuadEdge> toTinfourGraph(IIncrementalTin triangulation) {
		final SimpleGraph<Vertex, IQuadEdge> graph = new SimpleWeightedGraph<>(IQuadEdge.class);
		final boolean notConstrained = triangulation.getConstraints().isEmpty();
		triangulation.edges().forEach(e -> {
//			if (isEdgeOnPerimeter(e)) {
//				return; // skip to next triangle
//			}
			if ((notConstrained || e.isConstraintRegionMember())) {
				final IQuadEdge base = e.getBaseReference();
				graph.addVertex(base.getA());
				graph.addVertex(base.getB());
				graph.addEdge(base.getA(), base.getB(), base);
				graph.setEdgeWeight(base, base.getLength());
			}
		});
		return graph;
	}

	/**
	 * Finds the dual-graph of a triangulation.
	 * <p>
	 * A dual graph of a triangulation has a vertex for each constrained triangle of
	 * the input, and an edge connecting each pair of triangles that are adjacent.
	 * <p>
	 * Edges are weighted, where weights represent a quality measure ([0...1], where
	 * 1 is a perfect square) for the quadrilateral formed by collapsing the shared
	 * edge between the two triangles.
	 * 
	 * @param triangulation triangulation mesh
	 * @since 1.3.0
	 * @return a simple weighted dual graph
	 * @see #toTinfourGraph(IIncrementalTin)
	 */
	public static SimpleGraph<SimpleTriangle, DefaultEdge> toDualGraph(IIncrementalTin triangulation) {
		final SimpleGraph<SimpleTriangle, DefaultEdge> graph = new SimpleWeightedGraph<>(DefaultEdge.class);
		// vertex supplier required in some graph algorithms
		QuadEdge dummy = new QuadEdge(0);
		graph.setVertexSupplier(() -> new SimpleTriangle(null, dummy, dummy, dummy));

		final boolean notConstrained = triangulation.getConstraints().isEmpty();
		final Map<IQuadEdge, SimpleTriangle> edgeMap = new HashMap<>(triangulation.countTriangles().getCount() * 3);

		TriangleCollector.visitSimpleTriangles(triangulation, t -> {
			final IConstraint constraint = t.getContainingRegion();
			if (notConstrained || (constraint != null && constraint.definesConstrainedRegion())) {
				edgeMap.put(t.getEdgeA(), t);
				edgeMap.put(t.getEdgeB(), t);
				edgeMap.put(t.getEdgeC(), t);
				graph.addVertex(t); // triangle
			}
		});

		/*
		 * Set the edge weight. This provides a quality measure for the quadrilateral
		 * element formed by collapsing the shared edge between two triangles. This
		 * weight is used in Blossom Quad algorithm.
		 */
		graph.vertexSet().forEach(t -> {
			final SimpleTriangle t1 = edgeMap.get(t.getEdgeA().getDual());
			final SimpleTriangle t2 = edgeMap.get(t.getEdgeB().getDual());
			final SimpleTriangle t3 = edgeMap.get(t.getEdgeC().getDual());
			if (t1 != null) {
				var edge = graph.addEdge(t, t1);
				if (edge != null) {
					double weight = computeQuadQuality(t, t1, t.getEdgeA());
					graph.setEdgeWeight(edge, weight);
				}
			}
			if (t2 != null) {
				var edge = graph.addEdge(t, t2);
				if (edge != null) {
					double weight = computeQuadQuality(t, t2, t.getEdgeB());
					graph.setEdgeWeight(edge, weight);
				}
			}
			if (t3 != null) {
				var edge = graph.addEdge(t, t3);
				if (edge != null) {
					double weight = computeQuadQuality(t, t3, t.getEdgeC());
					graph.setEdgeWeight(edge, weight);
				}
			}
		});

		return graph;
	}

	static PVector toPVector(final Vertex v) {
		return new PVector((float) v.getX(), (float) v.getY());
	}

	static PEdge toPEdge(final IQuadEdge e) {
		return new PEdge(toPVector(e.getA()), toPVector(e.getB()));
	}

	/**
	 * Determines whether an edge or its dual is on the perimeter.
	 *
	 * @param edge a valid instance
	 * @return true if the edge is on the perimeter; otherwise, false.
	 */
	private static boolean isEdgeOnPerimeter(IQuadEdge edge) {
		/*
		 * The logic here is that each edge defines one side of a triangle with vertices
		 * A, B, and C. Vertices A and B are the first and second vertices of the edge,
		 * vertex C is the opposite one. Triangles lying outside the Delaunay
		 * Triangulation have a "ghost" vertex for vertex C. Tinfour represents a ghost
		 * vertex with a null reference. So we test both the edge and its dual to see if
		 * their vertex C reference is null. Also note that vertex C is the second
		 * vertex of the forward edge from our edge of interest. Thus the C =
		 * edge.getForward().getB().
		 */
		return edge.getForward().getB() == null || edge.getForwardFromDual().getB() == null;
	}

	/**
	 * Computes the centroid/barycentre of a triangle.
	 */
	private static Coordinate centroid(final SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new Coordinate(x, y);
	}

	/**
	 * Computes the quality measure of a quadrilateral formed by merging two
	 * adjacent triangles. The quality measure is based on the maximum deviation of
	 * the interior angles of the quadrilateral from 90 degrees (π/2 radians). A
	 * perfect rectangle or square will have a quality measure of 1, while a
	 * quadrilateral with an angle of 0 or π radians will have a quality measure of
	 * 0.
	 * 
	 * @param t1         The first triangle sharing the edge with the second
	 *                   triangle.
	 * @param t2         The second triangle sharing the edge with the first
	 *                   triangle.
	 * @param sharedEdge The edge shared by t1 and t2, belonging to t1. The dual of
	 *                   this edge (sharedEdge.getDual()) represents the same edge
	 *                   but belonging to t2.
	 * @return A quality measure between 0 and 1, where 1 indicates a perfect
	 *         quadrilateral (all interior angles are 90 degrees) and 0 indicates a
	 *         degenerate quadrilateral (at least one interior angle is 0 or π
	 *         radians).
	 */
	private static double computeQuadQuality(SimpleTriangle t1, SimpleTriangle t2, IQuadEdge sharedEdge) {
		PVector a = toPVector(sharedEdge.getReverse().getA());
		PVector b = toPVector(sharedEdge.getA());
		PVector c = toPVector(sharedEdge.getReverseFromDual().getA());
		PVector d = toPVector(sharedEdge.getB());

		List<PVector> vertices = Arrays.asList(a, b, c, d);

		PShape quad = PGS_Conversion.fromPVector(vertices);

		Map<PVector, Double> angles = PGS_ShapePredicates.interiorAngles(quad);

		// Calculate maximum deviation from 90 degrees
		double maxDeviation = 0.0;
		for (double angle : angles.values()) {
			double deviation = Math.abs(Math.PI / 2 - angle);
			if (deviation > maxDeviation) {
				maxDeviation = deviation;
			}
		}

		double cost = 1 - (2 / Math.PI) * maxDeviation;
		cost = Math.max(cost, 0); // Ensure cost is non-negative
		return cost;
	}

}
