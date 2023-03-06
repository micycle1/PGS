package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.operation.overlayng.CoverageUnion;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pgs.commons.PEdge;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Boolean set-operations for 2D shapes.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_ShapeBoolean {

	private PGS_ShapeBoolean() {
	}

	/**
	 * Computes a shape representing the area which is common to both input shapes
	 * (i.e. the shape formed by intersection of <code>a</code> and <code>b</code>).
	 * <p>
	 * Note: Intersecting a polygon with a path/linestring will crop the path to the
	 * polygon.
	 * <p>
	 * Note: The intersecting parts of faces of a mesh-like shape will be collapsed
	 * into a single area during intersection; to intersect such a shape and
	 * preserve how each face is intersected individually, use
	 * {@link #intersectMesh(PShape, PShape) intersectMesh()}.
	 * 
	 * @param a The first shape.
	 * @param b The second shape.
	 * @return A new shape representing the intersection of the two input shapes.
	 */
	public static PShape intersect(final PShape a, final PShape b) {
		return toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.INTERSECTION));
	}

	/**
	 * Intersects a mesh-like shape / polygonal coverage with a polygonal area,
	 * preserving individual faces/features of the mesh during the operation.
	 * <p>
	 * When a mesh-like shape / polygonal coverage is intersected with a whole
	 * polygon using {@link #intersect(PShape, PShape) intersect(a, b)}, the result
	 * is a <b>single</b> polygon comprising the combined/dissolved area of all
	 * intersecting mesh faces. Sometimes this behaviour is desired whereas other
	 * times it is not -- this method can be used to preserve how each face is
	 * intersected individually.
	 * <p>
	 * Using this method is faster than calling {@link #intersect(PShape, PShape)
	 * intersect(a, b)} repeatedly for every face of a mesh-like shape
	 * <code>a</code> against an area <code>b</code>.
	 * 
	 * @param mesh a mesh-like GROUP shape
	 * @param area a polygonal shape
	 * @return a GROUP shape, where each child shape is the union of one mesh face
	 *         and the area
	 * @since 1.3.0
	 */
	public static PShape intersectMesh(final PShape mesh, final PShape area) {
		final Geometry g = fromPShape(area);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(g);
		// @formatter:off
		List<Geometry> faces = PGS_Conversion.getChildren(mesh).parallelStream()
				.map(PGS_Conversion::fromPShape)
				.map(f -> cache.containsProperly(f) ? f : OverlayNG.overlay(f, g, OverlayNG.INTERSECTION))
				.collect(Collectors.toList());
		// @formatter:on
		return PGS_Conversion.toPShape(faces);
	}

	/**
	 * "Glues" shapes together so they become a single combined shape with the sum
	 * of its areas.
	 * 
	 * @return A∪B
	 * @see #union(PShape...)
	 */
	public static PShape union(final PShape a, final PShape b) {
		return toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.UNION));
	}

	/**
	 * Unions any variable number of shapes.
	 * 
	 * @param shapes
	 * @return
	 * @see #union(PShape, PShape)
	 * @see #union(PShape...)
	 */
	public static PShape union(final List<PShape> shapes) {
		Collection<Geometry> polygons = new ArrayList<>();
		shapes.forEach(s -> polygons.add(fromPShape(s)));
		return toPShape(UnaryUnionOp.union(polygons));
	}

	/**
	 * Unions any variable number of shapes.
	 * 
	 * @param shapes varArgs
	 * @return
	 * @see #union(PShape, PShape)
	 * @see #union(List)
	 */
	public static PShape union(PShape... shapes) {
		return union(Arrays.asList(shapes));
	}

	/**
	 * Unions/flattens/merges/dissolves a mesh-like PShape (that is, a GROUP PShape
	 * whose children represent faces that share edges) into a single shape that
	 * represents the boundary of the mesh. This method is optimised for meshes, and
	 * is accordingly much faster than unioning the mesh faces together using other
	 * methods.
	 * 
	 * @param mesh a GROUP pshape whose children shapes form a mesh (join/overlap at
	 *             edges). The mesh may contain holes.
	 * @return
	 * @since 1.2.0
	 */
	public static PShape unionMesh(final PShape mesh) {
		if (mesh.getChildCount() < 2 || mesh.getKind() != PConstants.GROUP) {
			if (mesh.getChildCount() == 1) {
				return mesh.getChild(0);
			}
			return mesh;
		}

		return unionMeshWithHoles(mesh);
	}
	
	private static PShape unionMeshWithHoles(final PShape mesh) {
		return toPShape(CoverageUnion.union(PGS_Conversion.fromPShape(mesh)));
	}

	/**
	 * Faster than JTS CoverageUnion but does not handle meshes with holes.
	 * 
	 * @deprecated
	 */
	private static PShape unionMeshWithoutHoles(final PShape mesh) {
		Map<PEdge, Integer> edges = new HashMap<>(mesh.getChildCount() * 3);

		final List<PEdge> allEdges;

		/*
		 * Compute set of unique edges belonging to the mesh (this set is equivalent to
		 * the boundary, assuming a holeless mesh).
		 */
		for (PShape child : mesh.getChildren()) {
			for (int i = 0; i < child.getVertexCount(); i++) {
				final PVector a = child.getVertex(i);
				final PVector b = child.getVertex((i + 1) % child.getVertexCount());
				if (!a.equals(b)) {
					PEdge edge = new PEdge(a, b);
					edges.merge(edge, 1, Integer::sum);
				}
			}
		}

		allEdges = edges.entrySet().stream().filter(e -> e.getValue() == 1).map(e -> e.getKey()).collect(Collectors.toList());

		/*
		 * Now find a sequential/winding order for the vertices of the boundary. The
		 * vertices output fromEdges() is not closed, so close it afterwards (assumes
		 * the input to unionMesh() was indeed closed and valid).
		 */
		final List<PVector> orderedVertices = PGS.fromEdges(allEdges);
		if (!orderedVertices.get(0).equals(orderedVertices.get(orderedVertices.size() - 1))) {
			orderedVertices.add(orderedVertices.get(0)); // close vertex list for fromPVector()
		}

		return PGS_Conversion.fromPVector(orderedVertices);
	}

	/**
	 * Subtracts one shape (b) from another shape (a) and returns the resulting
	 * shape. This procedure is also known as "difference".
	 * <p>
	 * Subtract is the opposite of {@link #union(PShape, PShape) union()}.
	 * 
	 * @param a The PShape from which the other PShape will be subtracted.
	 * @param b The PShape that will be subtracted from the first PShape.
	 * @return A new PShape representing the difference between the two input
	 *         shapes.
	 * @see #simpleSubtract(PShape, PShape)
	 */
	public static PShape subtract(final PShape a, final PShape b) {
		return toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.DIFFERENCE));
	}

	/**
	 * Subtracts <code>holes</code> from the <code>shell</code>, without geometric
	 * processing.
	 * <p>
	 * Rather than performing geometric overlay, this method simply appends the
	 * holes as contours to the shell. For this reason, it is faster than
	 * {@link #subtract(PShape, PShape) subtract()} but this method produces valid
	 * results only if <b>all holes lie inside the shell</b>.
	 * 
	 * @param shell polygonal shape
	 * @param holes single polygon, or GROUP shape, whose children are holes that
	 *              lie within the shell
	 * @since 1.3.1
	 * @see #subtract(PShape, PShape)
	 * @return
	 */
	public static PShape simpleSubtract(PShape shell, PShape holes) {
		List<PShape> children = PGS_Conversion.getChildren(holes);
		List<List<PVector>> holez;
		if (holes.getChildCount() == 0) {
			children.add(holes);
		}
		holez = children.stream().map(PGS_Conversion::toPVector).collect(Collectors.toList());

		return PGS_Conversion.fromPVector(PGS_Conversion.toPVector(shell), holez);
	}

	/**
	 * Subtracts a polygonal area from a mesh-like shape / polygonal coverage,
	 * preserving individual faces/features of the mesh during the operation.
	 * <p>
	 * When polygon is subtracted from a mesh-like shape / polygonal coverage using
	 * {@link #subtract(PShape, PShape) subtract(a, b)}, the result is a
	 * <b>single</b> polygon comprising the combined/dissolved area of all remaining
	 * mesh face parts. Sometimes this behaviour is desired whereas other times it
	 * is not -- this method can be used to preserve how each face is subtracted
	 * from individually.
	 * <p>
	 * Using this method is faster than calling {@link #subtract(PShape, PShape)
	 * subtract(a, b)} repeatedly on every face of a mesh-like shape <code>a</code>.
	 * 
	 * @param mesh a mesh-like GROUP shape
	 * @param area a polygonal shape
	 * @return a GROUP shape, where each child shape is the subtraction of the area
	 *         from one mesh face
	 * @since 1.3.0
	 */
	public static PShape subtractMesh(PShape mesh, PShape area) {
		final Geometry g = fromPShape(area);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(g);
		// @formatter:off
		List<Geometry> faces = PGS_Conversion.getChildren(mesh).parallelStream()
				.map(PGS_Conversion::fromPShape)
				.map(f -> cache.containsProperly(f) ? null : OverlayNG.overlay(f, g, OverlayNG.DIFFERENCE))
				.filter(Objects::nonNull)
				.collect(Collectors.toList());
		// @formatter:on
		return PGS_Conversion.toPShape(faces);
	}

	/**
	 * Computes the <i>symmetric difference</i> (the parts that the shapes do not
	 * have in common) between two shapes and returns the result as a new shape.
	 * 
	 * @param a The first shape.
	 * @param b The second shape.
	 * @return A new shape representing the symmetric difference between the two
	 *         input shapes.
	 */
	public static PShape symDifference(PShape a, PShape b) {
		return toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.SYMDIFFERENCE));
	}

	/**
	 * Computes the shape's complement (or inverse) against a rectangle of specified
	 * width and height.
	 * 
	 * @param shape  The PShape for which the complement is to be calculated.
	 * @param width  width of the rectangle plane to subtract shape from
	 * @param height height of the rectangle plane to subtract shape from
	 * @return A new PShape representing the complement of the input shape within
	 *         the specified rectangle.
	 */
	public static PShape complement(PShape shape, double width, double height) {
		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(4);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createRectangle().difference(fromPShape(shape)));
	}

}
