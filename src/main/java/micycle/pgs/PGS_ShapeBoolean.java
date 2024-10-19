package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.util.GeometryFixer;
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
	 * Calculates the intersection of two shapes, producing a new shape representing
	 * the shared area.
	 * <p>
	 * When intersecting a polygon with a path or linestring, the method will trim
	 * the path to the polygon's boundaries.
	 * <p>
	 * Note: When intersecting a mesh-like shape with a polygon, the remaining
	 * intersecting parts of individual faces will be collapsed into a single area.
	 * To preserve individual faces during intersection, use
	 * {@link #intersectMesh(PShape, PShape) intersectMesh()}.
	 * 
	 * @param a The first shape to be intersected.
	 * @param b The second shape to intersect with the first.
	 * @return A new shape representing the area of intersection between the two
	 *         input shapes. The resulting shape retains the style of the first
	 *         input shape, 'a'.
	 */
	public static PShape intersect(final PShape a, final PShape b) {
		Geometry shapeA = fromPShape(a);
		Geometry result = OverlayNG.overlay(shapeA, fromPShape(b), OverlayNG.INTERSECTION);
		result.setUserData(shapeA.getUserData()); // preserve shape style (if any)
		return toPShape(result);
	}

	/**
	 * Performs an intersection operation between a mesh-like shape (a polygonal
	 * coverage) and a polygonal area, while preserving the individual features of
	 * the mesh during the operation.
	 * <p>
	 * When a mesh-like shape or polygonal coverage is intersected with a polygon
	 * using the general {@link #intersect(PShape, PShape) intersect(a, b)} method,
	 * the result is a single polygon that consists of the unified or dissolved area
	 * of all intersecting mesh faces. Depending on the requirements of your task,
	 * this behavior may not be optimal. The current method addresses this issue by
	 * preserving the individual intersections of each mesh face.
	 * <p>
	 * This method performs faster than invoking the
	 * {@link #intersect(PShape, PShape) intersect(a, b)} method repeatedly for
	 * every face of a mesh-like shape <code>a</code> against an area
	 * <code>b</code>.
	 * 
	 * @param mesh A mesh-like GROUP shape that will be intersected with the
	 *             polygonal area.
	 * @param area A polygonal shape with which the mesh will be intersected.
	 * @return A GROUP shape where each child shape represents the union of one mesh
	 *         face with the area.
	 * @since 1.3.0
	 */
	public static PShape intersectMesh(final PShape mesh, final PShape area) {
		final Geometry g = fromPShape(area);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(g);

		List<Geometry> faces = PGS_Conversion.getChildren(mesh).parallelStream().map(s -> {
			final Geometry f = PGS_Conversion.fromPShape(s);
			if (cache.containsProperly(f)) {
				return f;
			} else {
				// preserve the fill etc of the PShape during intersection
				Geometry boundaryIntersect = OverlayNG.overlay(f, g, OverlayNG.INTERSECTION);
				boundaryIntersect.setUserData(f.getUserData());
				return boundaryIntersect;
			}
		}).collect(Collectors.toList());
		return PGS_Conversion.toPShape(faces);
	}

	/**
	 * Combines two shapes into a single new shape, representing the total area of
	 * both input shapes.
	 * <p>
	 * This method performs a geometric union operation, "gluing" the two input
	 * shapes together to create a new shape that encompasses the combined area of
	 * both inputs. If the input shapes overlap, the overlapping area is included
	 * only once in the resulting shape.
	 *
	 * @param a The first shape to be unified.
	 * @param b The second shape to be unified with the first.
	 * @return A new PShape representing the union of the two input shapes. The
	 *         resulting shape retains the style of the first input shape, 'a'.
	 * @see #union(PShape...) For union operations on multiple shapes.
	 */
	public static PShape union(final PShape a, final PShape b) {
		Geometry shapeA = fromPShape(a);
		Geometry result = OverlayNG.overlay(shapeA, fromPShape(b), OverlayNG.UNION);
		result.setUserData(shapeA.getUserData()); // preserve shape style (if any)
		return toPShape(result);
	}

	/**
	 * Performs a geometric union operation on a collection of shapes, merging them
	 * into a new shape that represents the total area of all the input shapes.
	 * Overlapping areas among the shapes are included only once in the resulting
	 * shape.
	 *
	 * @param shapes A list of PShapes to be unified.
	 * @return A new PShape object representing the union of the input shapes.
	 * @see #union(PShape, PShape) For a union operation on two shapes.
	 * @see #union(PShape...) For union operations on a variable number of shapes.
	 */
	public static PShape union(final Collection<PShape> shapes) {
		try {
			return toPShape(UnaryUnionOp.union(shapes.stream().map(s -> fromPShape(s)).toList()));
		} catch (Exception e) {
			return toPShape(UnaryUnionOp.union(shapes.stream().map(s -> GeometryFixer.fix(fromPShape(s))).toList()));
		}
	}

	/**
	 * Performs a geometric union operation on a variable number of shapes, merging
	 * them into a new shape that encompasses the total area of all input shapes.
	 * Overlapping areas among the shapes are included only once in the resulting
	 * shape.
	 * 
	 * @param shapes A variable number of PShape instances to be unified.
	 * @return A new PShape object representing the union of the input shapes.
	 * @see #union(PShape, PShape) For a union operation on two shapes.
	 * @see #union(List) For a union operation on a list of shapes.
	 */

	public static PShape union(PShape... shapes) {
		return union(Arrays.asList(shapes));
	}

	/**
	 * @see #unionMesh(PShape)
	 * @param faces collection of faces comprising a mesh
	 * @return A single PShape representing the boundary of the mesh.
	 */
	public static PShape unionMesh(final Collection<PShape> faces) {
		return unionMesh(PGS_Conversion.flatten(faces));
	}

	/**
	 * Merges a mesh-like PShape (i.e., a GROUP PShape whose children represent
	 * faces with shared edges) into a single shape that denotes the boundary of the
	 * entire mesh. This method is specifically optimized for meshes and
	 * significantly outperforms other methods of unioning the mesh faces.
	 * <p>
	 * The mesh can contain holes and these will be correctly processed and
	 * reflected in the final output.
	 *
	 * @param mesh A GROUP PShape, where each child shape forms part of a mesh,
	 *             meaning they join or overlap at edges.
	 * @return A single PShape representing the boundary of the mesh.
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
	 * Unifies a collection of mesh shapes without handling holes, providing a more
	 * faster approach than {@link #unionMesh(PShape)} if the input is known to have
	 * no holes.
	 * <p>
	 * This method calculates the set of unique edges belonging to the mesh, which
	 * is equivalent to the boundary, assuming a mesh without holes. It then
	 * determines a sequential/winding order for the vertices of the boundary.
	 * <p>
	 * Note: This method does not account for meshes with holes.
	 *
	 * @param mesh A collection of shapes representing a mesh.
	 * @return A new PShape representing the union of the mesh shapes.
	 * @deprecated This method is deprecated due to the lack of support for meshes
	 *             with holes.
	 */
	public static PShape unionMeshWithoutHoles(final Collection<PShape> mesh) {
		Map<PEdge, Integer> edges = new HashMap<>();

		final List<PEdge> allEdges;

		/*
		 * Compute set of unique edges belonging to the mesh (this set is equivalent to
		 * the boundary, assuming a holeless mesh).
		 */
		for (PShape child : mesh) {
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
	 *         shapes; the new shape will have the style of shape a.
	 * @see #simpleSubtract(PShape, PShape)
	 */
	public static PShape subtract(final PShape a, final PShape b) {
		Geometry shapeA = fromPShape(a);
		Geometry result = OverlayNG.overlay(shapeA, fromPShape(b), OverlayNG.DIFFERENCE);
		result.setUserData(shapeA.getUserData()); // preserve shape style (if any)
		return toPShape(result);
	}

	/**
	 * Subtracts <code>holes</code> from the <code>shell</code>, without geometric
	 * processing.
	 * <p>
	 * Rather than performing geometric overlay, this method simply appends the
	 * holes as contours to the shell. For this reason, it is faster than
	 * {@link #subtract(PShape, PShape) subtract()} but this method produces valid
	 * results only if <b>all holes lie inside the shell</b> and holes are <b>not
	 * nested</b>.
	 * 
	 * @param shell polygonal shape
	 * @param holes single polygon, or GROUP shape, whose children are holes that
	 *              lie within the shell
	 * @since 1.4.0
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
	 * Subtracts a polygonal area from a mesh-like shape or polygonal coverage,
	 * ensuring each individual face or feature of the mesh is preserved during the
	 * operation.
	 * <p>
	 * The behaviour of this method differs to {@link #subtract(PShape, PShape)
	 * subtract(a, b)} when <code>a</code> is a mesh; upon subtracting a polygon
	 * from a mesh, that method produces a single (unioned) polygon that represents
	 * the combined and dissolved area of all remaining mesh face parts. In
	 * contrast, this method preserves faces and how they are individually affected
	 * by the subtraction.
	 * <p>
	 * This method is more efficient than repeatedly calling
	 * {@link #subtract(PShape, PShape) subtract(a, b)} on each face of a mesh-like
	 * shape.
	 *
	 * @param mesh A GROUP PShape that represents a mesh-like shape.
	 * @param area A polygonal PShape from which the mesh shape is subtracted.
	 * @return A GROUP PShape, where each child shape is the result of the area
	 *         being subtracted from a single mesh face.
	 * @since 1.3.0
	 */
	public static PShape subtractMesh(PShape mesh, PShape area) {
		final Geometry g = fromPShape(area);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(g);

		List<Geometry> faces = PGS_Conversion.getChildren(mesh).parallelStream().map(s -> {
			final Geometry f = PGS_Conversion.fromPShape(s);
			if (cache.containsProperly(f)) {
				return null; // inside -- remove
			} else {
				if (cache.disjoint(f)) {
					return f; // outside -- keep
				}
				// preserve the fill etc of the PShape during subtraction
				Geometry boundarySubtract = OverlayNG.overlay(f, g, OverlayNG.DIFFERENCE);
				boundarySubtract.setUserData(f.getUserData());
				return boundarySubtract;
			}
		}).collect(Collectors.toList());

		return PGS_Conversion.toPShape(faces);
	}

	/**
	 * Calculates the symmetric difference between two shapes. The symmetric
	 * difference is the set of regions that exist in either of the two shapes but
	 * not in their intersection.
	 * 
	 * @param a The first shape.
	 * @param b The second shape.
	 * @return A new shape representing the symmetric difference between the two
	 *         input shapes; the new shape will have the style of shape a.
	 */
	public static PShape symDifference(PShape a, PShape b) {
		Geometry shapeA = fromPShape(a);
		Geometry result = OverlayNG.overlay(shapeA, fromPShape(b), OverlayNG.SYMDIFFERENCE);
		result.setUserData(shapeA.getUserData()); // preserve shape style (if any)
		return toPShape(result);
	}

	/**
	 * Calculates the complement (or inverse) of the provided shape within a
	 * rectangular boundary of specified width and height, anchored at (0, 0).
	 * <p>
	 * The resulting shape corresponds to the portion of the rectangle not covered
	 * by the input shape. The operation is essentially a subtraction of the input
	 * shape from the rectangle.
	 * 
	 * @param shape  The input shape for which the complement is to be determined.
	 * @param width  The width of the rectangular boundary.
	 * @param height The height of the rectangular boundary.
	 * @return A new PShape representing the inverse of the input shape within the
	 *         specified rectangular boundary.
	 */
	public static PShape complement(PShape shape, double width, double height) {
		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(4);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createRectangle().difference(fromPShape(shape)));
	}

}
