package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.operation.linemerge.LineMerger;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pgs.PGS.PEdge;
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
	 * 
	 * @return A∩B
	 */
	public static PShape intersect(PShape a, PShape b) {
		PShape out = toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.INTERSECTION));
//		PGS_Conversion.setAllFillColor(out, Blending.add(getPShapeFillColor(a), getPShapeFillColor(b)));
		return out;
	}

	/**
	 * "Glues" shapes together so they become a single combined shape with the sum
	 * of its areas.
	 * 
	 * @return A∪B
	 * @see #union(PShape...)
	 */
	public static PShape union(PShape a, PShape b) {
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
	public static PShape union(List<PShape> shapes) {
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
	 * Unions/flattens/merges a mesh-like PShape (that is, a GROUP PShape whose
	 * children represent faces that share edges) into a single shape that
	 * represents the boundary of the mesh. This method is optimised for meshes, and
	 * is accordingly much faster than unioning the mesh faces together using other
	 * methods.
	 * 
	 * @param mesh a GROUP pshape whose children shapes form a mesh (join/overlap at
	 *             edges)
	 * @return
	 * @since 1.2.0
	 */
	public static PShape unionMesh(PShape mesh) {

		if (mesh.getChildCount() == 0 || mesh.getKind() != PConstants.GROUP) {
			System.err.println("unionMesh Error: Input shape was not a GROUP shape, or had 0 children.");
			return new PShape();
		}

		final Set<PEdge> allEdges = new HashSet<>(mesh.getChildCount() * 3);
		final Set<PEdge> duplicateEdges = new HashSet<>(allEdges.size());

		/*
		 * Compute set of unique edges belonging to the mesh (this set is equivalent to
		 * the boundary); ignore edges if they are seen more than once.
		 */
		for (PShape child : mesh.getChildren()) {
			for (int i = 0; i < child.getVertexCount(); i++) {
				final PVector a = child.getVertex(i);
				final PVector b = child.getVertex((i + 1) % child.getVertexCount());
				if (!a.equals(b)) {
					PEdge edge = new PEdge(a, b);
					if (!allEdges.add(edge)) {
						duplicateEdges.add(edge);
					}
				}
			}
		}

		allEdges.removeAll(duplicateEdges); // allEdges now contains boundary edges only
		if (allEdges.isEmpty()) {
			return new PShape();
		}

		/*
		 * Now get the vertices belonging to boundary edges in sequential/winding order.
		 */
		final List<PVector> orderedVertices = computePolygonVertices(allEdges);
		return PGS_Conversion.fromPVector(orderedVertices);
	}

	/**
	 * Subtract is the opposite of Union. Subtract removes the area of a shape b
	 * from the base shape a. A.k.a "difference".
	 * 
	 * @return shape A - shape B
	 */
	public static PShape subtract(PShape a, PShape b) {
		return toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.DIFFERENCE));
	}

	/**
	 * Computes the parts that the shapes do not have in common.
	 * 
	 * @return A∪B - A∩B
	 */
	public static PShape symDifference(PShape a, PShape b) {
		return toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.SYMDIFFERENCE));
	}

	/**
	 * Computes the shape's complement (or inverse) against a plane.
	 * 
	 * @param shape
	 * @param width  width of the rectangle plane to subtract shape from
	 * @param height height of the rectangle plane to subtract shape from
	 * @return
	 */
	public static PShape complement(PShape shape, double width, double height) {
		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(4);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createRectangle().difference(fromPShape(shape)));
	}

	/**
	 * Computes an <b>ordered</b> list of <b>vertices</b> that make up the boundary
	 * of a polygon from an <b>unordered</b> collection of <b>edges</b>. The
	 * underlying approach is around ~10x faster than JTS .buffer(0) and ~3x faster
	 * than {@link LineMerger}.
	 * <p>
	 * For now, this method does not properly support multi-shapes.
	 * <p>
	 * Notably, unlike {@link LineMerger} this approach does not merge successive
	 * boundary segments that together form a straight line into a single longer
	 * segment.
	 * 
	 * @param edges unordered/random collection of edges (containing no duplicates),
	 *              that together constitute the boundary of a single polygon / a
	 *              closed ring
	 * @return sequential list of vertices belonging to the polygon that follow some
	 *         constant winding (may wind clockwise or anti-clockwise)
	 */
	private static List<PVector> computePolygonVertices(Collection<PEdge> edges) {
		final HashMap<PVector, HashSet<PEdge>> vertexEdges = new HashMap<>(); // map of vertex to the 2 edges that share it

		/*
		 * Build up map of vertex->edge to later find edges sharing a given vertex in
		 * O(1). When the input is valid (edges form a closed loop) every vertex is
		 * shared by 2 edges.
		 */
		for (PEdge e : edges) {
			if (vertexEdges.containsKey(e.a)) {
				vertexEdges.get(e.a).add(e);
			} else {
				HashSet<PEdge> h = new HashSet<>();
				h.add(e);
				vertexEdges.put(e.a, h);
			}
			if (vertexEdges.containsKey(e.b)) {
				vertexEdges.get(e.b).add(e);
			} else {
				HashSet<PEdge> h = new HashSet<>();
				h.add(e);
				vertexEdges.put(e.b, h);
			}
		}

		ArrayList<PVector> vertices = new ArrayList<>(edges.size() + 1); // boundary vertices

		// begin by choosing a random edge
		final PEdge startingEdge = edges.iterator().next();
		vertices.add(startingEdge.a);
		vertices.add(startingEdge.b);
		vertexEdges.get(startingEdge.a).remove(startingEdge);
		vertexEdges.get(startingEdge.b).remove(startingEdge);

		while (vertices.size() < edges.size()) {
			final PVector lastVertex = vertices.get(vertices.size() - 1);
			HashSet<PEdge> connectedEdges = vertexEdges.get(lastVertex);

			if (connectedEdges.isEmpty()) {
				/*
				 * This will be hit if the input is malformed (contains multiple disjoint shapes
				 * for example), and break when the first loop is closed. On valid inputs the
				 * while loop will break before this statement can be hit.
				 */
				break;
			}

			final PEdge nextEdge = connectedEdges.iterator().next();
			if (nextEdge.a.equals(lastVertex)) {
				vertices.add(nextEdge.b);
				vertexEdges.get(nextEdge.b).remove(nextEdge);
			} else {
				vertices.add(nextEdge.a);
				vertexEdges.get(nextEdge.a).remove(nextEdge);
			}
			connectedEdges.remove(nextEdge); // remove this edge from vertex mapping
			if (connectedEdges.isEmpty()) {
				vertexEdges.remove(lastVertex); // have used both edges connected to this vertex -- now remove!
			}
		}

		return vertices;
	}

}
