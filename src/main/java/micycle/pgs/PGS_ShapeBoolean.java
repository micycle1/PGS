package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.util.GeometricShapeFactory;

import processing.core.PShape;

/**
 * Boolean set-operations for 2D shapes.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_ShapeBoolean {

	private PGS_ShapeBoolean() {
	}

	/**
	 * Intersect creates a boolean group whose shape consists only of the
	 * overlapping parts of the given shapes.
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

}
