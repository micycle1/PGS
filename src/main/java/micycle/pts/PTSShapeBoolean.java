package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import static micycle.pts.PTS.getPShapeFillColor;

import java.util.ArrayList;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pts.color.Blending;
import processing.core.PShape;

/**
 * Provides Boolean Set-Operations on shapes.
 * 
 * @author Michael Carleton
 *
 */
public class PTSShapeBoolean {

	private PTSShapeBoolean() {

	}

	/**
	 * Intersect creates a boolean group whose shape consists only of the
	 * overlapping parts of the given shapes.
	 * 
	 * @return A∩B
	 */
	public static PShape intersect(PShape a, PShape b) {
		PShape out = toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.INTERSECTION));
		Conversion.setAllFillColor(out, Blending.add(getPShapeFillColor(a), getPShapeFillColor(b)));
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

	public static PShape union(PShape... shapes) {
		// same as flatten?
		ArrayList<Geometry> geoms = new ArrayList<>();
		for (int i = 0; i < shapes.length; i++) {
			geoms.add(fromPShape(shapes[i]));
		}
		return toPShape(UnaryUnionOp.union(geoms));
	}

	/**
	 * Subtract is the opposite of Union. Subtract removes the area of a shape b
	 * from the base shape. A.k.a difference.
	 * 
	 * @return shape A - shape B
	 */
	public static PShape subtract(PShape a, PShape b) {
		return toPShape(OverlayNG.overlay(fromPShape(a), fromPShape(b), OverlayNG.DIFFERENCE));
	}

	/**
	 * Symmetric difference: Compute the parts that the shapes do not have in
	 * common.
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
	 * @param width  width of the plane
	 * @param height height of the plane
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
