package micycle.pts;

import static micycle.pts.PTS.fromPShape;
import static micycle.pts.PTS.toPShape;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.util.AffineTransformation;

import processing.core.PShape;
import processing.core.PVector;

/**
 * Shape transformations
 * 
 * @author MCarleton
 *
 */
public class Transform {

	// TODO shear+scale

	/**
	 * Flip based on the centre point TODO move near rotate
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape flipHorizontal(PShape shape) {
		Geometry g = fromPShape(shape);
		Point c = g.getCentroid();
		AffineTransformation t = AffineTransformation.reflectionInstance(-1, c.getY(), 1, c.getY());
		return toPShape(t.transform(g));
	}

	/**
	 * Flip based on a line given by its Y location
	 * 
	 * @param shape
	 * @param y
	 * @return
	 */
	public static PShape flipHorizontal(PShape shape, float y) {
		AffineTransformation t = AffineTransformation.reflectionInstance(-1, y, 1, y);
		return toPShape(t.transform(fromPShape(shape)));
	}

	public static PShape flipVertical(PShape shape) {
		Geometry g = fromPShape(shape);
		Point c = g.getCentroid();
		AffineTransformation t = AffineTransformation.reflectionInstance(c.getX(), -1, c.getX(), 1);
		return toPShape(t.transform(g));
	}

	public static PShape flipVertical(PShape shape, float x) {
		AffineTransformation t = AffineTransformation.reflectionInstance(x, -1, x, 1);
		return toPShape(t.transform(fromPShape(shape)));
	}

	/**
	 * Rotate a shape around a given point.
	 * 
	 * @param shape
	 * @param point
	 * @param angle
	 * @return
	 */
	public static PShape rotate(PShape shape, PVector point, float angle) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.rotationInstance(angle, point.x, point.y);
		return toPShape(t.transform(g));
	}

	/**
	 * Rotate a shape around its centroid.
	 * 
	 * @param shape
	 * @param angle
	 * @return
	 */
	public static PShape rotateAroundCenter(PShape shape, float angle) {
		Geometry g = fromPShape(shape);
		Point center = g.getCentroid();
		AffineTransformation t = AffineTransformation.rotationInstance(angle, center.getX(), center.getY());
		return toPShape(t.transform(g));
	}

}
