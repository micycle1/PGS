package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
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

	public static PShape translate(PShape shape, float x, float y) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.translationInstance(x, y);
		return toPShape(t.transform(g));
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

	/**
	 * Calculate a Homothetic Transformation of the shape.
	 * 
	 * @param shape  shape input
	 * @param x0     X position of the center of the operation
	 * @param y0     Y position of the center of the operation
	 * @param scaleX X scale factor
	 * @param scaleY Y scale factor
	 */
	public static PShape homotheticTransformation(PShape shape, double x0, double y0, double scaleX, double scaleY) {
		Polygon geom = (Polygon) fromPShape(shape);

		// external contour
		Coordinate[] coord = geom.getExteriorRing().getCoordinates();
		Coordinate[] coord_ = new Coordinate[coord.length];
		for (int i = 0; i < coord.length; i++)
			coord_[i] = new Coordinate(x0 + scaleX * (coord[i].x - x0), y0 + scaleY * (coord[i].y - y0));
		LinearRing lr = geom.getFactory().createLinearRing(coord_);

		// the holes
		LinearRing[] holes = new LinearRing[geom.getNumInteriorRing()];
		for (int j = 0; j < geom.getNumInteriorRing(); j++) {
			Coordinate[] hole_coord = geom.getInteriorRingN(j).getCoordinates();
			Coordinate[] hole_coord_ = new Coordinate[hole_coord.length];
			for (int i = 0; i < hole_coord.length; i++)
				hole_coord_[i] = new Coordinate(x0 + scaleY * (hole_coord[i].x - x0),
						y0 + scaleY * (hole_coord[i].y - y0));
			holes[j] = geom.getFactory().createLinearRing(hole_coord_);
		}
		return toPShape(PTS.GEOM_FACTORY.createPolygon(lr, holes));
	}

}
