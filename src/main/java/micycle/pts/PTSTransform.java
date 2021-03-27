package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;

import org.locationtech.jts.algorithm.Distance;
import org.locationtech.jts.algorithm.distance.DistanceToPoint;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.operation.distance.DistanceOp;

import processing.core.PShape;
import processing.core.PVector;

/**
 * Geometric transformations for PShapes
 * 
 * @author Michael Carleton
 *
 */
public class PTSTransform {

	/**
	 * Scales the shape relative to its center
	 * 
	 * @param shape
	 * @param scale
	 * @return
	 */
	public static PShape scale(PShape shape, double scale) {
		Coordinate c = fromPShape(shape).getCentroid().getCoordinate();
		AffineTransformation t = AffineTransformation.scaleInstance(scale, scale, c.x, c.y);
		return toPShape(t.transform(fromPShape(shape)));
	}

	public static PShape scale(PShape shape, double scaleX, double scaleY) {
		AffineTransformation t = AffineTransformation.scaleInstance(scaleX, scaleY);
		return toPShape(t.transform(fromPShape(shape)));
	}

	/**
	 * Scales a shape (based on its centroid) so that it touches the boundary of
	 * another shape. Shape centroid must be outside the container.
	 * 
	 * @param shape     its centroid should be outside container
	 * @param container
	 * @param tolerance >=0
	 */
	public static PShape touchScale(PShape shape, PShape container, double tolerance) {
		tolerance = Math.max(tolerance, 1);
		Geometry scaleShape = fromPShape(shape);

		final Coordinate centroid = scaleShape.getCentroid().getCoordinate();

		double dist = 999999;
		final int maxIter = 75;
		int iter = 0;
		// NOTE uses DistanceOp.nearestPoints() which is n^2
		while (dist > tolerance && iter < maxIter) {
			Coordinate[] coords = DistanceOp.nearestPoints(scaleShape, fromPShape(container));
			dist = PTS.distance(coords[0], coords[1]);

			/**
			 * If dist == 0, then shape is either fully contained within the container or
			 * covers it. We attempt to first shrink the shape so that no longer covers the
			 * container. If dist remains zero after repeated shrinking we conclude the
			 * shape is inside the container.
			 */
			if (dist == 0) {
				int z = 7;
				while (z > 0) {
					AffineTransformation t = AffineTransformation.scaleInstance(0.5, 0.5, centroid.x, centroid.y);
					scaleShape = t.transform(scaleShape);
					coords = DistanceOp.nearestPoints(scaleShape, fromPShape(container));
					dist = PTS.distance(coords[0], coords[1]);
					if (dist > 0) {
						break;
					}
					z--;
				}
				if (dist == 0) { // still 0? probably contained inside, so just return shape
					return shape;
				}
			}
			double d1 = PTS.distance(centroid, coords[0]);
			double d2 = PTS.distance(centroid, coords[1]);
			AffineTransformation t = AffineTransformation.scaleInstance(d2 / d1, d2 / d1, centroid.x, centroid.y);
			scaleShape = t.transform(scaleShape);
			iter++;
		}
		return toPShape(scaleShape);
	}

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
	public static PShape flipHorizontal(PShape shape, double y) {
		AffineTransformation t = AffineTransformation.reflectionInstance(-1, y, 1, y);
		return toPShape(t.transform(fromPShape(shape)));
	}

	public static PShape flipVertical(PShape shape) {
		Geometry g = fromPShape(shape);
		Point c = g.getCentroid();
		AffineTransformation t = AffineTransformation.reflectionInstance(c.getX(), -1, c.getX(), 1);
		return toPShape(t.transform(g));
	}

	public static PShape flipVertical(PShape shape, double x) {
		AffineTransformation t = AffineTransformation.reflectionInstance(x, -1, x, 1);
		return toPShape(t.transform(fromPShape(shape)));
	}

	/**
	 * Objects are sheared around their relative position to the origin.
	 * 
	 * @param shape
	 * @param angleX radians
	 * @param angleY radians
	 * @return
	 */
	public static PShape shear(PShape shape, double angleX, double angleY) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.shearInstance(angleX, angleY);
		return toPShape(t.transform(g));
	}

	public static PShape translate(PShape shape, double x, double y) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.translationInstance(x, y);
		return toPShape(t.transform(g));
	}

	/**
	 * Translates the shape such that its centroid is equivalent to the given
	 * coordinates.
	 * 
	 * @param shape
	 * @param x     target centroid X
	 * @param y     target centroid Y
	 * @return translated shape
	 */
	public static PShape translateTo(PShape shape, double x, double y) {
		Geometry g = fromPShape(shape);
		Point c = g.getCentroid();
		double translateX = x - c.getX();
		double translateY = y - c.getY();
		AffineTransformation t = AffineTransformation.translationInstance(translateX, translateY);
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
	public static PShape rotate(PShape shape, PVector point, double angle) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.rotationInstance(angle, point.x, point.y);
		return toPShape(t.transform(g));
	}

	/**
	 * Rotate a shape around its centroid.
	 * 
	 * @param shape
	 * @param angle the rotation angle, in radians
	 * @return
	 */
	public static PShape rotateAroundCenter(PShape shape, double angle) {
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
	public static PShape homotheticTransformation(PShape shape, PVector center, double scaleX, double scaleY) {
		Polygon geom = (Polygon) fromPShape(shape);

		// external contour
		Coordinate[] coord = geom.getExteriorRing().getCoordinates();
		Coordinate[] coord_ = new Coordinate[coord.length];
		for (int i = 0; i < coord.length; i++)
			coord_[i] = new Coordinate(center.x + scaleX * (coord[i].x - center.x),
					center.y + scaleY * (coord[i].y - center.y));
		LinearRing lr = geom.getFactory().createLinearRing(coord_);

		// holes
		LinearRing[] holes = new LinearRing[geom.getNumInteriorRing()];
		for (int j = 0; j < geom.getNumInteriorRing(); j++) {
			Coordinate[] hole_coord = geom.getInteriorRingN(j).getCoordinates();
			Coordinate[] hole_coord_ = new Coordinate[hole_coord.length];
			for (int i = 0; i < hole_coord.length; i++)
				hole_coord_[i] = new Coordinate(center.x + scaleY * (hole_coord[i].x - center.x),
						center.y + scaleY * (hole_coord[i].y - center.y));
			holes[j] = geom.getFactory().createLinearRing(hole_coord_);
		}
		return toPShape(PTS.GEOM_FACTORY.createPolygon(lr, holes));
	}

}