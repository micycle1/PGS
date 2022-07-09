package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;

import processing.core.PShape;
import processing.core.PVector;

/**
 * Various geometric and affine transformations for PShapes that affect vertex
 * coordinates.
 * <p>
 * Notably, these transformation methods affect the vertex coordinates of
 * PShapes, unlike Processing's transform methods that affect the affine matrix
 * of shapes only (and thereby leave vertex coordinates in-tact).
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Transformation {

	private PGS_Transformation() {
	}

	/**
	 * Scales the shape relative to its center point.
	 * 
	 * @param shape
	 * @param scale X and Y axis scale factor
	 * @return
	 */
	public static PShape scale(PShape shape, double scale) {
		Geometry g = fromPShape(shape);
		Coordinate c = g.getCentroid().getCoordinate();
		AffineTransformation t = AffineTransformation.scaleInstance(scale, scale, c.x, c.y);
		return toPShape(t.transform(g));
	}

	/**
	 * Scales the shape relative to its center point.
	 * 
	 * @param shape
	 * @param scaleX X-axis scale factor
	 * @param scaleY Y-axis scale factor
	 * @return
	 */
	public static PShape scale(PShape shape, double scaleX, double scaleY) {
		Geometry g = fromPShape(shape);
		Point c = g.getCentroid();
		AffineTransformation t = AffineTransformation.scaleInstance(scaleX, scaleY, c.getX(), c.getY());
		return toPShape(t.transform(g));
	}

	/**
	 * Scales the shape relative to the origin (0,0).
	 * 
	 * @param shape
	 * @param scale cale factor
	 * @return
	 * @since 1.2.1
	 */
	public static PShape originScale(PShape shape, double scale) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.scaleInstance(scale, scale, 0, 0);
		return toPShape(t.transform(g));
	}

	/**
	 * Resizes a shape (based on its envelope) to the given dimensions, relative to
	 * its center point.
	 * 
	 * @param shape
	 * @param targetWidth  width of the output copy
	 * @param targetHeight height of the output copy
	 * @return resized copy of input shape
	 */
	public static PShape resize(PShape shape, double targetWidth, double targetHeight) {
		targetWidth = Math.max(targetWidth, 0.001);
		targetHeight = Math.max(targetHeight, 0.001);
		Geometry geometry = fromPShape(shape);
		Envelope e = geometry.getEnvelopeInternal();
		Point c = geometry.getCentroid();

		AffineTransformation t = AffineTransformation.scaleInstance(targetWidth / e.getWidth(), targetHeight / e.getHeight(), c.getX(),
				c.getY());
		return translateToOrigin(toPShape(t.transform(geometry)));
	}

	/**
	 * Resizes a shape (based on its envelope) to the given width relative to its
	 * center point; the height is resized accordingly to maintain the shape's
	 * aspect ratio.
	 * 
	 * @param shape       the shape to resize
	 * @param targetWidth width of the output
	 * @return resized copy of input shape
	 * @since 1.2.1
	 * @see #resizeByHeight(PShape, double)
	 */
	public static PShape resizeByWidth(PShape shape, double targetWidth) {
		targetWidth = Math.max(targetWidth, 1e-5);

		Geometry geometry = fromPShape(shape);
		Envelope e = geometry.getEnvelopeInternal();
		Point c = geometry.getCentroid();

		AffineTransformation t = AffineTransformation.scaleInstance(targetWidth / e.getWidth(), targetWidth / e.getWidth(), c.getX(),
				c.getY());
		return toPShape(t.transform(geometry));
	}

	/**
	 * Resizes a shape (based on its envelope) to the given height relative to its
	 * center point; the width is resized accordingly to maintain the shape's aspect
	 * ratio.
	 * 
	 * @param shape        the shape to resize
	 * @param targetHeight height of the output
	 * @return resized copy of input shape
	 * @since 1.2.1
	 * @see #resizeByWidth(PShape, double)
	 */
	public static PShape resizeByHeight(PShape shape, double targetHeight) {
		targetHeight = Math.max(targetHeight, 1e-5);

		Geometry geometry = fromPShape(shape);
		Envelope e = geometry.getEnvelopeInternal();
		Point c = geometry.getCentroid();

		AffineTransformation t = AffineTransformation.scaleInstance(targetHeight / e.getHeight(), targetHeight / e.getHeight(), c.getX(),
				c.getY());
		return toPShape(t.transform(geometry));
	}

	/**
	 * Scales a shape (based on its centroid) so that it touches the boundary of
	 * another shape. The scaling shape's centroid must lie outside the other shape.
	 * 
	 * @param shape     the shape to scale. its centroid should be outside container
	 * @param boundary
	 * @param tolerance >=0
	 */
	public static PShape touchScale(PShape shape, PShape boundary, double tolerance) {
		tolerance = Math.max(tolerance, 0.01);
		final IndexedFacetDistance index = new IndexedFacetDistance(fromPShape(boundary));
		Geometry scaleShape = fromPShape(shape);

		final Coordinate centroid = scaleShape.getCentroid().getCoordinate();

		double dist = Double.MAX_VALUE;
		final int maxIter = 75;
		int iter = 0;

		while (dist > tolerance && iter < maxIter) {
			Coordinate[] coords = index.nearestPoints(scaleShape);
			dist = coords[0].distance(coords[1]);

			/*
			 * When dist == 0, the shape is either fully contained within the boundary or
			 * partially covers it. We attempt to first shrink the shape so that no longer
			 * covers the boundary. If dist remains zero after repeated shrinking we
			 * conclude the shape is properly inside the boundary.
			 */
			if (dist == 0) {
				int z = 0;
				while (z++ < 7) {
					AffineTransformation t = AffineTransformation.scaleInstance(0.5, 0.5, centroid.x, centroid.y);
					scaleShape = t.transform(scaleShape);
					dist = index.distance(scaleShape);
					if (dist > 0) {
						break;
					}
				}
				if (dist == 0) {
					/*
					 * If dist still == 0, then shape's centroid lies on the perimeter of the
					 * boundary.
					 */
					return toPShape(scaleShape);
				}
			}
			double d1 = centroid.distance(coords[1]);
			double d2 = centroid.distance(coords[0]);
			AffineTransformation t = AffineTransformation.scaleInstance(d2 / d1, d2 / d1, centroid.x, centroid.y);
			scaleShape = t.transform(scaleShape);
			iter++;
		}
		return toPShape(scaleShape);
	}

	/**
	 * Translates a shape by the given coordinates.
	 * 
	 * @param shape
	 * @param x
	 * @param y
	 * @return translated copy
	 */
	public static PShape translate(PShape shape, double x, double y) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.translationInstance(x, y);
		return toPShape(t.transform(g));
	}

	/**
	 * Translates a shape such that its centroid is equivalent to the given
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
	 * Translates a shape such that the top-left corner of its bounding box is at
	 * (0, 0) (in Processing coordinates).
	 * 
	 * @param shape
	 * @return translated copy of input
	 */
	public static PShape translateToOrigin(PShape shape) {
		final Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();
		AffineTransformation t = AffineTransformation.translationInstance(-e.getMinX(), -e.getMinY());
		return toPShape(t.transform(g));
	}

	/**
	 * Calculate a Homothetic Transformation of the shape.
	 * 
	 * @param shape  shape input
	 * @param center coordinate of the center/origin position of the operation
	 * @param scaleX X scale factor
	 * @param scaleY Y scale factor
	 */
	public static PShape homotheticTransformation(PShape shape, PVector center, double scaleX, double scaleY) {
		Polygon geom = (Polygon) fromPShape(shape);

		// external contour
		Coordinate[] coord = geom.getExteriorRing().getCoordinates();
		Coordinate[] coord_ = new Coordinate[coord.length];
		for (int i = 0; i < coord.length; i++)
			coord_[i] = new Coordinate(center.x + scaleX * (coord[i].x - center.x), center.y + scaleY * (coord[i].y - center.y));
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
		return toPShape(PGS.GEOM_FACTORY.createPolygon(lr, holes));
	}

	/**
	 * Rotates a shape around a given point.
	 * 
	 * @param shape
	 * @param point
	 * @param angle
	 * @return
	 * @see #rotateAroundCenter(PShape, double)
	 */
	public static PShape rotate(PShape shape, PVector point, double angle) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.rotationInstance(angle, point.x, point.y);
		return toPShape(t.transform(g));
	}

	/**
	 * Rotates a shape around its centroid.
	 * 
	 * @param shape
	 * @param angle the rotation angle, in radians
	 * @return
	 * @see #rotate(PShape, PVector, double)
	 */
	public static PShape rotateAroundCenter(PShape shape, double angle) {
		Geometry g = fromPShape(shape);
		Point center = g.getCentroid();
		AffineTransformation t = AffineTransformation.rotationInstance(angle, center.getX(), center.getY());
		return toPShape(t.transform(g));
	}

	/**
	 * Flips the shape horizontally based on its centre point.
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
	 * Flips the shape horizontally based on a line given by its Y location.
	 * 
	 * @param shape
	 * @param y
	 * @return
	 */
	public static PShape flipHorizontal(PShape shape, double y) {
		AffineTransformation t = AffineTransformation.reflectionInstance(-1, y, 1, y);
		return toPShape(t.transform(fromPShape(shape)));
	}

	/**
	 * Flips the shape vertically based on its centre point.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape flipVertical(PShape shape) {
		Geometry g = fromPShape(shape);
		Point c = g.getCentroid();
		AffineTransformation t = AffineTransformation.reflectionInstance(c.getX(), -1, c.getX(), 1);
		return toPShape(t.transform(g));
	}

	/**
	 * Flips the shape vertically based on a line given by its X location.
	 * 
	 * @param shape
	 * @param x
	 * @return
	 */
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

}
