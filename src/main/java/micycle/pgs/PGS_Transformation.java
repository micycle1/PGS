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

import micycle.pgs.commons.ProcrustesAlignment;
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
	 * Scales the dimensions of the shape by a scaling factor relative to its
	 * centroid.
	 * 
	 * @param shape
	 * @param scale X and Y axis scale factor
	 */
	public static PShape scale(PShape shape, double scale) {
		Geometry g = fromPShape(shape);
		Coordinate c = g.getCentroid().getCoordinate();
		AffineTransformation t = AffineTransformation.scaleInstance(scale, scale, c.x, c.y);
		return toPShape(t.transform(g));
	}

	/**
	 * Scales the shape relative to its centroid.
	 * 
	 * @param shape
	 * @param scaleX X-axis scale factor
	 * @param scaleY Y-axis scale factor
	 */
	public static PShape scale(PShape shape, double scaleX, double scaleY) {
		Geometry g = fromPShape(shape);
		Point c = g.getCentroid();
		AffineTransformation t = AffineTransformation.scaleInstance(scaleX, scaleY, c.getX(), c.getY());
		return toPShape(t.transform(g));
	}

	/**
	 * Scale a shape around a point.
	 * 
	 * @since 2.0
	 */
	public static PShape scale(PShape shape, double scaleX, double scaleY, PVector point) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.scaleInstance(scaleX, scaleY, point.x, point.y);
		return toPShape(t.transform(g));
	}

	/**
	 * Scale a shape around a point.
	 * 
	 * @since 2.0
	 */
	public static PShape scale(PShape shape, double scale, double x, double y) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.scaleInstance(scale, scale, x, y);
		return toPShape(t.transform(g));
	}

	/**
	 * Scales a shape relative to the origin (0,0).
	 * 
	 * @param shape
	 * @param scale scale factor
	 * @since 1.3.0
	 */
	public static PShape originScale(PShape shape, double scale) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.scaleInstance(scale, scale, 0, 0);
		return toPShape(t.transform(g));
	}

	/**
	 * Scales the area of a given shape by a specified scale factor. The shape is
	 * scaled relative to its centroid.
	 * 
	 * @param shape The PShape to be scaled.
	 * @param scale The scale factor by which the area of the shape should be
	 *              scaled.
	 * @return A new PShape representing the scaled shape.
	 * @since 1.4.0
	 */
	public static PShape scaleArea(PShape shape, double scale) {
		Geometry geometry = fromPShape(shape);
		double scalingFactor = Math.sqrt(scale);
		Coordinate c = geometry.getCentroid().getCoordinate();
		AffineTransformation t = AffineTransformation.scaleInstance(scalingFactor, scalingFactor, c.x, c.y);
		return toPShape(t.transform(geometry));
	}

	/**
	 * Scales the given PShape to the target area, relative to its centroid.
	 * 
	 * @param shape      The PShape to be scaled.
	 * @param targetArea The target area for the shape.
	 * @return The scaled PShape (now having an area of <code>targetArea</code>).
	 * @since 1.4.0
	 */
	public static PShape scaleAreaTo(PShape shape, double targetArea) {
		Geometry geometry = fromPShape(shape);
		double area = geometry.getArea();
		double scalingFactor = Math.sqrt(targetArea / area);
		Coordinate c = geometry.getCentroid().getCoordinate();
		AffineTransformation t = AffineTransformation.scaleInstance(scalingFactor, scalingFactor, c.x, c.y);
		return toPShape(t.transform(geometry));
	}

	/**
	 * Resizes a shape (based on its envelope) to the given dimensions, relative to
	 * its centroid.
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

		AffineTransformation t = AffineTransformation.scaleInstance(targetWidth / e.getWidth(), targetHeight / e.getHeight(), c.getX(), c.getY());
		return translateToOrigin(toPShape(t.transform(geometry)));
	}

	/**
	 * Resizes a shape (based on its envelope) to the given width relative to its
	 * centroid; the height is resized accordingly to maintain the shape's aspect
	 * ratio.
	 * 
	 * @param shape       the shape to resize
	 * @param targetWidth width of the output
	 * @return resized copy of input shape
	 * @since 1.3.0
	 * @see #resizeByHeight(PShape, double)
	 * @see #resizeByMajorAxis(PShape, double)
	 */
	public static PShape resizeByWidth(PShape shape, double targetWidth) {
		targetWidth = Math.max(targetWidth, 1e-5);

		Geometry geometry = fromPShape(shape);
		Envelope e = geometry.getEnvelopeInternal();
		Point c = geometry.getCentroid();

		AffineTransformation t = AffineTransformation.scaleInstance(targetWidth / e.getWidth(), targetWidth / e.getWidth(), c.getX(), c.getY());
		return toPShape(t.transform(geometry));
	}

	/**
	 * Resizes a shape (based on its envelope) to the given height relative to its
	 * centroid; the width is resized accordingly to maintain the shape's aspect
	 * ratio.
	 * 
	 * @param shape        the shape to resize
	 * @param targetHeight height of the output
	 * @return resized copy of input shape
	 * @since 1.3.0 resized copy of input shape
	 * @see #resizeByWidth(PShape, double)
	 * @see #resizeByMajorAxis(PShape, double)
	 */
	public static PShape resizeByHeight(PShape shape, double targetHeight) {
		targetHeight = Math.max(targetHeight, 1e-5);

		Geometry geometry = fromPShape(shape);
		Envelope e = geometry.getEnvelopeInternal();
		Point c = geometry.getCentroid();

		AffineTransformation t = AffineTransformation.scaleInstance(targetHeight / e.getHeight(), targetHeight / e.getHeight(), c.getX(), c.getY());
		return toPShape(t.transform(geometry));
	}

	/**
	 * Resizes a shape (based on the longest axis of its envelope) to the given size
	 * relative to its centroid.
	 * <p>
	 * For example, if the shape's width is larger than its height, the width is set
	 * to <code>targetSize</code> and the height is resized to maintain the shape's
	 * aspect ratio.
	 * 
	 * @param shape      the shape to resize
	 * @param targetSize the new length of its longest axis
	 * @return resized copy of input shape
	 * @since 1.3.0
	 * @see #resizeByWidth(PShape, double)
	 * @see #resizeByHeight(PShape, double)
	 */
	public static PShape resizeByMajorAxis(PShape shape, double targetSize) {
		Envelope e = fromPShape(shape).getEnvelopeInternal();
		if (e.getWidth() > e.getHeight()) {
			return resizeByWidth(shape, targetSize);
		} else {
			return resizeByHeight(shape, targetSize);
		}
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
	 * @param shape shape to translate
	 * @param x     the value to translate by in the x direction
	 * @param y     the value to translate by in the y direction
	 * @return translated copy of input
	 */
	public static PShape translate(PShape shape, double x, double y) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.translationInstance(x, y);
		return toPShape(t.transform(g));
	}

	/**
	 * Translates a shape such that its <b>centroid</b> is equivalent to the given
	 * coordinates.
	 * 
	 * @param shape shape to translate
	 * @param x     target centroid X
	 * @param y     target centroid Y
	 * @return translated shape
	 * @deprecated Use {@link #translateCentroidTo(PShape, double, double)
	 *             translateCentroidTo()} instead.
	 */
	@Deprecated
	public static PShape translateTo(PShape shape, double x, double y) {
		return translateCentroidTo(shape, x, y);
	}

	/**
	 * Translates a shape such that its <b>centroid</b> aligns with the specified
	 * (x, y) coordinates. The centroid of a polygon is its "center of mass" or
	 * geometric center.
	 *
	 * @param shape The PShape instance to be translated.
	 * @param x     The horizontal coordinate to which the centroid of the shape's
	 *              bounding polygon should be aligned. Measured in pixels from the
	 *              left of the container.
	 * @param y     The vertical coordinate to which the centroid of the shape's
	 *              bounding polygon should be aligned. Measured in pixels from the
	 *              top of the container.
	 * @return A new PShape instance that is a translation of the input shape such
	 *         that the centroid of its bounding polygon aligns with the specified
	 *         coordinates.
	 * @since 1.3.0 (superceeds {@link #translateTo(PShape, double, double)
	 *        translateTo()}
	 */
	public static PShape translateCentroidTo(PShape shape, double x, double y) {
		Geometry g = fromPShape(shape);
		if (g.getNumPoints() == 0) {
			return shape;
		}
		Point c = g.getCentroid();
		double translateX = x - c.getX();
		double translateY = y - c.getY();
		AffineTransformation t = AffineTransformation.translationInstance(translateX, translateY);
		return toPShape(t.transform(g));
	}

	/**
	 * Translates a shape such that the <b>centroid</b> of its <b>bounding box</b>
	 * is equivalent to the given coordinates. This method differs a little from
	 * {@link #translateTo(PShape, double, double) translateTo()} because that uses
	 * center of shape's "mass", whereas this is visual center (however both are
	 * usually similar, if not identical).
	 * 
	 * @param shape shape to translate
	 * @param x     the x-coordinate of new position of the shape's bounding box
	 *              centroid point
	 * @param y     the y-coordinate of new position of the shape's bounding box
	 *              centroid point
	 * @return translated shape
	 * @since 1.3.0
	 */
	public static PShape translateEnvelopeTo(PShape shape, double x, double y) {
		Geometry g = fromPShape(shape);
		if (g.getNumPoints() == 0) {
			return shape;
		}
		Point c = g.getEnvelope().getCentroid();
		double translateX = x - c.getX();
		double translateY = y - c.getY();
		AffineTransformation t = AffineTransformation.translationInstance(translateX, translateY);
		return toPShape(t.transform(g));
	}

	/**
	 * Translates the given shape such that the upper-left corner of its bounding
	 * box aligns with the specified (x, y) coordinates.
	 *
	 * @param shape The PShape instance to be translated.
	 * @param x     The horizontal coordinate to which the upper-left corner of the
	 *              shape's bounding box is to be aligned. Measured in pixels from
	 *              the left of the container.
	 * @param y     The vertical coordinate to which the upper-left corner of the
	 *              shape's bounding box is to be aligned. Measured in pixels from
	 *              the top of the container.
	 * @return A new PShape instance that is a translation of the input shape such
	 *         that the upper-left corner of its bounding box aligns with the
	 *         specified coordinates.
	 * @since 1.3.0
	 */
	public static PShape translateCornerTo(PShape shape, double x, double y) {
		Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();
		AffineTransformation t = AffineTransformation.translationInstance(x - e.getMinX(), y - e.getMinY());
		return toPShape(t.transform(g));
	}

	/**
	 * Translates the given shape such that the upper-left corner of its bounding
	 * box aligns with the origin point (0, 0) of the Processing coordinate system.
	 *
	 * @param shape The PShape instance to be translated.
	 * @return A new PShape instance that is a translation of the input shape such
	 *         that the upper-left corner of its bounding box aligns with the origin
	 *         point (0, 0).
	 */
	public static PShape translateToOrigin(PShape shape) {
		final Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();
		AffineTransformation t = AffineTransformation.translationInstance(-e.getMinX(), -e.getMinY());
		return toPShape(t.transform(g));
	}

	/**
	 * Calculates a Homothetic Transformation of a shape.
	 * <p>
	 * A Homothetic Transformation is a special geometric transformation that
	 * enlarges or shrinks geometries by a scale factor that is the same in all
	 * directions according to a centric point.
	 * 
	 * @param shape  shape input
	 * @param center coordinate of the center/origin position of the operation
	 * @param scaleX X scale factor
	 * @param scaleY Y scale factor
	 * @deprecated
	 */
	@Deprecated
	public static PShape homotheticTransformation(PShape shape, PVector center, double scaleX, double scaleY) {
		Polygon geom = (Polygon) fromPShape(shape);

		// exterior contour
		Coordinate[] coord = geom.getExteriorRing().getCoordinates();
		Coordinate[] coord_ = new Coordinate[coord.length];
		for (int i = 0; i < coord.length; i++) {
			coord_[i] = new Coordinate(center.x + scaleX * (coord[i].x - center.x), center.y + scaleY * (coord[i].y - center.y));
		}
		LinearRing lr = geom.getFactory().createLinearRing(coord_);

		// holes
		LinearRing[] holes = new LinearRing[geom.getNumInteriorRing()];
		for (int j = 0; j < geom.getNumInteriorRing(); j++) {
			Coordinate[] hole_coord = geom.getInteriorRingN(j).getCoordinates();
			Coordinate[] hole_coord_ = new Coordinate[hole_coord.length];
			for (int i = 0; i < hole_coord.length; i++) {
				hole_coord_[i] = new Coordinate(center.x + scaleY * (hole_coord[i].x - center.x), center.y + scaleY * (hole_coord[i].y - center.y));
			}
			holes[j] = geom.getFactory().createLinearRing(hole_coord_);
		}
		return toPShape(PGS.GEOM_FACTORY.createPolygon(lr, holes));
	}

	/**
	 * Aligns one polygon shape to another, using Procrustes analysis to find the
	 * optimal transformation. The transformation includes translation, rotation and
	 * scaling to maximize overlap between the two shapes.
	 * 
	 * @param shapeToAlign   the polygon shape to be transformed and aligned to the
	 *                       other shape.
	 * @param referenceShape the shape that the other shape will be aligned to.
	 * @return a new PShape that is the transformed and aligned version of
	 *         sourceShape.
	 * @since 1.4.0
	 */
	public static PShape align(PShape shapeToAlign, PShape referenceShape) {
		return align(shapeToAlign, referenceShape, 1);
	}

	/**
	 * Aligns one polygon shape to another, using Procrustes analysis to find the
	 * optimal transformation. The transformation includes translation, rotation and
	 * scaling to maximize overlap between the two shapes.
	 * <p>
	 * This method signature aligns the shape according to a provided ratio,
	 * indicating how much alignment transformation to apply.
	 * 
	 * @param shapeToAlign   the polygon shape to be transformed and aligned to the
	 *                       reference shape.
	 * @param referenceShape the shape that the other shape will be aligned to.
	 * @param alignmentRatio a value between 0 and 1 indicating the degree of
	 *                       alignment transformation to apply.
	 * @return a new PShape that is the transformed and aligned version of
	 *         alignShape.
	 * @since 1.4.0
	 */
	public static PShape align(PShape shapeToAlign, PShape referenceShape, double alignmentRatio) {
		return align(shapeToAlign, referenceShape, alignmentRatio, true, true, true);
	}

	/**
	 * Aligns one polygon shape to another, allowing for control over the
	 * application of scaling, translation, and rotation transformations
	 * individually.
	 *
	 * @param shapeToAlign     the polygon shape to be aligned to the reference
	 *                         shape.
	 * @param referenceShape   the reference shape to which the
	 *                         <code>shapeToAlign</code> will be aligned.
	 * @param alignmentRatio   a value between 0 and 1 indicating the degree of
	 *                         alignment transformation to apply.
	 * @param applyScale       if true, applies scaling alignment
	 * @param applyTranslation if true, applies transformation alignment
	 * @param applyRotation    if true, applies rotation alignment
	 * @return a new PShape that is the transformed and aligned version of
	 *         <code>shapeToAlign</code>.
	 * @since 2.0
	 */
	public static PShape align(PShape shapeToAlign, PShape referenceShape, double alignmentRatio, boolean applyScale, boolean applyTranslation,
			boolean applyRotation) {
		final double[] params = getProcrustesParams(shapeToAlign, referenceShape);

		final Geometry g1 = fromPShape(shapeToAlign);
		Coordinate c = g1.getCentroid().getCoordinate();

		double scale = applyScale ? 1 + (params[2] - 1) * alignmentRatio : 1;
		double rotation = applyRotation ? params[3] * alignmentRatio : 0;
		double translateX = applyTranslation ? params[0] * alignmentRatio : 0;
		double translateY = applyTranslation ? params[1] * alignmentRatio : 0;

		AffineTransformation transform = AffineTransformation.scaleInstance(scale, scale, c.x, c.y).rotate(rotation, c.x, c.y).translate(translateX,
				translateY);

		Geometry aligned = transform.transform(g1);

		return toPShape(aligned);
	}

	/**
	 * @return an array having 4 values: the optimal translation (x, y), scale, and
	 *         rotation angle (radians, clockwise).
	 */
	private static double[] getProcrustesParams(PShape shapeToAlign, PShape referenceShape) {
		final Geometry g1 = fromPShape(shapeToAlign);
		final Geometry g2 = fromPShape(referenceShape);
		if (!g1.getGeometryType().equals(Geometry.TYPENAME_POLYGON) || !g2.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			throw new IllegalArgumentException("Inputs to align() must be polygons.");
		}
		if (((Polygon) g1).getNumInteriorRing() > 0 || ((Polygon) g2).getNumInteriorRing() > 0) {
			throw new IllegalArgumentException("Polygon inputs to align() must be holeless.");
		}

		// both shapes need same vertex quantity, simplify rather than densify
		final int vertices = Math.min(shapeToAlign.getVertexCount(), referenceShape.getVertexCount());
		PShape referenceShapeT = referenceShape;
		PShape shapeToAlignT = shapeToAlign;
		if (shapeToAlign.getVertexCount() > vertices) {
			shapeToAlignT = PGS_Morphology.simplifyDCE(shapeToAlign, (v, r, verticesRemaining) -> verticesRemaining <= vertices);
		}
		if (referenceShape.getVertexCount() > vertices) {
			referenceShapeT = PGS_Morphology.simplifyDCE(referenceShape, (v, r, verticesRemaining) -> verticesRemaining <= vertices);
		}

		return ProcrustesAlignment.transform((Polygon) fromPShape(referenceShapeT), (Polygon) fromPShape(shapeToAlignT));
	}

	/**
	 * Rotates a shape around a given point.
	 * 
	 * @param shape the shape to tranform/rotate
	 * @param point rotation point
	 * @param angle the rotation angle, in radians
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
	 * @param y     y-coordinate of horizontal reflection line
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
	 * @param x     x-coordinate of vertical reflection line
	 * @return
	 */
	public static PShape flipVertical(PShape shape, double x) {
		AffineTransformation t = AffineTransformation.reflectionInstance(x, -1, x, 1);
		return toPShape(t.transform(fromPShape(shape)));
	}

	/**
	 * Shears a given shape by specified angles along the x and y axis and returns
	 * the result as a new PShape. Shapes are sheared around their relative position
	 * to the origin.
	 * 
	 * @param shape  The shape to be sheared.
	 * @param angleX The angle by which the shape should be sheared along the
	 *               x-axis, in radians.
	 * @param angleY The angle by which the shape should be sheared along the
	 *               y-axis, in radians.
	 * @return A new shape representing the sheared shape.
	 */
	public static PShape shear(PShape shape, double angleX, double angleY) {
		Geometry g = fromPShape(shape);
		AffineTransformation t = AffineTransformation.shearInstance(angleX, angleY);
		return toPShape(t.transform(g));
	}

}
