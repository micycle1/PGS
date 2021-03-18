package micycle.pts;

import static micycle.pts.Conversion.fromPShape;

import org.locationtech.jts.algorithm.match.HausdorffSimilarityMeasure;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

import processing.core.PApplet;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Various shape metrics & predicates
 * 
 * @author Michael Carleton
 *
 */
public class PTSShapeMetrics {

	// https://doc.cgal.org/latest/Polygon/index.html#Chapter_2D_Polygons

	public static boolean contains(PShape outer, PShape inner) {
		return fromPShape(outer).contains(fromPShape(inner));
	}

	public static boolean containsPoint(PShape shape, PVector point) {
		return fromPShape(shape).contains(PTS.pointFromPVector(point));
	}

	/**
	 * Returns the minimum distance between two shapes.
	 * 
	 * @param a shape A
	 * @param b shape B
	 * @return
	 */
	public static float distance(PShape a, PShape b) {
		return (float) fromPShape(a).distance(fromPShape(b));
	}

	public static float area(PShape shape) {
		return (float) fromPShape(shape).getArea();
	}

	/**
	 * 
	 * @param shape
	 * @return null if point is empty (geometry empty)
	 */
	public static PVector centroid(PShape shape) {
		Point point = fromPShape(shape).getCentroid();
		if (!point.isEmpty()) {
			return new PVector((float) point.getX(), (float) point.getY());
		}
		return null;
	}

	/**
	 * Returns the length of a shape. Linear shapes return their length; areal
	 * shapes (polygons) return their perimeter.
	 */
	public static float length(PShape shape) {
		return (float) fromPShape(shape).getLength();
	}

	/**
	 * * Calculates the Miller circularity index for a polygon. This index, between
	 * 0 and 1, is equal to 1 if the polygon is perfectly circular and tends towards
	 * 0 for a segment.
	 * 
	 * @param shape
	 * @return
	 */
	public static float circularity(PShape shape) {
		Polygon poly = (Polygon) fromPShape(shape);
		return (float) (4 * PApplet.PI * poly.getArea()
				/ (poly.getBoundary().getLength() * poly.getBoundary().getLength()));
	}

	/**
	 * Measures the degree of similarity between two Geometrysusing the Hausdorff
	 * distance metric. The measure is normalized to lie in the range [0, 1]. Higher
	 * measures indicate a great degree of similarity.
	 * 
	 * @param a first shape
	 * @param b second shape
	 * @return the value of the similarity measure, in [0.0, 1.0]
	 */
	public static float similarity(PShape a, PShape b) {
		HausdorffSimilarityMeasure sm = new HausdorffSimilarityMeasure();
		return (float) sm.measure(fromPShape(a), fromPShape(b));
	}

	/**
	 * Checks whether a shape is simple. A shape is simple if it has no points of
	 * self-tangency, self-intersection or other anomalous points.
	 */
	public static boolean isSimple(PShape shape) {
		return fromPShape(shape).isSimple();
	}

	/**
	 * Determines whether a shape is convex. A shape is convex if its interior
	 * angles less than or equal to 180°.
	 */
	public static boolean isConvex(PShape shape) {
		final Geometry g = fromPShape(shape);
		final double area = g.getArea();
		return ((g.convexHull().getArea() - area) / area < 0.001);
	}

}
