package micycle.pts;

import static micycle.pts.Conversion.fromPShape;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.algorithm.match.HausdorffSimilarityMeasure;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Various shape metrics & predicates
 * 
 * @author Michael Carleton
 *
 */
public class PTSShapePredicates {

	// https://doc.cgal.org/latest/Polygon/index.html#Chapter_2D_Polygons

	/**
	 * A shape does contain itself.
	 * 
	 * @param outer
	 * @param inner
	 * @return
	 */
	public static boolean contains(PShape outer, PShape inner) {
		return fromPShape(outer).covers(fromPShape(inner));
	}

	public static boolean containsPoint(PShape shape, PVector point) {
		return fromPShape(shape).covers(PTS.pointFromPVector(point));
	}

	/**
	 * Faster than calling {@link #containsPoint(PShape, PVector)} repeatedly.
	 * Points that lie on the boundary of the shape considered to be contained.
	 * 
	 * @param shape
	 * @param points list of PVectors to check
	 * @return true if every point is contained within the shape
	 */
	public static boolean containsAllPoints(PShape shape, List<PVector> points) {
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(fromPShape(shape));
		for (PVector p : points) {
			if (pointLocator.locate(new Coordinate(p.x, p.y)) == Location.EXTERIOR) {
				return false;
			}
		}
		return true;
	}

	/**
	 * This method checks every point individually, returning a boolean for each
	 * point. Faster than calling {@link #containsPoint(PShape, PVector)}
	 * repeatedly. Points that lie on the boundary of the shape considered to be
	 * contained.
	 * 
	 * @param shape
	 * @param points list of PVectors to check
	 * @return a list of booleans corresponding to whether the shape contains the
	 *         point at same index
	 */
	public static List<Boolean> containsPoints(PShape shape, List<PVector> points) {
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(fromPShape(shape));
		ArrayList<Boolean> bools = new ArrayList<>(points.size());
		for (PVector p : points) {
			if (pointLocator.locate(new Coordinate(p.x, p.y)) == Location.EXTERIOR) {
				bools.add(false);
			} else {
				bools.add(true);
			}
		}
		return bools;
	}

	/**
	 * intersect/overlap: at least one point in common
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static boolean intersect(PShape a, PShape b) {
		return fromPShape(a).intersects(fromPShape(b));
	}

	/**
	 * Checks whether the have at least one point in common, but their interiors do
	 * not intersect.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static boolean touch(PShape a, PShape b) {
		return fromPShape(a).touches(fromPShape(b));
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
		return (float) (4 * PConstants.PI * poly.getArea()
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
