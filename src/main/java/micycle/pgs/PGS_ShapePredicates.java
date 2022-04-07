package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.algorithm.match.HausdorffSimilarityMeasure;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

import micycle.trapmap.TrapMap;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Various shape metrics &amp; predicates
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_ShapePredicates {

	private PGS_ShapePredicates() {
	}

	/**
	 * Determines whether the outer shape contains the inner shape (meaning every
	 * point of the inner shape is a point of the outer shape). A shape is
	 * considered to contain itself itself.
	 * 
	 * @param outer
	 * @param inner
	 * @return
	 */
	public static boolean contains(PShape outer, PShape inner) {
		return fromPShape(outer).covers(fromPShape(inner));
	}

	/**
	 * Determines whether a shape contains a point. Points that lie on the boundary
	 * of the shape are considered to be contained.
	 * 
	 * @param shape
	 * @param point
	 * @return
	 * @see #containsAllPoints(PShape, Collection)
	 * @see #containsPoints(PShape, Collection)
	 */
	public static boolean containsPoint(PShape shape, PVector point) {
		return fromPShape(shape).covers(PGS.pointFromPVector(point));
	}

	/**
	 * Determines whether a shape contains every point from a list of points. It is
	 * faster to use method rather than than calling
	 * {@link #containsPoint(PShape, PVector)} repeatedly. Points that lie on the
	 * boundary of the shape are considered to be contained.
	 * 
	 * @param shape
	 * @param points list of points to check
	 * @return true if every point is contained within the shape
	 */
	public static boolean containsAllPoints(PShape shape, Collection<PVector> points) {
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(fromPShape(shape));
		for (PVector p : points) {
			if (pointLocator.locate(new Coordinate(p.x, p.y)) == Location.EXTERIOR) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Measures for each point in the input whether it is contained in the given
	 * shape. This method checks every point individually, returning a boolean for
	 * each point. Using this method is faster than calling
	 * {@link #containsPoint(PShape, PVector) containsPoint()} repeatedly. Points
	 * that lie on the boundary of the shape are considered to be contained.
	 * 
	 * @param shape
	 * @param points list of points to check
	 * @return a list of booleans corresponding to whether the shape contains the
	 *         point at same index
	 */
	public static List<Boolean> containsPoints(PShape shape, Collection<PVector> points) {
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(fromPShape(shape));
		ArrayList<Boolean> bools = new ArrayList<>(points.size());
		for (PVector p : points) {
			bools.add(pointLocator.locate(new Coordinate(p.x, p.y)) != Location.EXTERIOR);
		}
		return bools;
	}

	/**
	 * Tests for each point in the input whether it is contained in/inside the given
	 * shape; if it is, then the point is included in the output list. This method
	 * does not mutate the input; it returns a filtered copy. Points that lie on the
	 * boundary of the shape are considered to be contained.
	 * <p>
	 * Using this method is faster than calling
	 * {@link #containsPoint(PShape, PVector) containsPoint()} repeatedly.
	 * 
	 * @param shape
	 * @param points list of points to check
	 * @return a filtered view of the input points
	 */
	public static List<PVector> findContainedPoints(PShape shape, Collection<PVector> points) {
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(fromPShape(shape));
		List<PVector> contained = new ArrayList<>();
		for (PVector p : points) {
			if (pointLocator.locate(new Coordinate(p.x, p.y)) != Location.EXTERIOR) {
				contained.add(p);
			}
		}
		return contained;
	}

	/**
	 * Finds the single child shape/cell (if any) that contains the query point from
	 * a GROUP shape input (a shape that has non-overlapping children).
	 * <p>
	 * This method locates the containing shape in log(n) time (after some
	 * pre-processing overhead).
	 * 
	 * @param groupShape a GROUP shape
	 * @param point      the query point
	 * @return the child shape that contains the query point, or null if no child
	 *         shape contains the point
	 * @since 1.2.1
	 */
	public static PShape findContainingShape(PShape groupShape, PVector point) {
		if (groupShape.getKind() != PConstants.GROUP) { // handle non-mesh shape
			if (containsPoint(groupShape, point)) {
				return groupShape;
			} else {
				return null;
			}
		}

		TrapMap map;
		try {
			map = new TrapMap(PGS_Conversion.getChildren(groupShape));
		} catch (Exception e) {
			/*
			 * Handle error thrown by TrapMap on degenerate/strange inputs. Generally
			 * shearing will fix the problem (ideally this would be done within TrapMap).
			 */
			try {
				map = new TrapMap(PGS_Conversion.getChildren(PGS_Transformation.shear(groupShape, .00001, 0)));
			} catch (Exception e2) {
				System.err.println(e.getLocalizedMessage());
				return new PShape();
			}
		}
		return map.findContainingPolygon(point.x, point.y);
	}

	/**
	 * Determines whether the shapes intersect/overlap (meaning that have at least
	 * one point in common).
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static boolean intersect(PShape a, PShape b) {
		return fromPShape(a).intersects(fromPShape(b));
	}

	/**
	 * Determines whether the have at least one point in common, but where their
	 * interiors do not intersect.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static boolean touch(PShape a, PShape b) {
		return fromPShape(a).touches(fromPShape(b));
	}

	/**
	 * Computes the minimum distance between two shapes.
	 * 
	 * @param a shape A
	 * @param b shape B
	 * @return
	 */
	public static float distance(PShape a, PShape b) {
		return (float) fromPShape(a).distance(fromPShape(b));
	}

	/**
	 * Computes the area of the given shape.
	 * 
	 * @param shape
	 * @return
	 */
	public static float area(PShape shape) {
		return (float) fromPShape(shape).getArea();
	}

	/**
	 * Computes the ratio (density) of the shape's area compared to the area of it's
	 * envelope.
	 * 
	 * @param shape
	 * @return Density value. A rectangular shape will have a value of 1.
	 */
	public static float density(PShape shape) {
		Geometry g = fromPShape(shape);
		return (float) (g.getArea() / g.getEnvelope().getArea());
	}

	/**
	 * Computes the centroid of a shape. A centroid is the center of mass of the
	 * shape.
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
	 * Computes the horizontal width of a shape (the width of its bounding-box).
	 */
	public static float width(PShape shape) {
		return (float) fromPShape(shape).getEnvelopeInternal().getWidth();
	}

	/**
	 * Computes the vertical height of a shape (the height of its bounding-box).
	 */
	public static float height(PShape shape) {
		return (float) fromPShape(shape).getEnvelopeInternal().getHeight();
	}

	/**
	 * Returns the length of a shape. Linear shapes return their length; areal
	 * shapes (polygons) return their perimeter.
	 */
	public static float length(PShape shape) {
		return (float) fromPShape(shape).getLength();
	}

	/**
	 * Returns the diameter of a shape. Diameter is the maximum distance between any
	 * 2 coordinates on the shape perimeter; this is equal to the diameter of the
	 * circumscribed circle.
	 * 
	 * @param shape
	 * @return
	 * @since 1.1.3
	 */
	public static float diameter(PShape shape) {
		List<PVector> farPoints = PGS_Optimisation.farthestPointPair(PGS_Conversion.toPVector(shape));
		return farPoints.get(0).dist(farPoints.get(1));
	}

	/**
	 * Calculates the Miller circularity index for a shape. This index, between 0
	 * and 1, is equal to 1 if the polygon is perfectly circular and tends towards 0
	 * for a segment.
	 * 
	 * @param shape
	 * @return
	 */
	public static float circularity(PShape shape) {
		final Polygon poly = (Polygon) fromPShape(shape);
		final double length = poly.getBoundary().getLength();
		return (float) (4 * PConstants.PI * poly.getArea() / (length * length));
	}

	/**
	 * Measures the degree of similarity between two shapes using the Hausdorff
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
	 * Measures the sphericity of a shape; the ratio of the maximum inscribed circle
	 * to the minimum bounding circle.
	 * 
	 * @param shape
	 * @return a value in [0, 1]
	 */
	public static double sphericity(final PShape shape) {
		Geometry g = fromPShape(shape);
		MinimumBoundingCircle circle1 = new MinimumBoundingCircle(g);
		final double rO = circle1.getRadius();
		MaximumInscribedCircle circle2 = new MaximumInscribedCircle(g, 1);
		Point center = circle2.getCenter();
		Point radiusPoint = circle2.getRadiusPoint();
		final double rI = center.distance(radiusPoint);
		return Math.min(1, rI / rO);
	}

	/**
	 * Measures the elongation of a shape; the ratio of the shapes bounding box's
	 * length to its width.
	 * 
	 * @param shape
	 * @return a value in [0, 1]
	 */
	public static double elongation(final PShape shape) {
		Geometry obb = MinimumDiameter.getMinimumRectangle(fromPShape(shape));
		Polygon rect = (Polygon) obb;
		Coordinate c0 = rect.getCoordinates()[0];
		Coordinate c1 = rect.getCoordinates()[1];
		Coordinate c2 = rect.getCoordinates()[2];
		double l = c0.distance(c1);
		double w = c1.distance(c2);
		if (l >= w) {
			return w / l;
		} else {
			return l / w;
		}
	}

	/**
	 * Computes the number of holes in a shape.
	 * 
	 * @return
	 */
	public static int holes(PShape shape) {
		return ((Polygon) fromPShape(shape)).getNumInteriorRing(); // NOTE assume a single polygon
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
	 * angles are less than or equal to 180°.
	 */
	public static boolean isConvex(PShape shape) {
		final Geometry g = fromPShape(shape);
		final double area = g.getArea();
		return ((g.convexHull().getArea() - area) / area < 0.001);
	}

}
