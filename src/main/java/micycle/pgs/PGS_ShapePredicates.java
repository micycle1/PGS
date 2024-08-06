package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point3d;
import javax.vecmath.Point4d;

import org.locationtech.jts.algorithm.Angle;
import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.algorithm.match.HausdorffSimilarityMeasure;
import org.locationtech.jts.coverage.CoverageUnion;
import org.locationtech.jts.coverage.CoverageValidator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.operation.valid.IsValidOp;

import micycle.pgs.commons.EllipticFourierDesc;
import micycle.pgs.commons.GeometricMedian;
import micycle.trapmap.TrapMap;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Various shape metrics, predicates and descriptors.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_ShapePredicates {

	private PGS_ShapePredicates() {
	}

	/**
	 * Determines whether the outer shape fully contains the inner shape. A shape is
	 * considered to contain itself. Points of the inner shape that lie on the
	 * boundary of the outer shape are considered to be contained.
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
	 * {@link #containsPoint(PShape, PVector) containsPoint()} repeatedly. Any
	 * points that lie on the boundary of the shape are considered to be contained.
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
	 * @since 1.3.0
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
				System.err.println(e.getMessage());
				return new PShape();
			}
		}
		return map.findContainingPolygon(point.x, point.y);
	}

	/**
	 * Determines whether the shapes intersect/overlap (meaning that have at least
	 * one point in common).
	 * <p>
	 * Note that the input shapes may be lineal (open path) or polygonal (closed
	 * path), and this affects the meaning of the method. The following intersection
	 * tests are performed based on the type combinations of the shapes:
	 * <ul>
	 * <li>Polygon-line intersection: the polygon <b>area</b> and line share at
	 * least one point in common. This means a path contained entirely inside a
	 * polygon will return true as it needn't intersect with the polygon's
	 * perimeter.</li>
	 * <li>Line-line intersection: the two lines intersect at least once.</li>
	 * <li>Polygon-polygon intersection: the two polygons share at least one point
	 * in common (from their area or perimeter).</li>
	 * </ul>
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
	public static double distance(PShape a, PShape b) {
		return fromPShape(a).distance(fromPShape(b));
	}

	/**
	 * Computes the area of the given shape.
	 * 
	 * @param shape
	 * @return
	 */
	public static double area(PShape shape) {
		return fromPShape(shape).getArea();
	}

	/**
	 * Computes the ratio (density) of the shape's area compared to the area of it's
	 * envelope.
	 * 
	 * @param shape
	 * @return Density value. A rectangular shape will have a value of 1.
	 */
	public static double density(PShape shape) {
		Geometry g = fromPShape(shape);
		return (g.getArea() / g.getEnvelope().getArea());
	}

	/**
	 * Computes the centroid of a shape. A centroid is the center of mass of the
	 * shape.
	 * <p>
	 * If the input is a polygon, the centroid will always lie inside it. For other
	 * shapes, this may not be the case.
	 * 
	 * @param shape the PShape object representing the shape for which the centroid
	 *              is to be computed
	 * @return a PVector representing the centroid of the shape, or null if the
	 *         shape's geometry is empty
	 */
	public static PVector centroid(PShape shape) {
		Point point = fromPShape(shape).getCentroid();
		if (!point.isEmpty()) {
			return new PVector((float) point.getX(), (float) point.getY());
		}
		return null;
	}

	/**
	 * Computes the center of the bounding box of a shape. The bounding box is the
	 * smallest rectangle that completely contains the shape.
	 *
	 * @param shape the PShape object representing the shape for which the bounding
	 *              box center is to be computed
	 * @return a PVector representing the center of the shape's bounding box, or
	 *         null if the shape's geometry is empty
	 * @since 2.0
	 */
	public static PVector boundsCenter(PShape shape) {
		Point point = fromPShape(shape).getEnvelope().getCentroid();
		if (!point.isEmpty()) {
			return new PVector((float) point.getX(), (float) point.getY());
		}
		return null;
	}

	/**
	 * Computes the geometric median location of a shape's vertices.
	 * <p>
	 * The median point is the point that minimises the sum of distances to the
	 * shape vertices. If the input is a concave polygon, the median may not lie
	 * inside it.
	 * 
	 * @param shape
	 * @return median point
	 * @since 1.4.0
	 */
	public static PVector median(PShape shape) {
		List<PVector> points = PGS_Conversion.toPVector(shape);
		Point4d[] wp = points.stream().map(p -> new Point4d(p.x, p.y, 0, 1)).toArray(Point4d[]::new);
		Point3d median = GeometricMedian.median(wp, 1e-3, 50);
		return new PVector((float) median.x, (float) median.y);
	}

	/**
	 * Computes the horizontal width of a shape (the width of its bounding-box).
	 */
	public static double width(PShape shape) {
		return fromPShape(shape).getEnvelopeInternal().getWidth();
	}

	/**
	 * Computes the vertical height of a shape (the height of its bounding-box).
	 */
	public static double height(PShape shape) {
		return fromPShape(shape).getEnvelopeInternal().getHeight();
	}

	/**
	 * Returns the length of a shape. Linear shapes return their length; areal
	 * shapes (polygons) return their perimeter.
	 */
	public static double length(PShape shape) {
		return fromPShape(shape).getLength();
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
	public static double diameter(PShape shape) {
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
	public static double circularity(PShape shape) {
		final Polygon poly = (Polygon) fromPShape(shape);
		final double length = poly.getBoundary().getLength();
		return (4 * PConstants.PI * poly.getArea() / (length * length));
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
	public static double similarity(PShape a, PShape b) {
		HausdorffSimilarityMeasure sm = new HausdorffSimilarityMeasure();
		return sm.measure(fromPShape(a), fromPShape(b));
	}

	/**
	 * Measures the degree of <b>mutual overlap</b> between two shapes.
	 * <p>
	 * This metric aggregates how much each shape is overlapped (fractional),
	 * weighted by its respective area.
	 * 
	 * @param a first shape
	 * @param b second shape
	 * @return overlap metric, in [0.0, 1.0]
	 * @since 1.3.0
	 */
	public static double overlap(PShape a, PShape b) {
		Geometry g1 = fromPShape(a);
		Geometry g2 = fromPShape(b);
		Geometry overlap = g1.intersection(g2);
		double aOverlap = overlap.getArea();
		if (aOverlap == 0) {
			return 0;
		}
		double a1 = g1.getArea();
		double a2 = g2.getArea();
		double total = a1 + a2;
		double w1 = a1 / total;
		double w2 = a2 / total;
		return w1 * (aOverlap / a1) + w2 * (aOverlap / a2);
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
	 * Measures the elongation of a shape; the ratio of a shape's bounding box
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
	 * Computes the convexity of a shape using a simple area-based measure of
	 * convexity.
	 * 
	 * @param shape
	 * @return a value in [0, 1]
	 * @since 1.4.0
	 */
	public static double convexity(PShape shape) {
		// also see 'A New Convexity Measure for Polygons'
		Geometry g = fromPShape(shape);
		return g.getArea() / g.convexHull().getArea();
	}

	/**
	 * Counts the number of holes in a shape.
	 * <p>
	 * If the shape forms a polygon coverage (a mesh), then this method will count
	 * holes from gaps within the mesh.
	 * 
	 * @param shape a polygonal shape (can be a GROUP shape having multiple
	 *              polygons)
	 * @return total number of holes in the shape
	 */
	public static int holes(PShape shape) {
		Geometry g = fromPShape(shape);

		/*
		 * Attempt to convert to a polygon coverage to identify and count mesh holes.
		 */
		if (g.getNumGeometries() > 2) { // only 3 or more faces can form a mesh hole
			Geometry[] geoms = new Geometry[g.getNumGeometries()];
			for (int i = 0; i < g.getNumGeometries(); i++) {
				geoms[i] = g.getGeometryN(i);
			}
			if (CoverageValidator.isValid(geoms)) {
				g = CoverageUnion.union(geoms);
			}
		}

		@SuppressWarnings("unchecked")
		List<Polygon> polygons = PolygonExtracter.getPolygons(g);

		int holes = 0;
		for (Polygon p : polygons) {
			holes += p.getNumInteriorRing();
		}
		return holes;
	}

	/**
	 * Computes the maximum/largest interior angle of a polygon.
	 * 
	 * @param shape simple polygonal shape
	 * @return an angle in the range [0, 2PI]
	 * @since 1.3.0
	 */
	public static double maximumInteriorAngle(PShape shape) {
		final Coordinate[] coordz = fromPShape(shape).getCoordinates();
		final CoordinateList coords = new CoordinateList(coordz);
		coords.remove(coords.size() - 1); // remove duplicate/closed coordinate
		if (Orientation.isCCW(coordz)) {
			Collections.reverse(coords); // CCW -> CW
		}
		double maxAngle = 0;
		for (int i = 0; i < coords.size(); i++) {
			Coordinate p0 = coords.get(i);
			Coordinate p1 = coords.get((i + 1) % coords.size());
			Coordinate p2 = coords.get((i + 2) % coords.size());
			maxAngle = Math.max(maxAngle, Angle.interiorAngle(p0, p1, p2));
		}
		return maxAngle;
	}

	/**
	 * Quantifies the similarity between two shapes, by using the pairwise euclidean
	 * distance between each shape's <i>Elliptic Fourier Descriptors</i> (EFD).
	 * <p>
	 * Smaller values indicate greater similarity or equivalence, and the measure is
	 * translation and rotation invariant.
	 * <p>
	 * This method can be useful in shape recognition tasks where it is necessary to
	 * quantify the difference or similarity between two shapes.
	 *
	 * @param a polygonal shape
	 * @param b polygonal shape
	 * @return The EFD distance between the two provided PShapes. Smaller values
	 *         indicate greater similarity or equivalence between the shapes.
	 * @since 1.4.0
	 */
	public static double efdSimilarity(PShape a, PShape b) {
		int n = Math.min(a.getVertexCount(), b.getVertexCount()) / 2;
		n = Math.min(n, 50); // max of 50 descriptors
		EllipticFourierDesc efdA = new EllipticFourierDesc(((Polygon) fromPShape(a)).getExteriorRing(), n);
		EllipticFourierDesc efdB = new EllipticFourierDesc(((Polygon) fromPShape(b)).getExteriorRing(), n);
		return EllipticFourierDesc.computeEFDDistance(efdA.getEFD(), efdB.getEFD());
	}

	/**
	 * Tests two shapes for <b>structural equality</b>. In simple terms, this means
	 * that they must have the same number of vertices, in the same locations, and
	 * in the same order.
	 * <p>
	 * Note: If two Polygons have matching vertices, but one is arranged clockwise
	 * while the other is counter-clockwise, then then this method will return
	 * false.
	 * 
	 * @param a shape a
	 * @param b shape b
	 * @return true if both shapes have identical structure and point values.
	 * @since 1.3.0
	 * @see #equalsNorm(PShape, PShape)
	 * @see #equalsTopo(PShape, PShape)
	 */
	public static boolean equalsExact(PShape a, PShape b) {
		return fromPShape(a).equalsExact(fromPShape(b));
	}

	/**
	 * Tests two shapes for <b>normalised structural equality</b>. In simple terms,
	 * this means that they must have the same number of vertices in the same
	 * locations. Unlike {@link #equalsExact(PShape, PShape)}, vertices do not need
	 * to be in the same order for the shapes to be considered equal.
	 * 
	 * @param a shape a
	 * @param b shape b
	 * @return true the shapes are exactly equal in their normalized form
	 * @since 1.3.0
	 * @see #equalsExact(PShape, PShape)
	 * @see #equalsTopo(PShape, PShape)
	 */
	public static boolean equalsNorm(PShape a, PShape b) {
		return fromPShape(a).equalsNorm(fromPShape(b));
	}

	/**
	 * Tests two shapes for <b>topological equality</b>. In simple terms, this is
	 * equivalent to drawing the two shapes and seeing if all of their component
	 * edges overlap. It is the most robust kind of comparison but also the most
	 * computationally expensive.
	 * 
	 * @param a shape a
	 * @param b shape b
	 * @return true if the two shapes are topologically equal
	 * @since 1.3.0
	 * @see #equalsExact(PShape, PShape)
	 * @see #equalsTopo(PShape, PShape)
	 */
	public static boolean equalsTopo(PShape a, PShape b) {
		return fromPShape(a).equalsTopo(fromPShape(b));
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
	 * angles are all less than or equal to 180°.
	 */
	public static boolean isConvex(PShape shape) {
		final Geometry g = fromPShape(shape);
		final double area = g.getArea();
		return ((g.convexHull().getArea() - area) / area < 0.001);
	}

	/**
	 * Determines whether a GROUP shape forms a conforming mesh / valid polygon
	 * coverage.
	 * <p>
	 * Conforming meshes comprise faces that do not intersect; any adjacent faces
	 * not only share edges, but every pair of shared edges are <b>identical</b>
	 * (having the same coordinates) (such as a triangulation).
	 * 
	 * @param mesh shape to test
	 * @return true if the shape is a conforming mesh
	 * @since 1.4.0
	 */
	public static boolean isConformingMesh(PShape mesh) {
		Geometry[] geoms = PGS_Conversion.getChildren(mesh).stream().map(f -> fromPShape(f)).toArray(Geometry[]::new);
		return CoverageValidator.isValid(geoms);
	}

	/**
	 * Checks if a PShape is valid, and reports the validation error if it is
	 * invalid.
	 * <p>
	 * An invalid shape is one that violates the rules of geometric validity. Some
	 * common reasons for a shape to be considered invalid include:
	 * <ul>
	 * <li>Self-intersection: The shape intersects itself at one or more points or
	 * segments, creating overlapping or self-crossing areas.
	 * <li>Invalid topology: The shape's topology is incorrect, such as having
	 * dangling edges, disconnected components, or invalid ring configurations in
	 * polygons.
	 * <li>Degenerate geometry: The shape has collapsed or degenerate components,
	 * such as zero-length lines, zero-area polygons, or overlapping vertices.
	 * </ul>
	 * 
	 * @param shape The PShape to validate.
	 * @return {@code true} if the shape is valid, {@code false} otherwise.
	 * @since 1.4.0
	 */
	public static boolean isValid(PShape shape) {
		IsValidOp validate = new IsValidOp(fromPShape(shape));
		if (validate.getValidationError() == null) {
			return true;
		} else {
			System.err.println(validate.getValidationError());
			return false;
		}
	}

}
