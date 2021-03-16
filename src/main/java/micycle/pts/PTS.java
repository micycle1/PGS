package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import org.geodelivery.jap.concavehull.SnapHull;
import org.geotools.geometry.jts.JTS;
import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.algorithm.match.HausdorffSimilarityMeasure;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.shape.random.RandomPointsBuilder;
import org.locationtech.jts.shape.random.RandomPointsInGridBuilder;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.simplify.VWSimplifier;
import org.locationtech.jts.util.GeometricShapeFactory;
import micycle.pts.color.Blending;
import micycle.pts.utility.PolygonDecomposition;
import micycle.pts.utility.RandomPolygon;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.concavehull.ConcaveHull;
import uk.osgb.algorithm.concavehull.TriCheckerChi;
import uk.osgb.algorithm.minkowski_sum.Minkowski_Sum;

/**
 * PTS | Processing Topology Suite
 * <p>
 * <ul>
 * <li>https://github.com/IGNF/CartAGen
 * <li>https://ignf.github.io/CartAGen/docs/algorithms/others/spinalize.html
 * <li>https://discourse.processing.org/t/straight-skeleton-or-how-to-draw-a-center-line-in-a-polygon-or-shape/17208/9
 * </p>
 * </ul>
 * 
 * TODO https://ignf.github.io/CartAGen/docs/algorithms.html TODO take into
 * account strokweight (buffer by stroke amount?) TODO maintain pshape fill, etc
 * on output
 * 
 * @author Michael Carleton
 */
public class PTS implements PConstants {

	// TODO check for getCoordinates() in loops (and replace) (if lots of child
	// geometries)

	// TODO replace line() in for-each with beginshape(LINES)...vertex...endShape()

	// TODO use LinearRingIterator when possible (refactor)

	/**
	 * Calling Polygon#union repeatedly is one way to union several Polygons
	 * together. But here’s a trick that can be significantly faster (seconds rather
	 * than minutes) – add the Polygons to a GeometryCollection, then apply a buffer
	 * with zero distance
	 */

	/**
	 * Defines number of vertex samples per bezier line (during PShape->JTS)
	 */
	protected static final int CURVE_SAMPLES = 20;

	public static GeometryFactory GEOM_FACTORY = new GeometryFactory(
			new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	static {
		Minkowski_Sum.setGeometryFactory(GEOM_FACTORY);
	}

	/**
	 * The Maximum Inscribed Circle is determined by a point in the interior of the
	 * area which has the farthest distance from the area boundary, along with a
	 * boundary point at that distance.
	 * 
	 * @param shape
	 * @param tolerance the distance tolerance for computing the centre point
	 */
	public static PShape maximumInscribedCircle(PShape shape, float tolerance) {
		MaximumInscribedCircle mic = new MaximumInscribedCircle(fromPShape(shape), tolerance);

		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(new Coordinate(mic.getCenter().getX(), mic.getCenter().getY()));
		shapeFactory.setWidth(mic.getRadiusLine().getLength() * 2); // r*2 for total width & height
		shapeFactory.setHeight(mic.getRadiusLine().getLength() * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());

	}

	/**
	 * Return the maximum circle (at a given centerpoint inside/outside the circle)
	 * 
	 * @param shape
	 * @param centerPoint
	 * @return
	 */
	public static PShape maximumInscribedCircle(PShape shape, PVector centerPoint) {
		Geometry g = fromPShape(shape);
		Point p = pointFromPVector(centerPoint);
		Coordinate closestEdgePoint = DistanceOp.nearestPoints(g.getBoundary(), p)[0];

		double radius = distance(GEOM_FACTORY.createPoint(closestEdgePoint), p);
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(p.getCoordinate());
		shapeFactory.setWidth(radius * 2); // r*2 for total width & height
		shapeFactory.setHeight(radius * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * Computes the Minimum Bounding Circle (MBC)for the points in a Geometry. The
	 * MBC is the smallest circle which covers all the vertices of the input shape
	 * (this is also known as the Smallest Enclosing Circle). This is equivalent to
	 * computing the Maximum Diameter of the input vertex set.
	 */
	public static PShape minimumBoundingCircle(PShape shape) {
		MinimumBoundingCircle mbc = new MinimumBoundingCircle(fromPShape(shape));
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(new Coordinate(mbc.getCentre().getX(), mbc.getCentre().getY()));
		shapeFactory.setWidth(mbc.getRadius() * 2); // r*2 for total width & height
		shapeFactory.setHeight(mbc.getRadius() * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * Gets the minimum rectangle enclosing a shape.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumBoundingRectangle(PShape shape) {
		Polygon md = (Polygon) MinimumDiameter.getMinimumRectangle(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * 
	 * 
	 * Computes the minimum diameter of a shape. The minimum diameter is defined to
	 * be the width of the smallest band that contains the shape, where a band is a
	 * strip of the plane defined by two parallel lines. This can be thought of as
	 * the smallest hole that the geometry can bemoved through, with a single
	 * rotation.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumDiameter(PShape shape) {
		LineString md = (LineString) MinimumDiameter.getMinimumDiameter(fromPShape(shape));
		return toPShape(md);
	}

//	public static PShape smallestSurroundingRectangle(PShape shape) {
//		return toPShape(
//				SmallestSurroundingRectangleComputation.getSSRPreservedArea(geometryFactory.buildGeometry(null)));
//	}

	/**
	 * 
	 * @param shape
	 * @param fit   0...1
	 * @return
	 */
	private static Polygon smooth(Polygon shape, float fit) {
		// http://lastresortsoftware.blogspot.com/2010/12/smooth-as.html
		return (Polygon) JTS.smooth(shape, fit);
	}

	/**
	 * @param shape
	 * @param fit   tightness of fit from 0 (loose) to 1 (tight)
	 * @return
	 */
	protected static Geometry smooth(Geometry shape, float fit) {
		return JTS.smooth(shape, fit);
	}

	/**
	 * @return shape A - shape B
	 */
	public static PShape difference(PShape a, PShape b) {
		return toPShape(fromPShape(a).difference(fromPShape(b)));
	}

	/**
	 * Compute the parts that the shapes do not have in common.
	 * 
	 * @return A∪B - A∩B
	 */
	public static PShape symDifference(PShape a, PShape b) {
		return toPShape(fromPShape(a).symDifference(fromPShape(b)));
	}

	/**
	 * @return A∩B
	 */
	public static PShape intersection(PShape a, PShape b) {
		PShape out = toPShape(fromPShape(a).intersection(fromPShape(b)));
//		a.draw(p.getGraphics());
		out.setFill(Blending.screen(getPShapeFillColor(a), getPShapeFillColor(b))); // TODO
		return out;
	}

	private static final int getPShapeFillColor(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("fillColor");
			f.setAccessible(true);
			return f.getInt(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * @return A∪B
	 * @see #union(PShape...)
	 */
	public static PShape union(PShape a, PShape b) {
//		OverlayNG.overlay(null, null, CURVE_SAMPLES) // TODO investigate
		return toPShape(fromPShape(a).union(fromPShape(b)));
	}

	public static PShape union(PShape... shapes) {
		ArrayList<Geometry> geoms = new ArrayList<>();
		for (int i = 0; i < shapes.length; i++) {
			geoms.add(fromPShape(shapes[i]));
		}
		return toPShape(UnaryUnionOp.union(geoms));
	}

	/**
	 * 
	 * @param shape
	 * @param buffer extent/width of the buffer (may be positive or negative)
	 * @return
	 */
	public static PShape buffer(PShape shape, float buffer) {
		// TODO read
		// https://locationtech.github.io/jts/javadoc/org/locationtech/jts/operation/buffer/BufferOp.html
		return toPShape(fromPShape(shape).buffer(buffer, 4));
	}

	/**
	 * small edges are removed, while the general structure of the shape is
	 * preserved. This process is known as "opening" in computer vision.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape erosionDilation(PShape shape, float buffer) {
		buffer = Math.abs(buffer);
		return toPShape(fromPShape(shape).buffer(-buffer).buffer(buffer));
	}

	/**
	 * Simplifies a PShape using the Douglas-Peucker algorithm.
	 * 
	 * @param shape
	 * @param distanceTolerance
	 * @return
	 * @see #topologySimplify(PShape, float)
	 */
	public static PShape simplify(PShape shape, float distanceTolerance) {
		return toPShape(DouglasPeuckerSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Simplifies a shape using the Visvalingam-Whyatt area-based algorithm.
	 * 
	 * @param shape
	 * @param distanceTolerance The simplification tolerance is specified as a
	 *                          distance.This is converted to an area tolerance by
	 *                          squaring it.
	 * @return
	 */
	public static PShape simplifyVW(PShape shape, float distanceTolerance) {
		return toPShape(VWSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Preserves topological structure (holes, etc.)
	 * 
	 * @param shape
	 * @param distanceTolerance
	 * @return
	 * @see #simplify(PShape, float)
	 */
	public static PShape topologySimplify(PShape shape, float distanceTolerance) {
		return toPShape(TopologyPreservingSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Densifies a Geometry by inserting extra vertices along the line segments
	 * contained in the geometry. Specify maximum length of segments
	 */
	public static PShape densify(PShape shape, float distanceTolerance) {
		Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(distanceTolerance);
		d.setValidate(false);
		return toPShape(d.getResultGeometry());
	}

	/**
	 * Smoothes a geometry
	 * 
	 * @param shape
	 * @param fit   tightness of fit from 0 (loose) to 1 (tight)
	 * @return
	 */
	public static PShape smooth(PShape shape, float fit) {
		return toPShape(smooth(fromPShape(shape), fit));
	}

	public static PShape flatten(PShape shape) {
		// TODO iterate over geometries shapes, then create group PShape
		Polygon poly = (Polygon) fromPShape(shape).union().getGeometryN(0);
		return toPShape(poly.getExteriorRing());
	}

	public static PShape convexHull(PShape... shapes) {
		Geometry g = fromPShape(shapes[0]);
		for (int i = 1; i < shapes.length; i++) {
			g = g.union(fromPShape(shapes[i]));
		}
		return toPShape(g.convexHull());
	}

	/**
	 * Uses Chi heurstic
	 * 
	 * @param points
	 * @param threshold euclidean distance threshold
	 * @return
	 */
	public static PShape concaveHull(ArrayList<PVector> points, float threshold) {

		// calls Ordnance Survey implementation

		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
			points.add(points.get(0)); // close geometry
		} else { // already closed
			coords = new Coordinate[points.size()];
		}

		for (int i = 0; i < coords.length; i++) {
			coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
		}

		Geometry g = GEOM_FACTORY.createPolygon(coords);
		ConcaveHull hull = new ConcaveHull(g);
		return toPShape(hull.getConcaveHullBFS(new TriCheckerChi(threshold), false, false).get(0));
	}

	/**
	 * Has a more "organic" structure compared to other concave method.
	 * 
	 * @param points
	 * @param threshold 0...1 (Normalized length parameter). Setting λP = 1 means
	 *                  that no edges will be removed from the Delaunay
	 *                  triangulation, so the resulting polygon will be the convex
	 *                  hull. Setting λP = 0 means that all edges that can be
	 *                  removed subject to the regularity constraint will be removed
	 *                  (however polygons that are eroded beyond the point where
	 *                  they provide a desirable characterization of the shape).
	 *                  Although the optimal parameter value varies for different
	 *                  shapes and point distributions, values of between 0.05–0.2
	 *                  typically produce optimal or near-optimal shape
	 *                  characterization across a wide range of point distributions.
	 * @return
	 */
	public static PShape concaveHull2(ArrayList<PVector> points, float threshold) {

		/**
		 * (from https://doi.org/10.1016/j.patcog.2008.03.023) It is more convenient to
		 * normalize the threshold parameter with respect to a particular set of points
		 * P by using the maximum and minimum edge lengths of the Delaunay triangulation
		 * of P. Increasing l beyond the maximum edge length of the Delaunay
		 * triangulation cannot reduce the number of edges that will be removed (which
		 * will be zero anyway). Decreasing l beyond the minimum edge length of the
		 * Delaunay triangulation cannot increase the number of edges that will be
		 * removed.
		 */

		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
			points.add(points.get(0)); // close geometry
		} else { // already closed
			coords = new Coordinate[points.size()];
		}

		for (int i = 0; i < coords.length; i++) {
			coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
		}

		Geometry g = GEOM_FACTORY.createPolygon(coords);

		// TODO test AVG threshold heuristic
		org.geodelivery.jap.concavehull.ConcaveHull hull = new org.geodelivery.jap.concavehull.ConcaveHull(threshold);

		return toPShape(hull.transform(g));
	}

	/**
	 * Adjust segment factor to change between
	 * 
	 * @param shape
	 * @param segmentFactor
	 * @return
	 */
	public static PShape snapHull(PShape shape, float segmentFactor) {
		return toPShape(SnapHull.snapHull(fromPShape(shape), segmentFactor));
	}

	/**
	 * Minkowski addition a.k.a dilation
	 * 
	 * @return
	 */
	public static PShape minkSum(PShape source, PShape addition) {
//		Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
		Geometry sum = Minkowski_Sum.compMinkSum(fromPShape(source), fromPShape(addition), false, false);
		return toPShape(sum);
	}

	/**
	 * TODO check, a.k.a erosion
	 * 
	 * @param source
	 * @param addition
	 * @return
	 */
	public static PShape minkDifference(PShape source, PShape addition) {
//			Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
		Geometry sum = Minkowski_Sum.compMinkDiff(fromPShape(source), fromPShape(addition), false, false);
		return toPShape(sum);
	}

	/**
	 * TODO: outline of exterior ring only
	 * 
	 * @param shape
	 * @param points
	 * @param offsetDistance the distance the point is offset from the
	 *                       segment(positive is to the left, negative is to the
	 *                       right)
	 * @return
	 * @see #equidistantOutlineByDistance(PShape, float, float)
	 */
	public static PVector[] equidistantOutline(PShape shape, int points, float offsetDistance) {
		Polygon p = (Polygon) fromPShape(shape); // TODO CHECK CAST (ITERATE OVER GROUP)

		LengthIndexedLine l = new LengthIndexedLine(p.getExteriorRing());
		PVector[] outlinePoints = new PVector[points * (p.getNumInteriorRing() + 1)];
		// exterior ring
		for (int i = 0; i < points; i++) {
			Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
			outlinePoints[i] = new PVector((float) q.x, (float) q.y);
		}

		for (int j = 0; j < p.getNumInteriorRing(); j++) {
			l = new LengthIndexedLine(p.getInteriorRingN(j));
			for (int i = 0; i < points; i++) {
				Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
				outlinePoints[(j + 1) * points + i] = new PVector((float) q.x, (float) q.y);
			}
		}

		return outlinePoints;
	}

	/**
	 * 
	 * @param shape
	 * @param interPointDistance Distance between each point on outline
	 * @param offsetDistance
	 * @return nearest distance such that every distance is equal (maybe be
	 *         different due to rounding)
	 */
	public static PVector[] equidistantOutlineByDistance(PShape shape, float interPointDistance, float offsetDistance) {
		Polygon p = (Polygon) fromPShape(shape); // TODO CHECK CAST (ITERATE OVER GROUP)

		LengthIndexedLine l = new LengthIndexedLine(p.getExteriorRing());
		if (interPointDistance > l.getEndIndex()) {
			System.err.println("Interpoint greater than shape length");
			return null;
		}
		int points = (int) Math.round(l.getEndIndex() / interPointDistance);
		ArrayList<PVector> outlinePoints = new ArrayList<>();
		// exterior ring
		for (int i = 0; i < points; i++) {
			Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
			outlinePoints.add(new PVector((float) q.x, (float) q.y));
		}

		for (int j = 0; j < p.getNumInteriorRing(); j++) {
			l = new LengthIndexedLine(p.getInteriorRingN(j));
			points = (int) Math.round(l.getEndIndex() / interPointDistance);
			for (int i = 0; i < points; i++) {
				Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
				outlinePoints.add(new PVector((float) q.x, (float) q.y));
			}
		}

		final PVector[] out = new PVector[outlinePoints.size()];
		Arrays.setAll(out, outlinePoints::get);
		return out;
	}

	/**
	 * 
	 * @param shape
	 * @param distance 0...1 around shape perimeter; or -1...0 (other direction)
	 * @return
	 */
	public static PVector pointOnOutline(PShape shape, float distance, float offsetDistance) {
		distance %= 1;
		// TODO CHECK CAST (ITERATE OVER GROUP)
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
		return new PVector((float) coord.x, (float) coord.y);
	}

	public static PVector getCentroid(PShape shape) {
		Point point = fromPShape(shape).getCentroid();
		return new PVector((float) point.getX(), (float) point.getY());
	}

	/**
	 * Returns a hole-less version of the shape, or boundary of group of shapes
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape boundary(PShape shape) {
		return toPShape(fromPShape(shape).getBoundary());
	}

	/**
	 * Get N random points contained within the PShape region. Points are
	 * distributed randomly.
	 * 
	 * @param shape
	 * @param points number of points to generate
	 * @return
	 */
	public static ArrayList<PVector> generateRandomPoints(PShape shape, int points) {
		RandomPointsBuilder r = new RandomPointsBuilder();
		r.setExtent(fromPShape(shape));
		r.setNumPoints(points);

		ArrayList<PVector> vertices = new ArrayList<>();

		for (Coordinate coord : r.getGeometry().getCoordinates()) {
			vertices.add(new PVector((float) coord.x, (float) coord.y));
		}
		return vertices;
	}

	/**
	 * Random points generated in a grid of cells (one point randomly located in
	 * each cell) from the envelope of the shape
	 * 
	 * @param shape
	 * @param maxPoints           number of points, if this shape was its own
	 *                            envelope
	 * @param constrainedToCircle
	 * 
	 *                            Sets whether generated points are constrained to
	 *                            liewithin a circle contained within each grid
	 *                            cell.This provides greater separation between
	 *                            pointsin adjacent cells.
	 * @param gutterFraction      Sets the fraction of the grid cell side which will
	 *                            be treated asa gutter, in which no points will be
	 *                            created.The provided value is clamped to the range
	 *                            [0.0, 1.0].
	 * 
	 * @return
	 */
	public static ArrayList<PVector> generateRandomGridPoints(PShape shape, int maxPoints, boolean constrainedToCircle,
			double gutterFraction) {
		Geometry g = fromPShape(shape);
		IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);

		RandomPointsInGridBuilder r = new RandomPointsInGridBuilder();
		r.setConstrainedToCircle(constrainedToCircle);
		r.setExtent(g.getEnvelopeInternal());
		r.setNumPoints(maxPoints);
		r.setGutterFraction(gutterFraction);

		ArrayList<PVector> vertices = new ArrayList<>();

		for (Coordinate coord : r.getGeometry().getCoordinates()) {
			if (pointLocator.locate(coord) != Location.EXTERIOR) {
				vertices.add(new PVector((float) coord.x, (float) coord.y));
			}
		}
		return vertices;
	}

	/**
	 * Returns a copy of the shape with its small holes (i.e. inner rings with area
	 * < given threshold) removed.
	 * 
	 * @param polygon
	 * @return
	 */
	public static PShape removeSmallHoles(PShape shape, float area) {
		Polygon polygon = (Polygon) fromPShape(shape);
		Polygon noHolePol = GEOM_FACTORY.createPolygon(polygon.getExteriorRing());
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			LinearRing hole = polygon.getInteriorRingN(i);
			if (hole.getArea() < area)
				continue;
			noHolePol = (Polygon) noHolePol.difference(hole);

		}
		return toPShape(noHolePol);
	}

	/**
	 * 
	 * @param n    number of vertices
	 * @param xMax
	 * @param yMax
	 * @return
	 */
	public static PShape randomPolygon(int n, float xMax, float yMax) {
		return pshapeFromPVector(RandomPolygon.generateRandomConvexPolygon(n, xMax, yMax));
	}

	/**
	 * Simple polygon from coordinate list
	 */
	private static PShape pshapeFromPVector(List<PVector> coords) {
		PShape shape = new PShape();
		shape.setFamily(PShape.GEOMETRY);
//		lines.setStrokeCap(ROUND);
		shape.setStroke(true);
		shape.setStrokeWeight(2);
		shape.setStroke(0);
		shape.setFill(-123712);
		shape.setFill(true);
		shape.beginShape();

		for (PVector v : coords) {
			shape.vertex(v.x, v.y);
		}

		shape.endShape(CLOSE);
		return shape;
	}

	/**
	 * TODO return list of points when shape is group
	 * 
	 * @param shape
	 * @param point
	 * @return
	 */
	public static PVector closestVertexToPoint(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		Coordinate coord = DistanceOp.nearestPoints(g, pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

	// TODO split these and similar methods into "Overlay" class?

	public static boolean contains(PShape outer, PShape inner) {
		return fromPShape(outer).contains(fromPShape(inner));
	}

	public static boolean containsPoint(PShape shape, PVector point) {
		return fromPShape(shape).contains(pointFromPVector(point));
	}

	public static float distance(PShape a, PShape b) {
		return (float) fromPShape(a).distance(fromPShape(b));
	}

	public static float area(PShape shape) {
		return (float) fromPShape(shape).getArea();
	}

	/**
	 * * Calculates the Miller circularity index for a polygon. This index, between
	 * 0 and 1, is equal to 1 if the polygon is perfectly circular and tends towards
	 * 0 for a segment.
	 * 
	 * @param shape
	 * @return
	 */
	public static float getCircularity(PShape shape) {
		Polygon poly = (Polygon) fromPShape(shape);
		return (float) (4 * PApplet.PI * poly.getArea()
				/ (poly.getBoundary().getLength() * poly.getBoundary().getLength()));
	}

	/**
	 * 
	 * 
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
	 * Same as getVertices for geometry PShapes; different (subdivides) for circles
	 * etc.
	 * 
	 * @param shape
	 * @return
	 */
	public static PVector[] vertices(PShape shape) {
		Coordinate[] coords = fromPShape(shape).getCoordinates();
		PVector[] vertices = new PVector[coords.length];
		for (int i = 0; i < coords.length; i++) {
			Coordinate coord = coords[i];
			vertices[i] = new PVector((float) coord.x, (float) coord.y);
		}
		return vertices;
	}

	/**
	 * aka envelope
	 * 
	 * @param shape
	 * @return float[] of [X,Y,W,H]
	 */
	public static float[] bound(PShape shape) {
		Envelope e = (Envelope) fromPShape(shape).getEnvelopeInternal();
		return new float[] { (float) e.getMinX(), (float) e.getMinY(), (float) e.getWidth(), (float) e.getHeight() };
	}

	/**
	 * aka envelope
	 * 
	 * @param shape
	 * @return float[] of [X1, Y1, X2, Y2]
	 */
	public static float[] boundCoords(PShape shape) {
		Envelope e = (Envelope) fromPShape(shape).getEnvelopeInternal();
		return new float[] { (float) e.getMinX(), (float) e.getMinY(), (float) e.getMaxX(), (float) e.getMaxY() };
	}

	/**
	 * 
	 * @param x      centre X
	 * @param y      center Y
	 * @param width  ?total? width
	 * @param height ?total? height
	 * @return
	 */
	public static PShape createSquircle(float x, float y, float width, float height) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(40);
		shapeFactory.setCentre(new Coordinate(x, y));
//		shapeFactory.setBase(new Coordinate(x, y));
//		shapeFactory.setRotation(1);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createSquircle());
	}

	public static PShape createArcPolygon(float x, float y, float width, float height) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(40);
		shapeFactory.setCentre(new Coordinate(x, y));
//		shapeFactory.setBase(new Coordinate(x, y));
//		shapeFactory.setRotation(1);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createArcPolygon(1, 2));
	}

	/**
	 * Splits a polygon into 4 equal quadrants
	 * 
	 * @param p
	 * @return
	 */
	public static ArrayList<PShape> split(PShape shape) {
		// https://stackoverflow.com/questions/64252638/how-to-split-a-jts-polygon
		Geometry p = fromPShape(shape);
		ArrayList<PShape> ret = new ArrayList<>();

		final Envelope envelope = p.getEnvelopeInternal();
		double minX = envelope.getMinX();
		double maxX = envelope.getMaxX();
		double midX = minX + (maxX - minX) / 2.0;
		double minY = envelope.getMinY();
		double maxY = envelope.getMaxY();
		double midY = minY + (maxY - minY) / 2.0;

		Envelope llEnv = new Envelope(minX, midX, minY, midY);
		Envelope lrEnv = new Envelope(midX, maxX, minY, midY);
		Envelope ulEnv = new Envelope(minX, midX, midY, maxY);
		Envelope urEnv = new Envelope(midX, maxX, midY, maxY);
		Geometry ll = JTS.toGeometry(llEnv).intersection(p);
		Geometry lr = JTS.toGeometry(lrEnv).intersection(p);
		Geometry ul = JTS.toGeometry(ulEnv).intersection(p);
		Geometry ur = JTS.toGeometry(urEnv).intersection(p);
		ret.add(toPShape(ll));
		ret.add(toPShape(lr));
		ret.add(toPShape(ul));
		ret.add(toPShape(ur));

		return ret;
	}

	/**
	 * Bayazit convex decomposition algorithm for simple polygons.
	 * 
	 * @param shape
	 * @return lsit of convex polygons comprising the original shape
	 */
	public static ArrayList<PShape> decompose(PShape shape) {
		// retry GreedyPolygonSplitter()?

		Geometry g = fromPShape(shape);

		ArrayList<PShape> out = new ArrayList<>();

		for (int i = 0; i < g.getNumGeometries(); i++) {
			Geometry child = g.getGeometryN(i);
			if (child.getGeometryType() == Geometry.TYPENAME_POLYGON) { // skip any linestrings etc
				List<Polygon> decomposed = PolygonDecomposition.decompose((Polygon) child);
				for (Polygon polygon : decomposed) {
					out.add(toPShape(polygon));
				}
			}
		}

		return out;
	}

	/**
	 * Slice a shape using a line given by its start and endpoints
	 * 
	 * @param shape
	 * @param p1    must be outside shape
	 * @param p2    must be outside shape
	 * @return
	 */
	public static ArrayList<PShape> slice(PShape shape, PVector p1, PVector p2) {
		// https://gis.stackexchange.com/questions/189976/
		Geometry poly = fromPShape(shape);
		LineString line = lineFromPVectors(p1, p2);
		Geometry nodedLinework = poly.getBoundary().union(line);
		Geometry polys = polygonize(nodedLinework);

		// Only keep polygons which are inside the input
		ArrayList<Polygon> slices = new ArrayList<>();
		for (int i = 0; i < polys.getNumGeometries(); i++) {
			Polygon candpoly = (Polygon) polys.getGeometryN(i);
			// TODO geometry cache for contains check?
			if (poly.contains(candpoly.getInteriorPoint())) {
				slices.add(candpoly);
			}
		}

		ArrayList<PShape> output = new ArrayList<>();

		for (Polygon polygon : slices) {
			output.add(toPShape(polygon));
		}
		return output;
	}

	/**
	 * Used by slice()
	 */
	@SuppressWarnings("unchecked")
	private static Geometry polygonize(Geometry geometry) {
		List<LineString> lines = LineStringExtracter.getLines(geometry);
		Polygonizer polygonizer = new Polygonizer();
		polygonizer.add(lines);
		Collection<Polygon> polys = polygonizer.getPolygons();
		Polygon[] polyArray = GeometryFactory.toPolygonArray(polys);
		return geometry.getFactory().createGeometryCollection(polyArray);
	}

	private static double distance(Coordinate a, Coordinate b) {
		double deltaX = a.y - b.y;
		double deltaY = a.x - b.x;
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	private static double distance(Point a, Point b) {
		double deltaX = a.getY() - b.getY();
		double deltaY = a.getX() - b.getX();
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	private static void removeCollinearVertices(Geometry g) {
		JTS.removeCollinearVertices(g);
	}

	public static Point createPoint(float x, float y) {
		return GEOM_FACTORY.createPoint(new Coordinate(x, y));
	}

	protected static Point pointFromPVector(PVector p) {
		return GEOM_FACTORY.createPoint(new Coordinate(p.x, p.y));
	}

	protected static Coordinate coordFromPoint(Point p) {
		return new Coordinate(p.getX(), p.getY());
	}

	/**
	 * cirumcircle center must lie inside triangle
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	private static double smallestSide(Coordinate a, Coordinate b, Coordinate c) {
		double ab = Math.sqrt((b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x));
		double bc = Math.sqrt((c.y - b.y) * (c.y - b.y) + (c.x - b.x) * (c.x - b.x));
		double ca = Math.sqrt((a.y - c.y) * (a.y - c.y) + (a.x - c.x) * (a.x - c.x));
		return Math.min(Math.min(ab, bc), ca);
	}

	private static LineString lineFromPVectors(PVector a, PVector b) {
		return GEOM_FACTORY.createLineString(new Coordinate[] { new Coordinate(a.x, a.y), new Coordinate(b.x, b.y) });
	}

	/**
	 * A more convenient way to iterate over both the exterior and linear rings (if
	 * any) of a JTS geometry.
	 * 
	 * @author Michael Carleton
	 *
	 * @param <LinearRing>
	 */
	protected static class LinearRingIterator implements Iterable<LinearRing> {

		private LinearRing[] array;
		private int size;

		public LinearRingIterator(Geometry g) {
			Polygon poly = (Polygon) g;
			this.size = 1 + poly.getNumInteriorRing();
			this.array = new LinearRing[size];
			array[0] = poly.getExteriorRing();
			for (int i = 0; i < poly.getNumInteriorRing(); i++) {
				array[i + 1] = poly.getInteriorRingN(i);
			}
		}

		@Override
		public Iterator<LinearRing> iterator() {
			Iterator<LinearRing> it = new Iterator<LinearRing>() {

				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < size;
				}

				@Override
				public LinearRing next() {
					return array[currentIndex++];
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
			return it;
		}
	}

}
