package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.geodelivery.jap.concavehull.SnapHull;
import org.geotools.geometry.jts.JTS;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.noding.MCIndexSegmentSetMutualIntersector;
import org.locationtech.jts.noding.SegmentIntersectionDetector;
import org.locationtech.jts.noding.SegmentIntersector;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.SegmentStringUtil;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.shape.random.RandomPointsBuilder;
import org.locationtech.jts.shape.random.RandomPointsInGridBuilder;

import micycle.balaban.BalabanSolver;
import micycle.balaban.Point;
import micycle.balaban.Segment;
import micycle.pgs.utility.PolygonDecomposition;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.concavehull.ConcaveHull;
import uk.osgb.algorithm.concavehull.TriCheckerChi;

/**
 * Geometry Processing.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Processing {

	private PGS_Processing() {
	}

	public static PShape envelope(PShape shape) {
		return toPShape(fromPShape(shape).getEnvelope());
	}

	/**
	 * Densifies a Geometry by inserting extra vertices along the line segments
	 * contained in the geometry. Specify maximum length of segments
	 */
	public static PShape densify(PShape shape, double distanceTolerance) {
		Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(distanceTolerance);
		d.setValidate(false);
		return toPShape(d.getResultGeometry());
	}

	/**
	 * Extracts a point from the perimeter (exterior) of the given shape.
	 * 
	 * @param shape
	 * @param distance       0...1 around shape perimeter; or -1...0 (other
	 *                       direction)
	 * @param offsetDistance perpendicular offset distance, where 0 is exactly on
	 *                       the shape exteriod. Positive values offset the point
	 *                       away from the shape (outwards); negative values offset
	 *                       the point inwards.
	 * @return
	 * @see #pointsOnExterior(PShape, int, double)
	 */
	public static PVector pointOnExterior(PShape shape, double distance, double offsetDistance) {
		distance %= 1;
		// TODO CHECK CAST (ITERATE OVER GROUP); apply to interior rings too?
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Extracts many points from the perimeter (faster than calling other method
	 * lots)
	 * 
	 * @param shape
	 * @param points         number of points to return; evenly distibuted around
	 *                       the perimeter of the shape
	 * @param offsetDistance offset distance along a line perpendicular to the
	 *                       perimeter
	 * @return
	 * @see #pointOnExterior(PShape, double, double)
	 * @see #pointsOnExterior(PShape, double, double)
	 */
	public static List<PVector> pointsOnExterior(PShape shape, int points, double offsetDistance) {
		// TODO another method that returns concave hull of returned points (when
		// offset)
		ArrayList<PVector> coords = new ArrayList<>(points);
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		final double increment = 1d / points;
		for (double distance = 0; distance < 1; distance += increment) {
			Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
			coords.add(new PVector((float) coord.x, (float) coord.y));
		}
		return coords;
	}

	/**
	 * Generates a list of points that lie on the exterior/perimeter of the given
	 * shape.
	 * 
	 * @param shape
	 * @param interPointDistance distance between each exterior point
	 * @param offsetDistance
	 * @return
	 */
	public static List<PVector> pointsOnExterior(PShape shape, double interPointDistance, double offsetDistance) {
		// TODO points on holes
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		if (interPointDistance > l.getEndIndex()) {
			System.err.println("Interpoint length greater than shape length");
			return new ArrayList<PVector>();
		}
		final int points = (int) Math.round(l.getEndIndex() / interPointDistance);

		ArrayList<PVector> coords = new ArrayList<>(points);

		final double increment = 1d / points;
		for (double distance = 0; distance < 1; distance += increment) {
			Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
			coords.add(new PVector((float) coord.x, (float) coord.y));
		}
		return coords;
	}

	/**
	 * Calculate all the points of intersection between two shapes.
	 * 
	 * @param a one shape
	 * @param b the other shape
	 * @return list of all intersecting points represented by PVectors
	 */
	public static List<PVector> shapeIntersections(PShape a, PShape b) {

		final HashSet<PVector> points = new HashSet<>();

		final MCIndexSegmentSetMutualIntersector mci = new MCIndexSegmentSetMutualIntersector(
				SegmentStringUtil.extractSegmentStrings(fromPShape(a)));
		final SegmentIntersectionDetector sid = new SegmentIntersectionDetector();

		mci.process(SegmentStringUtil.extractSegmentStrings(fromPShape(b)), new SegmentIntersector() {
			public void processIntersections(SegmentString e0, int segIndex0, SegmentString e1, int segIndex1) {
				sid.processIntersections(e0, segIndex0, e1, segIndex1);
				if (sid.getIntersection() != null) {
					points.add(new PVector((float) sid.getIntersection().x, (float) sid.getIntersection().y));
				}
			}
			public boolean isDone() {
				return false;
			}
		});
		return new ArrayList<>(points);
	}

	/**
	 * Computes all points of intersection between segment pairs from a set of
	 * segments. The input set is first processed to remove degenerate segments
	 * (does not mutate the input).
	 * 
	 * @param lineSegments a list of PVectors where each pair (couplet) of PVectors
	 *                     represent the start and end point of one segment
	 * @return A list of PVectors each representing the intersection point of a
	 *         segment pair
	 */
	public static List<PVector> lineSegmentIntersections(List<PVector> lineSegments) {
		final List<PVector> intersections = new ArrayList<PVector>();
		if (lineSegments.size() % 2 != 0) {
			System.err.println("Error: detected an odd number of line segment vertices.");
			return intersections;
		}

		Collection<Segment> segments = new ArrayList<>();
		for (int i = 0; i < lineSegments.size(); i += 2) { // iterate pairwise
			final PVector p1 = lineSegments.get(i);
			final PVector p2 = lineSegments.get(i + 1);
			segments.add(new Segment(p1.x, p1.y, p2.x, p2.y));
		}

		final BalabanSolver balabanSolver = new BalabanSolver((a, b) -> {
			final Point pX = a.getIntersection(b);
			intersections.add(new PVector((float) pX.x, (float) pX.y));
		});
		segments.removeAll(balabanSolver.findDegenerateSegments(segments));
		balabanSolver.computeIntersections(segments);

		return intersections;
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
	 * @param maxPoints           max number of points, if this shape was its own
	 *                            envelope
	 * @param constrainedToCircle
	 * 
	 *                            Sets whether generated points are constrained to
	 *                            liewithin a circle contained within each grid
	 *                            cell. This provides greater separation between
	 *                            points in adjacent cells.
	 * @param gutterFraction      Sets the fraction of the grid cell side which will
	 *                            be treated as a gutter, in which no points will be
	 *                            created. The provided value is clamped to the
	 *                            range [0.0, 1.0].
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
	 * Returns a copy of the shape with small holes (i.e. inner rings with area <
	 * given threshold) are removed.
	 * 
	 * @param polygon
	 * @return
	 */
	public static PShape removeSmallHoles(PShape shape, double area) {
		Polygon polygon = (Polygon) fromPShape(shape);
		Polygon noHolePol = PGS.GEOM_FACTORY.createPolygon(polygon.getExteriorRing());
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			LinearRing hole = polygon.getInteriorRingN(i);
			if (hole.getArea() < area) {
				continue;
			}
			noHolePol = (Polygon) noHolePol.difference(hole);
		}
		return toPShape(noHolePol);
	}

	/**
	 * Computes the convex hull of multiple PShapes.
	 * 
	 * @param shapes
	 * @return
	 */
	public static PShape convexHull(PShape... shapes) {
		Geometry g = fromPShape(shapes[0]);
		for (int i = 1; i < shapes.length; i++) {
			g = g.union(fromPShape(shapes[i])); // TODO slow for many
		}
		return toPShape(g.convexHull());
	}

	/**
	 * Uses Chi heurstic
	 * 
	 * @param points
	 * @param threshold euclidean distance threshold
	 * @return
	 * @see #concaveHull2(List, double)
	 */
	public static PShape concaveHull(List<PVector> points, double threshold) {

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

		Geometry g = PGS.GEOM_FACTORY.createPolygon(coords);
		ConcaveHull hull = new ConcaveHull(g);
		return toPShape(hull.getConcaveHullBFS(new TriCheckerChi(threshold), false, false).get(0));
	}

	/**
	 * Computes the concave hull of a shape using a different algorithm. Has a more
	 * "organic" structure compared to other concave method.
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
	 * @see #concaveHull(List, double)
	 */
	public static PShape concaveHull2(List<PVector> points, double threshold) {

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

		Geometry g = PGS.GEOM_FACTORY.createPolygon(coords);

		// TODO test AVG threshold heuristic
		org.geodelivery.jap.concavehull.ConcaveHull hull = new org.geodelivery.jap.concavehull.ConcaveHull(threshold);

		return toPShape(hull.transform(g));
	}

	/**
	 * Computes the "snap hull" for a shape, which is a convex hull that snaps to
	 * the shape. Adjust segment factor to change between
	 * 
	 * @param shape
	 * @param segmentFactor default = 4
	 * @return
	 */
	public static PShape snapHull(PShape shape, double segmentFactor) {
		return toPShape(SnapHull.snapHull(fromPShape(shape), segmentFactor));
	}

	/**
	 * Splits a shape into 4 equal quadrants
	 * 
	 * @param shape
	 * @return list containing the 4 split quadrants of the input shape
	 */
	public static List<PShape> split(PShape shape) {
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
		Geometry UL = JTS.toGeometry(llEnv).intersection(p);
		Geometry UR = JTS.toGeometry(lrEnv).intersection(p);
		Geometry LL = JTS.toGeometry(ulEnv).intersection(p);
		Geometry LR = JTS.toGeometry(urEnv).intersection(p);
		ret.add(toPShape(UL));
		ret.add(toPShape(UR));
		ret.add(toPShape(LL));
		ret.add(toPShape(LR));

		return ret;
	}

	/**
	 * Splits a shape into 4^(1+recursions) rectangular partitions
	 * 
	 * @param shape
	 * @param splitDepth
	 */
	public static List<PShape> split(final PShape shape, int splitDepth) {
		splitDepth = Math.max(0, splitDepth);
		ArrayDeque<Geometry> stack = new ArrayDeque<>();
		stack.add(fromPShape(shape));

		ArrayList<PShape> ret = new ArrayList<>(); // add to when recursion depth reached
		ArrayList<Geometry> next = new ArrayList<>(); // add to when recursion depth reached

		int depth = 0;
		while (depth < splitDepth) {
			while (!stack.isEmpty()) {
				final Geometry slice = stack.pop();
				final Envelope envelope = slice.getEnvelopeInternal();
				final double minX = envelope.getMinX();
				final double maxX = envelope.getMaxX();
				final double midX = minX + (maxX - minX) / 2.0;
				final double minY = envelope.getMinY();
				final double maxY = envelope.getMaxY();
				final double midY = minY + (maxY - minY) / 2.0;

				Envelope llEnv = new Envelope(minX, midX, minY, midY);
				Envelope lrEnv = new Envelope(midX, maxX, minY, midY);
				Envelope ulEnv = new Envelope(minX, midX, midY, maxY);
				Envelope urEnv = new Envelope(midX, maxX, midY, maxY);
				Geometry UL = JTS.toGeometry(llEnv).intersection(slice);
				Geometry UR = JTS.toGeometry(lrEnv).intersection(slice);
				Geometry LL = JTS.toGeometry(ulEnv).intersection(slice);
				Geometry LR = JTS.toGeometry(urEnv).intersection(slice);
				next.add(UL);
				next.add(UR);
				next.add(LL);
				next.add(LR);
			}
			depth++;
			stack.addAll(next);
			next.clear();
		}

		stack.forEach(g -> ret.add(toPShape(g)));
		return ret;
	}

	/**
	 * Partitions a shape into simple polygons using Mark Bayazit's algorithm.
	 * 
	 * @param shape
	 * @return list of convex (simple) polygons comprising the original shape
	 */
	public static List<PShape> partition(PShape shape) {
		// https://mpen.ca/406/bayazit
		// retry GreedyPolygonSplitter()?

		Geometry g = fromPShape(shape);

		ArrayList<PShape> out = new ArrayList<>();

		for (int i = 0; i < g.getNumGeometries(); i++) {
			Geometry child = g.getGeometryN(i);
			if (child.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) { // skip any linestrings etc
				List<Polygon> decomposed = PolygonDecomposition.decompose((Polygon) child);
				for (Polygon polygon : decomposed) {
					out.add(toPShape(polygon));
				}
			}
		}

		return out;
	}

	/**
	 * Slice a shape using a line given by its start and endpoints.
	 * 
	 * @param shape
	 * @param p1    must be outside shape
	 * @param p2    must be outside shape
	 * @return a list containg two PShapes
	 */
	public static List<PShape> slice(PShape shape, PVector p1, PVector p2) {
		// adapted from https://gis.stackexchange.com/questions/189976/
		final Geometry poly = fromPShape(shape);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(poly);
		final LineSegment ls = new LineSegment(p1.x, p1.y, p2.x, p2.y);
		final LineString line = ls.toGeometry(PGS.GEOM_FACTORY);
		final Geometry nodedLinework = poly.getBoundary().union(line);
		final Geometry polys = polygonize(nodedLinework);

		final ArrayList<Polygon> leftSlices = new ArrayList<>();
		final ArrayList<Polygon> rightSlices = new ArrayList<>();

		for (int i = 0; i < polys.getNumGeometries(); i++) {
			final Polygon candpoly = (Polygon) polys.getGeometryN(i);
			if (cache.contains(candpoly.getInteriorPoint())) {
				if (ls.orientationIndex(candpoly.getCentroid().getCoordinate()) == Orientation.LEFT) {
					leftSlices.add(candpoly);
				} else {
					rightSlices.add(candpoly);
				}
			}
		}

		ArrayList<PShape> output = new ArrayList<>();
		output.add(toPShape(UnaryUnionOp.union(leftSlices)));
		output.add(toPShape(UnaryUnionOp.union(rightSlices)));
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

}
