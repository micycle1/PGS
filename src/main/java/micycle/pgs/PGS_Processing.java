package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.SplittableRandom;
import java.util.stream.StreamSupport;

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
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.noding.IteratedNoder;
import org.locationtech.jts.noding.MCIndexSegmentSetMutualIntersector;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.Noder;
import org.locationtech.jts.noding.SegmentIntersectionDetector;
import org.locationtech.jts.noding.SegmentIntersector;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.SegmentStringUtil;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.CascadedPolygonUnion;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.shape.random.RandomPointsInGridBuilder;
import org.tinfour.common.IConstraint;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;

import micycle.balaban.BalabanSolver;
import micycle.balaban.Point;
import micycle.balaban.Segment;
import micycle.pgs.utility.PolygonDecomposition;
import micycle.pgs.utility.SeededRandomPointsInGridBuilder;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.concavehull.ConcaveHull;
import uk.osgb.algorithm.concavehull.TriCheckerChi;

/**
 * Geometry Processing -- methods that process a shape in some way: compute
 * hulls, partition, slice, etc.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Processing {

	private PGS_Processing() {
	}

	/**
	 * Densifies a shape by inserting extra vertices along the line segments
	 * contained in the shape.
	 * 
	 * @param shape
	 * @param distanceTolerance the densification tolerance to use. All line
	 *                          segments in the densified geometry will be no longer
	 *                          than the distance tolerance. The distance tolerance
	 *                          must be positive.
	 */
	public static PShape densify(PShape shape, double distanceTolerance) {
		Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(distanceTolerance);
		d.setValidate(false);
		return toPShape(d.getResultGeometry());
	}

	/**
	 * Extracts a point from the perimeter (exterior) of the shape at a given
	 * fraction around it perimeter.
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

		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing();
		}
		LengthIndexedLine l = new LengthIndexedLine(g);

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

		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing();
		}
		LengthIndexedLine l = new LengthIndexedLine(g);

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
		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing();
		}
		LengthIndexedLine l = new LengthIndexedLine(g);

		if (interPointDistance > l.getEndIndex()) {
			System.err.println("Interpoint length greater than shape length");
			return new ArrayList<>();
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
	 * Computes all points of intersection between the edges of two shapes.
	 * 
	 * @param a one shape
	 * @param b another shape
	 * @return list of all intersecting points represented by PVectors
	 */
	public static List<PVector> shapeIntersection(PShape a, PShape b) {

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
	 * Computes all points of intersection between segments in a set of line
	 * segments. The input set is first processed to remove degenerate segments
	 * (does not mutate the input).
	 * 
	 * @param lineSegments a list of PVectors where each pair (couplet) of PVectors
	 *                     represent the start and end point of one line segment
	 * @return A list of PVectors each representing the intersection point of a
	 *         segment pair
	 */
	public static List<PVector> lineSegmentsIntersection(List<PVector> lineSegments) {
		final List<PVector> intersections = new ArrayList<>();
		if (lineSegments.size() % 2 != 0) {
			System.err.println(
					"The input to lineSegmentsIntersection() contained an odd number of line segment vertices. The method expects successive pairs of vertices");
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
	 * Generates N random points that lie within the shape region.
	 * <p>
	 * This is a fast method but note that the underlying algorithm makes a minor
	 * trade-off for its speed: the resulting point set is slightly more uniformly
	 * distributed over the input shape compared to a purely random approach (this
	 * arises because the shape is first divided into triangles; each triangle is
	 * then sampled a <b>fixed</b> number of times according to its area relative to
	 * the whole).
	 * 
	 * @param shape  defines the region in which random points are generated
	 * @param points number of points to generate within the shape region
	 * @return
	 * @see #generateRandomPoints(PShape, int, long)
	 * @see #generateRandomGridPoints(PShape, int, boolean, double)
	 */
	public static List<PVector> generateRandomPoints(PShape shape, int points) {
		return generateRandomPoints(shape, points, System.currentTimeMillis());
	}

	/**
	 * Generates N random points that are contained within the shape region. Points
	 * are distributed completely randomly. This method accepts a seed for the RNG
	 * when identical sequences of random points are required.
	 * <p>
	 * This is a fast method but note that the underlying algorithm makes a minor
	 * trade-off for its speed: the resulting point set is slightly more uniformly
	 * distributed over the input shape compared to a purely random approach (this
	 * arises because the shape is first divided into triangles; each triangle is
	 * then sampled a <b>fixed</b> number of times according to its area relative to
	 * the whole).
	 * 
	 * @param shape  defines the region in which random points are generated
	 * @param points number of points to generate within the shape region
	 * @param seed   number used to initialize the underlying pseudorandom number
	 *               generator
	 * @return
	 * @since 1.1.0
	 * @see #generateRandomPoints(PShape, int)
	 * @see #generateRandomGridPoints(PShape, int, boolean, double)
	 */
	public static List<PVector> generateRandomPoints(PShape shape, int points, long seed) {
		final ArrayList<PVector> randomPoints = new ArrayList<>(points); // random points out

		final IncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape, null, true, 0, false);
		final double totalArea = StreamSupport.stream(tin.getConstraints().spliterator(), false)
				.mapToDouble(c -> ((PolygonConstraint) c).getArea()).sum();

		// use arrays to hold variables (to enable assignment during consumer)
		final SimpleTriangle[] largestTriangle = new SimpleTriangle[1];
		final double[] largestArea = new double[1];

		final SplittableRandom r = new SplittableRandom(seed);

		TriangleCollector.visitSimpleTriangles(tin, triangle -> {
			final IConstraint constraint = triangle.getContainingRegion();
			if (constraint != null && constraint.definesConstrainedRegion()) {
				final Vertex a = triangle.getVertexA();
				final Vertex b = triangle.getVertexB();
				final Vertex c = triangle.getVertexC();

				final double triangleArea = 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
				if (triangleArea > largestArea[0]) {
					largestTriangle[0] = triangle;
					largestArea[0] = triangleArea;
				}

				/*
				 * Rather than choose a random triangle for each sample, pre-determine the
				 * number of samples per triangle and sample this number of points in each
				 * triangle successively. I conjecture that this results in a slightly more
				 * uniform random distribution, the downside of which is the resulting
				 * distribution has less entropy.
				 */
				double areaWeight = (triangleArea / totalArea) * points;
				int samples = (int) Math.round(areaWeight);
				if (r.nextDouble() <= (areaWeight - samples)) {
					samples += 1;
				}
				for (int i = 0; i < samples; i++) {
					final double s = r.nextDouble();
					final double t = Math.sqrt(r.nextDouble());
					final double rX = (1 - t) * a.x + t * ((1 - s) * b.x + s * c.x);
					final double rY = (1 - t) * a.y + t * ((1 - s) * b.y + s * c.y);
					randomPoints.add(new PVector((float) rX, (float) rY));
				}
			}
		});

		final int remaining = points - randomPoints.size(); // due to rounding, may be a few above/below target number
		if (remaining > 0) {
			final Vertex a = largestTriangle[0].getVertexA();
			final Vertex b = largestTriangle[0].getVertexB();
			final Vertex c = largestTriangle[0].getVertexC();
			for (int i = 0; i < remaining; i++) {
				double s = r.nextDouble();
				double t = Math.sqrt(r.nextDouble());
				double rX = (1 - t) * a.x + t * ((1 - s) * b.x + s * c.x);
				double rY = (1 - t) * a.y + t * ((1 - s) * b.y + s * c.y);
				randomPoints.add(new PVector((float) rX, (float) rY));
			}
		} else if (remaining < 0) {
			Collections.shuffle(randomPoints, new Random(seed)); // shuffle so that points are removed from regions randomly
			return randomPoints.subList(0, points);
		}

		return randomPoints;
	}

	/**
	 * Generates up to <code>maxPoints</code> number of random points that are
	 * contained within the shape region. Points are distributed according to a grid
	 * of cells (one point randomly located in each cell), based on the envelope of
	 * the shape.
	 * 
	 * @param shape               defines the region in which random points are
	 *                            generated
	 * @param maxPoints           maximum number of points, if this shape was its
	 *                            own envelope
	 * @param constrainedToCircle Sets whether generated points are constrained to
	 *                            lie within a circle contained within each grid
	 *                            cell. This provides greater separation between
	 *                            points in adjacent cells.
	 * @param gutterFraction      Sets the fraction of the grid cell side which will
	 *                            be treated as a gutter, in which no points will be
	 *                            created. The provided value is clamped to the
	 *                            range [0.0, 1.0]. Higher values increase how
	 *                            "grid-like" the point distribution is.
	 * 
	 * @return a list of random points, distributed according to a grid, contained
	 *         within the given shape
	 * @see #generateRandomPoints(PShape, int)
	 * @see #generateRandomGridPoints(PShape, int, boolean, double, long)
	 */
	public static List<PVector> generateRandomGridPoints(PShape shape, int maxPoints, boolean constrainedToCircle, double gutterFraction) {
		return generateRandomGridPoints(shape, maxPoints, constrainedToCircle, gutterFraction, System.currentTimeMillis());
	}

	/**
	 * Generates up to <code>maxPoints</code> number of random points that are
	 * contained within the shape region. Points are distributed according to a grid
	 * of cells (one point randomly located in each cell), based on the envelope of
	 * the shape.
	 * <p>
	 * This method accepts a seed for the RNG when identical sequences of random
	 * points are required.
	 * 
	 * @param shape               defines the region in which random points are
	 *                            generated
	 * @param maxPoints           maximum number of points, if this shape was its
	 *                            own envelope
	 * @param constrainedToCircle Sets whether generated points are constrained to
	 *                            lie within a circle contained within each grid
	 *                            cell. This provides greater separation between
	 *                            points in adjacent cells.
	 * @param gutterFraction      Sets the fraction of the grid cell side which will
	 *                            be treated as a gutter, in which no points will be
	 *                            created. The provided value is clamped to the
	 *                            range [0.0, 1.0]. Higher values increase how
	 *                            "grid-like" the point distribution is.
	 * @param randomSeed
	 * @return a list of random points, distributed according to a grid, contained
	 *         within the given shape
	 * @see #generateRandomGridPoints(PShape, int, boolean, double)
	 */
	public static List<PVector> generateRandomGridPoints(PShape shape, int maxPoints, boolean constrainedToCircle, double gutterFraction,
			long randomSeed) {
		Geometry g = fromPShape(shape);
		IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);

		RandomPointsInGridBuilder r = new SeededRandomPointsInGridBuilder(randomSeed);
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
	 * Returns a copy of the shape where small holes (i.e. inner rings with area <
	 * given threshold) are removed.
	 * 
	 * @param shape
	 * @param areaThreshold
	 * @return
	 */
	public static PShape removeSmallHoles(PShape shape, double areaThreshold) {
		Polygon polygon = (Polygon) fromPShape(shape);
		Polygon noHolePol = PGS.GEOM_FACTORY.createPolygon(polygon.getExteriorRing());
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			LinearRing hole = polygon.getInteriorRingN(i);
			if (hole.getArea() < areaThreshold) {
				continue;
			}
			noHolePol = (Polygon) noHolePol.difference(hole);
		}
		return toPShape(noHolePol);
	}

	/**
	 * Computes the polygonal faces formed by a set of intersecting line segments.
	 * 
	 * @param lineSegmentVertices a list of PVectors where each pair (couplet) of
	 *                            PVectors represent the start and end point of one
	 *                            line segment
	 * @return a list of polygonal PShapes each representing a face / enclosed area
	 *         formed between intersecting lines
	 * @since 1.1.2
	 */
	@SuppressWarnings("unchecked")
	public static PShape polygonizeLines(List<PVector> lineSegmentVertices) {
		// TODO constructor for LINES PShape
		if (lineSegmentVertices.size() % 2 != 0) {
			System.err.println(
					"The input to polygonizeLines() contained an odd number of vertices. The method expects successive pairs of vertices.");
			return new PShape();
		}

		List<SegmentString> segmentStrings = new ArrayList<>(lineSegmentVertices.size() / 2);
		for (int i = 0; i < lineSegmentVertices.size(); i += 2) {
			final PVector v1 = lineSegmentVertices.get(i);
			final PVector v2 = lineSegmentVertices.get(i + 1);
			segmentStrings.add(new NodedSegmentString(new Coordinate[] { PGS.coordFromPVector(v1), PGS.coordFromPVector(v2) }, null));
		}

		final Polygonizer polygonizer = new Polygonizer();
		polygonizer.setCheckRingsValid(false);
		final Noder noder = new IteratedNoder(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));
		noder.computeNodes(segmentStrings);
		noder.getNodedSubstrings().forEach(s -> {
			SegmentString ss = (SegmentString) s;
			polygonizer.add(PGS.GEOM_FACTORY.createLineString(new Coordinate[] { ss.getCoordinate(0), ss.getCoordinate(1) }));
		});
		Collection<Geometry> polygons = polygonizer.getPolygons();

		final PShape out = new PShape(PConstants.GROUP);
		polygons.forEach(p -> out.addChild(toPShape(p)));
		return out;
	}

	/**
	 * Computes the convex hull of multiple shapes.
	 * 
	 * @param shapes
	 * @return
	 * @see #convexHull(PShape...)
	 */
	public static PShape convexHull(List<PShape> shapes) {
		Collection<Polygon> polygons = new ArrayList<>();
		shapes.forEach(s -> polygons.add((Polygon) fromPShape(s)));
		return toPShape(CascadedPolygonUnion.union(polygons).convexHull());
	}

	/**
	 * Computes the convex hull of multiple shapes.
	 * 
	 * @param shapes varArgs
	 * @return
	 * @see #convexHull(List)
	 */
	public static PShape convexHull(PShape... shapes) {
		return convexHull(Arrays.asList(shapes));
	}

	/**
	 * Computes the concave hull of a point set using a breadth-first method.
	 * 
	 * @param points
	 * @param threshold euclidean distance threshold
	 * @return
	 * @since 1.1.0
	 * @see #concaveHullDFS(List, double)
	 * @see #concaveHull2(List, double)
	 */
	public static PShape concaveHullBFS(List<PVector> points, double threshold) {
		ConcaveHull hull = new ConcaveHull(prepareConcaveGeometry(points));
		return toPShape(hull.getConcaveHullBFS(new TriCheckerChi(threshold), false, false).get(0));
	}

	/**
	 * Computes the concave hull of a point set using a depth-first method. In
	 * contrast to the BFS method, the depth-first approach produces shapes that are
	 * more contiguous/less branching and spiral-like.
	 * 
	 * @param points
	 * @param threshold euclidean distance threshold
	 * @return
	 * @since 1.1.0
	 * @see #concaveHullBFS(List, double)
	 * @see #concaveHull2(List, double)
	 */
	public static PShape concaveHullDFS(List<PVector> points, double threshold) {
		ConcaveHull hull = new ConcaveHull(prepareConcaveGeometry(points));
		return toPShape(hull.getConcaveHullDFS(new TriCheckerChi(threshold)));
	}

	/**
	 * Computes the concave hull of a point set using a different algorithm. This
	 * approach has a more "organic" structure compared to other concaveBFS method.
	 * 
	 * @param points
	 * @param threshold 0...1 (Normalized length parameter). Setting threshold=1
	 *                  means that no edges will be removed from the Delaunay
	 *                  triangulation, so the resulting polygon will be the convex
	 *                  hull. Setting threshold=0 means that all edges that can be
	 *                  removed subject to the regularity constraint will be removed
	 *                  (however polygons that are eroded beyond the point where
	 *                  they provide a desirable characterization of the shape).
	 *                  Although the optimal parameter value varies for different
	 *                  shapes and point distributions, values of between 0.05–0.2
	 *                  typically produce optimal or near-optimal shape
	 *                  characterization across a wide range of point distributions.
	 * @return
	 * @see #concaveHullBFS(List, double)
	 */
	public static PShape concaveHull2(List<PVector> points, double threshold) {

		/*
		 * (from https://doi.org/10.1016/j.patcog.2008.03.023) It is more convenient to
		 * normalize the threshold parameter with respect to a particular set of points
		 * P by using the maximum and minimum edge lengths of the Delaunay triangulation
		 * of P. Increasing l beyond the maximum edge length of the Delaunay
		 * triangulation cannot reduce the number of edges that will be removed (which
		 * will be zero anyway). Decreasing l beyond the minimum edge length of the
		 * Delaunay triangulation cannot increase the number of edges that will be
		 * removed.
		 */

		org.geodelivery.jap.concavehull.ConcaveHull hull = new org.geodelivery.jap.concavehull.ConcaveHull(threshold);

		return toPShape(hull.transform(prepareConcaveGeometry(points)));
	}

	/**
	 * Prepares a multipoint geometry from a list of PVectors.
	 */
	private static Geometry prepareConcaveGeometry(List<PVector> points) {
		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
		} else { // already closed
			coords = new Coordinate[points.size()];
		}

		for (int i = 0; i < coords.length; i++) {
			if (i >= points.size()) {
				coords[i] = new Coordinate(points.get(0).x, points.get(0).y); // close geometry
			} else {
				coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
			}
		}

		return PGS.GEOM_FACTORY.createMultiPointFromCoords(coords);
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
	 * @return a GROUP PShape, where each child shape is some quadrant partition of
	 *         the original shape
	 * @see #split(PShape, int)
	 */
	public static PShape split(PShape shape) {
		// https://stackoverflow.com/questions/64252638/how-to-split-a-jts-polygon
		Geometry p = fromPShape(shape);

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

		final PShape partitions = new PShape(PConstants.GROUP);
		partitions.addChild(toPShape(UL));
		partitions.addChild(toPShape(UR));
		partitions.addChild(toPShape(LL));
		partitions.addChild(toPShape(LR));

		return partitions;
	}

	/**
	 * Splits a shape into 4^(1+recursions) rectangular partitions
	 * 
	 * @param shape
	 * @param splitDepth
	 * @return a GROUP PShape, where each child shape is some quadrant partition of
	 *         the original shape
	 * @see #split(PShape)
	 */
	public static PShape split(final PShape shape, int splitDepth) {
		splitDepth = Math.max(0, splitDepth);
		ArrayDeque<Geometry> stack = new ArrayDeque<>();
		stack.add(fromPShape(shape));

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

		final PShape partitions = new PShape(PConstants.GROUP);
		stack.forEach(g -> partitions.addChild(toPShape(g)));
		return partitions;
	}

	/**
	 * Partitions a shape into simple polygons using Mark Bayazit's algorithm.
	 * 
	 * @param shape
	 * @return a GROUP PShape, where each child shape is some convex partition of
	 *         the original shape
	 */
	public static PShape partition(PShape shape) {
		// https://mpen.ca/406/bayazit
		final Geometry g = fromPShape(shape);

		final PShape partitions = new PShape(PConstants.GROUP);
		for (int i = 0; i < g.getNumGeometries(); i++) {
			Geometry child = g.getGeometryN(i);
			if (child.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) { // skip any linestrings etc
				List<Polygon> decomposed = PolygonDecomposition.decompose((Polygon) child);
				for (Polygon polygon : decomposed) {
					partitions.addChild(toPShape(polygon));
				}
			}
		}

		return partitions;
	}

	/**
	 * Slices a shape using a line given by its start and endpoints.
	 * 
	 * @param shape PShape to slice into two shapes
	 * @param p1    must be outside shape
	 * @param p2    must be outside shape
	 * @return a GROUP PShape with two children, where each child shape one of the
	 *         slices
	 */
	public static PShape slice(PShape shape, PVector p1, PVector p2) {
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

		final PShape slices = new PShape(PConstants.GROUP);
		slices.addChild(toPShape(UnaryUnionOp.union(leftSlices)));
		slices.addChild(toPShape(UnaryUnionOp.union(rightSlices)));
		return slices;
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
