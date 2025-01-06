package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.awt.geom.Rectangle2D;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.SplittableRandom;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.locationtech.jts.algorithm.Angle;
import org.locationtech.jts.algorithm.Area;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.hull.ConcaveHullOfPolygons;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.noding.MCIndexSegmentSetMutualIntersector;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.Noder;
import org.locationtech.jts.noding.SegmentIntersectionDetector;
import org.locationtech.jts.noding.SegmentIntersector;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.SegmentStringUtil;
import org.locationtech.jts.noding.snap.SnappingNoder;
import org.locationtech.jts.operation.overlay.snap.GeometrySnapper;
import org.locationtech.jts.operation.overlayng.MultiOperationOverlayNG;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.shape.random.RandomPointsInGridBuilder;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.utils.TriangleCollector;
import org.tinfour.voronoi.BoundedVoronoiBuildOptions;
import org.tinfour.voronoi.BoundedVoronoiDiagram;

import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import it.unimi.dsi.util.XoRoShiRo128PlusRandomGenerator;
import micycle.balaban.BalabanSolver;
import micycle.balaban.Point;
import micycle.balaban.Segment;
import micycle.pgs.color.ColorUtils;
import micycle.pgs.color.Colors;
import micycle.pgs.commons.PolygonDecomposition;
import micycle.pgs.commons.SeededRandomPointsInGridBuilder;
import micycle.trapmap.TrapMap;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Methods that process shape geometry: partitioning, slicing, cleaning, etc.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Processing {

	private PGS_Processing() {
	}

	/**
	 * Densifies a shape by inserting additional vertices along its line segments.
	 * This results in a more detailed representation of the shape, affecting its
	 * structure without altering its geometry.
	 * 
	 * @param shape             The shape to be densified. It should be a lineal or
	 *                          polygonal shape.
	 * @param distanceTolerance The densification tolerance used. All line segments
	 *                          in the densified geometry will have lengths no
	 *                          longer than the specified distance tolerance. The
	 *                          distance tolerance must be a positive value.
	 * @return A new PShape object representing the densified geometry.
	 */
	public static PShape densify(PShape shape, double distanceTolerance) {
		Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(distanceTolerance);
		d.setValidate(false);
		return toPShape(d.getResultGeometry());
	}

	/**
	 * Extracts a point from the perimeter (exterior) of the given shape at a
	 * specific position along its perimeter.
	 * 
	 * @param shape             A lineal or polygonal shape. If the input is a GROUP
	 *                          shape, a single point will be extracted from its
	 *                          first child shape.
	 * @param perimeterPosition A normalised position along the perimeter of a shape
	 *                          [0...1]. 0 corresponds to the starting point of the
	 *                          shape's perimeter, and 1 corresponds to the ending
	 *                          point of the perimeter; any value between 0 and 1
	 *                          represents a proportional distance along the shape's
	 *                          boundary.
	 * @param offsetDistance    A perpendicular offset distance from the shape's
	 *                          perimeter. A value of 0 corresponds to exactly on
	 *                          the shape's exterior/ Positive values offset the
	 *                          point away from the shape (outwards); negative
	 *                          values offset the point inwards towards its
	 *                          interior.
	 * @return
	 * @see #pointsOnExterior(PShape, int, double)
	 */
	public static PVector pointOnExterior(PShape shape, double perimeterPosition, double offsetDistance) {
		perimeterPosition %= 1;
		LengthIndexedLine l = makeIndexedLine(shape);

		Coordinate coord = l.extractPoint(perimeterPosition * l.getEndIndex(), offsetDistance);
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Extracts a point from the perimeter (exterior) of the given shape at some
	 * distance along its perimeter.
	 * 
	 * @param shape             A lineal or polygonal shape. If the input is a GROUP
	 *                          shape, a single point will be extracted from its
	 *                          first child shape.
	 * @param perimeterDistance Distance along shape perimeter to extract the point.
	 *                          0 corresponds to the first vertex of the shape's
	 *                          perimeter.
	 * @param offsetDistance    A perpendicular offset distance from the shape's
	 *                          perimeter. A value of 0 corresponds to exactly on
	 *                          the shape's exterior/ Positive values offset the
	 *                          point away from the shape (outwards); negative
	 *                          values offset the point inwards towards its
	 *                          interior.
	 * @return
	 * @since 1.4.0
	 */
	public static PVector pointOnExteriorByDistance(PShape shape, double perimeterDistance, double offsetDistance) {
		LengthIndexedLine l = makeIndexedLine(shape);
		Coordinate coord = l.extractPoint(perimeterDistance % l.getEndIndex(), offsetDistance);
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Extracts multiple points evenly distributed along the boundary of individual
	 * rings within a shape, including both exterior and interior rings (i.e.,
	 * holes).
	 * <p>
	 * This method enhances the performance of boundary sampling by directly
	 * extracting points from each linear ring of the shape's polygons. It is more
	 * efficient than multiple individual point extractions, and it supports
	 * sampling from holes (interior rings) as well.
	 * 
	 * @param shape          The shape from which to extract points. Should have
	 *                       polygonal members.
	 * @param points         The number of points to extract <b>per ring</b>, evenly
	 *                       distributed around each ring's boundary.
	 * @param offsetDistance The offset distance measured perpendicular to each
	 *                       point on the ring's boundary. Positive values offset
	 *                       outwards, while negative values offset inwards.
	 * @return A list of PVector objects, each representing a point on the perimeter
	 *         or interior rings of the shape.
	 * @see #pointOnExterior(PShape, double, double)
	 * @see #pointsOnExterior(PShape, double, double)
	 * @since 1.3.0
	 */
	public static List<PVector> pointsOnExterior(PShape shape, int points, double offsetDistance) {
		// TODO another method that returns concave hull of returned points (when
		// offset)
		List<Polygon> polygons = PGS.extractPolygons(fromPShape(shape));
		List<PVector> coords = new ArrayList<>();
		polygons.forEach(polygon -> {
			PGS.extractLinearRings(polygon).forEach(ring -> {
				if (Orientation.isCCW(ring.getCoordinates())) {
					ring = ring.reverse();
				}
				final LengthIndexedLine l = new LengthIndexedLine(ring);
				final double increment = 1d / points;
				for (double distance = 0; distance < 1; distance += increment) {
					final Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
					coords.add(PGS.toPVector(coord));
				}
			});
		});

		return coords;
	}

	/**
	 * Generates a list of evenly distributed points along the boundary of each ring
	 * within a given polygonal shape, which may include its exterior and any
	 * interior rings (holes).
	 * <p>
	 * This method is used to obtain a set of points that approximate the polygonal
	 * shape's outline and interior boundaries, with the points spaced at
	 * approximate equal intervals determined by the <code>interPointDistance</code>
	 * parameter. It supports complex shapes with interior rings (holes) by
	 * extracting points from all rings.
	 * 
	 * @param shape              The shape from which to generate the points. It
	 *                           should be a polygonal shape.
	 * @param interPointDistance The desired distance between consecutive points
	 *                           along each ring's boundary. This controls the
	 *                           spacing of points and the granularity of the
	 *                           representation.
	 * @param offsetDistance     The offset distance perpendicular to each point on
	 *                           the ring's boundary. Positive values offset points
	 *                           outwards, while negative values bring them towards
	 *                           the interior.
	 * @return A list of PVector objects representing the points along the exterior
	 *         and interior boundaries of the shape.
	 * @see #pointOnExterior(PShape, double, double)
	 * @see #densify(PShape, double)
	 * @since 1.3.0
	 */
	public static List<PVector> pointsOnExterior(PShape shape, double interPointDistance, double offsetDistance) {
		List<Polygon> polygons = PGS.extractPolygons(fromPShape(shape));
		List<PVector> coords = new ArrayList<>();
		polygons.forEach(polygon -> {
			PGS.extractLinearRings(polygon).forEach(ring -> {
				if (Orientation.isCCW(ring.getCoordinates())) {
					ring = ring.reverse();
				}
				final LengthIndexedLine l = new LengthIndexedLine(ring);
				final int points = (int) Math.round(l.getEndIndex() / interPointDistance);
				final double increment = 1d / points;
				for (double distance = 0; distance < 1; distance += increment) {
					final Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
					coords.add(PGS.toPVector(coord));
				}

			});
		});

		return coords;
	}

	/**
	 * Extracts evenly spaced dashed line segments along the perimeter of a shape.
	 * This method ensures that the segments are distributed uniformly along the
	 * shape's boundary, with the possibility of adjusting the start position of the
	 * first line based on an offset.
	 * 
	 * @param shape             The shape from which to extract the segments.
	 * @param lineLength        The length of each segment. Must be a positive
	 *                          number.
	 * @param interLineDistance The distance between the end of one segment and the
	 *                          start of the next. Must be non-negative.
	 * @param offset            The starting position offset (around the perimeter
	 *                          [0...1]) for the first line. Values > |1| loop
	 *                          around the shape. Positive values indicate a
	 *                          clockwise (CW) direction, and negative values
	 *                          indicate a counter-clockwise (CCW) direction.
	 * @return A GROUP PShape whose children are the extracted segments.
	 * @since 2.0
	 */
	public static PShape segmentsOnExterior(PShape shape, double lineLength, double interLineDistance, double offset) {
		LengthIndexedLine l = makeIndexedLine(shape);
		lineLength = Math.max(lineLength, 0.5);
		interLineDistance = Math.max(interLineDistance, 0); // ensure >= 0

		final double perimeter = l.getEndIndex();

		double totalSegmentLength = lineLength + interLineDistance;
		int numberOfLines = (int) Math.floor(perimeter / totalSegmentLength);

		double adjustmentFactor = perimeter / (numberOfLines * totalSegmentLength);
		lineLength *= adjustmentFactor;
		interLineDistance *= adjustmentFactor;
		totalSegmentLength = lineLength + interLineDistance;

		offset = -offset; // positive values should wrap CW
		double startingPosition = ((offset % 1.0) + 1.0) % 1.0 * perimeter;

		List<PShape> lines = new ArrayList<>(numberOfLines);

		for (int i = 0; i < numberOfLines; i++) {
			double lineStart = (startingPosition + i * totalSegmentLength) % perimeter;
			double lineEnd = (lineStart + lineLength) % perimeter;

			Geometry segment;
			PShape line;
			if (lineStart < lineEnd) {
				segment = l.extractLine(lineStart, lineEnd);
				line = PGS_Conversion.toPShape(segment);
				line.setName(String.valueOf(lineStart));
				lines.add(line);
			} else {
				// Handle case where line wraps around the end of the shape (straddles 0)
				// combine 2 segments into a single linestring
				Coordinate[] c1 = l.extractLine(lineStart, perimeter).getCoordinates();
				Coordinate[] c2 = l.extractLine(0, lineEnd).getCoordinates();
				segment = PGS.GEOM_FACTORY.createLineString(ArrayUtils.addAll(c1, c2));
				line = PGS_Conversion.toPShape(segment);
				line.setName(String.valueOf(lineStart));
				lines.add(line);
			}
		}

		return PGS_Conversion.flatten(lines);
	}

	/**
	 * Creates an CW-oriented length-indexed line from a given PShape.
	 */
	private static LengthIndexedLine makeIndexedLine(PShape shape) {
		Geometry g = fromPShape(shape);
		if (g instanceof Polygonal) {
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)) {
				g = g.getGeometryN(0);
			}
			LinearRing e = ((Polygon) g).getExteriorRing();
			if (Orientation.isCCW(e.getCoordinates())) {
				e = e.reverse();
			}
			g = e;
		}
		return new LengthIndexedLine(g);
	}

	/**
	 * Extracts a portion/subline of the perimeter of a shape between two locations
	 * on the perimeter.
	 * 
	 * @param shape the shape from which to extract a the perimeter
	 * @param from  the starting location of the perimeter extract, given by a
	 *              fraction (0...1) of the total perimeter length
	 * @param to    the end location of the perimeter extract, given by a fraction
	 *              (0...1) of the total perimeter length
	 * @return lineal shape
	 * @since 1.2.0
	 */
	public static PShape extractPerimeter(PShape shape, double from, double to) {
		from = floatMod(from, 1);
		if (to != 1) { // so that value of 1 is not moduloed to equal 0
			to = floatMod(to, 1);
		}
		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing();
		}
		final LengthIndexedLine l = new LengthIndexedLine(g);
		final double length = l.getEndIndex(); // perimeter length

		if (from > to) {
			final Geometry l1 = l.extractLine(length * from, length);
			final Geometry l2 = l.extractLine(0, length * to);
			return toPShape(GEOM_FACTORY
					.createLineString(Stream.concat(Arrays.stream(l1.getCoordinates()), Arrays.stream(l2.getCoordinates())).toArray(Coordinate[]::new)));
		}

		/*
		 * The PGS toPShape() method treats a closed linestring as polygonal (having a
		 * fill), which occurs when from==0 and to==1. We don't want the output to be
		 * filled in, so build the PATH shape here without closing it.
		 */
		LineString string = (LineString) l.extractLine(length * from, length * to);
		PShape perimeter = new PShape();
		perimeter.setFamily(PShape.PATH);
		perimeter.setStroke(true);
		perimeter.setStroke(micycle.pgs.color.Colors.PINK);
		perimeter.setStrokeWeight(4);

		perimeter.beginShape();
		Coordinate[] coords = string.getCoordinates();
		for (int i = 0; i < coords.length; i++) {
			perimeter.vertex((float) coords[i].x, (float) coords[i].y);
		}
		perimeter.endShape();

		return perimeter;
	}

	/**
	 * Calculates the angle of the line tangent to the shape at a specific point on
	 * its perimeter. The position of the point is determined by the normalized
	 * distance along the shape's perimeter.
	 * 
	 * @param shape          The shape for which to find the tangent angle.
	 * @param perimeterRatio A normalised position along the perimeter of a shape
	 *                       [0...1]. 0 corresponds to the starting point of the
	 *                       shape's perimeter, and 1 corresponds to the ending
	 *                       point of the perimeter; any value between 0 and 1
	 *                       represents a proportional distance along the shape's
	 *                       boundary.
	 * @return the normalized angle (in radians) that the tangent line at the
	 *         specified position makes with the positive x-axis (east), orientated
	 *         clockwise.
	 * @since 1.3.0
	 */
	public static double tangentAngle(PShape shape, double perimeterRatio) {
		perimeterRatio %= 1;

		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			LinearRing e = ((Polygon) g).getExteriorRing();
			if (Orientation.isCCW(e.getCoordinates())) {
				e = e.reverse();
			}
			g = e;
		}
		LengthIndexedLine l = new LengthIndexedLine(g);

		final double position = perimeterRatio * l.getEndIndex();
		final double d1 = position - 1e-5;
		final double d2 = position + 1e-5;

		Coordinate coordA = l.extractPoint(d1); // CCW - point behind first
		Coordinate coordB = l.extractPoint(d2); // CCW - point after second

		double angle = Angle.angle(coordA, coordB); // CCW
		angle += Math.PI * 2; // [1PI...3PI]
		angle %= (Math.PI * 2); // [0...2PI], CW

		return angle;
	}

	/**
	 * Computes all <b>points</b> of intersection between the <b>perimeters</b> of
	 * two shapes.
	 * <p>
	 * NOTE: This method shouldn't be confused with
	 * {@link micycle.pgs.PGS_ShapeBoolean#intersect(PShape, PShape)
	 * PGS_ShapeBoolean.intersect()}, which finds the shape made by the intersecting
	 * shape areas.
	 * 
	 * @param a one shape
	 * @param b another shape
	 * @return list of all intersecting points (as PVectors)
	 */
	public static List<PVector> shapeIntersection(PShape a, PShape b) {
		final HashSet<PVector> points = new HashSet<>();

		final Collection<?> segmentStrings = SegmentStringUtil.extractSegmentStrings(fromPShape(a));
		final MCIndexSegmentSetMutualIntersector mci = new MCIndexSegmentSetMutualIntersector(segmentStrings);
		final SegmentIntersectionDetector sid = new SegmentIntersectionDetector();

		mci.process(SegmentStringUtil.extractSegmentStrings(fromPShape(b)), new SegmentIntersector() {
			@Override
			public void processIntersections(SegmentString e0, int segIndex0, SegmentString e1, int segIndex1) {
				sid.processIntersections(e0, segIndex0, e1, segIndex1);
				if (sid.hasIntersection()) {
					points.add(new PVector((float) sid.getIntersection().x, (float) sid.getIntersection().y));
				}
			}

			@Override
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
		return generateRandomPoints(shape, points, System.nanoTime());
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

		final IIncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape);
		final boolean constrained = !tin.getConstraints().isEmpty();
		final double totalArea = StreamSupport.stream(tin.getConstraints().spliterator(), false).mapToDouble(c -> ((PolygonConstraint) c).getArea()).sum();

		// use arrays to hold variables (to enable assignment during consumer)
		final SimpleTriangle[] largestTriangle = new SimpleTriangle[1];
		final double[] largestArea = new double[1];

		final SplittableRandom r = new SplittableRandom(seed);
		TriangleCollector.visitSimpleTriangles(tin, triangle -> {
			final IConstraint constraint = triangle.getContainingRegion();
			if (!constrained || (constraint != null && constraint.definesConstrainedRegion())) {
				final Vertex a = triangle.getVertexA();
				final Vertex b = triangle.getVertexB();
				final Vertex c = triangle.getVertexC();

				// TODO more robust area (dense input produces slivers)
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
			Collections.shuffle(randomPoints, new XoRoShiRo128PlusRandom(seed)); // shuffle so that points are removed from regions randomly
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
		return generateRandomGridPoints(shape, maxPoints, constrainedToCircle, gutterFraction, System.nanoTime());
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
	public static List<PVector> generateRandomGridPoints(PShape shape, int maxPoints, boolean constrainedToCircle, double gutterFraction, long randomSeed) {
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
				vertices.add(PGS.toPVector(coord));
			}
		}
		return vertices;
	}

	/**
	 * Creates a nested shape having n levels of inner polygons; each inner polygon
	 * is obtained by joining the points at some fractional distance <code>r</code>
	 * along each side of the previous polygon.
	 * <p>
	 * This is also called a "derived polygon".
	 *
	 * @param shape a PShape representing the single input polygon for which the
	 *              nested shape will be created. Although all polygon types are
	 *              compatible, simpler/convex polygons make for a more effective
	 *              result.
	 * @param n     the number of nested levels to be created inside the input
	 *              shape.
	 * @param r     fractional distance between each nested level of polygons, where
	 *              <code>0.5</code> produces maximal nesting.
	 * @return A PShape representing the nested shape, including the input shape and
	 *         all the nested levels.
	 * @since 1.4.0
	 */
	public static PShape nest(PShape shape, int n, double r) {
		final Polygon polygon = (Polygon) fromPShape(shape);

		if (r != 1) {
			r %= 1;
		}
		final Polygon[] derivedPolygons = new Polygon[n + 1];
		derivedPolygons[0] = polygon;
		Polygon currentPolygon = polygon;

		for (int i = 0; i < n; i++) {
			Coordinate[] inputCoordinates = currentPolygon.getCoordinates();
			int numVertices = inputCoordinates.length - 1;

			Coordinate[] derivedCoordinates = new Coordinate[numVertices + 1];

			for (int k = 0; k < numVertices; k++) {
				double x = inputCoordinates[k].x * (1 - r) + inputCoordinates[(k + 1) % numVertices].x * r;
				double y = inputCoordinates[k].y * (1 - r) + inputCoordinates[(k + 1) % numVertices].y * r;
				derivedCoordinates[k] = new Coordinate(x, y);
			}
			derivedCoordinates[numVertices] = derivedCoordinates[0]; // close the ring

			Polygon derivedPolygon = PGS.GEOM_FACTORY.createPolygon(derivedCoordinates);

			derivedPolygons[i + 1] = derivedPolygon;
			currentPolygon = derivedPolygon;
		}

		return toPShape(GEOM_FACTORY.createMultiPolygon(derivedPolygons));
	}

	/**
	 * Removes overlap between polygons contained in a <code>GROUP</code> shape,
	 * preserving only visible line segments suitable for pen plotting and similar
	 * applications.
	 * <p>
	 * This method processes a <code>GROUP</code> shape consisting of lineal or
	 * polygonal child shapes, aiming to create linework that represents only the
	 * segments visible to a human, rather than a computer. The resulting linework
	 * is useful for pen plotters or other applications where only the visible paths
	 * are desired.
	 * <p>
	 * During the operation, any overlapping lines are also removed to ensure a
	 * clean and clear representation of the shapes. It's important to note that the
	 * order of shape layers in the input GROUP shape is significant. The method
	 * considers the last child shape of the input to be "on top" of all other
	 * shapes, as is the case visually.
	 * 
	 * @param shape A GROUP shape containing lineal or polygonal child shapes.
	 * @return The resulting linework of the overlapping input as a LINES PShape,
	 *         representing only visible line segments.
	 * @since 1.3.0
	 */
	public static PShape removeHiddenLines(PShape shape) {
		if (shape.getChildCount() == 0) {
			return shape;
		}

		List<PShape> layers = PGS_Conversion.getChildren(shape); // visual top last
		Collections.reverse(layers); // visual top first
		final List<Geometry> geometries = layers.stream().map(PGS_Conversion::fromPShape).collect(Collectors.toList());
		Geometry union = geometries.get(0); // start of cascading union

		List<Geometry> culledGeometries = new ArrayList<>(geometries.size());
		Iterator<Geometry> i = geometries.iterator();
		culledGeometries.add(i.next());

		// for each shape, subtract the union of shapes visually above it
		while (i.hasNext()) {
			final Geometry layer = i.next();
			MultiOperationOverlayNG overlay = new MultiOperationOverlayNG(layer, union);
			Geometry occulted = overlay.getResult(OverlayNG.DIFFERENCE); // occulted version of layer
			union = overlay.getResult(OverlayNG.UNION);

			culledGeometries.add(occulted);
		}

		Geometry dissolved = LineDissolver.dissolve(GEOM_FACTORY.createGeometryCollection(culledGeometries.toArray(new Geometry[0])));
		PShape out = toPShape(dissolved);
		PGS_Conversion.setAllStrokeColor(out, ColorUtils.setAlpha(Colors.PINK, 192), 4);

		return out;
	}

	/**
	 * Returns a copy of the shape where holes having an area <b>less than</b> the
	 * specified threshold are removed.
	 * 
	 * @param shape         a single polygonal shape or GROUP polygonal shape
	 * @param areaThreshold removes any holes with an area smaller than this value
	 * @return
	 */
	public static PShape removeSmallHoles(PShape shape, double areaThreshold) {
		final Geometry g = fromPShape(shape);
		@SuppressWarnings("unchecked")
		final List<Polygon> polygons = PolygonExtracter.getPolygons(g);
		return toPShape(polygons.stream().map(p -> removeSmallHoles(p, areaThreshold)).collect(Collectors.toList()));
	}

	private static Polygon removeSmallHoles(Polygon polygon, double areaThreshold) {
		LinearRing noHolePol = polygon.getExteriorRing();
		List<LinearRing> validHoles = new ArrayList<>();
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			final LinearRing hole = polygon.getInteriorRingN(i);
			double holeArea = Area.ofRing(polygon.getInteriorRingN(i).getCoordinates());
			if (holeArea < areaThreshold) {
				continue;
			}
			validHoles.add(hole);
		}
		return GEOM_FACTORY.createPolygon(noHolePol, GeometryFactory.toLinearRingArray(validHoles));
	}

	/**
	 * Extracts all the holes from a shape, returning them as if they are polygons.
	 * 
	 * @param shape the PShape to extract holes from
	 * @return a new PShape that represents the holes extracted from the input
	 *         shape. If the input had multiple holes, the output is a GROUP PShape
	 *         where each child is polygon of one hole.
	 * @since 1.4.0
	 */
	public static PShape extractHoles(PShape shape) {
		final Geometry g = fromPShape(shape);
		@SuppressWarnings("unchecked")
		final List<Polygon> polygons = PolygonExtracter.getPolygons(g);

		List<PShape> holes = new ArrayList<>(polygons.size());
		for (Polygon polygon : polygons) {
			for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
				final LinearRing hole = polygon.getInteriorRingN(i);
				holes.add(toPShape(hole));
			}
		}

		return PGS_Conversion.fromChildren(holes);
	}

	/**
	 * Finds the polygonal faces formed by a set of intersecting line segments.
	 * 
	 * @param lineSegmentVertices a list of PVectors where each pair (couplet) of
	 *                            PVectors represent the start and end point of one
	 *                            line segment
	 * @return a GROUP PShape where each child shape is a face / enclosed area
	 *         formed between intersecting lines
	 * @since 1.1.2
	 */
	public static PShape polygonizeLines(List<PVector> lineSegmentVertices) {
		// TODO constructor for LINES PShape
		if (lineSegmentVertices.size() % 2 != 0) {
			System.err.println("The input to polygonizeLines() contained an odd number of vertices. The method expects successive pairs of vertices.");
			return new PShape();
		}

		final List<SegmentString> segmentStrings = new ArrayList<>(lineSegmentVertices.size() / 2);
		for (int i = 0; i < lineSegmentVertices.size(); i += 2) {
			final PVector v1 = lineSegmentVertices.get(i);
			final PVector v2 = lineSegmentVertices.get(i + 1);
			if (!v1.equals(v2)) {
				segmentStrings.add(new NodedSegmentString(new Coordinate[] { PGS.coordFromPVector(v1), PGS.coordFromPVector(v2) }, null));
			}
		}

		return PGS.polygonizeSegments(segmentStrings, true);
	}

	/**
	 * Splits a shape into four equal quadrants (based on the envelope of the
	 * shape).
	 * 
	 * @param shape The shape to be split into quadrants.
	 * @return A GROUP PShape where each child shape represents one of the four
	 *         quadrants partitioned from the original shape.
	 * @see #split(PShape, int)
	 */
	public static PShape split(PShape shape) {
		return split(shape, 1);
	}

	/**
	 * Splits a shape into <code>4^(1+recursions)</code> rectangular partitions.
	 * <p>
	 * This method performs a recursive split operation on the given shape, dividing
	 * it into a finer grid of rectangular partitions. Each recursion further
	 * divides the envelopes of the parent shapes into four quadrants, resulting in
	 * a more refined partitioning of the original shape.
	 * <p>
	 * It's important to note that this operation differs from merely overlaying a
	 * regular square grid on the shape and then splitting it. The recursive process
	 * ensures that each subdivision is based on the envelope of the parent shape,
	 * which may be rectangular.
	 * 
	 * @param shape      The shape to be split into rectangular partitions.
	 * @param splitDepth The number of split recursions to perform, determining the
	 *                   level of partitioning and grid refinement.
	 * @return A GROUP PShape, where each child shape represents one of the quadrant
	 *         partitions of the original shape.
	 * @see #split(PShape)
	 */
	public static PShape split(final PShape shape, int splitDepth) {
		// https://stackoverflow.com/questions/64252638/how-to-split-a-jts-polygon
		splitDepth = Math.max(0, splitDepth);
		Deque<Geometry> stack = new ArrayDeque<>();
		stack.add(fromPShape(shape));

		List<Geometry> next = new ArrayList<>(); // add to when recursion depth reached

		Noder noder = new SnappingNoder(1e-11);
		int depth = 0;
		while (depth < splitDepth) {
			while (!stack.isEmpty()) {
				final Geometry slice = GeometryFixer.fix(stack.pop());
				final Envelope envelope = slice.getEnvelopeInternal();

				final double minX = envelope.getMinX();
				final double maxX = envelope.getMaxX();
				final double midX = minX + envelope.getWidth() / 2.0;
				final double minY = envelope.getMinY();
				final double maxY = envelope.getMaxY();
				final double midY = minY + envelope.getHeight() / 2.0;

				Envelope llEnv = new Envelope(minX, midX, minY, midY);
				Envelope lrEnv = new Envelope(midX, maxX, minY, midY);
				Envelope ulEnv = new Envelope(minX, midX, midY, maxY);
				Envelope urEnv = new Envelope(midX, maxX, midY, maxY);

				Geometry UL = OverlayNG.overlay(slice, toGeometry(ulEnv), OverlayNG.INTERSECTION, noder);
				Geometry UR = OverlayNG.overlay(slice, toGeometry(urEnv), OverlayNG.INTERSECTION, noder);
				Geometry LL = OverlayNG.overlay(slice, toGeometry(llEnv), OverlayNG.INTERSECTION, noder);
				Geometry LR = OverlayNG.overlay(slice, toGeometry(lrEnv), OverlayNG.INTERSECTION, noder);

				// Geometries may not be polygonal; in which case, do not include in output.
				if (UL instanceof Polygonal && !UL.isEmpty()) {
					next.add(UL);
				}
				if (UR instanceof Polygonal && !UR.isEmpty()) {
					next.add(UR);
				}
				if (LL instanceof Polygonal && !LL.isEmpty()) {
					next.add(LL);
				}
				if (LR instanceof Polygonal && !LR.isEmpty()) {
					next.add(LR);
				}

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
	 * Partitions shape(s) into convex (simple) polygons.
	 * 
	 * @param shape the shape to partition. can be a single polygon or a GROUP of
	 *              polygons
	 * @return a GROUP PShape, where each child shape is some convex partition of
	 *         the original shape
	 */
	public static PShape convexPartition(PShape shape) {
		// algorithm described in https://mpen.ca/406/bayazit
		final Geometry g = fromPShape(shape);

		final PShape polyPartitions = new PShape(PConstants.GROUP);
		@SuppressWarnings("unchecked")
		final List<Polygon> polygons = PolygonExtracter.getPolygons(g);
		polygons.forEach(p -> polyPartitions.addChild(toPShape(PolygonDecomposition.decompose(p))));

		if (polyPartitions.getChildCount() == 1) {
			return polyPartitions.getChild(0);
		} else {
			return polyPartitions;
		}
	}

	/**
	 * Randomly partitions a shape into N approximately equal-area polygonal cells.
	 * 
	 * @param shape a polygonal (non-group, no holes) shape to partition
	 * @param parts number of roughly equal area partitons to create
	 * @return a GROUP PShape, whose child shapes are partitions of the original
	 * @since 1.3.0
	 */
	public static PShape equalPartition(final PShape shape, final int parts) {
		return equalPartition(shape, parts, System.nanoTime());
	}

	/**
	 * Randomly (with a given seed) partitions a shape into N approximately
	 * equal-area polygonal cells.
	 * 
	 * @param shape a polygonal (non-group, no holes) shape to partition
	 * @param parts number of roughly equal area partitons to create
	 * @param seed  number used to initialize the underlying pseudorandom number
	 *              generator
	 * @return a GROUP PShape, whose child shapes are partitions of the original
	 * @since 1.4.0
	 */
	@SuppressWarnings("unchecked")
	public static PShape equalPartition(final PShape shape, final int parts, long seed) {
		final Geometry g = fromPShape(shape);
		if (g instanceof Polygonal) {
			final Envelope e = g.getEnvelopeInternal();

			int samples = (int) (e.getArea() / 100); // sample every ~10 units in x and y axes

			final List<PVector> samplePoints = PGS_Processing.generateRandomGridPoints(shape, samples, false, 0.8, seed);
			KMeansPlusPlusClusterer<Clusterable> clusterer = new KMeansPlusPlusClusterer<>(parts, 20, new EuclideanDistance(),
					new XoRoShiRo128PlusRandomGenerator(seed));
			List<? extends Clusterable> cl = samplePoints.stream().map(p -> (Clusterable) () -> new double[] { p.x, p.y }).collect(Collectors.toList());
			List<CentroidCluster<Clusterable>> clusters = clusterer.cluster((Collection<Clusterable>) cl);
			if (clusters.size() < 3) {
				// since voronoi needs 3+ vertices
				return shape;
			}
			List<Vertex> vertices = new ArrayList<>();
			clusters.forEach(c -> {
				double[] p = c.getCenter().getPoint();
				vertices.add(new Vertex(p[0], p[1], 0));
			});

			double x = e.getMinX();
			double y = e.getMinY();
			double w = e.getMaxX() - e.getMinX();
			double h = e.getMaxY() - e.getMinY();

			final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
			options.setBounds(new Rectangle2D.Double(x, y, w, h));

			final BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(vertices, options);
			List<Geometry> faces = v.getPolygons().stream().map(PGS_Voronoi::toPolygon).collect(Collectors.toList());

			faces = faces.parallelStream().map(f -> OverlayNG.overlay(f, g, OverlayNG.INTERSECTION)).collect(Collectors.toList());

			return PGS_Conversion.toPShape(faces);
		} else {
			System.err.println("equalPartition(): Input shape is not polygonal.");
			return shape;
		}
	}

	/**
	 * Decomposes/partitions a shape into axis-aligned (stip-like) trazepoids.
	 * <p>
	 * The output can contain some "degenerate" trapezoids that do indeed have 4
	 * vertices but look like triangles.
	 * 
	 * @param shape a polygonal or a GROUP shape
	 * @return a GROUP PShape comprising of trapezoid child shapes
	 * @since 1.3.0
	 */
	public static PShape trapezoidPartition(PShape shape) {
		final PShape trapezoids = new PShape(PConstants.GROUP);

		TrapMap map;
		try {
			map = new TrapMap(PGS_Conversion.getChildren(shape));
		} catch (Exception e) {
			/*
			 * Handle error thrown by TrapMap on degenerate/strange inputs. Generally
			 * shearing will fix the problem (ideally this would be done within TrapMap).
			 */
			try {
				map = new TrapMap(PGS_Conversion.getChildren(PGS_Transformation.shear(shape, .00001, 0)));
			} catch (Exception e2) {
				System.err.println(e.getLocalizedMessage());
				return trapezoids;
			}
		}

		final IndexedPointInAreaLocator locator = new IndexedPointInAreaLocator(PGS_Conversion.fromPShape(shape));
		map.getAllTrapezoids().forEach(t -> {
			if (t.getFace() != null) {
				/*
				 * Some spurious trapezoids that lie outside the shape boundary are not filtered
				 * out by 't.getFace() != null' (most are). For those that remain use a
				 * trapezoid's centroid to determine whether it makes up the shape.
				 */
				final PShape tz = t.getBoundaryPolygon();
				final PVector centroid = new PVector(0, 0);
				// don't need to use all 4 vertices
				centroid.add(tz.getVertex(0));
				centroid.add(tz.getVertex(1));
				centroid.add(tz.getVertex(2));
				centroid.div(3);
				if (locator.locate(PGS.coordFromPVector(centroid)) != Location.EXTERIOR) {
					trapezoids.addChild(tz); // include if not outside
				}
			}
		});

		PGS_Conversion.setAllFillColor(trapezoids, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(trapezoids, Colors.PINK, 1);
		return trapezoids;
	}

	/**
	 * Slices a shape using a line given by its start and endpoints.
	 * 
	 * @param shape PShape to slice into two shapes
	 * @param p1    must be outside shape
	 * @param p2    must be outside shape
	 * @return a GROUP PShape with two children, where each child shape is one of
	 *         the slices
	 */
	public static PShape slice(PShape shape, PVector p1, PVector p2) {
		// adapted from https://gis.stackexchange.com/questions/189976/
		final Geometry poly = fromPShape(shape);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(poly);
		final LineSegment ls = new LineSegment(p1.x, p1.y, p2.x, p2.y);
		final LineString line = ls.toGeometry(GEOM_FACTORY);
		final Geometry nodedLinework = poly.getBoundary().union(line);
		final Geometry polys = polygonize(nodedLinework);

		final List<Polygon> leftSlices = new ArrayList<>();
		final List<Polygon> rightSlices = new ArrayList<>();

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
	 * Removes narrow areas ("slivers") from a shape while preserving the geometry
	 * of the remaining parts.
	 * <p>
	 * This operation is similar to
	 * {@link PGS_Morphology#erosionDilation(PShape, double) erosionDilation()}, but
	 * better preserves the original geometry of remaining parts.
	 * <p>
	 * If the input is a single polygon and if when removing slivers, a multipolygon
	 * is produced, further processing occurs within the method to repair it back
	 * into a single polygon.
	 * 
	 * @param shape     polygonal shape
	 * @param threshold width threshold (probably no more than 10); parts narrower
	 *                  than this are eliminated
	 * @return a copy of the input shape having narrow areas/parts removed
	 * @since 1.3.0
	 */
	public static PShape eliminateSlivers(PShape shape, double threshold) {
		threshold = Math.max(threshold, 1e-5);
		Geometry g = fromPShape(shape);
		boolean multi = g.getNumGeometries() > 1;
		Geometry snapped = GeometrySnapper.snapToSelf(g, threshold, true);
		Geometry out = snapped;
		if (!multi && snapped.getNumGeometries() > 1) {
			ConcaveHullOfPolygons chp = new ConcaveHullOfPolygons(snapped);
			chp.setTight(true);
			chp.setHolesAllowed(false);

			double ratio = 0.02;
			try {
				do { // pick lowest ratio such that a single polygon is returned
					chp.setMaximumEdgeLengthRatio(ratio);
					out = chp.getHull();
					ratio += 0.01;
				} while (out.getNumGeometries() > 1);
			} catch (Exception e) { // catch "Unable to find a convex corner"
			}
		}
		return toPShape(out); // better on smaller thresholds
	}

	/**
	 * Attempts to fix shapes with invalid geometry, while preserving its original
	 * form and location as much as possible. See
	 * {@link org.locationtech.jts.geom.util.GeometryFixer GeometryFixer} for a full
	 * list of potential fixes.
	 * <p>
	 * Input shapes are always processed, so even valid inputs may have some minor
	 * alterations. The output is always a new geometry object.
	 * 
	 * @param shape The shape to be corrected.
	 * @return A modified version of the input shape that aligns with valid shape
	 *         geometry standards.
	 * @since 2.0
	 */
	public static PShape fix(PShape shape) {
		return toPShape(GeometryFixer.fix(fromPShape(shape)));
	}

	/**
	 * Filters out the children of a given PShape object based on a given Predicate
	 * function. Child shapes are filtered when the predicate is true: <i>"remove
	 * if..."</i>.
	 * <p>
	 * Example Lambda Code:
	 * 
	 * <pre>
	 * {@code
	 * PGS_Processing.filterChildren(myShape, child -> PGS_ShapePredicates.area(child) < 50); // discard small child shapes
	 * }
	 * </pre>
	 * 
	 * @param shape          The PShape whose children will be filtered.
	 * @param filterFunction A Predicate function that takes a PShape object as its
	 *                       argument and returns <code>true</code> if the child
	 *                       shape should be <b>filtered out/removed</b> from the
	 *                       shape. To retain child shapes, the Predicate function
	 *                       should return <code>false</code>. You can use a lambda
	 *                       expression or a method reference to implement the
	 *                       Predicate function.
	 * @since 1.4.0
	 * @return A new PShape that contains only the children shapes of the input
	 *         shape that satisfy the given Predicate function
	 *         (<code>==false</code>).
	 */
	public static PShape filterChildren(PShape shape, Predicate<PShape> filterFunction) {
		filterFunction = filterFunction.negate();
		List<PShape> filteredFaces = PGS_Conversion.getChildren(shape).stream().filter(filterFunction::test).collect(Collectors.toList());
		return PGS_Conversion.flatten(filteredFaces);
	}

	/**
	 * Applies a specified transformation function to each child of the given PShape
	 * and returns a new PShape containing the transformed children.
	 * <p>
	 * This method processes each child of the input shape using the provided
	 * function, which can modify, replace, or filter out shapes. The resulting
	 * transformed shapes are flattened into a new PShape.
	 * <p>
	 * The transformation function can:
	 * <ul>
	 * <li>Modify the shape in-place and return it</li>
	 * <li>Create and return a new shape</li>
	 * <li>Return null to remove the shape</li>
	 * </ul>
	 * <p>
	 * Note: This method creates a new PShape and does not modify the original shape
	 * or its children. The hierarchical structure of the original shape is not
	 * preserved in the result.
	 *
	 * @param shape    The PShape whose children will be transformed.
	 * @param function A UnaryOperator that takes a PShape as input and returns a
	 *                 transformed PShape. If the function returns null for a shape,
	 *                 that shape will be excluded from the result.
	 * @return A new PShape containing the transformed children, flattened into a
	 *         single level.
	 * @see #transformWithIndex(PShape, BiFunction)
	 * @since 2.0
	 */
	public static PShape transform(PShape shape, UnaryOperator<PShape> function) {
		return PGS_Conversion.flatten(PGS_Conversion.getChildren(shape).stream().map(function).filter(Objects::nonNull).toList());
	}

	/**
	 * Applies a specified transformation function to each child of the given
	 * PShape, providing the index of each child, and returns a new PShape
	 * containing the transformed children.
	 * <p>
	 * This method processes each child of the input shape alongside its index using
	 * the provided BiFunction. The function can modify, replace, or filter out
	 * shapes. The resulting transformed shapes are flattened into a new PShape.
	 * <p>
	 * The transformation function can:
	 * <ul>
	 * <li>Modify the shape in-place and return it</li>
	 * <li>Create and return a new shape</li>
	 * <li>Return null to exclude the shape from the final result</li>
	 * </ul>
	 * <p>
	 * Note: This method creates a new PShape and does not modify the original shape
	 * or its children. The hierarchical structure of the original shape is not
	 * preserved in the result.
	 * <p>
	 * Example usage:
	 * 
	 * <pre>{@code
	 * PShape newShape = PGS_Processing.transformWithIndex(originalShape, (index, child) -> {
	 * 	return (index % 2 == 0) ? child : null; // Include only even-indexed children
	 * });
	 * }</pre>
	 *
	 * @param shape    The PShape whose children will be transformed.
	 * @param function A BiFunction that takes an integer index and a PShape as
	 *                 input and returns a transformed PShape. If the function
	 *                 returns null for a shape, that shape will be excluded from
	 *                 the result.
	 * @return A new PShape containing the transformed children, flattened into a
	 *         single level.
	 * @see #transform(PShape, UnaryOperator)
	 * @since 2.0
	 */
	public static PShape transformWithIndex(PShape shape, BiFunction<Integer, PShape, PShape> function) {
		List<PShape> children = PGS_Conversion.getChildren(shape);
		return PGS_Conversion.flatten(IntStream.range(0, children.size()).mapToObj(i -> function.apply(i, children.get(i))).filter(Objects::nonNull).toList());
	}

	/**
	 * Applies a specified function to each child of the given PShape.
	 * <p>
	 * This method iterates over each child of the input PShape, applying the
	 * provided Consumer function to each. The function can perform any operation on
	 * the shapes, such as modifying their properties or applying effects, but does
	 * not inherently alter the structure of the PShape or its hierarchy.
	 * <p>
	 * The changes are made in place; hence, the original PShape is modified, and
	 * the same reference is returned for convenience in chaining or further use.
	 * <p>
	 * Example usage:
	 * 
	 * <pre>{@code
	 * PShape shape = PGS_Processing.apply(shape, child -> child.setFill(false));
	 * }</pre>
	 *
	 * @param shape         The PShape whose children will be processed.
	 * @param applyFunction A Consumer that takes a PShape as input and performs
	 *                      operations on it.
	 * @return The original PShape with the function applied to each child.
	 * @see #applyWithIndex(PShape, BiConsumer)
	 * @since 2.0
	 */
	public static PShape apply(PShape shape, Consumer<PShape> applyFunction) {
		for (PShape child : PGS_Conversion.getChildren(shape)) {
			applyFunction.accept(child);
		}
		return shape;
	}

	/**
	 * Applies a specified function to each child of the given PShape, providing the
	 * index of each child.
	 * <p>
	 * This method iterates over each child of the input PShape and applies the
	 * provided BiConsumer function to each child along with its index. The function
	 * can perform any operation on the shapes such as modifying properties or
	 * applying effects. The index can be used for operations that depend on the
	 * position of the child within its parent.
	 * <p>
	 * The changes are made in place; hence, the original PShape is modified. The
	 * method returns the same PShape reference for convenience in chaining or
	 * further use.
	 * <p>
	 * Example usage:
	 * 
	 * <pre>{@code
	 * PShape shape = PGS_Processing.applyWithIndex(shape, (index, child) -> {
	 * 	if (index % 2 == 0)
	 * 		child.setFill(true);
	 * });
	 * }</pre>
	 *
	 * @param shape         The PShape whose children will be processed.
	 * @param applyFunction A BiConsumer that takes an integer index and a PShape as
	 *                      input and performs operations on it.
	 * @return The original PShape with the function applied to each indexed child.
	 * @see #apply(PShape, Consumer)
	 * @since 2.0
	 */
	public static PShape applyWithIndex(PShape shape, BiConsumer<Integer, PShape> applyFunction) {
		List<PShape> children = PGS_Conversion.getChildren(shape);
		for (int i = 0; i < children.size(); i++) {
			applyFunction.accept(i, children.get(i));
		}
		return shape;
	}

	private static Polygon toGeometry(Envelope envelope) {
		return GEOM_FACTORY.createPolygon(GEOM_FACTORY.createLinearRing(new Coordinate[] { new Coordinate(envelope.getMinX(), envelope.getMinY()),
				new Coordinate(envelope.getMaxX(), envelope.getMinY()), new Coordinate(envelope.getMaxX(), envelope.getMaxY()),
				new Coordinate(envelope.getMinX(), envelope.getMaxY()), new Coordinate(envelope.getMinX(), envelope.getMinY()) }));
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

	private static double floatMod(double x, double y) {
		// x mod y behaving the same way as Math.floorMod but with doubles
		return (x - Math.floor(x / y) * y);
	}

}
