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
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

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
import org.tinfour.common.Vertex;
import org.tinfour.voronoi.BoundedVoronoiBuildOptions;
import org.tinfour.voronoi.BoundedVoronoiDiagram;

import com.github.micycle1.geoblitz.IndexedLengthIndexedLine;
import com.github.micycle1.geoblitz.YStripesPointInAreaLocator;

import it.unimi.dsi.util.XoRoShiRo128PlusRandomGenerator;
import micycle.balaban.BalabanSolver;
import micycle.balaban.Point;
import micycle.balaban.Segment;
import micycle.pgs.color.ColorUtils;
import micycle.pgs.color.Colors;
import micycle.pgs.commons.PolygonDecomposition;
import micycle.pgs.commons.SeededRandomPointsInGridBuilder;
import micycle.pgs.commons.ShapeRandomPointSampler;
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
	 * @see #pointsOnExterior(PShape, int, double)
	 */
	public static PVector pointOnExterior(PShape shape, double perimeterPosition, double offsetDistance) {
		perimeterPosition %= 1;
		var l = makeIndexedLine(shape);

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
	 * @since 1.4.0
	 */
	public static PVector pointOnExteriorByDistance(PShape shape, double perimeterDistance, double offsetDistance) {
		var l = makeIndexedLine(shape);
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
		return pointsOnExterior(shape, points, offsetDistance, 0);
	}

	/**
	 * Extract a fixed number of evenly spaced samples from every linear component.
	 *
	 * <p>
	 * Samples exactly {@code points} positions from each linear component (each
	 * LinearRing or LineString) found in {@code shape}. Sampling is done
	 * independently per component — i.e., {@code points} samples are produced for
	 * each ring or line. {@code startOffset} sets the fractional start position
	 * along each component (0..1). {@code offsetDistance} is applied perpendicular
	 * to the boundary when extracting each sample.
	 *
	 * <p>
	 * Orientation is normalized per component (reversed if necessary) so offsets
	 * are applied consistently. The method only collects points and does not modify
	 * the input geometry.
	 *
	 * @param shape          the input PShape containing rings or LineStrings to
	 *                       sample
	 * @param points         number of samples to take from each linear component
	 * @param offsetDistance perpendicular offset applied to each sampled point
	 * @param startOffset    fractional start position along each component (0..1)
	 * @return a list of PVector points sampled from every linear component; empty
	 *         if none produce samples
	 * @since 1.3.0
	 */
	public static List<PVector> pointsOnExterior(PShape shape, int points, double offsetDistance, double startOffset) {
		List<PVector> coords = new ArrayList<>();
		// normalise startOffset into [0,1)
		double startNorm = startOffset % 1.0;

		PGS.applyToLinealGeometries(shape, ring -> {
			// Normalise orientation so sampling/offset direction is consistent
			if (Orientation.isCCW(ring.getCoordinates())) {
				ring = ring.reverse();
			}
			final var l = new IndexedLengthIndexedLine(ring);
			if (l.getEndIndex() == 0) {
				return ring; // skip zero-length components
			}
			final double increment = 1.0 / points;
			double start = startNorm;
			for (int i = 0; i < points; i++) {
				double posFrac = (start + i * increment) % 1.0;
				final Coordinate coord = l.extractPoint(posFrac * l.getEndIndex(), offsetDistance);
				coords.add(PGS.toPVector(coord));
			}
			return ring; // we don't modify geometries here
		});

		return coords;
	}

	/**
	 * Sample points along every linear component of a shape.
	 *
	 * <p>
	 * Samples points independently for each lineal component (polygon exterior
	 * rings, interior rings/holes, and standalone LineStrings) found in
	 * {@code shape}. Each component is sampled along its length with approximately
	 * {@code interPointDistance} spacing. {@code offsetDistance} is applied
	 * perpendicular to the component when extracting each point.
	 *
	 * <p>
	 * Orientation is normalised per component before sampling so offsets are
	 * applied consistently. Components with zero length (or that yield zero
	 * samples) are skipped. The method collects points only and does not modify the
	 * input geometry.
	 *
	 * @param shape              the input PShape containing polygon rings or
	 *                           LineStrings to sample
	 * @param interPointDistance approximate spacing between consecutive samples on
	 *                           each linear component
	 * @param offsetDistance     perpendicular offset applied to each sampled point
	 * @return a list of PVector points sampled from every linear component; empty
	 *         if no samples are produced
	 * @see #applyToLinealGeometries(PShape, java.util.function.UnaryOperator)
	 * @since 1.3.0
	 */
	public static List<PVector> pointsOnExterior(PShape shape, double interPointDistance, double offsetDistance) {
		List<PVector> coords = new ArrayList<>();
		PGS.applyToLinealGeometries(shape, ring -> {
			// Normalise orientation so sampling/offset direction is consistent
			if (Orientation.isCCW(ring.getCoordinates())) {
				ring = ring.reverse();
			}
			final var l = new IndexedLengthIndexedLine(ring);
			if (l.getEndIndex() == 0) {
				return ring;
			}
			final int points = (int) Math.round(l.getEndIndex() / interPointDistance);
			final double increment = 1d / points;
			for (double distance = 0; distance < 1; distance += increment) {
				final Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
				coords.add(PGS.toPVector(coord));
			}
			return ring;
		});

		return coords;
	}

	/**
	 * Extracts evenly spaced dashed line segments along every boundary of a shape.
	 *
	 * <p>
	 * For each linear component found in {@code shape} (each polygon perimeter or
	 * path) this method places repeated line segments along that component's
	 * length. Sampling is performed independently per component: exactly the same
	 * spacing/length rules are applied to each component, but counts and positions
	 * are computed from that component's own perimeter.
	 *
	 * @param shape             input PShape containing polygons or paths (or a mix
	 *                          of both)
	 * @param lineLength        desired length of each segment
	 * @param interLineDistance gap between consecutive segments along a component
	 * @param offset            a fractional offset (a fraction of the component's
	 *                          perimeter) that sets where the first segment starts.
	 *                          Any real value is accepted and is normalised by
	 *                          modulo 1.0 (wrapped into [0,1)). Increasing offset
	 *                          values move the start point clockwise around the
	 *                          perimeter, while decreasing or negative values move
	 *                          it counter-clockwise. Varying {@code offset} over
	 *                          time will effectively "animate" the segments around
	 *                          the component.
	 * @return if the input contains only one linear element (i.e. a holeless
	 *         polygon), a GROUP shape of segments; otherwise a GROUP PShape whose
	 *         children are GROUPs of segment PShapes (one child group per linear
	 *         component).
	 * @since 2.0
	 */
	public static PShape segmentsOnExterior(PShape shape, double lineLength, double interLineDistance, double offset) {
		// Normalise parameters
		double lineLengthFinal = Math.max(lineLength, 0.1);
		double interLineDistanceFinal = Math.max(interLineDistance, 0.0);

		PShape topGroup = new PShape(PConstants.GROUP);

		// Process every linear component independently
		PGS.applyToLinealGeometries(shape, ring -> {
			// Normalise orientation so offsets are consistent
			if (Orientation.isCCW(ring.getCoordinates())) {
				ring = ring.reverse();
			}

			LengthIndexedLine l = new LengthIndexedLine(ring);
			double perimeter = l.getEndIndex();
			if (perimeter <= 0) {
				return ring; // skip empty components
			}

			double totalSegmentLength = lineLengthFinal + interLineDistanceFinal;
			int numberOfLines = (int) Math.floor(perimeter / totalSegmentLength);
			if (numberOfLines <= 0) {
				numberOfLines = 1; // ensure at least one placement to avoid division by zero
			}

			// Adjust lengths so segments + gaps tile the perimeter evenly
			double adjustmentFactor = perimeter / (numberOfLines * totalSegmentLength);
			double adjLineLength = lineLengthFinal * adjustmentFactor;
			double adjInterLineDistance = interLineDistanceFinal * adjustmentFactor;
			totalSegmentLength = adjLineLength + adjInterLineDistance;

			// offset convention: positive values should wrap CW
			double startingPosition = (((-offset) % 1.0) + 1.0) % 1.0 * perimeter;

			PShape compGroup = new PShape(PConstants.GROUP);
			for (int i = 0; i < numberOfLines; i++) {
				double lineStart = (startingPosition + i * totalSegmentLength) % perimeter;
				double lineEnd = (lineStart + adjLineLength) % perimeter;

				Geometry segmentGeom;
				if (lineStart < lineEnd) {
					segmentGeom = l.extractLine(lineStart, lineEnd);
				} else {
					// Wrap-around case, combine two pieces
					Coordinate[] c1 = l.extractLine(lineStart, perimeter).getCoordinates();
					Coordinate[] c2 = l.extractLine(0, lineEnd).getCoordinates();
					Coordinate[] combined = new Coordinate[c1.length + c2.length];
					System.arraycopy(c1, 0, combined, 0, c1.length);
					System.arraycopy(c2, 0, combined, c1.length, c2.length);
					segmentGeom = PGS.GEOM_FACTORY.createLineString(combined);
				}

				PShape segShape = PGS_Conversion.toPShape(segmentGeom);
				segShape.setName(String.format("i=%s@%s", i, lineStart));
				compGroup.addChild(segShape);
			}

			if (!PGS.isEmptyShape(compGroup)) {
				topGroup.addChild(compGroup);
			}

			return ring;
		});

		if (topGroup.getChildCount() == 1) {
			return topGroup.getChild(0);
		}

		return topGroup;
	}

	/**
	 * Creates an CW-oriented length-indexed line from a given PShape. NOTE extracts
	 * first ring from multipolygon
	 */
	private static IndexedLengthIndexedLine makeIndexedLine(PShape shape) {
		Geometry g = fromPShape(shape);
		if (g instanceof Polygonal) {
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)) {
				g = g.getGeometryN(0);
			}
			LinearRing e = ((Polygon) g).getExteriorRing();
			if (Orientation.isCCW(e.getCoordinates())) {
				e = (LinearRing) e.copy().reverse();
			}
			g = e;
		}
		return new IndexedLengthIndexedLine(g);
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
		for (Coordinate coord : coords) {
			perimeter.vertex((float) coord.x, (float) coord.y);
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
	 * Computes all <b>points</b> of intersection between the <b>linework</b> of two
	 * shapes.
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
		final Collection<?> segmentStringsA = SegmentStringUtil.extractSegmentStrings(fromPShape(a));
		final Collection<?> segmentStringsB = SegmentStringUtil.extractSegmentStrings(fromPShape(b));

		return intersections(segmentStringsA, segmentStringsB);
	}

	static List<PVector> intersections(Collection<?> segmentStringsA, Collection<?> segmentStringsB) {
		final Collection<?> larger, smaller;
		if (segmentStringsA.size() > segmentStringsB.size()) {
			larger = segmentStringsA;
			smaller = segmentStringsB;
		} else {
			larger = segmentStringsB;
			smaller = segmentStringsA;

		}

		final Set<PVector> points = new HashSet<>();
		// finds possibly overlapping bounding boxes
		final MCIndexSegmentSetMutualIntersector mci = new MCIndexSegmentSetMutualIntersector(larger);
		// checks if two segments actually intersect
		final SegmentIntersectionDetector sid = new SegmentIntersectionDetector();

		mci.process(smaller, new SegmentIntersector() {
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
	 * @since 1.1.0
	 * @see #generateRandomPoints(PShape, int)
	 * @see #generateRandomGridPoints(PShape, int, boolean, double)
	 */
	public static List<PVector> generateRandomPoints(PShape shape, int points, long seed) {
		var sampler = new ShapeRandomPointSampler(shape, seed);
		return sampler.getRandomPoints(points);
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
		YStripesPointInAreaLocator pointLocator = new YStripesPointInAreaLocator(g);

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
		final double rActual = r == 1 ? r : r % 1;

		PShape out = new PShape();
		PGS_Processing.apply(shape, child -> {
			final Polygon polygon = (Polygon) fromPShape(child);

			final Polygon[] derivedPolygons = new Polygon[n + 1];
			derivedPolygons[0] = polygon;
			Polygon currentPolygon = polygon;

			for (int i = 0; i < n; i++) {
				Coordinate[] inputCoordinates = currentPolygon.getCoordinates();
				int numVertices = inputCoordinates.length - 1;

				Coordinate[] derivedCoordinates = new Coordinate[numVertices + 1];

				for (int k = 0; k < numVertices; k++) {
					double x = inputCoordinates[k].x * (1 - r) + inputCoordinates[(k + 1) % numVertices].x * rActual;
					double y = inputCoordinates[k].y * (1 - r) + inputCoordinates[(k + 1) % numVertices].y * rActual;
					derivedCoordinates[k] = new Coordinate(x, y);
				}
				derivedCoordinates[numVertices] = derivedCoordinates[0]; // close the ring

				Polygon derivedPolygon = PGS.GEOM_FACTORY.createPolygon(derivedCoordinates);
				derivedPolygon.setUserData(polygon.getUserData()); // copy styling

				derivedPolygons[i + 1] = derivedPolygon;
				currentPolygon = derivedPolygon;
			}
			out.addChild(toPShape(GEOM_FACTORY.createMultiPolygon(derivedPolygons)));
		});

		return out.getChildCount() == 1 ? out.getChild(1) : out;
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
		// https://stackoverflow.com/questions/64252638/
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

				Geometry UL = OverlayNG.overlay(slice, GEOM_FACTORY.toGeometry(ulEnv), OverlayNG.INTERSECTION, noder);
				Geometry UR = OverlayNG.overlay(slice, GEOM_FACTORY.toGeometry(urEnv), OverlayNG.INTERSECTION, noder);
				Geometry LL = OverlayNG.overlay(slice, GEOM_FACTORY.toGeometry(llEnv), OverlayNG.INTERSECTION, noder);
				Geometry LR = OverlayNG.overlay(slice, GEOM_FACTORY.toGeometry(lrEnv), OverlayNG.INTERSECTION, noder);

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
	 * Splits the input shape into multiple wedge-shaped regions by connecting the
	 * centroid to each <b>vertex</b>.
	 * <p>
	 * This method computes the centroid of the given shape, and then draws a line
	 * from the centroid to every vertex on the shape’s exterior. The shape is then
	 * split along these lines, partitioning it into distinct regions—one for each
	 * side of the original shape. For polygons with {@code n} vertices, this
	 * typically creates {@code n} wedge-shaped regions. In the case of concave
	 * polygons, this operation may result in more than {@code n} regions due to the
	 * underlying geometry.
	 * </p>
	 * <p>
	 * The returned shape is a GROUP {@code PShape}, whose children are the
	 * resulting partitioned regions.
	 * </p>
	 *
	 * @param shape The shape to be partitioned. Typically a polygonal
	 *              {@code PShape}.
	 * @return A GROUP {@code PShape} which contains one child for each wedge-shaped
	 *         partition created by splitting from the centroid to each vertex.
	 * @since 2.1
	 * @see #centroidSplit(PShape, int, double)
	 */
	public static PShape centroidSplit(PShape shape) {
		var splits = PGS.prepareLinesPShape(null, null, null);
		var c = PGS_ShapePredicates.centroid(shape);
		PGS_Conversion.toPVector(shape).forEach(p -> {
			splits.vertex(p.x, p.y);
			splits.vertex(c.x, c.y);
		});
		splits.endShape();

		var croppedLines = toPShape(fromPShape(shape).intersection(fromPShape(splits)));
		var splitPolygons = PGS_ShapeBoolean.unionLines(shape, croppedLines);
		return PGS_Optimisation.radialSortFaces(splitPolygons, c, 0);
	}

	/**
	 * Splits the input shape into {@code n} wedge-shaped regions by connecting the
	 * centroid to points along the perimeter.
	 * <p>
	 * This method computes the centroid of the given shape, and then splits the
	 * shape by connecting the centroid to {@code n} points sampled evenly (with
	 * optional offset) around the outer ring (perimeter) of the shape. The split
	 * lines start at these points and end at the centroid. The operation results in
	 * approximately {@code n} regions for convex shapes, but may produce more than
	 * {@code n} partitions for highly concave input.
	 * </p>
	 * <p>
	 * The {@code offset} parameter determines the rotation of the sampling around
	 * the perimeter.
	 * </p>
	 * <p>
	 * The returned shape is a GROUP {@code PShape}, whose children are the
	 * resulting wedge-shaped partitions, radially sorted around the centroid.
	 * </p>
	 *
	 * @param shape  The shape to be split; typically a polygonal {@code PShape}.
	 * @param n      The number of splits (rays) to create from the centroid;
	 *               usually determines number of regions.
	 * @param offset The offset for the split lines, as a fraction of the perimeter,
	 *               in [0, 1).
	 * @return A GROUP {@code PShape}, with child shapes being the regions created
	 *         by the centroid splits, sorted radially.
	 * @see #centroidSplit(PShape)
	 */
	public static PShape centroidSplit(PShape shape, int n, double offset) {
		// 1) build your n radial “spoke” lines from the centroid out to the boundary
		PVector c = PGS_ShapePredicates.centroid(shape);
		PShape splits = PGS.prepareLinesPShape(null, null, null);
		var in = fromPShape(shape);
		if (in instanceof Polygon) {
			in = ((Polygon) in).getExteriorRing();
		}
		LengthIndexedLine l = new LengthIndexedLine(in);

		for (int i = 0; i < n; i++) {
			double f = (i / (double) n + offset) % 1.0;
			Coordinate coord = l.extractPoint(f * l.getEndIndex());
			PVector p = PGS.toPVector(coord);
			splits.vertex(p.x, p.y);
			splits.vertex(c.x, c.y);
		}
		splits.endShape();

		// 2) slice the shape by those lines
		PShape cropLines = toPShape(fromPShape(shape).intersection(fromPShape(splits)));
		PShape splitPolygons = PGS_ShapeBoolean.unionLines(shape, cropLines);

		// order faces so that face[0] is the wedge starting at offset
		return PGS_Optimisation.radialSortFaces(splitPolygons, c, offset);
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
	 * Dissolves the linear components of a shape (or group of shapes) into a set of
	 * maximal-length lines in which each unique segment appears once.
	 * <p>
	 * Example uses: avoid double-drawing shared polygon edges; merge contiguous
	 * segments for cleaner rendering/export; or extract non-redundant network edges
	 * for topology or routing.
	 * </p>
	 * <p>
	 * This method does not node the input lines. Crossing segments without a vertex
	 * at the intersection remain crossing in the output.
	 * </p>
	 *
	 * @param shape The {@code PShape} containing linear geometry.
	 * @return A GROUP {@code PShape} whose children are the dissolved
	 *         maximal-length lines.
	 * @since 2.1
	 */
	public static PShape dissolve(PShape shape) {
		var g = fromPShape(shape);
		return toPShape(LineDissolver.dissolve(g));
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
	 * Normalises a shape by standardising its vertex ordering and orientation:
	 * <ul>
	 * <li>The outer shell (exterior contour) is oriented clockwise (CW).</li>
	 * <li>All holes (interior contours) are oriented counterclockwise (CCW).</li>
	 * <li>Each contour (shell or hole) is rotated so its sequence starts at the
	 * vertex with the minimum coordinate (lexicographically by x, then y).</li>
	 * </ul>
	 *
	 * @param shape the {@code PShape} to normalise
	 * @return a new {@code PShape} instance with standardised vertex winding and
	 *         canonicalised vertex rotation
	 * @since 2.1
	 */
	public static PShape normalise(PShape shape) {
		var g = fromPShape(shape);
		g.normalize();
		return toPShape(g);
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
		List<PShape> filteredFaces = PGS_Conversion.getChildren(shape).stream().filter(filterFunction::test).toList();
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
	 * Applies a specified transformation function to each child of the given PShape
	 * and returns a list of results produced by the function.
	 * <p>
	 * This method processes each child of the input shape using the provided
	 * function, which can transform the shape into any desired type T. The function
	 * can return null to exclude a shape from the result list.
	 * <p>
	 * Unlike the {@link #transform(PShape, UnaryOperator)} method, this method does
	 * not flatten the results into a PShape. Instead, it returns a list of
	 * arbitrary objects (type T) produced by the transformation function. This
	 * makes it more flexible for use cases where the transformation does not
	 * necessarily produce PShape objects.
	 * <p>
	 * The transformation function can:
	 * <ul>
	 * <li>Transform the shape into a new object of type T</li>
	 * <li>Return null to exclude the shape from the result list</li>
	 * </ul>
	 * <p>
	 * Note: This method does not modify the original shape or its children. It only
	 * applies the transformation function to each child and collects the results.
	 *
	 * @param <T>      The type of the objects produced by the transformation
	 *                 function.
	 * @param shape    The PShape whose children will be transformed.
	 * @param function A Function that takes a PShape as input and returns an object
	 *                 of type T. If the function returns null for a shape, that
	 *                 shape will be excluded from the result list.
	 * @return A list of objects of type T produced by applying the transformation
	 *         function to each child of the input shape.
	 * @see #transform(PShape, UnaryOperator)
	 * @since 2.1
	 */
	public static <T> List<T> forEachShape(PShape shape, Function<PShape, T> function) {
		return PGS_Conversion.getChildren(shape).stream().map(function).filter(Objects::nonNull).collect(Collectors.toList());
	}

	/**
	 * Applies a specified transformation function to each child of the given PShape
	 * along with its index and returns a list of results produced by the function.
	 * <p>
	 * This method processes each child of the input shape using the provided
	 * function, which takes both the index of the child and the child itself as
	 * input. The function can transform the shape into any desired type T or return
	 * null to exclude the shape from the result list.
	 * <p>
	 * Unlike the {@link #transformWithIndex(PShape, BiFunction)} method, this
	 * method does not flatten the results into a PShape. Instead, it returns a list
	 * of arbitrary objects (type T) produced by the transformation function. This
	 * makes it more flexible for use cases where the transformation does not
	 * necessarily produce PShape objects.
	 * <p>
	 * The transformation function can:
	 * <ul>
	 * <li>Transform the shape into a new object of type T</li>
	 * <li>Return null to exclude the shape from the result list</li>
	 * </ul>
	 * <p>
	 * Note: This method does not modify the original shape or its children. It only
	 * applies the transformation function to each child and collects the results.
	 *
	 * @param <T>      The type of the objects produced by the transformation
	 *                 function.
	 * @param shape    The PShape whose children will be transformed.
	 * @param function A BiFunction that takes an integer index and a PShape as
	 *                 input and returns an object of type T. If the function
	 *                 returns null for a shape, that shape will be excluded from
	 *                 the result list.
	 * @return A list of objects of type T produced by applying the transformation
	 *         function to each child of the input shape along with its index.
	 * @see #transformWithIndex(PShape, BiFunction)
	 * @see #forEachShape(PShape, Function)
	 * @since 2.1
	 */
	public static <T> List<T> forEachShapeWithIndex(PShape shape, BiFunction<Integer, PShape, T> function) {
		List<PShape> children = PGS_Conversion.getChildren(shape);
		return IntStream.range(0, children.size()).mapToObj(i -> function.apply(i, children.get(i))).filter(Objects::nonNull).collect(Collectors.toList());
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
