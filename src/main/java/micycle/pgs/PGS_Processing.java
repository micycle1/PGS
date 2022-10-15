package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.SplittableRandom;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.locationtech.jts.algorithm.Angle;
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

import com.vividsolutions.jcs.conflate.coverage.CoverageCleaner;
import com.vividsolutions.jcs.conflate.coverage.CoverageCleaner.Parameters;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.feature.FeatureDatasetFactory;
import com.vividsolutions.jump.feature.FeatureUtil;
import com.vividsolutions.jump.task.DummyTaskMonitor;

import de.incentergy.geometry.impl.RandomPolygonSplitter;
import micycle.balaban.BalabanSolver;
import micycle.balaban.Point;
import micycle.balaban.Segment;
import micycle.pgs.PGS.GeometryIterator;
import micycle.pgs.color.RGB;
import micycle.pgs.commons.PolygonDecomposition;
import micycle.pgs.commons.SeededRandomPointsInGridBuilder;
import micycle.trapmap.TrapMap;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

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
	 * fraction around its perimeter.
	 * 
	 * @param shape          a lineal or polygonal shape. If the input is a GROUP
	 *                       shape, a single point will be extracted from its first
	 *                       child.
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
		if (g instanceof Polygonal) {
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)) {
				g = g.getGeometryN(0);
			}
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
	 * @param shape          a lineal or polygonal shape. If the input is a GROUP
	 *                       shape, a single point will be extracted from its first
	 *                       child.
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
		if (g instanceof Polygonal) {
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)) {
				g = g.getGeometryN(0);
			}
			g = ((Polygon) g).getExteriorRing();
		}
		LengthIndexedLine l = new LengthIndexedLine(g);

		final double increment = 1d / points;
		for (double distance = 0; distance < 1; distance += increment) {
			final Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
			coords.add(PGS.toPVector(coord));
		}
		return coords;
	}

	/**
	 * Generates a list of points that lie on the exterior/perimeter of the given
	 * shape.
	 * 
	 * @param shape              a lineal or polygonal shape. If the input is a
	 *                           GROUP shape, a single point will be extracted from
	 *                           its first child.
	 * @param interPointDistance distance between each exterior point
	 * @param offsetDistance
	 * @return
	 */
	public static List<PVector> pointsOnExterior(PShape shape, double interPointDistance, double offsetDistance) {
		// TODO points on holes
		Geometry g = fromPShape(shape);
		if (g instanceof Polygonal) {
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)) {
				g = g.getGeometryN(0);
			}
			g = ((Polygon) g).getExteriorRing();
		}
		LengthIndexedLine l = new LengthIndexedLine(g);

		final int points = (int) Math.round(l.getEndIndex() / interPointDistance);

		ArrayList<PVector> coords = new ArrayList<>(points);

		final double increment = 1d / points;
		for (double distance = 0; distance < 1; distance += increment) {
			final Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
			coords.add(PGS.toPVector(coord));
		}
		return coords;
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
		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing();
		}
		final LengthIndexedLine l = new LengthIndexedLine(g);
		final double length = l.getEndIndex(); // perimeter length

		if (from > to) {
			final Geometry l1 = l.extractLine(length * from, length);
			final Geometry l2 = l.extractLine(0, length * to);
			return toPShape(GEOM_FACTORY.createLineString(
					Stream.concat(Arrays.stream(l1.getCoordinates()), Arrays.stream(l2.getCoordinates())).toArray(Coordinate[]::new)));
		}

		return toPShape(l.extractLine(length * from, length * to));
	}

	/**
	 * Finds the angle of the line tangent to the shape at a certain point on its
	 * perimeter (given by the some fraction of the distance around the perimeter).
	 * <p>
	 * The tangent line is orientated clockwise with respect to the shape and the
	 * output angle is normalized to be in the range [ -PI, PI ].
	 * 
	 * @param shape
	 * @param distanceFraction the distance fraction around the perimeter [0...1]
	 * @return the normalized angle (in radians) that a line tangent to the
	 *         perimeter of the shape at the given position makes with the positive
	 *         x-axis, where 0 is north.
	 * @since 1.3.0
	 */
	public static double tangentAngle(PShape shape, double distanceFraction) {
		distanceFraction %= 1;

		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing();
		}
		LengthIndexedLine l = new LengthIndexedLine(g);

		double d1 = (distanceFraction * l.getEndIndex()) - 1e-7;
		double d2 = (distanceFraction * l.getEndIndex()) + 1e-7;

		Coordinate coordA = l.extractPoint(d1); // CCW - point behind first
		Coordinate coordB = l.extractPoint(d2); // CCW - point after second

		return Angle.angle(coordA, coordB);
	}

	/**
	 * Computes all <b>points</b> of intersection between the <b>perimeter</b> of
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
				if (sid.getIntersection() != null) {
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

		final IIncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape);
		final boolean constrained = !tin.getConstraints().isEmpty();
		final double totalArea = StreamSupport.stream(tin.getConstraints().spliterator(), false)
				.mapToDouble(c -> ((PolygonConstraint) c).getArea()).sum();

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
				vertices.add(PGS.toPVector(coord));
			}
		}
		return vertices;
	}

	/**
	 * Removes overlap between polygons contained in a GROUP shape, preserving only
	 * line segments that are visible to a human, rather than a computer (to use as
	 * input for a pen plotter, for example).
	 * <p>
	 * Any overlapping lines are also removed during the operation.
	 * <p>
	 * Order of shape layers is important: the method will consider the last child
	 * shape of the input to be "on top" of all other shapes (as is the case
	 * visually).
	 * 
	 * @param shape a GROUP shape consisting of lineal or polygonal child shapes
	 * @return linework of the overlapping input (LINES PShape)
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
		PGS_Conversion.setAllStrokeColor(out, RGB.setAlpha(RGB.PINK, 128), 8);

		return out;
	}

	/**
	 * Returns a copy of the shape where small holes (i.e. inner rings with area <
	 * given threshold) are removed.
	 * 
	 * @param shape         a single polygonal shape or GROUP polygonal shape
	 * @param areaThreshold remove any holes with an area smaller than this value
	 * @return
	 */
	public static PShape removeSmallHoles(PShape shape, double areaThreshold) {
		final ArrayList<Geometry> polygons = new ArrayList<>();
		final Geometry g = fromPShape(shape);
		for (final Geometry geom : new GeometryIterator(g)) {
			if (geom.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
				polygons.add(removeSmallHoles((Polygon) geom, areaThreshold));
			}
		}
		return toPShape(polygons);
	}

	private static Polygon removeSmallHoles(Polygon polygon, double areaThreshold) {
		Polygon noHolePol = GEOM_FACTORY.createPolygon(polygon.getExteriorRing());
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			final LinearRing hole = polygon.getInteriorRingN(i);
			if (hole.getArea() < areaThreshold) {
				continue;
			}
			noHolePol = (Polygon) noHolePol.difference(hole);
		}
		return noHolePol;
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
			System.err.println(
					"The input to polygonizeLines() contained an odd number of vertices. The method expects successive pairs of vertices.");
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
	 * Splits a shape into 4 equal (as measured be envelope area) quadrants.
	 * 
	 * @param shape the shape to split
	 * @return a GROUP PShape, where each child shape is some quadrant partition of
	 *         the original shape
	 * @see #split(PShape, int)
	 */
	public static PShape split(PShape shape) {
		return split(shape, 1);
	}

	/**
	 * Splits a shape into 4^(1+recursions) rectangular partitions.
	 * 
	 * @param shape
	 * @param splitDepth
	 * @return a GROUP PShape, where each child shape is some quadrant partition of
	 *         the original shape
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
		stack.forEach(g -> {
			partitions.addChild(toPShape(g));
		});
		return partitions;
	}

	/**
	 * Partitions a shape into simple convex polygons.
	 * 
	 * @param shape the shape to partition
	 * @return a GROUP PShape, where each child shape is some convex partition of
	 *         the original shape
	 */
	public static PShape convexPartition(PShape shape) {
		// algorithm described in https://mpen.ca/406/bayazit
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
	 * Partitions a shape into N approximately equal-area polygonal cells.
	 * <p>
	 * This method produces a voronoi-like output.
	 * 
	 * @param shape   a polygonal (non-group, no holes) shape
	 * @param parts   number of roughly equal area partitons to create
	 * @param precise whether to use a subroutine that partitions the shape into
	 *                more precisely equal partitions. The tradeoff here is
	 *                computation time vs partition quality
	 * @return a GROUP PShape, whose child shapes are partitions of the original
	 * @since 1.3.0
	 */
	public static PShape equalPartition(final PShape shape, final int parts, boolean precise) {
		final Geometry g = fromPShape(shape);
		if (g.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			RandomPolygonSplitter splitter = new RandomPolygonSplitter();
			List<? extends Geometry> partitions;
			if (precise) {
				partitions = splitter.split((Polygon) g, parts, 10000, 3);
			} else {
				partitions = splitter.split((Polygon) g, parts, 2000, 5);
			}
			return toPShape(partitions);
		} else {
			System.err.println("equalPartition(): Input shape is not a polygon.");
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

		PGS_Conversion.setAllFillColor(trapezoids, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(trapezoids, RGB.PINK, 1);
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
	 * Removes gaps and overlaps from meshes/polygon collections that are intended
	 * to satisfy the following conditions:
	 * <ul>
	 * <li>Vector-clean - edges between neighbouring polygons must either be
	 * identical or intersect only at endpoints.</li>
	 * <li>Non-overlapping - No two polygons may overlap. Equivalently, polygons
	 * must be interior-disjoint.</li>
	 * </ul>
	 * <p>
	 * It may not always be possible to perfectly clean the input.
	 * <p>
	 * While this method is intended to be used to fix malformed coverages, it can
	 * be used to snap collections of disparate polygons together.
	 * 
	 * @param coverage          a GROUP shape, consisting of the polygonal faces to
	 *                          clean
	 * @param distanceTolerance the distance below which segments and vertices are
	 *                          considered to match
	 * @param angleTolerance    the maximum angle difference between matching
	 *                          segments, in degrees
	 * @return GROUP shape whose child polygons satisfy a (hopefully) valid coverage
	 * @since 1.3.0
	 */
	public static PShape cleanCoverage(PShape coverage, double distanceTolerance, double angleTolerance) {
		final List<Geometry> geometries = PGS_Conversion.getChildren(coverage).stream().map(PGS_Conversion::fromPShape)
				.collect(Collectors.toList());
		final FeatureCollection features = FeatureDatasetFactory.createFromGeometry(geometries);

		final CoverageCleaner cc = new CoverageCleaner(features, new DummyTaskMonitor());
		cc.process(new Parameters(distanceTolerance, angleTolerance));

		final List<Geometry> cleanedGeometries = FeatureUtil.toGeometries(cc.getUpdatedFeatures().getFeatures());
		final PShape out = PGS_Conversion.toPShape(cleanedGeometries);
		PGS_Conversion.setAllStrokeColor(out, RGB.PINK, 2);
		return out;
	}

	private static Polygon toGeometry(Envelope envelope) {
		return GEOM_FACTORY.createPolygon(
				GEOM_FACTORY.createLinearRing(new Coordinate[] { new Coordinate(envelope.getMinX(), envelope.getMinY()),
						new Coordinate(envelope.getMaxX(), envelope.getMinY()), new Coordinate(envelope.getMaxX(), envelope.getMaxY()),
						new Coordinate(envelope.getMinX(), envelope.getMaxY()), new Coordinate(envelope.getMinX(), envelope.getMinY()) }),
				null);
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
