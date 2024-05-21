package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.SplittableRandom;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;
import org.locationtech.jts.util.GeometricShapeFactory;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinspin.index.Index.PointEntryKnn;
import org.tinspin.index.PointDistance;
import org.tinspin.index.PointMap;
import org.tinspin.index.covertree.CoverTree;

import micycle.pgs.commons.FrontChainPacker;
import micycle.pgs.commons.LargestEmptyCircles;
import micycle.pgs.commons.RepulsionCirclePack;
import micycle.pgs.commons.TangencyPack;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Circle packings of shapes, subject to varying constraints and patterns of
 * tangencies.
 * <p>
 * Each method produces a circle packing with different characteristics using a
 * different technique; for this reason input arguments vary across the methods.
 * <p>
 * The output of each method is a list of PVectors, each representing one
 * circle: (.x, .y) represent the center point and .z represents radius.
 * <p>
 * Where applicable, packings will include circles that overlap with the shape,
 * rather than only including those circles whose center point lies inside the
 * shape.
 * 
 * @author Michael Carleton
 * @since 1.1.0
 *
 */
public final class PGS_CirclePacking {

	/*-
	 * Roadmap (see/implement): 'A LINEARIZED CIRCLE PACKING ALGORITHM'? 
	 * 'A note on circle packing' Young Joon AHN.
	 */

	private PGS_CirclePacking() {
	}

	/**
	 * Packs circles of varying radii within a given shape, whilst respecting
	 * pointal obstacles using the Largest Empty Circle (LEC) algorithm. The method
	 * continues to generate circles until the sum of the areas of the circles
	 * exceeds a specified proportion of the area of the given shape.
	 * 
	 * @param shape          The shape within which circles will be packed. The
	 *                       shape should be in the form of PShape.
	 * @param pointObstacles A collection of PVector points representing obstacles,
	 *                       around which circles are packed. Only points contained
	 *                       within the shape are relevant.
	 * @param areaCoverRatio The target ratio of the total area of the circles to
	 *                       the area of the shape. This parameter should be a
	 *                       double between 0 and 1. Circle generation will stop
	 *                       when this ratio is reached.
	 * @return A list of PVectors, where each PVector represents a circle. The x and
	 *         y components of the PVector represent the center of the circle, and
	 *         the z component represents the radius of the circle.
	 * @since 1.4.0
	 */
	public static List<PVector> obstaclePack(PShape shape, Collection<PVector> pointObstacles, double areaCoverRatio) {
		final Geometry geometry = fromPShape(shape);

		LargestEmptyCircles lec = new LargestEmptyCircles(fromPShape(PGS_Conversion.toPointsPShape(pointObstacles)), geometry,
				areaCoverRatio > 0.95 ? 0.5 : 1);

		final double shapeArea = geometry.getArea();
		double circlesArea = 0;
		List<PVector> circles = new ArrayList<>();

		while (circlesArea / shapeArea < areaCoverRatio) {
			double[] currentLEC = lec.findNextLEC();
			circles.add(new PVector((float) currentLEC[0], (float) currentLEC[1], (float) currentLEC[2]));
			circlesArea += Math.PI * currentLEC[2] * currentLEC[2];
			if (currentLEC[2] < 0.5) {
				break;
			}
		}
		return circles;
	}

	/**
	 * Generates a circle packing of the input shape, using the inscribed circles
	 * (or incircles) of triangles from a triangulation of the shape.
	 * <p>
	 * Circles in this packing do not overlap and are contained entirely within the
	 * shape. However, not every circle is necessarily tangent to others.
	 * 
	 * @param shape       the shape from which to generate a circle packing
	 * @param points      the number of random points to insert into the
	 *                    triangulation as steiner points. Larger values lead to
	 *                    more circles that are generally smaller.
	 * @param refinements number of times to refine the underlying triangulation.
	 *                    Larger values lead to more circles that are more regularly
	 *                    spaced and sized. 0...3 is a suitable range for this
	 *                    parameter
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 */
	public static List<PVector> trinscribedPack(PShape shape, int points, int refinements) {
		final List<PVector> steinerPoints = PGS_Processing.generateRandomPoints(shape, points);
		final IIncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape, steinerPoints, true, refinements, true);
		return StreamSupport.stream(tin.triangles().spliterator(), false).filter(filterBorderTriangles).map(t -> inCircle(t))
				.collect(Collectors.toList());
	}

	/**
	 * Generates a random circle packing of the input shape by generating random
	 * points one-by-one and calculating the maximum radius a circle at each point
	 * can have (such that it's tangent to its nearest circle or a shape vertex).
	 * <p>
	 * Notably, the {@code points} argument defines the number of random point
	 * attempts (or circle attempts), and not the number of circles in the final
	 * packing output, since a point is rejected if it lies in an existing circle or
	 * whose nearest circle is less than minRadius distance away. In other words,
	 * {@code points} defines the maximum number of circles the packing can have; in
	 * practice, the packing will contain somewhat fewer circles.
	 * <p>
	 * Circles in this packing do not overlap and are contained entirely within the
	 * shape. However, not every circle is necessarily tangent to other circles (in
	 * which case, such a circle will be tangent to a shape vertex).
	 * 
	 * @param shape             the shape from which to generate a circle packing
	 * @param points            number of random points to generate (this is not the
	 *                          number of circles in the packing).
	 * @param minRadius         filter (however not simply applied at the end, so
	 *                          affects how the packing operates during packing)
	 * @param triangulatePoints when true, triangulates an initial random point set
	 *                          and uses triangle centroids as the random point set
	 *                          instead; this results in a packing that covers the
	 *                          shape more evenly (particularly when points is
	 *                          small), which is sometimes desirable
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 */
	public static List<PVector> stochasticPack(final PShape shape, final int points, final double minRadius, boolean triangulatePoints) {
		return stochasticPack(shape, points, minRadius, triangulatePoints, System.nanoTime());
	}

	/**
	 * 
	 * Generates a seeded random circle packing within the input shape. Circles are
	 * created one-by-one by generating random points and calculating the maximum
	 * possible radius for a circle at each point, ensuring it is tangent to its
	 * nearest circle or shape vertex.
	 * <p>
	 * The {@code points} parameter defines the number of random point attempts (or
	 * circle attempts), but not the exact number of circles in the final packing
	 * output. An attempted point is rejected if it lies within an existing circle
	 * or if its nearest circle is less than {@code minRadius} distance away.
	 * <p>
	 * Thus, {@code points} defines the maximum number of circles the packing can
	 * have; in practice, the packing will contain somewhat fewer circles.
	 * <p>
	 * The generated circles in the packing do not overlap and are contained
	 * entirely within the shape. However, not every circle is necessarily tangent
	 * to other circles; such circles will be tangent to a shape vertex.
	 * 
	 * @param shape             the shape within which to generate the circle
	 *                          packing
	 * @param points            number of random points to generate (not necessarily
	 *                          equal to the number of circles in the packing)
	 * @param minRadius         minimum allowed radius for circles in the packing
	 *                          (affects packing generation, not applied as a filter
	 *                          at the end)
	 * @param triangulatePoints when true, triangulates the initial random point set
	 *                          and uses triangle centroids as the random point set
	 *                          instead; results in a more evenly distributed
	 *                          packing (particularly when the number of points is
	 *                          small), which may be desirable
	 * @param seed              random seed used to initialize the underlying random
	 *                          number generator
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point, and .z represents the radius.
	 */
	public static List<PVector> stochasticPack(final PShape shape, final int points, final double minRadius, boolean triangulatePoints,
			long seed) {

		final PointMap<PVector> tree = CoverTree.create(3, 2, circleDistanceMetric);

		List<PVector> steinerPoints = PGS_Processing.generateRandomPoints(shape, points, seed);
		if (triangulatePoints) {
			final IIncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape, steinerPoints, true, 1, true);
			steinerPoints = StreamSupport.stream(tin.triangles().spliterator(), false).filter(filterBorderTriangles)
					.map(PGS_CirclePacking::centroid).collect(Collectors.toList());
		}

		IndexedFacetDistance indexedFacetDistance = new IndexedFacetDistance(fromPShape(shape));
		var alpha = shape.getVertex(0); // seed the tree
		tree.insert(new double[] { alpha.x, alpha.y, 0 }, alpha);

		/*
		 * "To find the circle nearest to a center (x, y), do a proximity search at (x,
		 * y, R), where R is greater than or equal to the maximum radius of a circle."
		 */
		double largestR = 0; // the radius of the largest circle in the tree
		final List<PVector> out = new ArrayList<>();
		for (PVector p : steinerPoints) {
			final PointEntryKnn<PVector> nn = tree.query1nn(new double[] { p.x, p.y, largestR }); // find nearest-neighbour circle

			/*
			 * nn.dist() does not return the radius (since it's a distance metric used to
			 * find nearest circle), so now calculate maximum radius for candidate circle
			 * using 2d euclidean distance between center points minus radius of nearest
			 * circle.
			 */
			final float dx = p.x - nn.value().x;
			final float dy = p.y - nn.value().y;
			final double distanceToNN = (Math.sqrt(dx * dx + dy * dy) - nn.value().z);
			final double distanceToBoundary = indexedFacetDistance.distance(PGS.pointFromPVector(p));
			final double packedDistance = Math.min(distanceToNN, distanceToBoundary);
			if (packedDistance > minRadius) {
				largestR = (packedDistance >= largestR) ? packedDistance : largestR;
				p.z = (float) packedDistance;
				tree.insert(new double[] { p.x, p.y, packedDistance }, p); // insert circle into tree
				out.add(p);
			}
		}
		return out;
	}

	/**
	 * Generates a random circle packing of tangential circles with varying radii
	 * that overlap the given shape. The method name references the packing
	 * algorithm used (Front Chain Packing), rather than any particular
	 * characteristic of the circle packing.
	 * <p>
	 * You can set <code>radiusMin</code> equal to <code>radiusMax</code> for a
	 * packing of equal-sized circles using this approach.
	 *
	 * @param shape     the shape within which to generate the circle packing
	 * @param radiusMin minimum radius of circles in the packing
	 * @param radiusMax maximum radius of circles in the packing
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 */
	public static List<PVector> frontChainPack(PShape shape, double radiusMin, double radiusMax) {
		return frontChainPack(shape, radiusMin, radiusMax, System.nanoTime());
	}

	/**
	 * Generates a random circle packing of tangential circles with varying radii
	 * that overlap the given shape. The method name references the packing
	 * algorithm used (Front Chain Packing), rather than any particular
	 * characteristic of the circle packing.
	 * <p>
	 * You can set <code>radiusMin</code> equal to <code>radiusMax</code> for a
	 * packing of equal-sized circles using this approach.
	 *
	 * @param shape     the shape within which to generate the circle packing
	 * @param radiusMin minimum radius of circles in the packing
	 * @param radiusMax maximum radius of circles in the packing
	 * @param seed      random seed used to initialize the underlying random number
	 *                  generator
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 */
	public static List<PVector> frontChainPack(PShape shape, double radiusMin, double radiusMax, long seed) {
		radiusMin = Math.max(1f, Math.min(radiusMin, radiusMax)); // choose min and constrain
		radiusMax = Math.max(1f, Math.max(radiusMin, radiusMax)); // choose max and constrain
		final Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();
		IndexedPointInAreaLocator pointLocator;

		final FrontChainPacker packer = new FrontChainPacker((float) e.getWidth(), (float) e.getHeight(), (float) radiusMin,
				(float) radiusMax, (float) e.getMinX(), (float) e.getMinY());

		if (radiusMin == radiusMax) {
			// if every circle same radius, use faster contains check
			pointLocator = new IndexedPointInAreaLocator(g.buffer(radiusMax));
			packer.getCircles().removeIf(p -> pointLocator.locate(PGS.coordFromPVector(p)) == Location.EXTERIOR);
		} else {
			pointLocator = new IndexedPointInAreaLocator(g);
			final PreparedGeometry cache = PreparedGeometryFactory.prepare(g);
			final GeometricShapeFactory circleFactory = new GeometricShapeFactory();
			circleFactory.setNumPoints(8); // approximate circles using octagon for intersects check
			packer.getCircles().removeIf(p -> {
				// first test whether shape contains circle center point (somewhat faster)
				if (pointLocator.locate(PGS.coordFromPVector(p)) != Location.EXTERIOR) {
					return false;
				}

				// if center point not in circle, check whether circle overlaps with shape using
				// intersects() (somewhat slower)
				circleFactory.setCentre(PGS.coordFromPVector(p));
				circleFactory.setSize(p.z * 2); // set diameter
				return !cache.intersects(circleFactory.createCircle());
			});
		}

		return packer.getCircles();
	}

	/**
	 * Packs a specified number of maximum inscribed circles within the given shape
	 * using the Largest Empty Circle (LEC) algorithm.
	 * <p>
	 * This method finds and returns the maximum inscribed circles up to the
	 * specified number (n), starting with the largest circle. It uses a tolerance
	 * value to control the accuracy of the LEC algorithm.
	 *
	 * @param shape     The input shape to pack maximum inscribed circles within.
	 * @param n         The number of maximum inscribed circles to find and pack.
	 * @param tolerance The tolerance value to control the LEC algorithm's accuracy.
	 *                  Higher values yield faster results but lower accuracy. A
	 *                  value of a 1 is good staring point.
	 * @return A list of PVector objects representing the centers (.x, .y) and radii
	 *         (.z) of the maximum inscribed circles.
	 */
	public static List<PVector> maximumInscribedPack(PShape shape, int n, double tolerance) {
		tolerance = Math.max(0.01, tolerance);
		LargestEmptyCircles mics = new LargestEmptyCircles(fromPShape(shape), null, tolerance);

		final List<PVector> out = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			double[] c = mics.findNextLEC();
			out.add(new PVector((float) c[0], (float) c[1], (float) c[2]));
		}

		return out;
	}

	/**
	 * Packs maximum inscribed circles within the given shape using the Largest
	 * Empty Circle (LEC) algorithm.
	 * <p>
	 * This method finds and returns the maximum inscribed circles with a radius
	 * equal to or larger than the specified minimum radius. It uses a tolerance
	 * value to control the accuracy of the LEC algorithm.
	 *
	 * @param shape     The input shape to pack maximum inscribed circles within.
	 * @param minRadius The minimum allowed radius for the inscribed circles.
	 * @param tolerance The tolerance value to control the LEC algorithm's accuracy.
	 *                  Higher values yield faster results but lower accuracy. A
	 *                  value of a 1 is good staring point.
	 * @return A list of PVector objects representing the centers (.x, .y) and radii
	 *         (.z) of the maximum inscribed circles.
	 */
	public static List<PVector> maximumInscribedPack(PShape shape, double minRadius, double tolerance) {
		tolerance = Math.max(0.01, tolerance);
		minRadius = Math.max(0.01, minRadius);
		LargestEmptyCircles mics = new LargestEmptyCircles(fromPShape(shape), null, tolerance);

		final List<PVector> out = new ArrayList<>();
		double[] currentLEC;
		do {
			currentLEC = mics.findNextLEC();
			if (currentLEC[2] >= minRadius) {
				out.add(new PVector((float) currentLEC[0], (float) currentLEC[1], (float) currentLEC[2]));
			}
		} while (currentLEC[2] >= minRadius);

		return out;
	}

	/**
	 * Generates a circle packing having a pattern of tangencies specified by a
	 * triangulation.
	 * 
	 * <p>
	 * This is an implementation of 'A circle packing algorithm' by Charles R.
	 * Collins & Kenneth Stephenson.
	 * 
	 * @param triangulation represents the pattern of tangencies; vertices connected
	 *                      by an edge inthe triangulation represent tangent circles
	 *                      in thepacking
	 * @param boundaryRadii radius of every circle associated with the
	 *                      boundary/perimeter vertices of the triangulation
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 */
	public static List<PVector> tangencyPack(IIncrementalTin triangulation, double boundaryRadii) {
		TangencyPack pack = new TangencyPack(triangulation, boundaryRadii);
		return pack.pack();
	}

	/**
	 * Generates a circle packing having a pattern of tangencies specified by a
	 * triangulation.
	 * <p>
	 * This is an implementation of 'A circle packing algorithm' by Charles R.
	 * Collins & Kenneth Stephenson.
	 * 
	 * @param triangulation represents the pattern of tangencies; vertices connected
	 *                      by an edge inthe triangulation represent tangent circles
	 *                      in the packing
	 * @param boundaryRadii list of radii of circles associated with the
	 *                      boundary/perimeter vertices of the triangulation. The
	 *                      list may have fewer radii than the number of boundary
	 *                      vertices; in this case, boundary radii will wrap around
	 *                      the list
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius. The packing is centered on
	 *         (0, 0) by default.
	 * @since 1.3.0
	 */
	public static List<PVector> tangencyPack(IIncrementalTin triangulation, double[] boundaryRadii) {
		TangencyPack pack = new TangencyPack(triangulation, boundaryRadii);
		return pack.pack();
	}

	/**
	 * Generates a random circle packing of circles with varying radii that overlap
	 * the given shape.
	 * <p>
	 * Repulsion-packing involves iterative pair-repulsion, in which overlapping
	 * circles move away from each other until there is no overlap. A packing is
	 * first computed for the envelope of the shape, and then any circles which do
	 * not overlap with the shape are discarded.
	 * 
	 * @param shape     the shape from which to generate a circle packing
	 * @param radiusMin minimum radius of circles in the packing. the radii
	 *                  parameters can be the same.
	 * @param radiusMax maximum radius of circles in the packing. the radii
	 *                  parameters can be the same.
	 * @param seed      for initial circle positions and radii
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 * @since 1.3.0
	 * @see #repulsionPack(PShape, List)
	 */
	public static List<PVector> repulsionPack(PShape shape, double radiusMin, double radiusMax, long seed) {
		final double rMinA = Math.max(1f, Math.min(radiusMin, radiusMax)); // actual min
		final double rMaxA = Math.max(1f, Math.max(radiusMin, radiusMax)); // actual max
		final Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();

		/*
		 * We want spawn N circles, such that there are enough to (theoretically) cover
		 * the envelope exactly without any overlap, assuming a packing efficiency of
		 * ~85% (close to optimum).
		 */
		double totalArea = e.getArea() * 0.85;
		/*
		 * Average area is not a simple mean since circle area is quadratic with regards
		 * to radius. The actual average area of circles with radii a...b is an integral
		 * of: pi*r^2 dr from r=a to b.
		 */
		double avgCircleArea = ((rMaxA * rMaxA * rMaxA) - (rMinA * rMinA * rMinA));
		avgCircleArea *= (Math.PI / (3 * (rMaxA - rMinA)));
		int n = (int) (totalArea / avgCircleArea);

		List<PVector> points = PGS_PointSet.poissonN(e.getMinX() + rMaxA, e.getMinY() + rMaxA, e.getMaxX() - rMaxA, e.getMaxY() - rMaxA, n,
				seed);
		SplittableRandom r = new SplittableRandom(seed);
		points.forEach(p -> p.z = rMaxA == rMinA ? (float) rMaxA : (float) r.nextDouble(rMinA, rMaxA));

		return repulsionPack(shape, points);
	}

	/**
	 * Generates a circle packing of a shape using a given collection of
	 * (overlapping) circles.
	 * <p>
	 * Circles in the input should be already bounded by the shape (since repulsion
	 * does not push lonely circles towards the shape, but only repulses overlapping
	 * circles); the intended input is one having circles with lots of overlap
	 * (perhaps seeded within a small rectangle), where they may be repulsed from
	 * each other to fill the shape.
	 * <p>
	 * Repulsion-packing involves iterative pair-repulsion, in which overlapping
	 * circles move away from each other until there is no overlap. A packing is
	 * first computed for the envelope of the shape, and then any circles which do
	 * not overlap with the shape are discarded.
	 * 
	 * @param shape   the shape from which to generate a circle packing
	 * @param circles the collection of circles to pack the shape with, specified as
	 *                PVectors, where .z is the radius (>=1) for each circle
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 * @since 1.3.0
	 * @see #repulsionPack(PShape, double, double, long)
	 */
	public static List<PVector> repulsionPack(PShape shape, List<PVector> circles) {
		final Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();

		float radiusMin = Float.MAX_VALUE;
		float radiusMax = Float.MIN_VALUE;
		for (PVector circle : circles) {
			radiusMax = Math.max(1f, Math.max(radiusMax, circle.z));
			radiusMin = Math.max(1f, Math.min(radiusMin, circle.z));
		}

		final RepulsionCirclePack packer = new RepulsionCirclePack(circles, e.getMinX() + radiusMin, e.getMaxX() - radiusMin,
				e.getMinY() + radiusMin, e.getMaxY() - radiusMin, false);

		final List<PVector> packing = packer.getPacking(); // packing result

		IndexedPointInAreaLocator pointLocator;
		if (radiusMin == radiusMax) {
			// if every circle same radius, use faster contains check
			pointLocator = new IndexedPointInAreaLocator(g.buffer(radiusMax));
			packing.removeIf(p -> pointLocator.locate(PGS.coordFromPVector(p)) == Location.EXTERIOR);
		} else {
			pointLocator = new IndexedPointInAreaLocator(g);
			IndexedFacetDistance distIndex = new IndexedFacetDistance(g);
			packing.removeIf(p -> {
				// first test whether shape contains circle center point (somewhat faster)
				if (pointLocator.locate(PGS.coordFromPVector(p)) != Location.EXTERIOR) {
					return false;
				}
				return distIndex.distance(PGS.pointFromPVector(p)) > p.z * 0.5;
			});
		}

		return packing;
	}

	/**
	 * Generates a tiled circle packing consisting of equal-sized circles arranged
	 * in a square lattice (or grid) bounded by the input shape.
	 * <p>
	 * Circles are included in the packing if they overlap with the given shape.
	 * 
	 * @param shape    the shape from which to generate a circle packing
	 * @param diameter diameter of every circle in the packing
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 * @see #hexLatticePack(PShape, double)
	 */
	public static List<PVector> squareLatticePack(PShape shape, double diameter) {
		diameter = Math.max(diameter, 0.1);
		final double radius = diameter / 2;

		final Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();
		// buffer the geometry to use InAreaLocator to test circles for overlap (this
		// works because all circles have the same diameter)
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g.buffer(radius * 0.95));
		final double w = e.getWidth() + diameter + e.getMinX();
		final double h = e.getHeight() + diameter + e.getMinY();

		final List<PVector> out = new ArrayList<>();

		for (double x = e.getMinX(); x < w; x += diameter) {
			for (double y = e.getMinY(); y < h; y += diameter) {
				if (pointLocator.locate(new Coordinate(x, y)) != Location.EXTERIOR) {
					out.add(new PVector((float) x, (float) y, (float) radius));
				}
			}
		}
		return out;
	}

	/**
	 * Generates a tiled circle packing consisting of equal-sized circles arranged
	 * in a hexagonal lattice bounded by the input shape.
	 * <p>
	 * Circles are included in the packing if they overlap with the given shape.
	 * 
	 * @param shape    the shape from which to generate a circle packing
	 * @param diameter diameter of every circle in the packing
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 * @see #squareLatticePack(PShape, double)
	 */
	public static List<PVector> hexLatticePack(PShape shape, double diameter) {
		diameter = Math.max(diameter, 0.1);
		final double radius = diameter / 2d;

		final Geometry g = fromPShape(shape);
		final Envelope e = g.getEnvelopeInternal();
		/*
		 * Buffer the geometry to use InAreaLocator to test circles for overlap (this
		 * works because all circles have the same diameter).
		 */
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g.buffer(radius * 0.95));
		final double w = e.getWidth() + diameter + e.getMinX();
		final double h = e.getHeight() + diameter + e.getMinY();

		final List<PVector> out = new ArrayList<>();

		final double z = radius * Math.sqrt(3); // hex distance between successive columns
		double offset = 0;
		for (double x = e.getMinX(); x < w; x += z) {
			offset = (offset == radius) ? 0 : radius;
			for (double y = e.getMinY() - offset; y < h; y += diameter) {
				if (pointLocator.locate(new Coordinate(x, y)) != Location.EXTERIOR) {
					out.add(new PVector((float) x, (float) y, (float) radius));
				}
			}
		}
		return out;
	}

	/**
	 * Computes the incircle of a triangle; the largest circle contained in a given
	 * triangle.
	 * 
	 * @param t triangle
	 * @return PVector, where x & y represent incenter coordinates, and z represents
	 *         incircle radius.
	 */
	private static PVector inCircle(SimpleTriangle t) {
		final double a = t.getEdgeA().getLength();
		final double b = t.getEdgeB().getLength();
		final double c = t.getEdgeC().getLength();

		double inCenterX = t.getVertexA().x * a + t.getVertexB().x * b + t.getVertexC().x * c;
		inCenterX /= (a + b + c);
		double inCenterY = t.getVertexA().y * a + t.getVertexB().y * b + t.getVertexC().y * c;
		inCenterY /= (a + b + c);

		final double s = (a + b + c) / 2; // semiPerimeter

		final double r = Math.sqrt(((s - a) * (s - b) * (s - c)) / s);

		return new PVector((float) inCenterX, (float) inCenterY, (float) r);
	}

	private static PVector centroid(SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new PVector((float) x, (float) y);
	}

	/**
	 * Calculate the distance between two points in 3D space, where each point
	 * represents a circle with (x, y, r) coordinates. This custom metric considers
	 * both the Euclidean distance between the centers of the circles and the
	 * absolute difference of their radii.
	 * <p>
	 * The metric is defined as follows: Given two points A and B, representing
	 * circles centered at (x1, y1) and (x2, y2) with radii r1 and r2 respectively,
	 * the distance is calculated as sqrt((x1 - x2)^2 + (y1 - y2)^2) + |r1 - r2|.
	 * <p>
	 * This metric can be used to find the nearest circle to a given center (x, y)
	 * in a proximity search. To perform the search, use a point (x, y, R) where R
	 * is greater than or equal to the maximum radius of a circle in the proximity
	 * structure.
	 *
	 * @param p1 3D point representing the first circle (x1, y1, r1)
	 * @param p2 3D point representing the second circle (x2, y2, r2)
	 * @return the distance between the two points based on the custom metric
	 */
	private static final PointDistance circleDistanceMetric = (p1, p2) -> {
		// from https://stackoverflow.com/a/21975136/
		final double dx = p1[0] - p2[0];
		final double dy = p1[1] - p2[1];
		final double dz = p1[2] - p2[2];

		double euclideanDistance = Math.sqrt(dx * dx + dy * dy);
		double absZDifference = Math.abs(dz);
		return euclideanDistance + absZDifference; // negative if inside
	};

	/**
	 * A streams filter to remove triangulation triangles that share at least one
	 * edge with the shape edge.
	 */
	private static final Predicate<SimpleTriangle> filterBorderTriangles = t -> t.getContainingRegion() != null
			&& !t.getEdgeA().isConstrainedRegionBorder() && !t.getEdgeB().isConstrainedRegionBorder()
			&& !t.getEdgeC().isConstrainedRegionBorder();

}
