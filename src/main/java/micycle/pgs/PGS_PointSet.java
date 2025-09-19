package micycle.pgs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.SplittableRandom;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.vecmath.Point3d;
import javax.vecmath.Point4d;

import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.Clusterer;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.alg.interfaces.SpanningTreeAlgorithm;
import org.jgrapht.alg.spanning.PrimMinimumSpanningTree;
import org.jgrapht.graph.SimpleGraph;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.Vertex;
import org.tinspin.index.IndexConfig;
import org.tinspin.index.PointMap;

import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import it.unimi.dsi.util.XoRoShiRo128PlusRandomGenerator;
import micycle.pgs.commons.GeometricMedian;
import micycle.pgs.commons.GreedyTSP;
import micycle.pgs.commons.PEdge;
import micycle.pgs.commons.PoissonDistributionJRUS;
import micycle.pgs.commons.ThomasPointProcess;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Generation of random sets of 2D points having a variety of different
 * distributions and constraints (and associated functions).
 * 
 * <p>
 * <strong>Note on Floating-Point Values and Collisions:</strong>
 * </p>
 * <p>
 * When generating many random points, collisions in coordinate values are
 * expected due to the limited precision of floating-point numbers. For example:
 * <ul>
 * <li>A {@code float} has ~24 bits of precision (23 explicit mantissa bits + 1
 * implicit leading bit), which allows for ~16.8 million unique values in the
 * range [0, 1).</li>
 * <li>When generating a large number of points (e.g., 100,000), the probability
 * of collisions follows the birthday paradox formula: <code>n² / (2m)</code>,
 * where <code>n</code> is the number of samples and <code>m</code> is the
 * number of possible unique values.</li>
 * <li>For 100,000 points, this results in ~300 expected collisions in the
 * x-coordinate, even when using a high-quality random number generator.</li>
 * </ul>
 * </p>
 * 
 * @author Michael Carleton
 * @since 1.2.0
 */
public final class PGS_PointSet {

	private static final float SQRT_3 = (float) Math.sqrt(3);
	/** Golden angle (in radians) */
	private static final float GOLDEN_ANGLE = (float) (Math.PI * (3 - Math.sqrt(5)));

	private PGS_PointSet() {
	}

	/**
	 * Returns a filtered copy of the input, containing no points that are within
	 * the <code>distanceTolerance</code> of each other.
	 * <p>
	 * This method can be used to convert a random point set into a blue-noise-like
	 * (poisson) point set.
	 * 
	 * @param points            list of points to filter
	 * @param distanceTolerance a point that is within this distance of a previously
	 *                          included point is not included in the output
	 * @return
	 */
	public static List<PVector> prunePointsWithinDistance(List<PVector> points, double distanceTolerance) {
		final PointMap<Object> tree = PointMap.Factory.createKdTree(IndexConfig.create(2).setDefensiveKeyCopy(false));
		final List<PVector> newPoints = new ArrayList<>();
		for (PVector p : points) {
			final double[] coords = new double[] { p.x, p.y };
			if (tree.size() == 0 || tree.query1nn(coords).dist() > distanceTolerance) {
				tree.insert(coords, p);
				newPoints.add(p);
			}
		}
		return newPoints;
	}

	/**
	 * Prunes a list of points by removing points that are considered not
	 * sufficiently dense.
	 * <p>
	 * A point's density is assessed based on its distance to its nearest neighbor;
	 * if the nearest neighbor of a point is farther away than the specified
	 * distance tolerance, the point is considered sparse and removed. In other
	 * words, only points that have at least one neighbor within the distance
	 * tolerance are kept.
	 *
	 * @param points            A List of {@code PVector} points to be analysed and
	 *                          pruned.
	 * @param distanceTolerance The maximum allowable distance for a point to be
	 *                          considered dense. Points with their nearest neighbor
	 *                          distance greater than the distance tolerance are
	 *                          pruned.
	 * @return A List of {@code PVector} points where each point has at least one
	 *         neighbor within the distance tolerance, i.e., the list of points
	 *         after sparse points have been removed.
	 * @since 2.0
	 */
	public static List<PVector> pruneSparsePoints(Collection<PVector> points, double distanceTolerance) {
		var tin = PGS_Triangulation.delaunayTriangulationMesh(points);
		var vertexDistanceMap = new Object2DoubleOpenHashMap<Vertex>();
		final double toleranceSquared = distanceTolerance * distanceTolerance;

		// Iterate over TIN edges to compute & store minimum distances
		for (var edge : tin.edges()) {
			var A = edge.getA();
			var B = edge.getB();
			double distance = A.getDistanceSq(B);
			// Update the minimum distances using merge
			vertexDistanceMap.merge(A, distance, Math::min);
			vertexDistanceMap.merge(B, distance, Math::min);
		}

		// Collect vertices within the distance tolerance
		return vertexDistanceMap.object2DoubleEntrySet().stream().filter(entry -> entry.getDoubleValue() <= toleranceSquared)
				.map(entry -> new PVector((float) entry.getKey().getX(), (float) entry.getKey().getY())).toList();
	}

	/**
	 * Sorts a list of points according to the Hilbert space-filling curve to ensure
	 * a high-degree of spatial locality in the sequence of points.
	 * 
	 * @param points list of points to sort. a list requires at least 24 points to
	 *               be sorted.
	 * @return a sorted <b>copy</b> of the input list, having a different order
	 *         according to points' Hilbert ranking of their (x, y) coordinate
	 * @since 1.3.0
	 */
	public static List<PVector> hilbertSort(List<PVector> points) {
		if (points.isEmpty() || points.size() < 24) {
			return points;
		}
		return hilbertSortRaw(points).stream().map(Pair::getSecond).collect(Collectors.toList());
	}

	/**
	 * Clusters points into N groups, using k-means clustering.
	 * <p>
	 * K-means finds the N cluster centers and assigns points to the nearest cluster
	 * center, such that the squared (euclidean) distances from the cluster are
	 * minimised.
	 * 
	 * @param points list of points to cluster
	 * @param groups desired number of clustered groups
	 * @since 1.4.0
	 * @see #cluster(Collection, int, long)
	 * @return list of groups, where each group is a list of PVectors
	 */
	public static List<List<PVector>> cluster(Collection<PVector> points, int groups) {
		return cluster(points, groups, System.nanoTime());
	}

	/**
	 * Clusters points into N groups, using k-means clustering.
	 * <p>
	 * K-means finds the N cluster centers and assigns points to the nearest cluster
	 * center, such that the squared (euclidean) distances from the cluster are
	 * minimised.
	 * 
	 * @param points list of points to cluster
	 * @param groups desired number of clustered groups
	 * @param seed   random seed
	 * @since 1.4.0
	 * @return list of groups, where each group is a list of PVectors
	 * @see #cluster(Collection, int)
	 */
	public static List<List<PVector>> cluster(Collection<PVector> points, int groups, long seed) {
		RandomGenerator r = new XoRoShiRo128PlusRandomGenerator(seed);
		Clusterer<CPVector> kmeans = new KMeansPlusPlusClusterer<>(groups, 25, new EuclideanDistance(), r);
		List<CPVector> pointz = points.stream().map(p -> new CPVector(p)).collect(Collectors.toList());

		List<List<PVector>> clusters = new ArrayList<>(groups);
		kmeans.cluster(pointz).forEach(cluster -> {
			clusters.add(cluster.getPoints().stream().map(p -> p.p).collect(Collectors.toList()));
		});

		return clusters;
	}

	/**
	 * Finds the geometric median point of a set of weighted sample points.
	 * <p>
	 * The median point is the point that minimises the sum of (weighted) distances
	 * to the sample points.
	 * <p>
	 * Points are expressed as PVectors; the z coordinate is used as the weight for
	 * each point. Weights must be positive. If every point has a weight of 0 (z=0),
	 * the function returns the median as if each point had an equal non-zero weight
	 * (set to 1).
	 * 
	 * @param points list of points, where the z coordinate is point weight
	 * @since 1.4.0
	 * @return 2D median point
	 */
	public static PVector weightedMedian(Collection<PVector> points) {
		boolean allZero = points.stream().allMatch(p -> p.z == 0);
		Point4d[] wp = points.stream().map(p -> new Point4d(p.x, p.y, 0, allZero ? 1 : p.z)).toArray(Point4d[]::new);
		Point3d median = GeometricMedian.median(wp, 1e-3, 50);
		return new PVector((float) median.x, (float) median.y);
	}

	/**
	 * Generates a set of random (uniform) points that lie within a bounding
	 * rectangle.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @return
	 * @see #random(double, double, double, double, int, long) seeded random()
	 */
	public static List<PVector> random(double xMin, double yMin, double xMax, double yMax, int n) {
		return random(xMin, yMin, xMax, yMax, n, System.nanoTime());
	}

	/**
	 * Generates a set of random (uniform) points that lie within a bounding
	 * rectangle, using the specified seed.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @param seed number used to initialize the underlying pseudorandom number
	 *             generator
	 * @return
	 * @see #random(double, double, double, double, int) non-seeded random()
	 */
	public static List<PVector> random(double xMin, double yMin, double xMax, double yMax, int n, long seed) {
		final SplittableRandom random = new SplittableRandom(seed);
		final List<PVector> points = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			final float x = (float) random.nextDouble(xMin, xMax);
			final float y = (float) random.nextDouble(yMin, yMax);
			points.add(new PVector(x, y));
		}
		return points;
	}

	/**
	 * Generates a set of random points having a gaussian/normal distribution. The
	 * point set is centered around the given center, given by mean coordinates.
	 * 
	 * @param centerX x coordinate of the center/mean of the point set
	 * @param centerY x coordinate of the center/mean of the point set
	 * @param sd      standard deviation, which specifies how much the values can
	 *                vary from the mean. 68% of point samples have a value within
	 *                one standard deviation of the mean; three standard deviations
	 *                account for 99.7% of the sample population
	 * @param n       number of points to generate
	 * @see #gaussian(double, double, double, int, long) seeded gaussian()
	 */
	public static List<PVector> gaussian(double centerX, double centerY, double sd, int n) {
		return gaussian(centerX, centerY, sd, n, System.nanoTime());
	}

	/**
	 * Generates a set of random points having a gaussian/normal distribution, using
	 * the specified seed. The point set is centered around the given center, given
	 * by mean coordinates.
	 * 
	 * @param centerX x coordinate of the center/mean of the point set
	 * @param centerY x coordinate of the center/mean of the point set
	 * @param sd      standard deviation, which specifies how much the values can
	 *                vary from the mean. 68% of point samples have a value within
	 *                one standard deviation of the mean; three standard deviations
	 *                account for 99.7% of the sample population
	 * @param n       number of points to generate
	 * @param seed    number used to initialize the underlying pseudorandom number
	 *                generator
	 * @return
	 * @see #gaussian(double, double, double, int) non-seeded gaussian()
	 */
	public static List<PVector> gaussian(double centerX, double centerY, double sd, int n, long seed) {
		final RandomGenerator random = new XoRoShiRo128PlusRandomGenerator(seed);
		final List<PVector> points = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			final float x = (float) (sd * random.nextGaussian() + centerX);
			final float y = (float) (sd * random.nextGaussian() + centerY);
			points.add(new PVector(x, y));
		}
		return points;
	}

	/**
	 * Generates a square grid/lattice of points that lie within a bounding
	 * rectangle.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum (inclusive)
	 * @param yMax y-coordinate of boundary maximum (inclusive)
	 */
	public static List<PVector> squareGrid(final double xMin, final double yMin, final double xMax, final double yMax, final double pointDistance) {
		final double width = xMax - xMin;
		final double height = yMax - yMin;

		final List<PVector> points = new ArrayList<>();

		for (double x = 0; x <= width; x += pointDistance) {
			for (double y = 0; y <= height; y += pointDistance) {
				points.add(new PVector((float) (x + xMin), (float) (y + yMin)));
			}
		}
		return points;
	}

	/**
	 * Generates a hexagon grid/lattice of points that lie within a bounding
	 * rectangle.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @return
	 * @see #hexGrid(double, double, double, double, double) hexGrid() where
	 *      inter-point distance is specified
	 */
	public static List<PVector> hexGrid(final double xMin, final double yMin, final double xMax, final double yMax, final int n) {
		final double width = xMax - xMin;
		final double height = yMax - yMin;

		final float h = (float) Math.sqrt((width * height * (Math.sqrt(5) / 2)) / n);
		final float v = (float) (h * (2 / Math.sqrt(5)));
		final List<PVector> points = new ArrayList<>(n);

		for (int i = 0; i < width / h; i++) {
			for (int j = 0; j < height / v; j++) {
				points.add(new PVector((i - (j % 2) / 2f) * h + (float) xMin, j * v + (float) yMin));
			}
		}
		return points;
	}

	/**
	 * Generates a hexagon grid of points that lie within a bounding rectangle.
	 * 
	 * @param xMin          x-coordinate of boundary minimum
	 * @param yMin          y-coordinate of boundary minimum
	 * @param xMax          x-coordinate of boundary maximum
	 * @param yMax          y-coordinate of boundary maximum
	 * @param pointDistance inter-point distance
	 * @return
	 * @see #hexGrid(double, double, double, double, int) hexGrid() where number of
	 *      points is specified
	 */
	public static List<PVector> hexGrid(final double xMin, final double yMin, final double xMax, final double yMax, final double pointDistance) {
		final double width = xMax - xMin;
		final double height = yMax - yMin;

		final List<PVector> points = new ArrayList<>();

		for (int i = 0; i < width / pointDistance; i++) {
			for (int j = 0; j < height / pointDistance; j++) {
				points.add(new PVector((float) ((i - (j % 2) / 2f) * pointDistance + xMin), (float) (j * pointDistance + yMin)));
			}
		}
		return points;
	}

	/**
	 * Generates a hexagonal grid of points <b>arranged in a hexagon pattern</b>.
	 * 
	 * @param centerX  x coordinate of the hexagon center point
	 * @param centerY  y coordinate of the hexagon center point
	 * @param length   layers/no. of points on each hexagon side
	 * @param distance inter-point distance
	 */
	public static List<PVector> hexagon(double centerX, double centerY, int length, double distance) {
		final float xOffset = (float) centerX;
		final float yOffset = (float) centerY;
		final float d = (float) distance;

		final List<PVector> points = new ArrayList<>();

		/*
		 * PVector .z is set to length so hexagon layer can be easily identified.
		 */
		for (int i = 0; i <= (length - 1); i++) {
			float y = (SQRT_3 * i * d) / 2.0f;
			for (int j = 0; j < (2 * length - 1 - i); j++) {
				float x = (-(2 * length - i - 2) * d) / 2.0f + j * d;
				points.add(new PVector(x + xOffset, y + yOffset, length));
				if (y != 0) {
					points.add(new PVector(x + xOffset, -y + yOffset, length));
				}
			}
		}
		return points;
	}

	/**
	 * Generates a set of n points that are randomly distributed on a ring
	 * (annulus).
	 * 
	 * @param centerX     x coordinate of the center/mean of the ring
	 * @param centerY     x coordinate of the center/mean of the ring
	 * @param innerRadius radius of the ring's hole
	 * @param outerRadius outer radius of the ring
	 * @param maxAngle    sweep angle of the ring (in radians). Can be negative
	 * @param n           the number of random points to generate
	 * @return a list of PVector objects representing the (x, y) coordinates of the
	 *         random points
	 * @see #ring(double, double, double, double, double, int, long) seeded ring()
	 */
	public static List<PVector> ring(double centerX, double centerY, double innerRadius, double outerRadius, double maxAngle, int n) {
		return ring(centerX, centerY, innerRadius, outerRadius, maxAngle, n, System.nanoTime());
	}

	/**
	 * Generates a set of points that are randomly distributed on a ring (annulus).
	 * 
	 * @param centerX     x coordinate of the center/mean of the ring
	 * @param centerY     x coordinate of the center/mean of the ring
	 * @param innerRadius radius of the ring's hole
	 * @param outerRadius outer radius of the ring
	 * @param maxAngle    angle of the ring (in radians). Can be negative
	 * @param n           the number of random points to generate
	 * @param seed        number used to initialize the underlying pseudorandom
	 *                    number generator
	 * @return a list of PVector objects representing the (x, y) coordinates of the
	 *         random points
	 * @see #ring(double, double, double, double, double, int) non-seeded ring()
	 */
	public static List<PVector> ring(double centerX, double centerY, double innerRadius, double outerRadius, double maxAngle, int n, long seed) {
		final SplittableRandom random = new SplittableRandom(seed);
		final List<PVector> points = new ArrayList<>(n);
		if (maxAngle == 0) {
			maxAngle = Double.MIN_VALUE;
		}
		for (int i = 0; i < n; i++) {
			double randomAngle = (maxAngle < 0 ? -1 : 1) * random.nextDouble(Math.abs(maxAngle));
			double randomRadius = random.nextDouble(innerRadius, outerRadius);
			double x = -Math.sin(randomAngle) * randomRadius;
			double y = Math.cos(randomAngle) * randomRadius;

			points.add(new PVector((float) (x + centerX), (float) (y + centerY)));
		}
		return points;
	}

	/**
	 * Generates a set of random points (constrained within a rectangular region)
	 * via Poisson Disk Sampling.
	 * <p>
	 * Poisson-disc sampling produces points that are tightly-packed, but no closer
	 * to each other than a specified minimum distance, resulting in a more natural
	 * and desirable pattern for many applications. This distribution is also
	 * described as blue noise.
	 * 
	 * @param xMin    x-coordinate of boundary minimum
	 * @param yMin    y-coordinate of boundary minimum
	 * @param xMax    x-coordinate of boundary maximum
	 * @param yMax    y-coordinate of boundary maximum
	 * @param minDist minimum euclidean distance between any two points
	 * @return
	 * @see #poisson(double, double, double, double, double, long) seeded poisson()
	 */
	public static List<PVector> poisson(double xMin, double yMin, double xMax, double yMax, double minDist) {
		return poisson(xMin, yMin, xMax, yMax, minDist, System.nanoTime());
	}

	/**
	 * Generates a set of random points (constrained within a rectangular region)
	 * via Poisson Disk Sampling, using the specified seed.
	 * <p>
	 * Poisson-disc sampling produces points that are tightly-packed, but no closer
	 * to each other than a specified minimum distance, resulting in a more natural
	 * and desirable pattern for many applications. This distribution is also
	 * described as blue noise.
	 * 
	 * @param xMin    x-coordinate of boundary minimum
	 * @param yMin    y-coordinate of boundary minimum
	 * @param xMax    x-coordinate of boundary maximum
	 * @param yMax    y-coordinate of boundary maximum
	 * @param minDist minimum euclidean distance between any two points
	 * @param seed    number used to initialize the underlying pseudorandom number
	 *                generator
	 * @return
	 * @see #poisson(double, double, double, double, double) non-seeded poisson()
	 */
	public static List<PVector> poisson(double xMin, double yMin, double xMax, double yMax, double minDist, long seed) {
		final PoissonDistributionJRUS pd = new PoissonDistributionJRUS(seed);
		return pd.generate(xMin, yMin, xMax, yMax, minDist);
	}

	/**
	 * Generates a poisson point set having N points constrained within a
	 * rectangular region using a random seed.
	 * <p>
	 * Poisson-disc sampling produces points that are tightly-packed, but no closer
	 * to each other than a specified minimum distance, resulting in a more natural
	 * and desirable pattern for many applications. This distribution is also
	 * described as blue noise.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    target size of poisson point set
	 * @return
	 * @see #poissonN(double, double, double, double, int, long)
	 */
	public static List<PVector> poissonN(double xMin, double yMin, double xMax, double yMax, int n) {
		return poissonN(xMin, yMin, xMax, yMax, n, System.nanoTime());
	}

	/**
	 * Generates a poisson point set having N points constrained within a
	 * rectangular region.
	 * <p>
	 * Poisson-disc sampling produces points that are tightly-packed, but no closer
	 * to each other than a specified minimum distance, resulting in a more natural
	 * and desirable pattern for many applications. This distribution is also
	 * described as blue noise.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    target size of poisson point set
	 * @param seed number used to initialize the underlying pseudorandom number
	 *             generator
	 * @return
	 */
	public static List<PVector> poissonN(double xMin, double yMin, double xMax, double yMax, int n, long seed) {
		final PoissonDistributionJRUS pd = new PoissonDistributionJRUS(seed);
		return pd.generate(xMin, yMin, xMax, yMax, n);
	}

	/**
	 * Generates random points having clustered properties using the Thomas Point
	 * Process.
	 * <p>
	 * Each cluster consists of child points normally distributed around a parent
	 * point.
	 * 
	 * @param xMin            the minimum x-coordinate of the boundary.
	 * @param yMin            the minimum y-coordinate of the boundary.
	 * @param xMax            the maximum x-coordinate of the boundary.
	 * @param yMax            the maximum y-coordinate of the boundary.
	 * @param parentsDensity  the density of parent points per unit area (scaled by
	 *                        a factor of 75x75 units).
	 * @param meanChildPoints the average number of child points generated per
	 *                        parent point (the actual values are gaussian
	 *                        distributed).
	 * @param childSpread     the first standard deviation of the distance between
	 *                        each parent point and its children.
	 * @param seed            number used to initialise the underlying pseudorandom
	 *                        number generator.
	 * @return a list of PVector objects representing the (x, y) coordinates of the
	 *         Thomas cluster points
	 * @since 2.0
	 */
	public static List<PVector> thomasClusters(double xMin, double yMin, double xMax, double yMax, double parentsDensity, double meanChildPoints,
			double childSpread, long seed) {
		ThomasPointProcess tpp = new ThomasPointProcess(seed);
		return tpp.sample(xMin, yMin, xMax, yMax, parentsDensity, meanChildPoints, childSpread);
	}

	/**
	 * Generates a set of points arranged in a phyllotaxis pattern (an arrangement
	 * similar to the florets in the head of a sunflower), using the golden ratio
	 * (the most irrational number) to position points with the least possible
	 * aliasing (which is arguably the "best" arrangement).
	 * 
	 * @param centerX x coordinate of the center of the point set
	 * @param centerY y coordinate of the center of the point set
	 * @param n       number of points to generate
	 * @param radius  radius of circular phyllotaxis extent (max distance of a point
	 *                from the center position)
	 * @return
	 */
	public static List<PVector> phyllotaxis(double centerX, double centerY, int n, double radius) {
		return phyllotaxis(centerX, centerY, n, radius, 2 * Math.PI - GOLDEN_ANGLE);
	}

	/**
	 * Generates a set of points arranged in a phyllotaxis pattern (an arrangement
	 * similar to the florets in the head of a sunflower), using a user-defined
	 * theta.
	 * 
	 * @param centerX x coordinate of the center of the point set
	 * @param centerY y coordinate of the center of the point set
	 * @param n       number of points to generate
	 * @param radius  radius of circular phyllotaxis extent (max distance of a point
	 *                from the center position)
	 * @param theta   angle (in radians) to turn after each point placement
	 * @return
	 */
	public static List<PVector> phyllotaxis(double centerX, double centerY, int n, double radius, double theta) {
		final double fillArea = radius * radius * Math.PI; // calculate area to be filled
		final double circleSpace = (fillArea / n); // area per circle
		final double fudge = 0.7; // Fudge factor: breathing space between circles
		final float circleRadius = (float) (Math.sqrt(circleSpace / Math.PI) * fudge);

		float cumArea = 0; // cumulative circle area

		final List<PVector> outList = new ArrayList<>();
		for (int i = 1; i <= n; ++i) {
			final double angle = i * theta; // rotation per circle
			cumArea += circleSpace; // add sm_area to cum_area every loop

			final double spiralR = Math.sqrt(cumArea / Math.PI); // expansion of spiral (distance of circle) per loop

			float pX = (float) (centerX + Math.cos(angle) * spiralR); // spiral rotation of golden angle per loop on X
			float pY = (float) (centerY + Math.sin(angle) * spiralR); // spiral rotation of golden angle per loop on Y

			outList.add(new PVector(pX, pY, circleRadius));
		}
		return outList;
	}

	/**
	 * Generates a set of deterministic stratified points (bounded by a rectangle)
	 * from a low discrepancy sequence (LDS) based on an irrational number (the
	 * plastic constant).
	 * <p>
	 * The <i>plastic LDS</i> has been <a href=
	 * "http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/">shown</a>
	 * to have superior low discrepancy properties amongst the quasirandom
	 * sequences, and is therefore recommended.
	 * <p>
	 * Low discrepancy sequences are deterministic (not randomized) number sequences
	 * that are low discrepancy - meaning the points tend not to clump together and
	 * leave holes; the resulting point set is more evenly spaced than a simple
	 * random distribution but less regular than lattices.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @return
	 */
	public static List<PVector> plasticLDS(double xMin, double yMin, double xMax, double yMax, int n) {
		// https://github.com/Atrix256/SampleZoo/blob/master/src/families/_2d/samples/irrational_numbers/irrational_numbers.cpp
		final double w = xMax - xMin;
		final double h = yMax - yMin;
		final double p = 1.32471795724474602596; // plastic constant
		final double a1 = 1.0 / p; // inverse of plastic number
		final double a2 = 1.0 / (p * p);

		final List<PVector> points = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			final float x = (float) (((0.5 + a1 * i) % 1) * w + xMin);
			final float y = (float) (((0.5 + a2 * i) % 1) * h + yMin);
			points.add(new PVector(x, y));
		}
		return points;
	}

	/**
	 * Generates a set of deterministic stratified points (bounded by a rectangle)
	 * from a low discrepancy sequence (LDS) based on an irrational number. In this
	 * method, a random jitter is added to points to give the point set
	 * blue-noise-like properties.
	 * <p>
	 * Low discrepancy sequences are deterministic (not randomized) number sequences
	 * that are low discrepancy - meaning the points tend not to clump together and
	 * leave holes; the resulting point set is more evenly spaced than a simple
	 * random distribution but less regular than lattices.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @return
	 * @see #plasticJitteredLDS(double, double, double, double, int, long) seeded
	 *      irrationalJitteredLDS()
	 */
	public static List<PVector> plasticJitteredLDS(double xMin, double yMin, double xMax, double yMax, int n) {
		return plasticJitteredLDS(xMin, yMin, xMax, yMax, n, System.nanoTime());
	}

	/**
	 * Generates a set of deterministic stratified points (bounded by a rectangle)
	 * from a low discrepancy sequence (LDS) based on an irrational number. In this
	 * method, a random jitter is added to points to give the point set
	 * blue-noise-like properties.
	 * <p>
	 * Low discrepancy sequences are deterministic (not randomized) number sequences
	 * that are low discrepancy - meaning the points tend not to clump together and
	 * leave holes; the resulting point set is more evenly spaced than a simple
	 * random distribution but less regular than lattices.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @param seed number used to initialize the underlying pseudorandom number
	 *             generator
	 * @return
	 * @see #plasticJitteredLDS(double, double, double, double, int) non-seeded
	 *      irrationalJitteredLDS()
	 */
	public static List<PVector> plasticJitteredLDS(double xMin, double yMin, double xMax, double yMax, int n, long seed) {
		// https://github.com/Atrix256/SampleZoo/blob/master/src/families/_2d/samples/irrational_numbers/irrational_numbers.cpp
		final double w = xMax - xMin;
		final double h = yMax - yMin;

		final SplittableRandom random = new SplittableRandom(seed);
		final double p = 1.32471795724474602596; // plastic constant
		final double a1 = 1.0 / p;
		final double a2 = 1.0 / (p * p);
		final double c_magicNumber = 0.732;

		final List<PVector> points = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			final float x = (float) (((random.nextDouble() * c_magicNumber / Math.sqrt(i + 1d) + a1 * i) % 1) * w + xMin);
			final float y = (float) (((random.nextDouble() * c_magicNumber / Math.sqrt(i + 1d) + a2 * i) % 1) * h + yMin);
			points.add(new PVector(x, y));
		}
		return points;
	}

	/**
	 * Generates a set of deterministic stratified points (bounded by a rectangle)
	 * from a low discrepancy sequence (LDS) based on a Halton sequence.
	 * <p>
	 * Low discrepancy sequences are deterministic (not randomized) number sequences
	 * that are low discrepancy - meaning the points tend not to clump together and
	 * leave holes; the resulting point set is more evenly spaced than a simple
	 * random distribution but less regular than lattices.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @return
	 */
	public static List<PVector> haltonLDS(double xMin, double yMin, double xMax, double yMax, int n) {
		final double w = xMax - xMin;
		final double h = yMax - yMin;
		float[][] values = new float[n][2];
		vanDerCorput(values, 2, 0, true, 0);
		vanDerCorput(values, 3, 1, true, 0);

		final List<PVector> points = new ArrayList<>(n);
		for (float[] point : values) {
			points.add(new PVector((float) (point[0] * w + xMin), (float) (point[1] * h + yMin)));
		}
		return points;
	}

	/**
	 * Generates a set of deterministic stratified points (bounded by a rectangle)
	 * from a low discrepancy sequence (LDS) based on a Hammersley sequence.
	 * <p>
	 * The Hammersley sequence in 2D is just the 1d Van Der Corput sequence on one
	 * axis, and regular sampling on the other axis.
	 * <p>
	 * Low discrepancy sequences are deterministic (not randomized) number sequences
	 * that are low discrepancy - meaning the points tend not to clump together and
	 * leave holes; the resulting point set is more evenly spaced than a simple
	 * random distribution but less regular than lattices.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @return
	 */
	public static List<PVector> hammersleyLDS(double xMin, double yMin, double xMax, double yMax, int n) {
		final double w = xMax - xMin;
		final double h = yMax - yMin;

		float[][] values = new float[n][2];
		vanDerCorput(values, 2, 0, false, 0);

		final float offset = 1.0f / (n * 2);
		for (int i = 0; i < n; ++i) {
			values[i][1] = offset + (i / (float) n);
		}

		final List<PVector> points = new ArrayList<>(n);
		for (float[] point : values) {
			points.add(new PVector((float) (point[0] * w + xMin), (float) (point[1] * h + yMin)));
		}
		return points;
	}

	/**
	 * Generates a set of random stratified points (bounded by a rectangle) based on
	 * the "N-Rooks" sampling pattern.
	 * <p>
	 * N-Rooks is a sampling pattern where you treat the boundary as if it were a
	 * chess board. Every sampling position is a rook that could move horizontally
	 * or vertically, and should be placed such that none of these rooks could
	 * capture/"see" any of the other rooks. In other words, every column has a
	 * single sample point in it, and every row has a single sample point in it.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @return
	 * @see #nRooksLDS(double, double, double, double, int, long)
	 */
	public static List<PVector> nRooksLDS(double xMin, double yMin, double xMax, double yMax, int n) {
		return nRooksLDS(xMin, yMin, xMax, yMax, n, System.nanoTime());
	}

	/**
	 * Generates a set of random stratified points (bounded by a rectangle) based on
	 * the "N-Rooks" sampling pattern.
	 * <p>
	 * N-Rooks is a sampling pattern where you treat the boundary as if it were a
	 * chess board. Every sampling position is a rook that could move horizontally
	 * or vertically, and should be placed such that none of these rooks could
	 * capture/"see" any of the other rooks. In other words, every column has a
	 * single sample point in it, and every row has a single sample point in it.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @param seed number used to initialize the underlying pseudorandom number
	 *             generator
	 * @return
	 * @see #nRooksLDS(double, double, double, double, int)
	 */
	public static List<PVector> nRooksLDS(double xMin, double yMin, double xMax, double yMax, int n, long seed) {
		final double w = xMax - xMin;
		final double h = yMax - yMin;

		final List<Integer> rookPositions = IntStream.range(0, n).boxed().collect(Collectors.toList());
		Collections.shuffle(rookPositions, new XoRoShiRo128PlusRandom(seed));

		final float offset = 1.0f / (n * 2);

		final List<PVector> points = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			float x = offset + (rookPositions.get(i) / (float) n);
			x *= w;
			x += xMin;
			float y = offset + i / (float) n;
			y *= h;
			y += yMin;
			points.add(new PVector(x, y));
		}

		return points;
	}

	/**
	 * Generates a 2D set of deterministic stratified points (bounded by a
	 * rectangle) from the Sobol low discrepancy sequence (LDS).
	 * <p>
	 * A Sobol sequence is a low-discrepancy sequence with the property that for all
	 * values of N,its subsequence (x1, ... xN) has a low discrepancy. It can be
	 * used to generate pseudo-randompoints in a space S, which are
	 * equi-distributed.
	 * 
	 * @param xMin x-coordinate of boundary minimum
	 * @param yMin y-coordinate of boundary minimum
	 * @param xMax x-coordinate of boundary maximum
	 * @param yMax y-coordinate of boundary maximum
	 * @param n    number of points to generate
	 * @since 1.4.0
	 * @return
	 */
	public static List<PVector> sobolLDS(double xMin, double yMin, double xMax, double yMax, int n) {
		final double w = xMax - xMin;
		final double h = yMax - yMin;
		final int dimension = 2;
		final int BITS = 52;
		final double SCALE = FastMath.pow(2, BITS);
		final long[][] direction = new long[dimension][BITS + 1];
		final long[] x = new long[dimension];
		final int[] m = new int[] { 0, 1 };
		final int a = 0;
		final int s = m.length - 1;

		for (int i = 1; i <= BITS; i++) {
			direction[0][i] = 1l << (BITS - i);
		}

		// init direction vector
		final int d = 1;
		for (int i = 1; i <= s; i++) {
			direction[d][i] = ((long) m[i]) << (BITS - i);
		}
		for (int i = s + 1; i <= BITS; i++) {
			direction[d][i] = direction[d][i - s] ^ (direction[d][i - s] >> s);
			for (int k = 1; k <= s - 1; k++) {
				direction[d][i] ^= ((a >> (s - 1 - k)) & 1) * direction[d][i - k];
			}
		}

		List<PVector> output = new ArrayList<>(n);
		for (int i = 1; i < n; i++) {

			// find the index c of the rightmost 0
			int c = 1;
			int value = i - 1;
			while ((value & 1) == 1) {
				value >>= 1;
				c++;
			}

			x[0] ^= direction[0][c];
			x[1] ^= direction[1][c];
			double vX = x[0] / SCALE;
			vX *= w;
			vX += xMin;
			double vY = x[1] / SCALE;
			vY *= h;
			vY += yMin;
			output.add(new PVector((float) vX, (float) vY));
		}

		return output;
	}

	/**
	 * Computes the <i>Euclidean minimum spanning tree</i> (EMST) of a set of
	 * points.
	 * <p>
	 * The EMST is a system of line segments, having only the given points as their
	 * endpoints, whose union includes all of the points in a connected set, and
	 * which has the minimum possible total length of any such system.
	 * 
	 * @param points the set of points over which to compute the EMST
	 * @return a LINES PShape
	 * @since 1.3.0
	 */
	public static PShape minimumSpanningTree(List<PVector> points) {
		/*
		 * The Euclidean minimum spanning tree in a plane is a subgraph of the Delaunay
		 * triangulation.
		 */
		IIncrementalTin triangulation = PGS_Triangulation.delaunayTriangulationMesh(points);
		SimpleGraph<PVector, PEdge> graph = PGS_Triangulation.toGraph(triangulation);
		SpanningTreeAlgorithm<PEdge> st = new PrimMinimumSpanningTree<>(graph); // faster than kruskal algorithm
		return PGS_SegmentSet.toPShape(st.getSpanningTree().getEdges());
	}

	/**
	 * Computes an <i>approximate</i> Traveling Salesman path for the set of points
	 * provided. Utilises a heuristic based TSP solver, followed by 2-opt heuristic
	 * improvements for further tour optimisation.
	 * <p>
	 * Note {@link PGS_Hull#concaveHullBFS(List, double) concaveHullBFS()} produces
	 * a similar result (somewhat longer tours, i.e. 10%) but is <b>much</b> more
	 * performant.
	 * 
	 * @param points the list of points for which to compute the approximate
	 *               shortest tour
	 * @return A closed polygon whose perimeter traces the shortest possible route
	 *         that visits each point exactly once (and returns to the path's
	 *         starting point).
	 * @since 2.0
	 */
	public static PShape findShortestTour(List<PVector> points) {
		var tour = new GreedyTSP<>(points, (a, b) -> a.dist(b));
		return PGS_Conversion.fromPVector(tour.getTour());
	}

	/**
	 * Applies random weights within a specified range to a list of points. The
	 * weights are assigned to the z-coordinate of each point using a random number
	 * generator with a random seed.
	 *
	 * @param points    The list of points to which random weights will be applied.
	 *                  Each point's x and y coordinates remain unchanged while the
	 *                  z coordinate is set to the random weight.
	 * @param minWeight The minimum weight value (inclusive) that can be assigned to
	 *                  a point
	 * @param maxWeight The maximum weight value (exclusive) that can be assigned to
	 *                  a point
	 * @return A new list of points where each point is a copy of the input point
	 *         with a random weight assigned to its z-coordinate
	 * @since 2.0
	 */
	public static List<PVector> applyRandomWeights(List<PVector> points, double minWeight, double maxWeight) {
		return applyRandomWeights(points, minWeight, maxWeight, System.nanoTime());
	}

	/**
	 * Applies random weights within a specified range to a list of points. The
	 * weights are assigned to the z-coordinate of each point using a random number
	 * generator with a given fixed seed.
	 *
	 * @param points    The list of points to which random weights will be applied.
	 *                  Each point's x and y coordinates remain unchanged while the
	 *                  z coordinate is set to the random weight.
	 * @param minWeight The minimum weight value (inclusive) that can be assigned to
	 *                  a point
	 * @param maxWeight The maximum weight value (exclusive) that can be assigned to
	 *                  a point
	 * @param seed      seed for the underlying random number generator
	 * @return A new list of points where each point is a copy of the input point
	 *         with a random weight assigned to its z-coordinate
	 * @since 2.0
	 */
	public static List<PVector> applyRandomWeights(List<PVector> points, double minWeight, double maxWeight, long seed) {
		final SplittableRandom random = new SplittableRandom(seed);
		return points.stream().map(p -> {
			p = p.copy();
			p.z = (float) random.nextDouble(minWeight, maxWeight);
			return p;
		}).collect(Collectors.toList());
	}

	private static List<Pair<Integer, PVector>> hilbertSortRaw(List<PVector> points) {
		double xMin, xMax, yMin, yMax;

		// find bounds
		PVector v = points.get(0);
		xMin = v.x;
		xMax = v.x;
		yMin = v.y;
		yMax = v.y;

		for (PVector PVector : points) {
			if (PVector.x < xMin) {
				xMin = PVector.x;
			} else if (PVector.x > xMax) {
				xMax = PVector.x;
			}
			if (PVector.y < yMin) {
				yMin = PVector.y;
			} else if (PVector.y > yMax) {
				yMax = PVector.y;
			}
		}

		double xDelta = xMax - xMin;
		double yDelta = yMax - yMin;
		// if (xDelta == 0 || yDelta == 0) {
		// return points;
		// }

		double hn = Math.log(points.size()) / 0.693147180559945 / 2.0;
		int nHilbert = (int) Math.floor(hn + 0.5);
		if (nHilbert < 4) {
			nHilbert = 4;
		}

		// could also use SortedMap<index -> point>
		List<Pair<Integer, PVector>> ranks = new ArrayList<>(points.size());
		double hScale = (1 << nHilbert) - 1.0;
		// scale coordinates to 2^n - 1
		for (PVector vh : points) {
			int ix = (int) (hScale * (vh.x - xMin) / xDelta);
			int iy = (int) (hScale * (vh.y - yMin) / yDelta);
			ranks.add(new Pair<>(xy2Hilbert(ix, iy, nHilbert), vh));
		}

		ranks.sort((a, b) -> Integer.compare(a.getFirst(), b.getFirst()));
		return ranks;
	}

	/**
	 * Computes the hilbert index of a coordinate on a hilbert curve of order n.
	 */
	private static int xy2Hilbert(final int px, final int py, final int n) {
		int i, xi, yi;
		int s, temp;

		int x = px;
		int y = py;
		s = 0; // Initialize.
		for (i = n - 1; i >= 0; i--) {
			xi = (x >> i) & 1; // Get bit i of x.
			yi = (y >> i) & 1; // Get bit i of y.

			if (yi == 0) {
				temp = x; // Swap x and y and,
				x = y ^ (-xi); // if xi = 1,
				y = temp ^ (-xi); // complement them.
			}
			s = 4 * s + 2 * xi + (xi ^ yi); // Append two bits to s.
		}
		return s;
	}

	/**
	 * @param values
	 * @param base
	 * @param axis         0 = x axis; 1 = y axis
	 * @param skipZero
	 * @param truncateBits
	 */
	private static void vanDerCorput(float[][] values, int base, int axis, boolean skipZero, int truncateBits) {
		// https://blog.demofox.org/2017/05/29/when-random-numbers-are-too-random-low-discrepancy-sequences/
		// https://github.com/Atrix256/SampleZoo/blob/master/src/families/_2d/samples/lds/LDS.cpp
		// figure out how many bits we are working in.
		int n = values.length;
		int value = 1;
		int numBits = 0;
		while (value < n) {
			value *= 2;
			++numBits;
		}
		int numBitsPreserved = numBits - truncateBits;
		int bitsPreservedMask = numBitsPreserved > 0 ? (1 << numBitsPreserved) - 1 : 0;

		for (int i = 0; i < n; ++i) {
			values[i][axis] = 0.0f;
			float denominator = base;
			int q = i + (skipZero ? 1 : 0);
			q &= bitsPreservedMask;
			while (q > 0) {
				int multiplier = q % base;
				values[i][axis] += multiplier / denominator;
				q = q / base;
				denominator *= base;
			}
		}
	}

	private static class CPVector implements Clusterable {
		final PVector p;
		final double[] point;

		CPVector(PVector p) {
			this.p = p;
			point = new double[] { p.x, p.y };
		}

		@Override
		public double[] getPoint() {
			return point;
		}
	}

}
