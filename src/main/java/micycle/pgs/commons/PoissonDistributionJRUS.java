package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

import processing.core.PVector;

/**
 * Generates sets of random points via <i>Poisson Disk Sampling</i>.
 * Poisson-disc sampling produces points that are tightly-packed, but no closer
 * to each other than a specified minimum distance, resulting in a natural and
 * desirable pattern for many applications. This distribution is also described
 * as blue noise.
 * <p>
 * The algorithm in this class is a Fork of <i>Martin Roberts’s</i> tweak to
 * <i>Bridson's Algorithm</i> for Poisson Disk sampling. This approach is faster
 * and better than the Bridson Algorithm, and balances performance with
 * distribution quality compared to Robert's tweak.
 * 
 * <p>
 * For more, see this <a href=
 * "https://observablehq.com/@fil/poisson-distribution-generators">anaylsis</a>
 * of different Poisson disk sampling functions.
 * 
 * @author Jacob Rus
 * @author Java port by Michael Carleton
 */
public final class PoissonDistributionJRUS {

	/*-
	 * Implements:
	 * 		https://observablehq.com/@jrus/bridson-fork/2
	 * 		https://observablehq.com/@fil/poisson-distribution-generators#C
	 */

	private double[] grid;
	private double cellSize;
	private int gridWidth;
	private List<double[]> queue;
	private SplittableRandom random;
	private float xOffset, yOffset;

	private List<PVector> points;

	public PoissonDistributionJRUS() {
		this(System.currentTimeMillis());
	}

	public PoissonDistributionJRUS(final long seed) {
		random = new SplittableRandom(seed);
		points = new ArrayList<>();
	}

	/**
	 * Returns the point set generated by most recent call to
	 * {@link #generate(double, double, double, double, double, int) generate()}.
	 * 
	 * @return
	 */
	public List<PVector> getPoints() {
		return points;
	}

	/**
	 * Generates a random point set, having a poisson/blue noise distribution.
	 * 
	 * @param xmin           x-coordinate of boundary minimum
	 * @param ymin           y-coordinate of boundary minimum
	 * @param xmax           x-coordinate of boundary maximum
	 * @param ymax           y-coordinate of boundary maximum
	 * @param minDist        minimum euclidean distance between any two points
	 * @param rejectionLimit the limit on the number of attempts to generate a
	 *                       random valid point around the previous point. Generally
	 *                       6 is sufficient.
	 * @return a set of random points
	 */
	public List<PVector> generate(double xmin, double ymin, double xmax, double ymax, double minDist, int rejectionLimit) {
		xOffset = (float) xmin;
		yOffset = (float) ymin;
		return generate(xmax - xmin, ymax - ymin, minDist, rejectionLimit);
	}

	/**
	 * Generates a random point set, having a poisson/blue noise distribution.
	 * 
	 * @param xmin    x-coordinate of boundary minimum
	 * @param ymin    y-coordinate of boundary minimum
	 * @param xmax    x-coordinate of boundary maximum
	 * @param ymax    y-coordinate of boundary maximum
	 * @param minDist minimum euclidean distance between any two points
	 * @return a set of random points
	 */
	public List<PVector> generate(double xmin, double ymin, double xmax, double ymax, double minDist) {
		return generate(xmin, ymin, xmax, ymax, minDist, 11);
	}

	private List<PVector> generate(double width, double height, double radius, int k) {
		int m = 1 + k * 2; // a number mutually prime to k
		cellSize = 1 / (radius * Math.sqrt(0.5));
		final double minDistSquared = radius * radius;
		/*
		 * Pad the grid on the sides to eliminate the need for special-case code near
		 * edges.
		 */
		gridWidth = (int) Math.ceil(width * cellSize) + 4;
		int gridHeight = (int) Math.ceil(height * cellSize) + 4;
		grid = new double[2 * gridWidth * gridHeight];
		queue = new ArrayList<>();
		final double rotx = Math.cos((2 * Math.PI * m) / k);
		final double roty = Math.sin((2 * Math.PI * m) / k);

		points.clear();

		sample(width * random(0.45, 0.55), height * random(0.45, 0.55));

		pick: while (!queue.isEmpty()) {
			final int i = random.nextInt(queue.size()); // Choose a point randomly from the active list, x
			double[] parent = queue.get(i); // parent

			final double epsilon = 1e-6;
			final double t = tanpi2(2 * random.nextDouble() - 1);
			final double q = 1d / (1 + t * t);
			double dw;
			double dx = q != 0 ? (1 - t * t) * q : -1;
			double dy = q != 0 ? 2 * t * q : 0;

			for (int j = 0; j < k; ++j) {
				dw = dx * rotx - dy * roty; // temporary name for dx
				dy = dx * roty + dy * rotx;
				dx = dw;

				final double r = radius * (1 + epsilon + 0.65 * random.nextDouble() * random.nextDouble());
				final double x = parent[0] + r * dx;
				final double y = parent[1] + r * dy;

				// Accept candidates that are inside the allowed extent
				// and farther than 2 * radius to all existing samples.
				if ((0 <= x) && (x < width) && (0 <= y) && (y < height) && far(x, y, minDistSquared)) {
					sample(x, y);
//					continue pick; // NOTE no continue is faster and negligibly worse
				}
			}
			queue.remove(i);
		}

		return new ArrayList<>(points);
	}

	/**
	 * Fast approximation of tan(πx/2).
	 */
	private static double tanpi2(double a) {
		double b = (1 - a * a);
		return a * (-0.0187108 * b + 0.31583526 + 1.27365776 / b);
	}

	private void sample(double x, double y) {
		// offset grid cell by 2 in each direction to account for border
		int i = (int) Math.floor(x * cellSize + 2);
		int j = (int) Math.floor(y * cellSize + 2);
		int index = 2 * (gridWidth * j + i);
		grid[index] = x;
		grid[index + 1] = y;
		queue.add(new double[] { x, y });
		points.add(new PVector((float) x + xOffset, (float) y + yOffset));
	}

	private boolean far(double x, double y, double minDistSquared) {
		int j0 = (int) Math.floor(y * cellSize);
		int i0 = (int) Math.floor(x * cellSize);
		for (int j = j0; j < j0 + 5; ++j) {
			int index0 = 2 * (j * gridWidth + i0);
			for (int i = index0; i < index0 + 10; i += 2) {
				double dx = grid[i] - x;
				double dy = grid[i + 1] - y;
				if (dx * dx + dy * dy < minDistSquared)
					return false;
			}
		}
		return true;
	}

	/**
	 * @param min - The minimum.
	 * @param max - The maximum.
	 * @return A random double between these numbers (inclusive the minimum and
	 *         maximum).
	 */
	private double random(double min, double max) {
		return (min + (max - min) * random.nextDouble());
	}

}