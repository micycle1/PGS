package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import it.unimi.dsi.util.XoRoShiRo128PlusRandomGenerator;
import processing.core.PVector;

/**
 * Implements the Thomas Point Process, a stochastic process used for simulating
 * clusters of points around parent locations in spatial analysis. Each cluster
 * consists of child points normally distributed around a parent point.
 * <p>
 * This Java implementation is a port of the simple Python implementation
 * available <a href=
 * "https://hpaulkeeler.com/simulating-a-thomas-cluster-point-process/">here</a>,
 * wherein parent points are seeded uniformly.
 * 
 * @author Michael Carleton
 */
public class ThomasPointProcess {

	// also https://github.com/For-a-few-DPPs-more/structure-factor/

	// even better impl here?
	// https://github.com/spatstat/spatstat.random/blob/main/src/rthomas.h

	private XoRoShiRo128PlusRandomGenerator random;
	private final long seed;

	/**
	 * Initalise with a random seed.
	 */
	public ThomasPointProcess() {
		this(System.nanoTime());
	}

	/**
	 * Initialise with a known seed.
	 * 
	 * @param seed
	 */
	public ThomasPointProcess(long seed) {
		this.seed = seed;
	}

	/**
	 * Generates a list of sample points within a specified rectangular boundary
	 * with clustering behavior. The method simulates a clustering pattern where
	 * each 'parent' point can spawn multiple 'child' points.
	 * 
	 * @param xMin            the minimum x-coordinate of the boundary.
	 * @param yMin            the minimum y-coordinate of the boundary.
	 * @param xMax            the maximum x-coordinate of the boundary.
	 * @param yMax            the maximum y-coordinate of the boundary.
	 * @param parentsDensity  the density of parent points per unit area (scaled by
	 *                        a factor of 75x75 units).
	 * @param meanChildPoints the average number of child points generated per
	 *                        parent point.
	 * @param childSpread     the first standard deviation of the distance between
	 *                        each parent point and its children.
	 * @return a list of PVector objects representing the generated points within
	 *         the boundary.
	 */
	public List<PVector> sample(double xMin, double yMin, double xMax, double yMax, double parentsDensity, double meanChildPoints,
			double childSpread) {
		random = new XoRoShiRo128PlusRandomGenerator(seed);
		final double boundaryBuffer = 0; // 6 * childSpread;
		final double xMinExt = xMin - boundaryBuffer;
		final double xMaxExt = xMax + boundaryBuffer;
		final double yMinExt = yMin - boundaryBuffer;
		final double yMaxExt = yMax + boundaryBuffer;

		// Rectangle dimensions
		final double xDeltaExt = xMaxExt - xMinExt;
		final double yDeltaExt = yMaxExt - yMinExt;
		final double areaTotalExt = xDeltaExt * yDeltaExt; // area of extended rectangle

		// seed a parent every 75 units on average
		// choosen as suitable density for val=1
		int numParentPoints = (int) ((areaTotalExt / (75 * 75)) * parentsDensity);

		// x and y coordinates of Poisson points for the parent
		final double[] xxParent = new double[numParentPoints];
		final double[] yyParent = new double[numParentPoints];
		for (int i = 0; i < numParentPoints; i++) {
			xxParent[i] = xMinExt + xDeltaExt * random.nextDouble();
			yyParent[i] = yMinExt + yDeltaExt * random.nextDouble();
		}

		/*
		 * Reinitalise random generator so that changing the density (more points) does
		 * not affect the locations of children.
		 */
		random = new XoRoShiRo128PlusRandomGenerator(Long.MAX_VALUE - seed);

		Set<PVector> points = new HashSet<>((int) (numParentPoints * meanChildPoints));

		for (int i = 0; i < numParentPoints; i++) {
			int numChildren = (int) gaussian(meanChildPoints, Math.sqrt(meanChildPoints));
			for (int j = 0; j < numChildren; j++) {
				double x = xxParent[i] + gaussian(0, childSpread);
				double y = yyParent[i] + gaussian(0, childSpread);
				if (x >= xMin && x <= xMax && y >= yMin && y <= yMax) {
					points.add(new PVector((float) x, (float) y));
				}
			}
		}

		return new ArrayList<>(points);
	}

	private final double gaussian(final double mean, final double stdDev) {
		return random.nextGaussian() * stdDev + mean;
	}

	private final double gaussianFast(final double mean, final double stdDev) {
		final double DELTA = 1.0 / 4294967296.0; // (1 / 2^32)
		long u = random.nextLong();

		// Split into 2 x 32 bits
		long major = (u >> 32); // Upper 32 bits
		long minor = (u & 0xFFFFFFFFL); // Lower 32 bits

		double x = Integer.bitCount((int) major); // x = random binomially distributed integer from 0 to 32
		x += minor * DELTA; // Linearly fill the gaps between integers
		x -= 16.5; // Re-center around 0 (the mean should be 16 + 0.5)
		x *= 0.3535534; // Scale to ~1 standard deviation

		return mean + stdDev * x;
	}
}