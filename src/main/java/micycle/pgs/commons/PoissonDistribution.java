package micycle.pgs.commons;

import static processing.core.PApplet.ceil;
import static processing.core.PApplet.floor;
import static processing.core.PApplet.max;
import static processing.core.PApplet.min;
import static processing.core.PApplet.sqrt;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

import processing.core.PVector;

/**
 * Generates sets of random points via <i>Poisson Disk Sampling</i>.
 * Poisson-disc sampling produces points that are tightly-packed, but no closer
 * to each other than a specified minimum distance, resulting in a more natural
 * and desirable pattern for many applications. This distribution is also
 * described as blue noise.
 * <p>
 * Implements <i>Martin Roberts’s</i> tweak of "Fast Poisson Disk Sampling in
 * Arbitrary Dimensions" by <i>Robert Bridson</i>.
 * 
 * @author Michael Carleton
 * @deprecated in favour of {@link micycle.pgs.commons.PoissonDistributionJRUS}
 */
@Deprecated
public final class PoissonDistribution {

	private ArrayList<ArrayList<PVector>> grid;
	private float cellSize;
	private int gridWidth, gridHeight;
	private float xmin, xmax, ymin, ymax;
	private ArrayList<PVector> points;
	private SplittableRandom random;

	public PoissonDistribution() {
		this(System.currentTimeMillis());
	}

	public PoissonDistribution(final long seed) {
		random = new SplittableRandom(seed);
		points = new ArrayList<>();
	}

	/**
	 * Returns the point set generated by most recent call to
	 * {@link #generate(double, double, double, double, double, int)}.
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
		return generate((float) xmin, (float) ymin, (float) xmax, (float) ymax, (float) minDist, rejectionLimit);
	}

	private List<PVector> generate(float xmin, float ymin, float xmax, float ymax, float minDist, int rejectionLimit) {
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
		final float minDistSquared = minDist * minDist;
		cellSize = (minDist) / sqrt(2); // Rmax/sqrt(d) where d is dimensions
		gridWidth = ceil((xmax - xmin) / cellSize);
		gridHeight = ceil((ymax - ymin) / cellSize);
		int s = gridWidth * gridHeight;
		grid = new ArrayList<>();
		for (int i = 0; i < s; i++) {
			grid.add(new ArrayList<>());
		}

		points.clear();
		final ArrayList<PVector> processList = new ArrayList<>();

		PVector p = new PVector(random(xmin, xmax), random(ymin, ymax));
		processList.add(p);
		points.add(p);
		addToGrid(p);

		while (!processList.isEmpty()) {
			final int i = random.nextInt(processList.size()); // Choose a point randomly from the active list, x
			p = processList.remove(i); // parent

			/*
			 * Create up to k new points uniformly at random in the spherical annulus
			 * between radii r and 2r centered on the parent point.
			 */
			final double seed = random.nextDouble();
			for (int k = 0; k < rejectionLimit; k++) {
				final PVector n = sampleAnnulus(p, minDist, k, rejectionLimit, seed);
				if (insideBoundaries(n) && testGrid(n, minDist, minDistSquared)) {
					processList.add(n);
					points.add(n);
					addToGrid(n);
				}
			}
		}

		return new ArrayList<>(points);
	}

	private boolean insideBoundaries(PVector p) {
		return (p.x >= xmin && p.x < xmax && p.y >= ymin && p.y < ymax); // keep center points in bounds
	}

	@Deprecated
	/**
	 * Create points randomly on a spherical annulus.
	 * 
	 * @param p
	 * @param minDist inner radius
	 * @param maxDist outer radius
	 * @deprecated in favour of Martin Roberts' non-uniform sampling
	 **/
	private PVector sampleAnnulus(PVector p, float minDist, float maxDist) {
		final float theta = random(0, (2 * Math.PI));
		final float r = random(minDist, maxDist);
		return new PVector(p.x + r * (float) Math.cos(theta), p.y + r * (float) Math.sin(theta));
	}

	/**
	 * Generates a point randomly from an annulus centered on p.
	 * 
	 * <p>
	 * This method uses an improved approach that is not only much faster than the
	 * naive approach, but it produces higher quality point distributions, as it
	 * allows for more tightly packed and consistent point distributions.
	 * <p>
	 * The improvement comes directly from the premise that we do not need to
	 * uniformly sample from the annulus. Rather, we would prefer to select points
	 * closer to the inner radius as this is results in neighboring points closer to
	 * the annulus center.
	 * 
	 * @param p annulus center
	 * @return
	 * @author Martin Roberts (extremelearning.com.au)
	 */
	private PVector sampleAnnulus(PVector p, float minDist, int k, int rejectionLimit, double seed) {
		double theta = 2 * Math.PI * (seed + 1.0 * k / rejectionLimit);
		float r = minDist + 0.0001f;
		return new PVector(p.x + r * (float) Math.cos(theta), p.y + r * (float) Math.sin(theta));
	}

	/**
	 * 
	 * @param p
	 * @param minDist
	 * @return true if there are no points inside the circle of minDist radius
	 *         around p
	 */
	private boolean testGrid(PVector p, float minDist, float minDistSquared) {
		int minX = floor(max(0, (p.x - minDist - xmin) / cellSize));
		int maxX = ceil(min(gridWidth - 1f, (p.x + minDist - xmin) / cellSize));
		int minY = floor(max(0, (p.y - minDist - ymin) / cellSize));
		int maxY = ceil(min(gridHeight - 1f, (p.y + minDist - ymin) / cellSize));

		for (int y = minY; y <= maxY; y++) {
			for (int x = minX; x <= maxX; x++) {
				ArrayList<PVector> cell = grid.get(y * gridWidth + x);
				for (PVector t : cell) {
					float dx = p.x - t.x;
					float dy = p.y - t.y;
					if (dx * dx + dy * dy <= minDistSquared) {
						return false;
					}
				}
			}
		}

		return true;
	}

	private void addToGrid(PVector p) {
		grid.get(index(p.x, p.y)).add(p);
	}

	private int index(float x, float y) {
		int gx = floor((x - xmin) / cellSize);
		int gy = floor((y - ymin) / cellSize);
		return gy * gridWidth + gx;
	}

	/**
	 * @param min - The minimum.
	 * @param max - The maximum.
	 * @return A random double between these numbers (inclusive the minimum and
	 *         maximum).
	 */
	private float random(double min, double max) {
		return (float) (min + (max - min) * random.nextDouble());
	}

}