package micycle.pgs.utility;

import static processing.core.PApplet.sin;
import static processing.core.PApplet.cos;
import static processing.core.PApplet.min;
import static processing.core.PApplet.max;
import static processing.core.PApplet.floor;
import static processing.core.PApplet.ceil;
import static processing.core.PApplet.sqrt;
import static processing.core.PApplet.dist;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import processing.core.PApplet;
import processing.core.PVector;

/**
 * From "Fast Poisson Disk Sampling in Arbitrary Dimensions" by Robert Bridson
 * http://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
 * 
 * Fast Poisson-Disc Sampling with constant r
 * 
 * a.ka blue noise (generated via disk sampling "dart throwing" but with
 * optimisation)
 *
 */
public final class PoissonDistribution {

	// TODO Phyllotaxis distrib

	private ArrayList<ArrayList<PVector>> grid;
	private float cellSize;
	private int gridWidth, gridHeight;
	private float xmin, xmax, ymin, ymax;
	private ArrayList<PVector> points;
	private Random random;

	public PoissonDistribution() {
		this(System.currentTimeMillis());
	}

	public PoissonDistribution(long seed) {
		random = new Random(seed);
		points = new ArrayList<>();
	}

	public List<PVector> getPoints() {
		return points;
	}

	public List<PVector> generate(double xmin, double ymin, double xmax, double ymax, double minDist,
			int rejectionLimit) {
		return generate((float) xmin, (float) ymin, (float) xmax, (float) ymax, (float) minDist, rejectionLimit);
	}

	/**
	 * rejection limit sometimes known as k
	 **/
	public List<PVector> generate(float xmin, float ymin, float xmax, float ymax, float minDist,
			int rejectionLimit) {
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
		cellSize = minDist / sqrt(2); // Rmax/sqrt(d) where d is dimensions
		gridWidth = ceil((xmax - xmin) / cellSize);
		gridHeight = ceil((ymax - ymin) / cellSize);
		int s = gridWidth * gridHeight;
		grid = new ArrayList<ArrayList<PVector>>();
		for (int i = 0; i < s; i++)
			grid.add(new ArrayList<>());

		points.clear();
		LinkedList<PVector> processList = new LinkedList<PVector>(); // active list

		PVector p = new PVector(random(xmin, xmax), random(ymin, ymax)); // ADJACENT SQUARES NEED A
																			// COMMON VALUE
		processList.add(p);
		points.add(p);
		addToGrid(p);

		while (processList.size() > 0) {
			int i = floor(random(0, processList.size())); // Choose a point randomly from the active list, x
			p = processList.get(i);
			processList.remove(i);
			for (i = 0; i < rejectionLimit; i++) // Create k new points uniformly at random in the spherical annulus
													// between radii r and 2r centered on xi
			{
				PVector n = createRandomPointAround(p, minDist, minDist * 2, i);
				if (insideBoundaries(n, minDist / 2) && testGrid(n, minDist)) {
					processList.add(n);
					points.add(n);
					addToGrid(n);
				}
			}
		}

		return new ArrayList<>(points);
	}

	private boolean insideBoundaries(PVector p, float border) {

		return (p.x >= xmin && p.x < xmax && p.y >= ymin && p.y < ymax); // keep center points in bounds
		// return (p.x >= _xmin+border && p.x < _xmax-border && p.y >= _ymin+border &&
		// p.y < _ymax-border); // keep circles in bounds
	}

	/**
	 * Create points randomly on a spherical annulus.
	 **/
	private PVector createRandomPointAround(PVector p, float minDist, float maxDist, int iteration) {
		float a = random(0, (float) (2 * Math.PI));
		float r = random(minDist, maxDist);
		// float a = rand((long)(p.x*iteration),(long)(p.y*iteration))*2*PI;
		// float r = minDist +
		// rand((long)(p.x*iteration),(long)(p.y*iteration))*minDist;
		return new PVector(p.x + r * cos(a), p.y + r * sin(a));
	}

	/**
	 * 
	 * @param p
	 * @param minDist
	 * @return true if there are no points inside the circle of minDist radius
	 *         around p
	 */
	private boolean testGrid(PVector p, float minDist) {
		int minX = floor(max(0, (p.x - minDist - xmin) / cellSize));
		int maxX = ceil(min(gridWidth - 1, (p.x + minDist - xmin) / cellSize));
		int minY = floor(max(0, (p.y - minDist - ymin) / cellSize));
		int maxY = ceil(min(gridHeight - 1, (p.y + minDist - ymin) / cellSize));

		for (int y = minY; y <= maxY; y++) {
			for (int x = minX; x <= maxX; x++) {
				ArrayList<PVector> cell = grid.get(y * gridWidth + x);
				for (PVector t : cell)
					if (dist(p.x, p.y, t.x, t.y) <= minDist) {
						return false;
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

//	private static float random(double min, double max) {
//		return (float) ThreadLocalRandom.current().nextDouble(min, max);
//	}

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