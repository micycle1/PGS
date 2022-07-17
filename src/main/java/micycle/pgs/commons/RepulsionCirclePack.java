package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import processing.core.PVector;

/**
 * Repulsion Packing attempts to arrange a set of circles of specified radii
 * within a rectangle such that there is no-overlap between circles.
 * 
 * <p>
 * It involves iterative pair-repulsion, in which overlapping circles move away
 * from each other. The distance moved by each circle is proportional to the
 * radius of the other to approximate inertia (very loosely), so that when a
 * small circle is overlapped by a large circle, the small circle moves
 * furthest. This process is repeated iteratively until no more movement takes
 * place (acceptable layout) or a maximum number of iterations is reached
 * (layout failure). To avoid edge effects, the bounding rectangle is treated as
 * a toroid. Each circle's centre is constrained to lie within the rectangle but
 * its edges are allowed to extend outside.
 * <p>
 * This forms a Java implementation of the algorithm from the <i>packcircles</i>
 * R package.
 * 
 * @author Michael Carleton
 *
 */
public class RepulsionCirclePack {

	// based on
	// https://github.com/mbedward/packcircles/blob/main/src/packcircles.cpp

	/*-
	 * TODO:
	 * R-tree (or similar) to do intersection checks,
	 * PShape version: place circles with weight=1 on perimeter, then pack circles inside shape
	 * https://web.archive.org/web/20080517011454/http://en.wiki.mcneel.com/default.aspx/McNeel/2DCirclePacking
	 */

	private static boolean USE_GRID = true;
	private static double almostZero_TOL = 1e-2;
	private static double MAX_ITER = 1000;

	private List<PVector> circles;
	/**
	 * Multiplicative weights for the distance a circle will move with pair
	 * repulsion [0...1]
	 */
	private List<Double> weights;

	double cellSize;

	double xmin;
	double xmax;
	double ymin;
	double ymax;

	int gridWidth, gridHeight;
	List<List<PVector>> grid;
	Map<PVector, float[]> updates;

	/**
	 * 
	 * @param circles overlapping circles, circle positions and sizes (x, y, radius)
	 * @param xmin
	 * @param xmax
	 * @param ymin
	 * @param ymax
	 * @param wrap
	 */
	public RepulsionCirclePack(List<PVector> xyr, double xmin, double xmax, double ymin, double ymax, final boolean wrap) {
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
		this.circles = xyr.stream().map(v -> v.copy()).collect(Collectors.toList()); // TODO hilbert sort?
//		this.circles = PGS_PointSet.hilbertSort(circles);
//		circles.sort((a,b) -> Float.compare(a.x+a.y, b.x+b.y));

		if (USE_GRID) {

			float maxR = circles.stream().max((a, b) -> Float.compare(a.z, b.z)).orElse(new PVector(0,0,1)).z;
			cellSize = (maxR * 2 + almostZero_TOL);
			gridWidth = (int) Math.ceil((xmax - xmin) / cellSize) + 4;
			gridHeight = (int) Math.ceil((ymax - ymin) / cellSize) + 4;
			int s = gridWidth * gridHeight;
			grid = new ArrayList<>(s);
			for (int i = 0; i < s; i++) {
				grid.add(new ArrayList<>());
			}
			circles.forEach(this::addToGrid);
//		System.out.println(grid == null);
			updates = new HashMap<>();
		}

		iterateLayout(xmin, xmax, ymin, ymax, wrap);
//		System.out.println(minR);
//		System.out.println(gridHeight*gridWidth);
	}

	private void addToGrid(PVector p) {
		grid.get(index(p.x, p.y)).add(p);
	}

	private int index(double x, double y) {
//offset grid cell by 2 in each direction to account for border
		int gx = (int) Math.floor((x - xmin) / cellSize) + 2;
		int gy = (int) Math.floor((y - ymin) / cellSize) + 2;
		return gy * gridWidth + gx;
	}

	private int indexY(double x, double y) {
		return (int) Math.floor((y - ymin) / cellSize) + 2;
	}

	private int indexX(double x, double y) {
		return (int) Math.floor((x - xmin) / cellSize) + 2;
	}

	private List<PVector> atIndex(int indexX, int indexY) {
		return grid.get(indexY * gridWidth + indexX);
	}

	void buildGrid() {
	}

	public List<PVector> getPacking() {
		return circles;
	}

	/**
	 * Given an input matrix of circle positions and sizes, attempts to position
	 * them without overlap by iterating the pair-repulsion algorithm.
	 * 
	 * @param circles 3 column matrix of
	 * @param weights vector of double values between 0 and 1, used as
	 *                multiplicative weights for the distance a circle will move
	 *                with pair repulsion
	 * @param xmin    bounds min X
	 * @param xmax    bounds max X
	 * @param ymin    bounds min Y
	 * @param ymax    bounds max Y
	 * @param maxiter maximum number of iterations
	 * @param wrap    allow coordinate wrapping across opposite bounds
	 * @return the number of iterations performed
	 */
	private int iterateLayout(double xmin, double xmax, double ymin, double ymax, final boolean wrap) {
		final int N = circles.size();
		int iter = 0;

		while (iter++ < MAX_ITER) {
			boolean moved = false;

			if (!USE_GRID) {
				for (int i = 0; i < N - 1; ++i) {
					for (int j = i + 1; j < N; ++j) {
						if (doRepulsion(circles.get(i), circles.get(j), xmin, xmax, ymin, ymax, wrap, false)) {
							moved = true;
						}
					}
				}
			} else {
				for (PVector circle : circles) {
//				int index = index(circle.x, circle.y);
					int indexX = indexX(circle.x, circle.y);
					int indexY = indexY(circle.x, circle.y);
					// MARK WHETHER we've already done the repulsion in the loop
					// since we do everything twice
					for (int i = indexX - 1; i < indexX + 1; i++) {
						for (int j = indexY - 1; j < indexY + 2; j++) { // NOTE +2
							for (PVector neighbor : atIndex(i, j)) {
								if (neighbor != circle && doRepulsion(circle, neighbor, xmin, xmax, ymin, ymax, wrap, true)) {
									moved = true;
								}
							}
							updateGrid();
						}
					}
				}
			}

			if (!moved) {
				break;
			}

		}

//		System.out.println(iter);
		return iter;
	}

	void updateGrid() {
		updates.forEach((c, pos) -> {
			int index = index(c.x, c.y); // original
			grid.get(index).remove(c);
			c.x = pos[0];
			c.y = pos[1];
			addToGrid(c); // reinsert at correct position
		});
		updates.clear();
	}

	/**
	 * Checks if two circles overlap excessively and, if so, moves them apart. The
	 * distance moved by each circle is proportional to the radius of the other to
	 * give some semblance of intertia.
	 * 
	 * @param c0   index of first circle
	 * @param c1   index of second circle
	 * @param xmin bounds min X
	 * @param xmax bounds max X
	 * @param ymin bounds min Y
	 * @param ymax bounds max Y
	 * @param wrap allow coordinate wrapping across opposite bounds
	 */
	private boolean doRepulsion(PVector c0, PVector c1, double xmin, double xmax, double ymin, double ymax, final boolean wrap,
			boolean useGrid) {

		// TODO every iteration, place into triangulation, repulse only connected
		// circles; find based on pinwheel ?

		// if both weights are zero, return zero to indicate
		// no movement
//		if (almostZero(weights.get(c0)) && almostZero(weights.get(c1))) {
//			return 0;
//		}

		double dx = c1.x - c0.x;
		double dy = c1.y - c0.y;
		double d = Math.sqrt(dx * dx + dy * dy);
		double r = c1.z + c0.z;
		double p;
		double w0;
		double w1;

		if (gtZero(r - d)) {
			if (almostZero(d)) {
				// The two centres are coincident or almost so.
				// Arbitrarily move along x-axis
				p = 1.0;
				dx = r - d;
			} else {
				p = (r - d) / d;
			}

			w0 = c1.z / r;
			w1 = c0.z / r;
			// w0*=weights.get(c0); // NOTE circle weights disabled
			// w1*=weights.get(c1); // NOTE circle weights disabled

			float c1XNew = ordinate(c1.x + p * dx * w1, xmin, xmax, wrap);
			float c1YNew = ordinate(c1.y + p * dy * w1, ymin, ymax, wrap);
			int c1Index = index(c1.x, c1.y); // original
			if (useGrid && c1Index != index(c1XNew, c1YNew)) {
				// TODO ideally update position immediately
				updates.put(c1, new float[] { c1XNew, c1YNew });
			} else {
				// update coords but not grid reference
				c1.x = c1XNew;
				c1.y = c1YNew;
			}

			float c0XNew = ordinate(c0.x - p * dx * w0, xmin, xmax, wrap);
			float c0YNew = ordinate(c0.y - p * dy * w0, ymin, ymax, wrap);

			int c0Index = index(c0.x, c0.y); // original
			if (useGrid && c0Index != index(c0XNew, c0YNew)) {
				updates.put(c0, new float[] { c0XNew, c0YNew });
			} else {
				// update coords but not grid reference
				c0.x = c0XNew;
				c0.y = c0YNew;
			}

			return true;
		}

		return false;
	}

	class Circle {
		float x, y;
		float r;
		int index; // used to reference to the object (since x,y can change) and update its grid
					// reference
		int provisionalGridIndex;
		
		@Override
		public int hashCode() {
			return index;
		}
	}

	private static boolean almostZero(double x) {
		return Math.abs(x) < almostZero_TOL;
	}

	private static boolean gtZero(double x) {
		return !almostZero(x) && (x > 0.0);
	}

	/**
	 * Adjust an X or Y ordinate to the given bounds by either wrapping (if `wrap`
	 * is true) or clamping (if `wrap` is false).
	 */

	private static float ordinate(double x, double lo, double hi, final boolean wrap) {
		if (wrap) {
			return (float) wrapOrdinate(x, lo, hi);
		} else {
			return (float) Math.max(lo, Math.min(hi, x));
		}
	}

	/**
	 * Map an X or Y ordinate to the toroidal interval [lo, hi).
	 *
	 * x - X or Y ordinate to be adjusted lo - lower coordinate bound hi - upper
	 * coordinate bound
	 */
	private static double wrapOrdinate(double x, double lo, double hi) {
		double w = hi - lo;
		while (x < lo) {
			x += w;
		}
		while (x >= hi) {
			x -= w;
		}
		return x;
	}

}
