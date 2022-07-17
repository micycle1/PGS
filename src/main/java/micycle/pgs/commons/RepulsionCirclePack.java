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
 * This Java code is based on an implementation of the algorithm from the
 * <i>packcircles</i> R package, but adds grid-based indexing to speed up the
 * packing convergence.
 * 
 * @author Michael Carleton
 *
 */
public class RepulsionCirclePack {

	// based on
	// https://github.com/mbedward/packcircles/blob/main/src/packcircles.cpp

	// consider
	// https://web.archive.org/web/20080517011454/http://en.wiki.mcneel.com/default.aspx/McNeel/2DCirclePacking

	private static boolean USE_GRID = true;
	private static double almostZero_TOL = 1e-2;
	private static double MAX_ITER = 1000;

	private List<PVector> circles;
	/**
	 * Multiplicative weights for the distance a circle will move with pair
	 * repulsion [0...1]
	 */
	private List<Double> weights;

	private double cellSize;
	private final float xmin;
	private final float xmax;
	private final float ymin;
	private final float ymax;
	private final boolean wrap;

	private int gridWidth, gridHeight;
	private List<List<PVector>> grid;
	private Map<PVector, float[]> updates;

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
		this.circles = xyr.stream().map(PVector::copy).collect(Collectors.toList());
		this.xmin = (float) xmin;
		this.xmax = (float) xmax;
		this.ymin = (float) ymin;
		this.ymax = (float) ymax;
		this.wrap = wrap;

		if (USE_GRID) {
			float maxR = circles.stream().max((a, b) -> Float.compare(a.z, b.z)).orElse(new PVector(0, 0, 2)).z;
			cellSize = (maxR * 2 + almostZero_TOL);
			gridWidth = (int) Math.ceil((xmax - xmin) / cellSize) + 4;
			gridHeight = (int) Math.ceil((ymax - ymin) / cellSize) + 4;
			int s = gridWidth * gridHeight;
			grid = new ArrayList<>(s);
			for (int i = 0; i < s; i++) {
				grid.add(new ArrayList<>());
			}
			circles.forEach(this::addToGrid);
			updates = new HashMap<>();
		}

		iterateLayout();
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
	 * @param maxiter maximum number of iterations
	 * 
	 * @return the number of iterations performed
	 */
	private int iterateLayout() {
		final int N = circles.size();
		int iter = 0;

		while (iter++ < MAX_ITER) {
			boolean moved = false;

			if (!USE_GRID) {
				for (int i = 0; i < N - 1; ++i) {
					for (int j = i + 1; j < N; ++j) {
						if (doRepulsion(circles.get(i), circles.get(j), false)) {
							moved = true;
						}
					}
				}
			} else {
				for (PVector circle : circles) {
					int indexX = indexX(circle.x);
					int indexY = indexY(circle.y);
					// MARK WHETHER we've already done the repulsion in the loop
					// since we do everything twice
					for (int i = indexX - 1; i < indexX + 1; i++) {
						for (int j = indexY - 1; j < indexY + 2; j++) { // NOTE +2
							for (PVector neighbor : atIndex(i, j)) {
								if (neighbor != circle && doRepulsion(circle, neighbor, true)) {
									moved = true;
								}
							}
							updateGrid();
						}
					}
				}
			}

			if (!moved) { // didn't move any -- no circles overlapping
				break;
			}

		}

		return iter;
	}

	/**
	 * Checks if two circles overlap excessively and, if so, moves them apart. The
	 * distance moved by each circle is proportional to the radius of the other to
	 * give some semblance of intertia.
	 * 
	 * @param c0 first circle
	 * @param c1 second circle
	 */
	private boolean doRepulsion(final PVector c0, final PVector c1, final boolean useGrid) {
		// if both weights are zero, return zero to indicate no movement
//		if (almostZero(weights.get(c0)) && almostZero(weights.get(c1))) {
//			return 0;
//		}

		float dx = c1.x - c0.x;
		float dy = c1.y - c0.y;
		float d = (float) Math.sqrt(dx * dx + dy * dy);
		float r = c1.z + c0.z;
		float p;
		float w0;
		float w1;

		if (gtZero(r - d)) {
			if (almostZero(d)) {
				// The two centres are coincident or almost so.
				// Arbitrarily move along x-axis
				p = 1.0f;
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
			int c1Index = index(c1.x, c1.y);
			if (useGrid && c1Index != index(c1XNew, c1YNew)) {
				/*
				 * Can't change since we're currently iterating one of the grid cells
				 */
				updates.put(c1, new float[] { c1XNew, c1YNew });
			} else { // update coords but not grid reference
				c1.x = c1XNew;
				c1.y = c1YNew;
			}

			float c0XNew = ordinate(c0.x - p * dx * w0, xmin, xmax, wrap);
			float c0YNew = ordinate(c0.y - p * dy * w0, ymin, ymax, wrap);

			int c0Index = index(c0.x, c0.y);
			if (useGrid && c0Index != index(c0XNew, c0YNew)) {
				updates.put(c0, new float[] { c0XNew, c0YNew });
			} else { // update coords but not grid reference
				c0.x = c0XNew;
				c0.y = c0YNew;
			}

			return true;
		}

		return false;
	}

	private void addToGrid(PVector p) {
		grid.get(index(p.x, p.y)).add(p);
	}

	private void updateGrid() {
		updates.forEach((c, pos) -> {
			int index = index(c.x, c.y); // original
			grid.get(index).remove(c);
			c.x = pos[0];
			c.y = pos[1];
			addToGrid(c); // reinsert at correct position
		});
		updates.clear();
	}

	private int index(double x, double y) {
		// offset grid cell by 2 in each direction to account for border
		int gx = (int) Math.floor((x - xmin) / cellSize) + 2;
		int gy = (int) Math.floor((y - ymin) / cellSize) + 2;
		return gy * gridWidth + gx;
	}

	private int indexY(double y) {
		return (int) Math.floor((y - ymin) / cellSize) + 2;
	}

	private int indexX(double x) {
		return (int) Math.floor((x - xmin) / cellSize) + 2;
	}

	private List<PVector> atIndex(int indexX, int indexY) {
		return grid.get(indexY * gridWidth + indexX);
	}

	private static boolean almostZero(float x) {
		return Math.abs(x) < almostZero_TOL;
	}

	private static boolean gtZero(float x) {
		return !almostZero(x) && (x > 0.0);
	}

	/**
	 * Adjust an X or Y ordinate to the given bounds by either wrapping (if `wrap`
	 * is true) or clamping (if `wrap` is false).
	 */
	private static float ordinate(float x, float lo, float hi, final boolean wrap) {
		if (wrap) {
			return wrapOrdinate(x, lo, hi);
		} else {
			return Math.max(lo, Math.min(hi, x));
		}
	}

	/**
	 * Map an X or Y ordinate to the toroidal interval [lo, hi].
	 *
	 * x - X or Y ordinate to be adjusted lo - lower coordinate bound hi - upper
	 * coordinate bound
	 */
	private static float wrapOrdinate(float x, float lo, float hi) {
		float w = hi - lo;
		while (x < lo) {
			x += w;
		}
		while (x >= hi) {
			x -= w;
		}
		return x;
	}

}
