package micycle.pgs.commons;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Deque;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygonal;

import com.github.micycle1.geoblitz.PointDistanceIndex;

/**
 * Computes a sequence of <em>largest empty circles</em> (LECs) whose centers
 * are constrained to lie within a polygonal {@code boundary} (holes respected).
 * <p>
 * The "emptiness" constraint is defined against:
 * <ul>
 * <li>the boundary rings (outer shell and holes), and</li>
 * <li>optional {@code obstacles}.</li>
 * </ul>
 *
 * <h2>Obstacles</h2> The optional {@code obstacles} geometry may contain a
 * mixture of:
 * <ul>
 * <li><b>Polygonal</b> components: treated as excluded regions (additional
 * holes). They affect the sign of the distance (points inside an obstacle
 * polygon are considered outside the feasible region).</li>
 * <li><b>Linear</b> components (lineal): contribute to the distance target
 * (circles must not cross them).</li>
 * <li><b>Puntal</b> components (pointal): contribute to the distance target
 * (circles must not cover them).</li>
 * </ul>
 *
 * <h2>Iteration and reuse</h2> Repeated calls to {@link #findNextLEC()} reuse
 * and refine a cached set of candidate cells, making successive extractions
 * faster than recomputing from scratch.
 */
public class LargestEmptyCircles {

	private final Geometry boundary; // polygonal
	private final Geometry obstacles; // nullable; may be any Geometry
	private final double tolerance;

	private PointDistanceIndex boundaryDistance;

	private Envelope gridEnv;
	private Cell farthestCell;

	private final Deque<Cell> cellStack = new ArrayDeque<>();
	private final List<Cell> nextIterCells = new ArrayList<>(4096);

	// Primitive circle store (x,y,r)
	private double[] cx = new double[64];
	private double[] cy = new double[64];
	private double[] cr = new double[64];
	private int circleCount = 0;

	private final Coordinate tmp = new Coordinate(); // scratch

	/**
	 * Creates an instance constrained only by a polygonal boundary.
	 *
	 * @param boundary  polygonal constraint region (shell and holes are respected)
	 * @param tolerance accuracy tolerance (> 0)
	 * @throws IllegalArgumentException if {@code boundary} is null/empty,
	 *                                  non-polygonal, or if {@code tolerance <= 0}
	 */
	public LargestEmptyCircles(Geometry boundary, double tolerance) {
		this(boundary, null, tolerance);
	}

	/**
	 * Creates an instance constrained by a polygonal boundary and optional
	 * obstacles.
	 *
	 * @param boundary  polygonal constraint region (shell and holes are respected)
	 * @param obstacles optional constraints geometry (may be any {@link Geometry}):
	 *                  <ul>
	 *                  <li>polygonal components are excluded regions (additional
	 *                  holes)</li>
	 *                  <li>lineal components contribute distance constraints</li>
	 *                  <li>puntal components contribute distance constraints</li>
	 *                  </ul>
	 *                  May be {@code null} or empty.
	 * @param tolerance accuracy tolerance (> 0). Smaller values increase work and
	 *                  accuracy.
	 * @throws IllegalArgumentException if {@code boundary} is null/empty,
	 *                                  non-polygonal, or if {@code tolerance <= 0}
	 */
	public LargestEmptyCircles(Geometry boundary, Geometry obstacles, double tolerance) {
		if (boundary == null || boundary.isEmpty()) {
			throw new IllegalArgumentException("Boundary geometry is null or empty.");
		}
		if (!(boundary instanceof Polygonal)) {
			throw new IllegalArgumentException("Boundary must be polygonal.");
		}
		if (tolerance <= 0) {
			throw new IllegalArgumentException("Accuracy tolerance is non-positive: " + tolerance);
		}
		this.boundary = boundary;
		this.obstacles = obstacles;
		this.tolerance = tolerance;
	}

	private void initBoundary() {
		gridEnv = boundary.getEnvelopeInternal();

		// Index distance to boundary rings AND obstacle linework.
		// Sign is determined by boundary (modified by polygonal obstacles).
		boundaryDistance = new PointDistanceIndex(boundary, obstacles);

		createInitialGrid(gridEnv, cellStack);
	}

	/**
	 * Computes signed distance to the constraint set at the given coordinate.
	 * <p>
	 * The returned value is:
	 * <ul>
	 * <li><b>positive</b> if the point lies inside the feasible region (inside
	 * {@code boundary} and outside any polygonal obstacles),</li>
	 * <li><b>negative</b> if the point lies outside the feasible region,</li>
	 * </ul>
	 * and the magnitude is the distance to the nearest constraining feature, which
	 * includes: boundary rings, obstacle lineal components, and obstacle puntal
	 * components.
	 *
	 * @param x x-ordinate
	 * @param y y-ordinate
	 * @return signed distance to constraints
	 */
	private double distanceToConstraints(double x, double y) {
		tmp.x = x;
		tmp.y = y;
		return boundaryDistance.distance(tmp); // signed distance
	}

	/**
	 * Computes the next {@code n} largest empty circles by repeatedly calling
	 * {@link #findNextLEC()}.
	 *
	 * @param n number of circles to compute
	 * @return array of circles as {@code [x, y, r]} triples (length {@code n})
	 */
	public double[][] findLECs(int n) {
		double[][] out = new double[n][3];
		for (int i = 0; i < n; i++) {
			out[i] = findNextLEC();
		}
		return out;
	}

	/**
	 * Finds the next largest empty circle and caches it internally.
	 * <p>
	 * On the first call, initialises the search structure. On subsequent calls,
	 * reuses candidate cells and updates them against the most recently found
	 * circle to avoid recomputing from scratch.
	 *
	 * @return the next circle as {@code [x, y, r]} where {@code (x,y)} is the
	 *         center and {@code r} is the radius (signed distance at the selected
	 *         center; typically {@code r > 0})
	 */
	public double[] findNextLEC() {
		double farthestD;

		if (gridEnv == null) { // first iteration
			initBoundary();

			farthestCell = createCentroidCell(boundary);
			farthestD = farthestCell.getDistance();
			for (Cell c : cellStack) {
				double d = c.getDistance();
				if (d > farthestD) {
					farthestD = d;
					farthestCell = c;
				}
			}
		} else {
			// update remaining candidates with newest circle only
			final double lastX = cx[circleCount - 1];
			final double lastY = cy[circleCount - 1];
			final double lastR = cr[circleCount - 1];

			for (Cell nextIterCell : nextIterCells) {
				nextIterCell.updateDistance(lastX, lastY, lastR);
			}

			cellStack.clear();
			cellStack.addAll(nextIterCells);
			nextIterCells.clear();

			farthestD = Double.NEGATIVE_INFINITY;
			for (Cell c : cellStack) {
				double d = c.getDistance();
				if (d > farthestD) {
					farthestD = d;
					farthestCell = c;
				}
			}
		}

		// Branch-and-bound
		while (!cellStack.isEmpty()) {
			Cell cell = cellStack.removeLast(); // DFS-like

			double d = cell.getDistance();
			if (d > farthestD) {
				farthestD = d;
				farthestCell = cell;
			}

			if (cell.isFullyOutside()) {
				continue;
			}

			if (cell.isOutside()) {
				if (cell.getMaxDistance() > tolerance) {
					enqueueChildren(cell, farthestD);
				}
			} else {
				if (cell.getMaxDistance() - farthestD > tolerance) {
					enqueueChildren(cell, farthestD);
				} else {
					nextIterCells.add(cell);
				}
			}
		}

		double x = farthestCell.getX();
		double y = farthestCell.getY();
		double r = farthestCell.getDistance();

		addCircle(x, y, r);
		return new double[] { x, y, r };
	}

	private void addCircle(double x, double y, double r) {
		if (circleCount == cx.length) {
			int n = cx.length << 1;
			cx = Arrays.copyOf(cx, n);
			cy = Arrays.copyOf(cy, n);
			cr = Arrays.copyOf(cr, n);
		}
		cx[circleCount] = x;
		cy[circleCount] = y;
		cr[circleCount] = r;
		circleCount++;
	}

	private void enqueueChildren(final Cell cell, final double farthestD) {
		final double h2 = cell.getHSide() / 2.0;

		// optimistic bound for any child of this cell
		final double maxChildPotential = cell.getDistance() + 2.0 * h2 * Cell.SQRT2;
		if (maxChildPotential <= farthestD + tolerance) {
			nextIterCells.add(cell);
			return;
		}

		// Create 4 kids, push all (no sorting)
		Cell c1 = createCellIfUseful(cell.x - h2, cell.y - h2, h2, farthestD);
		Cell c2 = createCellIfUseful(cell.x + h2, cell.y - h2, h2, farthestD);
		Cell c3 = createCellIfUseful(cell.x - h2, cell.y + h2, h2, farthestD);
		Cell c4 = createCellIfUseful(cell.x + h2, cell.y + h2, h2, farthestD);

		if (c1 != null) {
			cellStack.addLast(c1);
		}
		if (c2 != null) {
			cellStack.addLast(c2);
		}
		if (c3 != null) {
			cellStack.addLast(c3);
		}
		if (c4 != null) {
			cellStack.addLast(c4);
		}
	}

	private Cell createCellIfUseful(final double x, final double y, final double h, final double farthestD) {
		Cell c = createCell(x, y, h);

		if (c.getMaxDistance() > farthestD + tolerance) {
			return c;
		}

		if (!c.isFullyOutside()) {
			nextIterCells.add(c);
		}
		return null;
	}

	private void createInitialGrid(Envelope env, Collection<Cell> target) {
		double minX = env.getMinX(), maxX = env.getMaxX();
		double minY = env.getMinY(), maxY = env.getMaxY();
		double cellSize = Math.min(env.getWidth(), env.getHeight());
		double hSize = cellSize / 2.0;

		for (double x = minX; x < maxX; x += cellSize) {
			for (double y = minY; y < maxY; y += cellSize) {
				target.add(createCell(x + hSize, y + hSize, hSize));
			}
		}
	}

	private Cell createCell(final double x, final double y, final double h) {
		Cell c = new Cell(x, y, h, distanceToConstraints(x, y));
		c.updateDistanceAll(cx, cy, cr, circleCount);
		return c;
	}

	private Cell createCentroidCell(Geometry geom) {
		Point p = geom.getCentroid();
		Cell c = new Cell(p.getX(), p.getY(), 0, distanceToConstraints(p.getX(), p.getY()));
		c.updateDistanceAll(cx, cy, cr, circleCount);
		return c;
	}

	private static final class Cell {
		static final double SQRT2 = Math.sqrt(2);

		private final double x, y, hSide;
		private double distance; // signed
		private double maxDist;

		Cell(double x, double y, double hSide, double dist) {
			this.x = x;
			this.y = y;
			this.hSide = hSide;
			this.distance = dist;
			this.maxDist = dist + hSide * SQRT2;
		}

		void updateDistance(double cX, double cY, double cR) {
			final double dx = x - cX;
			final double dy = y - cY;
			final double dsq = dx * dx + dy * dy;

			double D = distance;
			double t = D + cR;
			if (t > 0) {
				double tsq = t * t;
				if (dsq < tsq) {
					double d = Math.sqrt(dsq) - cR;
					if (d < D) {
						distance = d;
						maxDist = d + hSide * SQRT2;
					}
				}
			}
		}

		void updateDistanceAll(double[] cx, double[] cy, double[] cr, int n) {
			double D = distance;
			for (int i = 0; i < n; i++) {
				final double r = cr[i];
				double t = D + r;
				if (t <= 0) {
					continue;
				}

				final double dx = x - cx[i];
				final double dy = y - cy[i];
				final double dsq = dx * dx + dy * dy;

				final double tsq = t * t;
				if (dsq < tsq) {
					double d = Math.sqrt(dsq) - r;
					if (d < D) {
						D = d;
					}
				}
			}
			if (D < distance) {
				distance = D;
				maxDist = D + hSide * SQRT2;
			}
		}

		boolean isFullyOutside() {
			return maxDist < 0;
		}

		boolean isOutside() {
			return distance < 0;
		}

		double getMaxDistance() {
			return maxDist;
		}

		double getDistance() {
			return distance;
		}

		double getHSide() {
			return hSide;
		}

		double getX() {
			return x;
		}

		double getY() {
			return y;
		}
	}
}