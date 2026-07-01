package micycle.pgs.commons;

import java.util.Arrays;
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
 * <h2>Algorithm</h2> Best-first branch-and-bound over a quadtree of cells,
 * ordered by a <b>lazy max-heap</b> keyed on each cell's optimistic upper
 * bound ({@code dist + halfSide·√2}):
 * <ul>
 * <li><b>Best-first</b>: the cell with the greatest upper bound is always
 * subdivided next, so the search only ever expands cells whose bound exceeds
 * the final answer — the provably minimal set for this bound function — and
 * terminates the moment {@code top.maxDist − incumbent ≤ tolerance}, with no
 * queue draining.</li>
 * <li><b>Persistent</b>: the heap survives across {@link #findNextLEC()}
 * calls. Cells left in the heap (bounds below the previous answer) are the
 * candidate set for the next extraction.</li>
 * <li><b>Lazy</b>: each cell carries a <em>stamp</em> — the number of circles
 * that had been found when its distance was last computed. Since new circles
 * only ever <em>decrease</em> a cell's distance (keys only shrink), a stale
 * cell's stored key is an upper bound on its true key, so the heap top always
 * dominates every true key. Cells are refreshed against newer circles only
 * when they surface at the top; the (typical) majority that never surface are
 * never touched.</li>
 * </ul>
 * Cells are stored in flat parallel primitive arrays (no per-cell objects, no
 * GC pressure, cache-friendly sifting).
 */
public class LargestEmptyCircles {

	private static final double SQRT2 = Math.sqrt(2);

	private final Geometry boundary; // polygonal
	private final Geometry obstacles; // nullable; may be any Geometry
	private final double tolerance;

	private PointDistanceIndex boundaryDistance;
	private boolean initialized = false;

	/*
	 * Lazy max-heap of cells, struct-of-arrays. Keyed on hmax (optimistic upper
	 * bound = hd + hh·√2). hstamp[i] is the circle count when hd[i] was last
	 * computed; hd only decreases as circles are added, so stale keys
	 * over-estimate — safe for a max-heap with refresh-at-top.
	 */
	private double[] hx = new double[4096];
	private double[] hy = new double[4096];
	private double[] hh = new double[4096]; // half side length
	private double[] hd = new double[4096]; // signed distance at center
	private double[] hmax = new double[4096]; // hd + hh·√2 (heap key)
	private int[] hstamp = new int[4096];
	private int heapSize = 0;

	// Incumbent (best evaluated center) of the current extraction.
	private double bestX, bestY, bestD;

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
		this(boundary, obstacles, null, tolerance);
	}

	/**
	 * Creates an instance seeded with an existing circle packing, so the search
	 * "fills in" the gaps of that packing rather than starting from an empty
	 * region.
	 * <p>
	 * Each seed circle behaves exactly like a circle previously returned by
	 * {@link #findNextLEC()}: subsequent circles will not overlap it (they treat
	 * its rim as a distance constraint). Seed circles are <em>not</em> returned
	 * by {@link #findNextLEC()} / {@link #findLECs(int)}; only newly-found
	 * circles are.
	 *
	 * @param boundary        polygonal constraint region (shell and holes are
	 *                        respected)
	 * @param obstacles       optional constraints geometry (see
	 *                        {@link #LargestEmptyCircles(Geometry, Geometry, double)});
	 *                        may be {@code null} or empty
	 * @param existingCircles optional list of pre-existing circles, one per
	 *                        {@link Coordinate}, where {@code x,y} is the circle
	 *                        center and {@code z} is its radius. A {@code NaN}
	 *                        radius is treated as {@code 0} (a point constraint).
	 *                        May be {@code null} or empty.
	 * @param tolerance       accuracy tolerance (> 0). Smaller values increase
	 *                        work and accuracy.
	 * @throws IllegalArgumentException if {@code boundary} is null/empty,
	 *                                  non-polygonal, if {@code tolerance <= 0},
	 *                                  or if a seed circle has a non-finite
	 *                                  center or negative radius
	 */
	public LargestEmptyCircles(Geometry boundary, Geometry obstacles, List<Coordinate> existingCircles, double tolerance) {
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

		if (existingCircles != null) {
			for (Coordinate c : existingCircles) {
				if (c == null) {
					continue;
				}
				final double r = Double.isNaN(c.getZ()) ? 0 : c.getZ();
				if (!Double.isFinite(c.x) || !Double.isFinite(c.y) || !Double.isFinite(r) || r < 0) {
					throw new IllegalArgumentException("Invalid seed circle: " + c);
				}
				addCircle(c.x, c.y, r);
			}
		}
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
	 * the persistent heap of candidate cells is reused; cells are lazily updated
	 * against circles found since they were last evaluated, only when they reach
	 * the top of the heap.
	 *
	 * @return the next circle as {@code [x, y, r]} where {@code (x,y)} is the
	 *         center and {@code r} is the radius (signed distance at the selected
	 *         center; typically {@code r > 0})
	 */
	public double[] findNextLEC() {
		if (!initialized) {
			init();
		} else {
			bestD = Double.NEGATIVE_INFINITY;
		}

		final double tol = tolerance;

		while (true) {
			refreshTop(); // ensure heap top (if any) is current w.r.t. all circles
			if (heapSize == 0) {
				break;
			}

			// top is fresh: its key is the global maximum of all true upper bounds
			final double d = hd[0];
			if (d > bestD) {
				bestD = d;
				bestX = hx[0];
				bestY = hy[0];
			}

			// Optimality gap: no cell anywhere can beat the incumbent by > tol.
			if (hmax[0] - bestD <= tol) {
				break;
			}

			// Pop the top and subdivide it into 4 children.
			final double x = hx[0];
			final double y = hy[0];
			final double h2 = hh[0] * 0.5;
			popTop();

			final double reach = h2 * SQRT2;
			evalChild(x - h2, y - h2, h2, reach);
			evalChild(x + h2, y - h2, h2, reach);
			evalChild(x - h2, y + h2, h2, reach);
			evalChild(x + h2, y + h2, h2, reach);
		}

		final double x = bestX, y = bestY, r = bestD;
		addCircle(x, y, r);
		return new double[] { x, y, r };
	}

	private void init() {
		final Envelope env = boundary.getEnvelopeInternal();

		// Index distance to boundary rings AND obstacle linework.
		// Sign is determined by boundary (modified by polygonal obstacles).
		boundaryDistance = new PointDistanceIndex(boundary, obstacles);

		// Seed the incumbent with the centroid (clamped by any seed circles).
		final Point p = boundary.getCentroid();
		bestX = p.getX();
		bestY = p.getY();
		bestD = clampedDistance(bestX, bestY);

		// Initial square grid over the envelope.
		final double minX = env.getMinX(), maxX = env.getMaxX();
		final double minY = env.getMinY(), maxY = env.getMaxY();
		final double cellSize = Math.min(env.getWidth(), env.getHeight());
		final double hSize = cellSize / 2.0;
		final double reach = hSize * SQRT2;

		for (double gx = minX; gx < maxX; gx += cellSize) {
			for (double gy = minY; gy < maxY; gy += cellSize) {
				evalChild(gx + hSize, gy + hSize, hSize, reach);
			}
		}
		initialized = true;
	}

	/**
	 * Evaluates a child cell against all constraints and circles, updates the
	 * incumbent, and pushes it onto the heap unless it is entirely outside the
	 * feasible region ({@code maxDist < 0}).
	 */
	private void evalChild(final double px, final double py, final double h, final double reach) {
		final double d = clampedDistance(px, py);
		if (d > bestD) {
			bestD = d;
			bestX = px;
			bestY = py;
		}
		final double max = d + reach;
		if (max >= 0) {
			push(px, py, h, d, max, circleCount);
		}
	}

	/** Raw signed distance to boundary/obstacle constraints at {@code (x,y)}. */
	private double signedDistance(final double x, final double y) {
		tmp.x = x;
		tmp.y = y;
		return boundaryDistance.distance(tmp);
	}

	/**
	 * Signed distance to the constraint set at {@code (x,y)}, clamped by all
	 * circles found so far (distance to a previous circle's rim caps the value so
	 * new circles remain empty of old ones).
	 */
	private double clampedDistance(final double x, final double y) {
		double D = signedDistance(x, y);
		final double[] pcx = cx, pcy = cy, pcr = cr;
		final int n = circleCount;
		for (int i = 0; i < n; i++) {
			final double r = pcr[i];
			final double t = D + r;
			if (t <= 0) {
				continue; // circle i cannot reduce D
			}
			final double dx = x - pcx[i];
			final double dy = y - pcy[i];
			final double dsq = dx * dx + dy * dy;
			if (dsq < t * t) { // sqrt(dsq) - r < D, so it improves
				final double d = Math.sqrt(dsq) - r;
				if (d < D) {
					D = d;
				}
			}
		}
		return D;
	}

	// ------------------------------------------------------------------
	// Lazy max-heap (keyed on hmax)
	// ------------------------------------------------------------------

	/**
	 * Ensures the heap top is current with respect to all circles found so far,
	 * discarding cells that become fully infeasible. Loops because a sift-down
	 * after refreshing may surface another stale cell.
	 * <p>
	 * Correctness of laziness: distances only decrease as circles are added, so a
	 * stale key over-estimates its true value; the heap top's stored key
	 * therefore dominates every true key in the heap, and once the top is fresh
	 * its key is the true global maximum.
	 */
	private void refreshTop() {
		final int cc = circleCount;
		while (heapSize > 0 && hstamp[0] != cc) {
			double D = hd[0];
			final double x = hx[0], y = hy[0];
			// apply only the circles added since this cell was last touched
			for (int i = hstamp[0]; i < cc; i++) {
				final double r = cr[i];
				final double t = D + r;
				if (t <= 0) {
					continue;
				}
				final double dx = x - cx[i];
				final double dy = y - cy[i];
				final double dsq = dx * dx + dy * dy;
				if (dsq < t * t) {
					final double d = Math.sqrt(dsq) - r;
					if (d < D) {
						D = d;
					}
				}
			}
			hstamp[0] = cc;
			if (D < hd[0]) {
				hd[0] = D;
				final double max = D + hh[0] * SQRT2;
				if (max < 0) { // now fully outside; discard
					popTop();
					continue;
				}
				hmax[0] = max;
				siftDown(0); // key decreased
			}
			// if key unchanged, top is fresh and still the max → loop exits
		}
	}

	private void push(final double x, final double y, final double h, final double d, final double max,
			final int stamp) {
		if (heapSize == hx.length) {
			grow();
		}
		int i = heapSize++;
		// sift up with a hole (no swaps)
		while (i > 0) {
			final int parent = (i - 1) >>> 1;
			if (hmax[parent] >= max) {
				break;
			}
			copyCell(parent, i);
			i = parent;
		}
		hx[i] = x;
		hy[i] = y;
		hh[i] = h;
		hd[i] = d;
		hmax[i] = max;
		hstamp[i] = stamp;
	}

	private void popTop() {
		final int last = --heapSize;
		if (last > 0) {
			copyCell(last, 0);
			siftDown(0);
		}
	}

	private void siftDown(int i) {
		final int n = heapSize;
		final double x = hx[i], y = hy[i], h = hh[i], d = hd[i], max = hmax[i];
		final int stamp = hstamp[i];
		final int half = n >>> 1; // nodes >= half are leaves
		while (i < half) {
			int child = (i << 1) + 1;
			final int right = child + 1;
			if (right < n && hmax[right] > hmax[child]) {
				child = right;
			}
			if (hmax[child] <= max) {
				break;
			}
			copyCell(child, i);
			i = child;
		}
		hx[i] = x;
		hy[i] = y;
		hh[i] = h;
		hd[i] = d;
		hmax[i] = max;
		hstamp[i] = stamp;
	}

	private void copyCell(final int from, final int to) {
		hx[to] = hx[from];
		hy[to] = hy[from];
		hh[to] = hh[from];
		hd[to] = hd[from];
		hmax[to] = hmax[from];
		hstamp[to] = hstamp[from];
	}

	private void grow() {
		final int n = hx.length << 1;
		hx = Arrays.copyOf(hx, n);
		hy = Arrays.copyOf(hy, n);
		hh = Arrays.copyOf(hh, n);
		hd = Arrays.copyOf(hd, n);
		hmax = Arrays.copyOf(hmax, n);
		hstamp = Arrays.copyOf(hstamp, n);
	}

	private void addCircle(final double x, final double y, final double r) {
		if (circleCount == cx.length) {
			final int n = cx.length << 1;
			cx = Arrays.copyOf(cx, n);
			cy = Arrays.copyOf(cy, n);
			cr = Arrays.copyOf(cr, n);
		}
		cx[circleCount] = x;
		cy[circleCount] = y;
		cr[circleCount] = r;
		circleCount++;
	}
}