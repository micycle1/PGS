package micycle.pgs.commons;

import org.locationtech.jts.geom.*;
import java.util.*;

/**
 * Generates a random domino tiling of an <b>Aztec diamond</b> as JTS geometry.
 *
 * <h2>What is an Aztec diamond?</h2> The <i>Aztec diamond of order n</i>,
 * denoted A(n), is a well-known region on the square grid. Informally, it is
 * the “diamond-shaped” union of unit grid squares whose centers satisfy
 * {@code |x| + |y| < n} (with a suitable half-integer centering). Its boundary
 * looks like a rotated square with “stair-step” edges.
 *
 * <p>
 * The region A(n) contains exactly {@code 2*n*(n+1)} unit squares. A <b>domino
 * tiling</b> of A(n) is a perfect cover of that region by {@code n*(n+1)}
 * dominoes, where each domino is a 1×2 or 2×1 rectangle made of two adjacent
 * unit squares, with no overlaps and no gaps.
 * 
 * <h2>Why are Aztec diamonds studied?</h2> Aztec diamonds are a “toy model” in
 * combinatorics and statistical physics: they are simple to define, but their
 * random domino tilings exhibit striking large-scale patterns. In particular, a
 * uniformly random tiling of a large Aztec diamond typically develops a frozen
 * outer region and a disordered inner region separated by an increasingly sharp
 * boundary that approaches a circle (the “arctic circle theorem”). Because the
 * model is exactly solvable, Aztec diamonds are used to study randomness, phase
 * transitions, and limit shapes, and they connect to other structures such as
 * non-intersecting lattice paths, determinantal point processes, and
 * alternating sign matrices.
 *
 * <h2>What does this library produce?</h2> This class implements a standard
 * <b>domino shuffling</b> / <b>local growth</b> procedure that produces a
 * random domino tiling of A(n). The result is exported as a JTS
 * {@link MultiPolygon} where <b>each domino is represented as a rectangle
 * {@link Polygon}</b>.
 *
 * <h2>Coordinate system</h2> This implementation uses an integer cell grid with
 * unit cell size, and outputs rectangles in that coordinate space:
 * <ul>
 * <li>X corresponds to column index; Y corresponds to row index.</li>
 * <li>Y increases downward (screen-like). For a typical Cartesian Y-up system,
 * negate Y during export or apply an affine transform.</li>
 * </ul>
 *
 * <h2>Randomness and reproducibility</h2> Random choices occur during the “fill
 * 2×2 blocks” step. Providing a seeded {@link Random} yields reproducible
 * tilings.
 * 
 * @author Michael Carleton
 */
public final class AztecDiamond {

	/**
	 * Domino "orientation" as used by the shuffling dynamics.
	 *
	 * <p>
	 * Important: this is not simply “north/south means vertical”. In this
	 * algorithm, the orientation determines how a domino moves during the “move
	 * tiles” step, and it also determines which adjacent grid cell is paired with
	 * its upper-left cell.
	 * </p>
	 */
	public enum Orientation {
		N, S, E, W
	}

	private static int dr(Orientation o) {
		return switch (o) {
			case N -> -1;
			case S -> 1;
			case E, W -> 0;
		};
	}

	private static int dc(Orientation o) {
		return switch (o) {
			case E -> 1;
			case W -> -1;
			case N, S -> 0;
		};
	}

	private static Orientation conflict(Orientation o) {
		return switch (o) {
			case N -> Orientation.S;
			case S -> Orientation.N;
			case E -> Orientation.W;
			case W -> Orientation.E;
		};
	}

	private static final class Domino {
		int r; // upper-left cell row, in centered coordinates
		int c; // upper-left cell col, in centered coordinates
		final Orientation o;

		Domino(int r, int c, Orientation o) {
			this.r = r;
			this.c = c;
			this.o = Objects.requireNonNull(o);
		}

		void step() {
			r += dr(o);
			c += dc(o);
		}

		/**
		 * @return true if this domino occupies two cells stacked vertically (1×2);
		 *         false if it occupies two cells horizontally (2×1).
		 */
		boolean isVerticalShape() {
			// Mirrors the Python implementation:
			// E/W occupy (r,c) and (r+1,c); N/S occupy (r,c) and (r,c+1)
			return o == Orientation.E || o == Orientation.W;
		}
	}

	private final GeometryFactory gf;
	private final Random rnd;

	private int order; // current order during growth
	private boolean[][] mask; // inside-diamond indicator for unit cells
	private Domino[][] grid; // per-cell reference to the occupying domino
	private final ArrayList<Domino> tiles = new ArrayList<>();

	/**
	 * Constructs a generator and immediately produces a random tiling of the Aztec
	 * diamond A(order).
	 *
	 * <p>
	 * The generation starts from A(1) and repeatedly applies one domino-shuffling
	 * "growth step" until the requested order is reached.
	 * </p>
	 *
	 * @param order the target order n for A(n); must be {@code > 0}.
	 * @param gf    JTS geometry factory used to create polygons; must not be null.
	 * @param rnd   source of randomness used by the algorithm; must not be null.
	 *              Supply a seeded instance for reproducible output.
	 * @throws IllegalArgumentException if {@code order <= 0}.
	 * @throws NullPointerException     if {@code gf} or {@code rnd} is null.
	 */
	public AztecDiamond(int order, GeometryFactory gf, Random rnd) {
		if (order <= 0)
			throw new IllegalArgumentException("order must be > 0");
		this.gf = Objects.requireNonNull(gf, "gf");
		this.rnd = Objects.requireNonNull(rnd, "rnd");

		// Start at order 1, then grow to requested order.
		this.order = 1;
		rebuildMaskAndGrid();
		fillTwoByTwos(); // fill A(1)
		while (this.order < order) {
			stepTileGeneration();
		}
	}

	/**
	 * Convenience factory that generates a random tiling of A(order) using a
	 * default {@link GeometryFactory} and a deterministic seed.
	 *
	 * <p>
	 * This is useful for one-liners and tests.
	 * </p>
	 *
	 * @param order the target order n for A(n); must be {@code > 0}.
	 * @param seed  seed for the pseudo-random number generator.
	 * @return a {@link MultiPolygon} where each component polygon is one domino
	 *         rectangle.
	 */
	public static MultiPolygon generate(int order, long seed) {
		GeometryFactory gf = new GeometryFactory();
		AztecDiamond gen = new AztecDiamond(order, gf, new Random(seed));
		return gen.toMultiPolygon(1.0, 0.0, 0.0);
	}

	/**
	 * Exports the current tiling as a JTS {@link MultiPolygon}.
	 *
	 * <p>
	 * The output contains one rectangle {@link Polygon} per domino. Each rectangle
	 * is axis-aligned in the internal grid coordinate system.
	 * </p>
	 *
	 * <h3>Per-domino “color/class” stored in {@code userData}</h3>
	 * <p>
	 * Each returned domino polygon has an {@code Integer} stored in
	 * {@link Geometry#getUserData()} (set via
	 * {@link Geometry#setUserData(Object)}). This integer encodes the domino’s
	 * <b>type</b> as used by the underlying Aztec-diamond tiling algorithm.
	 * </p>
	 *
	 * <p>
	 * In plain terms: the algorithm imagines that every domino belongs to one of
	 * four “classes”. The class says <b>which way that domino would like to
	 * slide</b> during the shuffling step (up, down, left, or right), and it also
	 * determines which neighboring domino it can cancel with. This is the same kind
	 * of 4-way labeling that is often visualized by coloring dominoes in different
	 * colors.
	 * </p>
	 *
	 * <p>
	 * Encoding (matching the original Python reference):
	 * </p>
	 * <ul>
	 * <li>{@code 0} = N (would move “north” / up)</li>
	 * <li>{@code 1} = S (would move “south” / down)</li>
	 * <li>{@code 2} = E (would move “east” / right)</li>
	 * <li>{@code 3} = W (would move “west” / left)</li>
	 * </ul>
	 *
	 * <p>
	 * Note: this “class” is a property of the shuffling dynamics. It is closely
	 * related to the common “checkerboard coloring” classification (vertical vs
	 * horizontal, and which checkerboard color is on the top/left half of the
	 * domino), but the exact correspondence depends on the coordinate and
	 * checkerboard conventions used.
	 * </p>
	 *
	 * @param cellSize scale factor applied to each unit grid cell (e.g., 10.0 makes
	 *                 each cell 10×10).
	 * @param originX  translation applied to exported X coordinates.
	 * @param originY  translation applied to exported Y coordinates.
	 * @return a {@link MultiPolygon} consisting of {@code n*(n+1)} domino polygons
	 *         for order n.
	 */
	public MultiPolygon toMultiPolygon(double cellSize, double originX, double originY) {
		Polygon[] polys = new Polygon[tiles.size()];
		int i = 0;
		for (Domino d : tiles) {
			polys[i++] = dominoPolygon(d, cellSize, originX, originY);
		}
		return gf.createMultiPolygon(polys);
	}

	/**
	 * Performs one full domino-shuffling growth step:
	 * <ol>
	 * <li>increase the order (A(k) → A(k+1))</li>
	 * <li>remove ("annihilate") adjacent dominoes that would move into each
	 * other</li>
	 * <li>move remaining dominoes one unit in their orientation direction</li>
	 * <li>fill all newly created 2×2 holes randomly with two dominoes</li>
	 * </ol>
	 *
	 * <p>
	 * This method is internal; the constructor repeatedly calls it until the target
	 * order is reached.
	 * </p>
	 */
	private void stepTileGeneration() {
		increaseOrder();
		cancelOpposingMovers();
		moveTiles();
		fillTwoByTwos();
	}

	/**
	 * Rebuilds the Aztec diamond mask and empties the occupancy grid for the
	 * current {@link #order}.
	 *
	 * <p>
	 * The mask marks which unit cells belong to A(order). We use the standard
	 * characterization in terms of cell centers around the point
	 * ({@code order-0.5}, {@code order-0.5}).
	 * </p>
	 */
	private void rebuildMaskAndGrid() {
		int size = 2 * order;
		mask = new boolean[size][size];
		grid = new Domino[size][size];

		// Aztec diamond membership test on unit cells using cell centers.
		double center = order - 0.5;
		for (int r = 0; r < size; r++) {
			for (int c = 0; c < size; c++) {
				double x = r - center;
				double y = c - center;
				mask[r][c] = (Math.abs(x) + Math.abs(y) < order);
			}
		}
	}

	/**
	 * Increases the order by one and embeds the previous occupancy grid into the
	 * new one.
	 *
	 * <p>
	 * This corresponds to padding the old 2n×2n grid into the center of the new
	 * 2(n+1)×2(n+1) grid. It mirrors the Python operation
	 * {@code new[1:-1,1:-1] = old}.
	 * </p>
	 */
	private void increaseOrder() {
		int oldOrder = order;
		Domino[][] oldGrid = grid;

		order = oldOrder + 1;
		rebuildMaskAndGrid();

		int oldSize = 2 * oldOrder;
		for (int r = 0; r < oldSize; r++) {
			System.arraycopy(oldGrid[r], 0, grid[r + 1], 1, oldSize);
		}
	}

	/**
	 * Cancels pairs of dominoes that would move into each other (opposing movers).
	 *
	 * <p>
	 * In domino shuffling, dominoes carry a "movement direction". If a domino would
	 * move into a neighboring cell that is occupied by a domino moving in the exact
	 * opposite direction, the pair is removed. This creates holes that are later
	 * re-filled in 2×2 blocks.
	 * </p>
	 */
	private void cancelOpposingMovers() {
		HashSet<Domino> removed = new HashSet<>();

		int size = 2 * order;
		for (int r = 0; r < size; r++) {
			for (int c = 0; c < size; c++) {
				if (!mask[r][c])
					continue;
				Domino d = grid[r][c];
				if (d == null || removed.contains(d))
					continue;

				int r2 = r + dr(d.o);
				int c2 = c + dc(d.o);
				if (r2 < 0 || r2 >= size || c2 < 0 || c2 >= size)
					continue;

				Domino d2 = grid[r2][c2];
				if (d2 == null || removed.contains(d2))
					continue;

				if (d2.o == conflict(d.o)) {
					removed.add(d);
					removed.add(d2);
					clearDominoFromGrid(d);
					clearDominoFromGrid(d2);
				}
			}
		}

		if (!removed.isEmpty()) {
			tiles.removeIf(removed::contains);
		}
	}

	/**
	 * Removes references to a domino from the current occupancy grid.
	 *
	 * @param d the domino to clear.
	 */
	private void clearDominoFromGrid(Domino d) {
		int br = d.r + order;
		int bc = d.c + order;
		if (0 <= br && br < grid.length && 0 <= bc && bc < grid.length && grid[br][bc] == d) {
			grid[br][bc] = null;
		}

		if (d.isVerticalShape()) {
			int r2 = br + 1, c2 = bc;
			if (0 <= r2 && r2 < grid.length && grid[r2][c2] == d)
				grid[r2][c2] = null;
		} else {
			int r2 = br, c2 = bc + 1;
			if (0 <= c2 && c2 < grid.length && grid[r2][c2] == d)
				grid[r2][c2] = null;
		}
	}

	/**
	 * Moves every domino one unit in its movement direction and rebuilds the
	 * occupancy grid.
	 *
	 * <p>
	 * This is the "shuffling" part: after cancellations, remaining dominoes are
	 * shifted. The algorithm guarantees that after cancellation, these shifts can
	 * be applied without creating overlaps inside the diamond mask.
	 * </p>
	 */
	private void moveTiles() {
		Domino[][] newGrid = new Domino[2 * order][2 * order];

		for (Domino d : tiles) {
			d.step();
			putDominoInGrid(d, newGrid);
		}

		grid = newGrid;
	}

	/**
	 * Places a domino into an occupancy grid by marking its two unit cells.
	 *
	 * @param d domino to place.
	 * @param g target grid.
	 */
	private void putDominoInGrid(Domino d, Domino[][] g) {
		int br = d.r + order;
		int bc = d.c + order;
		g[br][bc] = d;

		if (d.isVerticalShape()) {
			g[br + 1][bc] = d;
		} else {
			g[br][bc + 1] = d;
		}
	}

	/**
	 * Fills all empty cells (holes) inside the Aztec diamond mask by repeatedly
	 * selecting an empty cell and filling the surrounding 2×2 block with two
	 * dominoes in one of two random patterns.
	 *
	 * <p>
	 * This is the only randomized step in the growth iteration. Each 2×2 hole is
	 * filled either as two vertical dominoes side-by-side or as two horizontal
	 * dominoes stacked.
	 * </p>
	 */
	private void fillTwoByTwos() {
		int size = 2 * order;

		while (true) {
			int rr = -1, cc = -1;

			// Find any empty cell inside the mask.
			outer: for (int r = 0; r < size; r++) {
				for (int c = 0; c < size; c++) {
					if (mask[r][c] && grid[r][c] == null) {
						rr = r;
						cc = c;
						break outer;
					}
				}
			}
			if (rr < 0)
				break;

			// The shuffling algorithm ensures holes come in 2x2 blocks. Keep a safety
			// check.
			if (rr + 1 >= size || cc + 1 >= size) {
				throw new IllegalStateException("Unexpected hole too close to boundary at (" + rr + "," + cc + ")");
			}

			if (rnd.nextBoolean()) {
				// Two vertical-shape dominoes side-by-side
				Domino a = new Domino(rr - order, cc - order, Orientation.W);
				Domino b = new Domino(rr - order, (cc - order) + 1, Orientation.E);

				tiles.add(a);
				tiles.add(b);

				putDominoInGrid(a, grid);
				putDominoInGrid(b, grid);
			} else {
				// Two horizontal-shape dominoes stacked
				Domino a = new Domino(rr - order, cc - order, Orientation.N);
				Domino b = new Domino((rr - order) + 1, cc - order, Orientation.S);

				tiles.add(a);
				tiles.add(b);

				putDominoInGrid(a, grid);
				putDominoInGrid(b, grid);
			}
		}
	}

	private static int classCode(Orientation o) {
		// Match Python: ORIENTATIONS = N, S, E, W = range(4)
		return switch (o) {
			case N -> 0;
			case S -> 1;
			case E -> 2;
			case W -> 3;
		};
	}

	/**
	 * Converts one domino into a rectangle polygon.
	 *
	 * @param d        domino to export.
	 * @param cellSize scale factor for unit cells.
	 * @param originX  translation in X.
	 * @param originY  translation in Y.
	 * @return an axis-aligned rectangle polygon covering the domino area.
	 */
	private Polygon dominoPolygon(Domino d, double cellSize, double originX, double originY) {
		boolean vertical = d.isVerticalShape();
		int wCells = vertical ? 1 : 2;
		int hCells = vertical ? 2 : 1;

		double x1 = originX + d.c * cellSize;
		double y1 = originY + d.r * cellSize;
		double x2 = x1 + wCells * cellSize;
		double y2 = y1 + hCells * cellSize;

		Coordinate[] ring = new Coordinate[] { new Coordinate(x1, y1), new Coordinate(x2, y1), new Coordinate(x2, y2), new Coordinate(x1, y2),
				new Coordinate(x1, y1) };

		Polygon p = gf.createPolygon(ring);

		// Encode the “color/class” as an integer in userData:
		p.setUserData(classCode(d.o)); // Integer 0..3

		return p;
	}

	/**
	 * Example usage: prints WKT for a random tiling.
	 *
	 * @param args ignored
	 */
	public static void main(String[] args) {
		int n = 20;
		MultiPolygon mp = AztecDiamond.generate(n, 12345L);
		System.out.println("Dominoes: " + mp.getNumGeometries()); // should be n*(n+1)
		System.out.println(mp);
	}
}