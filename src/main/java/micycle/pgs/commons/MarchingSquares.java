package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.DoubleBinaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import processing.core.PConstants;
import processing.core.PShape;

/**
 * Fast isolines from a regular grid using Marching Squares + contour tracing.
 */
public final class MarchingSquares {

	private MarchingSquares() {
	}

	/**
	 * Builds isolines by sampling a regular grid over a rectangle and evaluating a
	 * user-provided height function z=f(x,y) at each sample.
	 *
	 * @param bounds               [xmin, ymin, xmax, ymax] of the sampling area
	 * @param sampleSpacing        grid spacing in pixels (smaller = more detail,
	 *                             slower)
	 * @param intervalValueSpacing contour interval spacing in "height units"
	 * @param fn                   (x,y) -> height
	 */
	public static Map<PShape, Float> isolines(double[] bounds, double sampleSpacing, double intervalValueSpacing, DoubleBinaryOperator fn) {
		return isolines(bounds, sampleSpacing, intervalValueSpacing, Double.NaN, Double.NaN, fn);
	}

	/**
	 * Builds isolines by sampling a regular grid over a rectangle and evaluating a
	 * user-provided height function z=f(x,y) at each sample.
	 *
	 * @param bounds               [xmin, ymin, xmax, ymax] of the sampling area
	 * @param sampleSpacing        grid spacing in pixels (smaller = more detail,
	 *                             slower)
	 * @param intervalValueSpacing contour interval spacing in "height units"
	 * @param isolineMin           minimum contour value (inclusive). Pass
	 *                             Double.NaN to auto-detect from data.
	 * @param isolineMax           maximum contour value (inclusive). Pass
	 *                             Double.NaN to auto-detect from data.
	 * @param fn                   (x,y) -> height
	 */
	public static Map<PShape, Float> isolines(double[] bounds, double sampleSpacing, double intervalValueSpacing, double isolineMin, double isolineMax,
			DoubleBinaryOperator fn) {
		if (sampleSpacing <= 0) {
			throw new IllegalArgumentException("sampleSpacing must be > 0");
		}
		if (intervalValueSpacing <= 0) {
			throw new IllegalArgumentException("intervalValueSpacing must be > 0");
		}
		if (!Double.isNaN(isolineMax) && !Double.isNaN(isolineMin) && isolineMax < isolineMin) {
			return Collections.emptyMap();
		}

		if (bounds.length < 4) {
			throw new IllegalArgumentException("bounds must be double[4] {xmin, ymin, xmax, ymax}");
		}
		double x = bounds[0];
		double y = bounds[1];
		double w = bounds[2] - x;
		double h = bounds[3] - y;

		// Regular grid counts
		final int nx = (int) Math.floor(w / sampleSpacing) + 1;
		final int ny = (int) Math.floor(h / sampleSpacing) + 1;

		final float x0 = (float) x;
		final float y0 = (float) y;
		final float dx = (float) sampleSpacing;
		final float dy = (float) sampleSpacing;

		final float[] z = new float[nx * ny];
		float minZ = Float.POSITIVE_INFINITY;
		float maxZ = Float.NEGATIVE_INFINITY;

		// Row-major fill: y outer, x inner => index = iy*nx + ix
		int idx = 0;
		for (int iy = 0; iy < ny; iy++) {
			final double yy = y + iy * sampleSpacing;
			for (int ix = 0; ix < nx; ix++, idx++) {
				final double xx = x + ix * sampleSpacing;
				float val = (float) fn.applyAsDouble(xx, yy);
				z[idx] = val;
				if (val < minZ) {
					minZ = val;
				}
				if (val > maxZ) {
					maxZ = val;
				}
			}
		}

		return isolinesRegularGrid(z, nx, ny, x0, y0, dx, dy, intervalValueSpacing, isolineMin, isolineMax, minZ, maxZ);
	}

	/**
	 * Computes isolines for a regular grid of data.
	 * 
	 * @param z                    row-major array of grid values
	 * @param nx                   number of points in x direction
	 * @param ny                   number of points in y direction
	 * @param bounds               [xmin, ymin, xmax, ymax] of the sampling area
	 * @param intervalValueSpacing vertical distance between contour levels
	 */
	public static Map<PShape, Float> isolines(float[] z, int nx, int ny, double[] bounds, double intervalValueSpacing) {
		return isolines(z, nx, ny, bounds, intervalValueSpacing, Double.NaN, Double.NaN);
	}

	/**
	 * Computes isolines for a regular grid of data.
	 * 
	 * @param z                    row-major array of grid values
	 * @param nx                   number of points in x direction
	 * @param ny                   number of points in y direction
	 * @param bounds               [xmin, ymin, xmax, ymax] of the sampling area
	 * @param intervalValueSpacing vertical distance between contour levels
	 * @param isolineMin           minimum value to contour. Pass Double.NaN to
	 *                             auto-detect from data.
	 * @param isolineMax           maximum value to contour. Pass Double.NaN to
	 *                             auto-detect from data.
	 */
	public static Map<PShape, Float> isolines(float[] z, int nx, int ny, double[] bounds, double intervalValueSpacing, double isolineMin, double isolineMax) {
		if (bounds.length < 4) {
			throw new IllegalArgumentException("bounds must be double[4] {xmin, ymin, xmax, ymax}");
		}
		double x = bounds[0];
		double y = bounds[1];
		double w = bounds[2] - x;
		double h = bounds[3] - y;

		final float x0 = (float) x;
		final float y0 = (float) y;
		final float dx = (float) (w / (nx - 1));
		final float dy = (float) (h / (ny - 1));

		// Compute min/max for auto-ranging and optimization
		float minZ = Float.POSITIVE_INFINITY, maxZ = Float.NEGATIVE_INFINITY;
		for (float val : z) {
			minZ = Math.min(minZ, val);
			maxZ = Math.max(maxZ, val);
		}

		return isolinesRegularGrid(z, nx, ny, x0, y0, dx, dy, intervalValueSpacing, isolineMin, isolineMax, minZ, maxZ);
	}

	/**
	 * Traces only the zero-contour (fn = 0) using marching squares. Useful for
	 * implicit curves like Voronoi edges where you want f(x,y)=0 only.
	 *
	 * @param bounds        [xmin, ymin, xmax, ymax] of the sampling area
	 * @param sampleSpacing grid spacing in pixels (smaller = more detail, slower)
	 * @param smoothing     unused
	 * @param fn            (x,y) -> value whose zero-set is traced
	 * @return map of PShape -> level (always 0f)
	 */
	public static Map<PShape, Float> isolineZero(double[] bounds, double sampleSpacing, int smoothing, DoubleBinaryOperator fn) {

		if (sampleSpacing <= 0) {
			throw new IllegalArgumentException("sampleSpacing must be > 0");
		}

		if (bounds.length < 4) {
			throw new IllegalArgumentException("bounds must be double[4] {xmin, ymin, xmax, ymax}");
		}
		double x = bounds[0];
		double y = bounds[1];
		double w = bounds[2] - x;
		double h = bounds[3] - y;

		// Regular grid counts (stable, no epsilon loops)
		final int nx = (int) Math.floor(w / sampleSpacing) + 1;
		final int ny = (int) Math.floor(h / sampleSpacing) + 1;

		if (nx < 2 || ny < 2) {
			return Collections.emptyMap();
		}

		final float x0 = (float) x;
		final float y0 = (float) y;
		final float dx = (float) sampleSpacing;
		final float dy = (float) sampleSpacing;

		final float[] z = new float[nx * ny];

		// Row-major fill: y outer, x inner => index = iy*nx + ix
		int idx = 0;
		for (int iy = 0; iy < ny; iy++) {
			final double yy = y + iy * sampleSpacing;
			for (int ix = 0; ix < nx; ix++, idx++) {
				final double xx = x + ix * sampleSpacing;
				z[idx] = (float) fn.applyAsDouble(xx, yy);
			}
		}

		final int cellsX = nx - 1, cellsY = ny - 1;
		final int cellCount = cellsX * cellsY;

		// Single level: 0
		return processLevel(0f, z, nx, ny, cellsX, cellsY, cellCount, x0, y0, dx, dy);
	}

	/**
	 * Computes isolines for a regular grid of data.
	 *
	 * @param z                    row-major array of grid values
	 * @param nx                   number of points in x direction
	 * @param ny                   number of points in y direction
	 * @param x0                   minimum x coordinate
	 * @param y0                   minimum y coordinate
	 * @param dx                   sampling spacing in x
	 * @param dy                   sampling spacing in y
	 * @param intervalValueSpacing vertical distance between contour levels
	 * @param isolineMin           minimum value to contour (can be NaN)
	 * @param isolineMax           maximum value to contour (can be NaN)
	 * @param dataMin              min value present in z
	 * @param dataMax              max value present in z
	 * @return a map of PShapes to their corresponding levels
	 */
	private static Map<PShape, Float> isolinesRegularGrid(float[] z, int nx, int ny, float x0, float y0, float dx, float dy, double intervalValueSpacing,
			double isolineMin, double isolineMax, float dataMin, float dataMax) {

		final int cellsX = nx - 1, cellsY = ny - 1;
		final int cellCount = cellsX * cellsY;

		double start = Double.isNaN(isolineMin) ? dataMin : isolineMin;
		double end = Double.isNaN(isolineMax) ? dataMax : isolineMax;

		final int levelCount = (int) Math.floor((end - start) / intervalValueSpacing) + 1;

		if (levelCount <= 0) {
			return Collections.emptyMap();
		}

		// Use pre-computed min/max for early exit optimization
		final float finalMinZ = dataMin;
		final float finalMaxZ = dataMax;

		// @formatter:off
	    return IntStream.range(0, levelCount)
	        .parallel()
	        .mapToObj(li -> {
	            final float level = (float) (start + li * intervalValueSpacing);
	            // Early exit if level out of bounds
	            if (level < finalMinZ || level > finalMaxZ) {
	                return Collections.<PShape, Float>emptyMap();
	            }
	            return processLevel(level, z, nx, ny, cellsX, cellsY, cellCount, x0, y0, dx, dy);
	        })
	        .flatMap(map -> map.entrySet().stream())
	        .collect(Collectors.toMap(
	            Map.Entry::getKey, 
	            Map.Entry::getValue,
	            (v1, v2) -> v1, // merge function (shouldn't happen, but required)
	            () -> new LinkedHashMap<>(Math.max(16, levelCount * 8))
	        ));
	    // @formatter:on
	}

	/**
	 * Process a single contour level and return all shapes at that level.
	 */
	private static Map<PShape, Float> processLevel(float level, float[] z, int nx, int ny, int cellsX, int cellsY, int cellCount, float x0, float y0, float dx,
			float dy) {

		// Thread-local arrays (each parallel stream gets its own)
		final byte[] codes = new byte[cellCount];
		final byte[] amb = new byte[cellCount];
		final byte[] visited = new byte[cellCount];

		// Build marching squares codes
		buildMarchingSquaresCodes(codes, amb, z, nx, cellsX, cellsY, level);

		// Trace all contours at this level
		return traceContoursAtLevel(level, codes, amb, visited, z, nx, ny, cellsX, cellsY, cellCount, x0, y0, dx, dy);
	}

	/**
	 * Build marching squares codes for all cells at given level.
	 */
	private static void buildMarchingSquaresCodes(byte[] codes, byte[] amb, float[] z, int nx, int cellsX, int cellsY, float level) {
		int c = 0;
		for (int iy = 0; iy < cellsY; iy++) {
			final int row0 = iy * nx;
			final int row1 = (iy + 1) * nx;

			for (int ix = 0; ix < cellsX; ix++, c++) {
				final float v0 = z[row0 + ix];
				final float v1 = z[row0 + ix + 1];
				final float v2 = z[row1 + ix + 1];
				final float v3 = z[row1 + ix];

				int code = ((v0 > level ? 1 : 0) << 0) | ((v1 > level ? 1 : 0) << 1) | ((v2 > level ? 1 : 0) << 2) | ((v3 > level ? 1 : 0) << 3);
				codes[c] = (byte) code;

				// Asymptotic decider (Nielson & Hamann) better?
				if (code == 5 || code == 10) {
					float center = 0.25f * (v0 + v1 + v2 + v3);
					amb[c] = (byte) (center > level ? 1 : 0);
				} else {
					amb[c] = 0;
				}
			}
		}
	}

	/**
	 * Trace all contours at a given level.
	 */
	private static Map<PShape, Float> traceContoursAtLevel(float level, byte[] codes, byte[] amb, byte[] visited, float[] z, int nx, int ny, int cellsX,
			int cellsY, int cellCount, float x0, float y0, float dx, float dy) {

		// Collect raw polylines first
		List<FloatPath> paths = new ArrayList<>();

		for (int cell = 0; cell < cellCount; cell++) {
			int code = codes[cell] & 0xFF;
			if (code == 0 || code == 15)
				continue;

			int edgesMask = edgesUsedMask(code);
			int unvisited = edgesMask & (~visited[cell] & 0x0F);

			while (unvisited != 0) {
				int startEdge = Integer.numberOfTrailingZeros(unvisited);

				FloatPath path = traceOne(level, cell, startEdge, codes, amb, visited, z, nx, ny, x0, y0, dx, dy);
				if (path != null && path.sizePairs() >= 2) {
					paths.add(path);
				}

				unvisited = edgesMask & (~visited[cell] & 0x0F);
			}
		}

		stitchPathsDirect(paths);

		Map<PShape, Float> result = new LinkedHashMap<>(paths.size() * 2);
		for (FloatPath p : paths) {
			p.snapClosed(1e-3f);

			PShape s = toPShape(p);
			if (s.getVertexCount() > 1) {
				result.put(s, level);
			}
		}
		return result;
	}

	private static void stitchPathsDirect(List<FloatPath> paths) {
		boolean merged;
		do {
			merged = false;

			outer: for (int i = 0; i < paths.size(); i++) {
				FloatPath a = paths.get(i);

				for (int j = i + 1; j < paths.size(); j++) {
					FloatPath b = paths.get(j);

					FloatPath res = tryMergeDirect(a, b);
					if (res != null) {
						// res is the merged path; it is either 'a' or 'b'
						if (res == a) {
							paths.set(i, a);
							paths.remove(j);
						} else { // res == b
							paths.set(i, b);
							paths.remove(j);
						}
						merged = true;
						break outer; // restart scanning after any merge
					}
				}
			}
		} while (merged);
	}

	private static boolean sameXY(float ax, float ay, float bx, float by) {
		return ax == bx && ay == by;
	}

	/**
	 * Try to merge two paths if any endpoints match exactly. Returns the merged
	 * path (either a or b) or null if no merge.
	 */
	private static FloatPath tryMergeDirect(FloatPath a, FloatPath b) {
		float aFx = a.firstX(), aFy = a.firstY();
		float aLx = a.lastX(), aLy = a.lastY();
		float bFx = b.firstX(), bFy = b.firstY();
		float bLx = b.lastX(), bLy = b.lastY();

		// A.end == B.start => A += B
		if (sameXY(aLx, aLy, bFx, bFy)) {
			a.appendPath(b, true);
			return a;
		}

		// A.end == B.end => reverse B, A += B
		if (sameXY(aLx, aLy, bLx, bLy)) {
			b.reversePairs();
			a.appendPath(b, true);
			return a;
		}

		// A.start == B.end => B += A (result is B)
		if (sameXY(aFx, aFy, bLx, bLy)) {
			b.appendPath(a, true);
			return b;
		}

		// A.start == B.start => reverse B, B += A (result is B)
		if (sameXY(aFx, aFy, bFx, bFy)) {
			b.reversePairs();
			b.appendPath(a, true);
			return b;
		}

		return null;
	}

	// Edge ids for a cell (i,j) with corners:
	// v0 = (i, j) bottom-left
	// v1 = (i+1, j) bottom-right
	// v2 = (i+1, j+1) top-right
	// v3 = (i, j+1) top-left
	// edges:
	// 0 = bottom (v0-v1)
	// 1 = right (v1-v2)
	// 2 = top (v3-v2)
	// 3 = left (v0-v3)

	private static final int E_BOTTOM = 0;
	private static final int E_RIGHT = 1;
	private static final int E_TOP = 2;
	private static final int E_LEFT = 3;

	/**
	 * Returns the edge opposite to the given edge in a square cell.
	 * 
	 * @param e edge index [0..3]
	 * @return opposite edge index
	 */
	private static int oppositeEdge(int e) {
		return switch (e) {
			case E_BOTTOM -> E_TOP;
			case E_TOP -> E_BOTTOM;
			case E_LEFT -> E_RIGHT;
			case E_RIGHT -> E_LEFT;
			default -> throw new IllegalArgumentException("bad edge " + e);
		};
	}

	/**
	 * Returns a bitmask representing the given edge.
	 */
	private static int edgeBit(int e) {
		return 1 << e;
	}

	/** Which edges are intersected for a given marching squares case (0..15). */
	private static int edgesUsedMask(int code) {
		return switch (code) {
			case 0, 15 -> 0;
			case 1, 14 -> edgeBit(E_LEFT) | edgeBit(E_BOTTOM);
			case 2, 13 -> edgeBit(E_BOTTOM) | edgeBit(E_RIGHT);
			case 3, 12 -> edgeBit(E_LEFT) | edgeBit(E_RIGHT);
			case 4, 11 -> edgeBit(E_RIGHT) | edgeBit(E_TOP);
			case 6, 9 -> edgeBit(E_BOTTOM) | edgeBit(E_TOP);
			case 7, 8 -> edgeBit(E_LEFT) | edgeBit(E_TOP);
			case 5, 10 -> edgeBit(E_BOTTOM) | edgeBit(E_RIGHT) | edgeBit(E_TOP) | edgeBit(E_LEFT);
			default -> 0;
		};
	}

	/**
	 * Given a cell and an entry edge that is intersected, return the exit edge (the
	 * other end of the segment inside that cell). For ambiguous cases 5/10 we
	 * consult amb[cell].
	 */
	private static int partnerEdge(int code, int ambMode, int entryEdge) {
		return switch (code) {
			case 1, 14 -> (entryEdge == E_LEFT) ? E_BOTTOM : E_LEFT;
			case 2, 13 -> (entryEdge == E_BOTTOM) ? E_RIGHT : E_BOTTOM;
			case 3, 12 -> (entryEdge == E_LEFT) ? E_RIGHT : E_LEFT;
			case 4, 11 -> (entryEdge == E_RIGHT) ? E_TOP : E_RIGHT;
			case 6, 9 -> (entryEdge == E_BOTTOM) ? E_TOP : E_BOTTOM;
			case 7, 8 -> (entryEdge == E_LEFT) ? E_TOP : E_LEFT;

			case 5 -> {
				// case 5: v0 and v2 are "high" corners (diagonal)
				// ambMode==1 connect high corners => pairs (0-3) and (1-2)
				// ambMode==0 connect other => pairs (0-1) and (2-3)
				if (ambMode == 1) {
					yield switch (entryEdge) {
						case E_BOTTOM -> E_LEFT;
						case E_LEFT -> E_BOTTOM;
						case E_RIGHT -> E_TOP;
						case E_TOP -> E_RIGHT;
						default -> -1;
					};
				} else {
					yield switch (entryEdge) {
						case E_BOTTOM -> E_RIGHT;
						case E_RIGHT -> E_BOTTOM;
						case E_TOP -> E_LEFT;
						case E_LEFT -> E_TOP;
						default -> -1;
					};
				}
			}

			case 10 -> {
				// case 10: v1 and v3 are "high" corners (other diagonal)
				// ambMode==1 connect high corners => pairs (0-1) and (2-3)
				// ambMode==0 connect other => pairs (0-3) and (1-2)
				if (ambMode == 1) {
					yield switch (entryEdge) {
						case E_BOTTOM -> E_RIGHT;
						case E_RIGHT -> E_BOTTOM;
						case E_TOP -> E_LEFT;
						case E_LEFT -> E_TOP;
						default -> -1;
					};
				} else {
					yield switch (entryEdge) {
						case E_BOTTOM -> E_LEFT;
						case E_LEFT -> E_BOTTOM;
						case E_RIGHT -> E_TOP;
						case E_TOP -> E_RIGHT;
						default -> -1;
					};
				}
			}

			default -> -1;
		};
	}

	/**
	 * Trace one polyline starting from (startCell,startEdge) at a given level.
	 * Marks visited edges along the way.
	 */
	private static FloatPath traceOne(float level, int startCell, int startEdge, byte[] codes, byte[] amb, byte[] visited, float[] z, int nx, int ny, float x0,
			float y0, float dx, float dy) {
		final int cellsX = nx - 1;
		final int cellsY = ny - 1;

		int cell = startCell;
		int edge = startEdge;

		FloatPath path = new FloatPath(128);

		final float[] tmp0 = new float[2];
		final float[] tmp1 = new float[2];

		edgePoint(level, cell, edge, z, nx, x0, y0, dx, dy, tmp0);
		path.add(tmp0[0], tmp0[1]);

		final int maxSteps = cellsX * cellsY * 4;

		for (int steps = 0; steps < maxSteps; steps++) {
			final int code = codes[cell] & 0xFF;
			final int exit = partnerEdge(code, amb[cell] & 0xFF, edge);
			if (exit < 0) {
				break;
			}

			visited[cell] = (byte) (visited[cell] | edgeBit(edge) | edgeBit(exit));

			edgePoint(level, cell, exit, z, nx, x0, y0, dx, dy, tmp1);
			path.addIfDifferent(tmp1[0], tmp1[1]);

			final int cx = cell % cellsX;
			final int cy = cell / cellsX;

			int ncx = cx, ncy = cy;
			switch (exit) {
				case E_BOTTOM -> ncy = cy - 1;
				case E_TOP -> ncy = cy + 1;
				case E_LEFT -> ncx = cx - 1;
				case E_RIGHT -> ncx = cx + 1;
			}

			// leaving grid => open contour
			if (ncx < 0 || ncx >= cellsX || ncy < 0 || ncy >= cellsY) {
				break;
			}

			final int nextCell = ncy * cellsX + ncx;
			final int nextEdge = oppositeEdge(exit);

			// closed loop
			if (nextCell == startCell && nextEdge == startEdge) {
				break;
			}

			// stop if next cell doesn't contain this entry edge, or already traced there
			final int nextCode = codes[nextCell] & 0xFF;
			if ((edgesUsedMask(nextCode) & edgeBit(nextEdge)) == 0) {
				break;
			}
			if ((visited[nextCell] & edgeBit(nextEdge)) != 0) {
				break;
			}

			cell = nextCell;
			edge = nextEdge;
		}

		return path.sizePairs() >= 2 ? path : null;
	}

	/** Compute interpolated (x,y) point where contour level crosses a cell edge. */
	private static void edgePoint(float level, int cell, int edge, float[] z, int nx, float x0, float y0, float dx, float dy, float[] outXY) {
		final int cellsX = nx - 1;
		final int ix = cell % cellsX;
		final int iy = cell / cellsX;

		final int row0 = iy * nx;
		final int row1 = (iy + 1) * nx;

		final float v0 = z[row0 + ix];
		final float v1 = z[row0 + ix + 1];
		final float v2 = z[row1 + ix + 1];
		final float v3 = z[row1 + ix];

		final float X0 = x0 + ix * dx;
		final float X1 = X0 + dx;
		final float Y0 = y0 + iy * dy;
		final float Y1 = Y0 + dy;

		float t, x, y;

		switch (edge) {
			case E_BOTTOM -> {
				t = interp(level, v0, v1);
				x = X0 + t * (X1 - X0);
				y = Y0;
			}
			case E_RIGHT -> {
				t = interp(level, v1, v2);
				x = X1;
				y = Y0 + t * (Y1 - Y0);
			}
			case E_TOP -> {
				t = interp(level, v3, v2);
				x = X0 + t * (X1 - X0);
				y = Y1;
			}
			case E_LEFT -> {
				t = interp(level, v0, v3);
				x = X0;
				y = Y0 + t * (Y1 - Y0);
			}
			default -> throw new IllegalArgumentException("bad edge " + edge);
		}

		outXY[0] = x;
		outXY[1] = y;
	}

	/**
	 * Linear interpolation between two values to find where a target level falls.
	 * 
	 * @param level target iso-value
	 * @param a     first value
	 * @param b     second value
	 * @return interpolation factor [0..1]
	 */
	private static float interp(float level, float a, float b) {
		float d = (b - a);
		if (d == 0f) {
			// Deterministic choice so adjacent cells compute identical points.
			// If the whole edge is exactly on the contour, pick the first endpoint.
			return (level == a) ? 0f : 0.5f;
		}
		return (level - a) / d;
	}

	/**
	 * Converts a FloatPath into a Processing PShape PATH.
	 */
	private static PShape toPShape(FloatPath path) {
		PShape s = new PShape();
		s.setFamily(PShape.PATH);
		s.setStroke(true);
		s.setStroke(128);
		s.setStrokeWeight(2);
		s.setFill(false);

		boolean closed = path.isClosed(1e-4f);
		int vertexCount = closed ? path.sizePairs() - 1 : path.sizePairs();

		s.beginShape();
		for (int i = 0; i < vertexCount; i++) {
			float x = path.x(i);
			float y = path.y(i);
			s.vertex(x, y);
		}
		s.endShape(closed ? PConstants.CLOSE : PConstants.OPEN);

		return s;
	}

	/**
	 * A lightweight list-of-floats structure for storing (x,y) coordinate pairs of
	 * a path. Avoids the overhead of many PVector objects.
	 */
	private static final class FloatPath {
		private float[] data;
		private int size; // number of floats used (even)

		/**
		 * @param initialPairsCapacity initial number of (x,y) pairs to allocate space
		 *                             for
		 */
		FloatPath(int initialPairsCapacity) {
			data = new float[Math.max(16, initialPairsCapacity * 2)];
			size = 0;
		}

		int sizePairs() {
			return size >> 1;
		}

		float x(int i) {
			return data[i << 1];
		}

		float y(int i) {
			return data[(i << 1) + 1];
		}

		void add(float x, float y) {
			int ns = size + 2;
			if (ns > data.length) {
				data = Arrays.copyOf(data, data.length << 1);
			}
			data[size] = x;
			data[size + 1] = y;
			size = ns;
		}

		void addIfDifferent(float x, float y) {
			if (size >= 2) {
				float lx = data[size - 2];
				float ly = data[size - 1];
				if (lx == x && ly == y) {
					return;
				}
			}
			add(x, y);
		}

		boolean isClosed(float eps) {
			if (sizePairs() < 3) {
				return false;
			}
			float x0 = data[0], y0 = data[1];
			float xn = data[size - 2], yn = data[size - 1];
			return (Math.abs(x0 - xn) <= eps) && (Math.abs(y0 - yn) <= eps);
		}

		float firstX() {
			return data[0];
		}

		float firstY() {
			return data[1];
		}

		float lastX() {
			return data[size - 2];
		}

		float lastY() {
			return data[size - 1];
		}

		void reversePairs() {
			int n = sizePairs();
			for (int i = 0, j = n - 1; i < j; i++, j--) {
				int ia = i << 1;
				int ja = j << 1;
				float tx = data[ia], ty = data[ia + 1];
				data[ia] = data[ja];
				data[ia + 1] = data[ja + 1];
				data[ja] = tx;
				data[ja + 1] = ty;
			}
		}

		void appendPath(FloatPath other, boolean skipFirstPair) {
			int startPair = skipFirstPair ? 1 : 0;
			int otherPairs = other.sizePairs();
			for (int i = startPair; i < otherPairs; i++) {
				add(other.x(i), other.y(i));
			}
		}

		void snapClosed(float eps) {
			if (sizePairs() < 3)
				return;
			float x0 = firstX(), y0 = firstY();
			float xn = lastX(), yn = lastY();
			if (Math.abs(x0 - xn) <= eps && Math.abs(y0 - yn) <= eps) {
				data[size - 2] = x0;
				data[size - 1] = y0;
			}
		}
	}
}