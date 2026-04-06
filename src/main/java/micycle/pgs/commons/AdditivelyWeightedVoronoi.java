package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.stream.IntStream;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;

import net.jafama.FastMath;

/**
 * Approximates 2D additively-weighted Voronoi cells (Apollonius diagram) for
 * input sites encoded as {@link Coordinate}s: (x,y) is the site center and z is
 * the additive weight (radius).
 *
 * <p>
 * Each cell boundary is represented in polar form around the site center as a
 * "lifted radius" function t(θ). For each direction θ, competitors induce a
 * candidate radius; the boundary is the lower envelope (minimum) over all
 * competitors. The envelope is approximated by angular sampling with adaptive
 * refinement. This idea is based on: Hans-Martin Will, <i>“Practical and
 * efficient computation of additively weighted Voronoi cells for applications
 * in molecular biology”, 1998</i>.
 *
 * <p>
 * Cells are bounded to a rectangular clip envelope. For each sampled direction,
 * the radius is upper-bounded by the ray's exit distance from the rectangle.
 * This produces a polygon that is already confined to the clip rectangle.
 *
 * <p>
 * Dominated (empty) cells are omitted from the returned list.
 *
 * <p>
 * Notes:
 * <ul>
 * <li>For best results, the clip envelope should contain all sites.</li>
 * <li>This implementation returns polygon shells only (no holes).</li>
 * </ul>
 * 
 * @author Michael Carleton
 */
public class AdditivelyWeightedVoronoi {

	private static class AV2DOptions {

		/**
		 * Adaptive sampling error tolerance (squared). Smaller values yield smoother
		 * curves at the cost of more vertices.
		 */
		public double errorToleranceSq = 0.01 * 0.01;

		/** Maximum recursion depth for adaptive sampling. */
		public int maxRecursionDepth = 10;

		/** Scale factor for the default clip envelope when none is provided. */
		public double boundsScale = 1.5;

		/**
		 * Small epsilon for denominator/degeneracy checks. This is an absolute epsilon;
		 * for very large coordinate magnitudes, consider increasing it.
		 */
		public double eps = 1e-12;
	}

	private final GeometryFactory gf;
	private final AV2DOptions opt;

	/**
	 * Creates an instance using the given geometry factory and an adaptive sampling
	 * tolerance.
	 *
	 * @param gf             geometry factory used to create output polygons
	 * @param errorTolerance maximum allowed deviation (in coordinate units) used by
	 *                       the adaptive sampler; smaller values produce smoother
	 *                       cell boundaries at the cost of more vertices
	 */
	public AdditivelyWeightedVoronoi(GeometryFactory gf, double errorTolerance) {
		this.gf = gf;
		this.opt = new AV2DOptions();
		this.opt.errorToleranceSq = errorTolerance * errorTolerance;
	}

	/**
	 * Computes clipped polygonal approximations of additively-weighted Voronoi
	 * cells.
	 *
	 * <p>
	 * Input sites are provided as {@link Coordinate}s where (x,y) is the site
	 * location and z is its additive weight. Cells are bounded to the supplied
	 * {@code bounds}, or to a default envelope derived from the data if
	 * {@code bounds} is {@code null}.
	 *
	 * <p>
	 * Dominated (empty) cells are omitted from the returned list. Each returned
	 * polygon has its {@linkplain Polygon#getUserData() userData} set to the
	 * corresponding input site index.
	 *
	 * @param sites  input sites; {@code (x,y)} is the site center and {@code z} is
	 *               the additive weight (non-finite z is treated as 0)
	 * @param bounds optional clip envelope; if {@code null}, a default is computed
	 *               from the input extent and weight magnitudes
	 * @return list of polygons for non-empty cells; order is unspecified due to
	 *         parallel execution
	 */
	public List<Polygon> computeCells(List<Coordinate> sites, Envelope bounds) {
		if (sites == null || sites.isEmpty()) {
			return List.of();
		}

		final int n = sites.size();
		final double[] x = new double[n];
		final double[] y = new double[n];
		final double[] w = new double[n];

		final Envelope dataEnv = new Envelope();
		double maxAbsW = 0.0;

		for (int i = 0; i < n; i++) {
			Coordinate c = sites.get(i);
			x[i] = c.getX();
			y[i] = c.getY();
			double zi = c.getZ();
			w[i] = Double.isFinite(zi) ? zi : 0.0;
			maxAbsW = Math.max(maxAbsW, Math.abs(w[i]));
			dataEnv.expandToInclude(x[i], y[i]);
		}

		final Envelope clipEnv = bounds != null ? new Envelope(bounds) : defaultClipEnvelope(dataEnv, maxAbsW);

		final double minX = clipEnv.getMinX(), maxX = clipEnv.getMaxX();
		final double minY = clipEnv.getMinY(), maxY = clipEnv.getMaxY();

		final Precompute pc = precomputePairs(x, y, w);

		List<Polygon> out = IntStream.range(0, n).parallel().mapToObj(i -> {
			if (pc.empty[i]) {
				return null;
			}

			// per-cell prune bound: farthest corner distance from site to clip rectangle
			final double Rsite = radiusToCoverEnvelopeFromSite(clipEnv, x[i], y[i]);

			Polygon approx = approximateCellPolygon(i, x, y, w, pc.candidates[i], Rsite, minX, maxX, minY, maxY);
			approx.setUserData(i);
			return approx;
		}).filter(Objects::nonNull).toList();

		return out;
	}

	private static double radiusToCoverEnvelopeFromSite(Envelope env, double xi, double yi) {
		double minX = env.getMinX(), maxX = env.getMaxX();
		double minY = env.getMinY(), maxY = env.getMaxY();

		double d1 = FastMath.hypot(minX - xi, minY - yi);
		double d2 = FastMath.hypot(minX - xi, maxY - yi);
		double d3 = FastMath.hypot(maxX - xi, minY - yi);
		double d4 = FastMath.hypot(maxX - xi, maxY - yi);

		double R = Math.max(Math.max(d1, d2), Math.max(d3, d4));
		return R * 1.0000001 + 1e-9;
	}

	private Polygon approximateCellPolygon(int i, double[] x, double[] y, double[] w, int[] candidateIndices, double Rsite, double minX, double maxX,
			double minY, double maxY) {
		final double xi = x[i], yi = y[i], wi = w[i];

		// no competitors => within the clip region the cell is the whole rectangle
		if (candidateIndices == null || candidateIndices.length == 0) {
			return createClipRectanglePolygon(minX, maxX, minY, maxY);
		}

		// sort candidates by a metric monotone with the lower bound on bisector
		// distance
		Candidate[] sortedCands = new Candidate[candidateIndices.length];
		for (int k = 0; k < candidateIndices.length; k++) {
			int j = candidateIndices[k];
			double d = FastMath.hypot(x[j] - xi, y[j] - yi);
			sortedCands[k] = new Candidate(j, d, d - w[j]);
		}
		Arrays.sort(sortedCands);

		// build compact primitive arrays and prune competitors that cannot matter
		// inside the clip
		int cN = sortedCands.length;
		double[] cxTmp = new double[cN];
		double[] cyTmp = new double[cN];
		double[] rTmp = new double[cN];
		double[] numTmp = new double[cN]; // |c|^2 - r^2
		double[] minBDTmp = new double[cN]; // lower bound on distance to bisector

		int m = 0;
		for (int t = 0; t < cN; t++) {
			int j = sortedCands[t].index;

			double dx = x[j] - xi;
			double dy = y[j] - yi;
			double d = sortedCands[t].dist;
			double rr = w[j] - wi;

			double minBisectorDist = (d - rr) * 0.5;

			if (minBisectorDist >= Rsite) {
				continue;
			}

			cxTmp[m] = dx;
			cyTmp[m] = dy;
			rTmp[m] = rr;
			numTmp[m] = (dx * dx + dy * dy) - (rr * rr);
			minBDTmp[m] = minBisectorDist;
			m++;
		}

		if (m == 0) {
			return createClipRectanglePolygon(minX, maxX, minY, maxY);
		}

		double[] cx = Arrays.copyOf(cxTmp, m);
		double[] cy = Arrays.copyOf(cyTmp, m);
		double[] r = Arrays.copyOf(rTmp, m);
		double[] num = Arrays.copyOf(numTmp, m);
		double[] minBisectorDist = Arrays.copyOf(minBDTmp, m);

		// adaptive sampling with carry-forward sector endpoints
		List<Coordinate> pts = new ArrayList<>();

		int initialSectors = 8;
		double dt = 2.0 * Math.PI / initialSectors;

		double tPrev = 0.0;
		double rPrev = computeRadius(cx, cy, r, num, minBisectorDist, xi, yi, minX, maxX, minY, maxY, tPrev);

		{
			double x1 = rPrev * FastMath.cos(tPrev);
			double y1 = rPrev * FastMath.sin(tPrev);
			pts.add(new Coordinate(xi + x1, yi + y1));
		}

		for (int k = 0; k < initialSectors; k++) {
			double tNext = (k + 1) * dt;
			double rNext = computeRadius(cx, cy, r, num, minBisectorDist, xi, yi, minX, maxX, minY, maxY, tNext);

			double x1 = rPrev * FastMath.cos(tPrev);
			double y1 = rPrev * FastMath.sin(tPrev);
			double x2 = rNext * FastMath.cos(tNext);
			double y2 = rNext * FastMath.sin(tNext);

			sampleRecursively(xi, yi, cx, cy, r, num, minBisectorDist, minX, maxX, minY, maxY, tPrev, rPrev, x1, y1, tNext, rNext, x2, y2, pts, 0);

			pts.add(new Coordinate(xi + x2, yi + y2));

			tPrev = tNext;
			rPrev = rNext;
		}

		pts = removeNearDuplicates(pts, 1e-6);
		if (pts.size() < 3) {
			return createClipRectanglePolygon(minX, maxX, minY, maxY);
		}

		if (!pts.get(0).equals2D(pts.get(pts.size() - 1))) {
			pts.add(new Coordinate(pts.get(0)));
		}

		return gf.createPolygon(pts.toArray(new Coordinate[0]));
	}

	private Polygon createClipRectanglePolygon(double minX, double maxX, double minY, double maxY) {
		Coordinate[] pts = new Coordinate[] { new Coordinate(minX, minY), new Coordinate(maxX, minY), new Coordinate(maxX, maxY), new Coordinate(minX, maxY),
				new Coordinate(minX, minY) };
		return gf.createPolygon(gf.createLinearRing(pts));
	}

	private void sampleRecursively(double xi, double yi, double[] cx, double[] cy, double[] r, double[] num, double[] minBisectorDist, double minX, double maxX,
			double minY, double maxY, double t1, double r1, double x1, double y1, double t2, double r2, double x2, double y2, List<Coordinate> resultPts,
			int depth) {

		if (depth >= opt.maxRecursionDepth) {
			return;
		}

		final double tMid = 0.5 * (t1 + t2);
		final double rMidActual = computeRadius(cx, cy, r, num, minBisectorDist, xi, yi, minX, maxX, minY, maxY, tMid);

		final double chordMx = 0.5 * (x1 + x2);
		final double chordMy = 0.5 * (y1 + y2);

		final double cosMid = FastMath.cos(tMid);
		final double sinMid = FastMath.sin(tMid);
		final double surfMx = rMidActual * cosMid;
		final double surfMy = rMidActual * sinMid;

		final double dx = chordMx - surfMx;
		final double dy = chordMy - surfMy;
		final double distSq = dx * dx + dy * dy;

		if (distSq > opt.errorToleranceSq) {
			sampleRecursively(xi, yi, cx, cy, r, num, minBisectorDist, minX, maxX, minY, maxY, t1, r1, x1, y1, tMid, rMidActual, surfMx, surfMy, resultPts,
					depth + 1);

			resultPts.add(new Coordinate(xi + surfMx, yi + surfMy));

			sampleRecursively(xi, yi, cx, cy, r, num, minBisectorDist, minX, maxX, minY, maxY, tMid, rMidActual, surfMx, surfMy, t2, r2, x2, y2, resultPts,
					depth + 1);
		}
	}

	private double computeRadius(double[] cx, double[] cy, double[] r, double[] num, double[] minBisectorDist, double xi, double yi, double minX, double maxX,
			double minY, double maxY, double theta) {

		final double px = FastMath.cos(theta);
		final double py = FastMath.sin(theta);

		// initialize with the ray's exit distance from the clip rectangle (tight upper
		// bound)
		double tMin = rayBoxExitDistance(xi, yi, px, py, minX, maxX, minY, maxY);
		if (!(tMin > 0.0) || !Double.isFinite(tMin)) {
			tMin = Double.POSITIVE_INFINITY;
		}

		for (int t = 0; t < cx.length; t++) {
			// candidates are sorted by a lower bound; stop once the bound exceeds the
			// current best
			if (minBisectorDist[t] > tMin) {
				break;
			}

			final double denomInner = r[t] + (px * cx[t] + py * cy[t]);
			if (denomInner <= opt.eps) {
				continue;
			}

			final double tt = num[t] / (2.0 * denomInner);
			if (tt > 0.0 && tt < tMin) {
				tMin = tt;
			}
		}

		return tMin;
	}

	private static double rayBoxExitDistance(double ox, double oy, double dx, double dy, double minX, double maxX, double minY, double maxY) {
		double tEnter = 0.0;
		double tExit = Double.POSITIVE_INFINITY;

		if (dx == 0.0) {
			if (ox < minX || ox > maxX) {
				return Double.NaN;
			}
		} else {
			double tx1 = (minX - ox) / dx;
			double tx2 = (maxX - ox) / dx;
			if (tx1 > tx2) {
				double tmp = tx1;
				tx1 = tx2;
				tx2 = tmp;
			}
			tEnter = Math.max(tEnter, tx1);
			tExit = Math.min(tExit, tx2);
			if (tExit < tEnter) {
				return Double.NaN;
			}
		}

		if (dy == 0.0) {
			if (oy < minY || oy > maxY) {
				return Double.NaN;
			}
		} else {
			double ty1 = (minY - oy) / dy;
			double ty2 = (maxY - oy) / dy;
			if (ty1 > ty2) {
				double tmp = ty1;
				ty1 = ty2;
				ty2 = tmp;
			}
			tEnter = Math.max(tEnter, ty1);
			tExit = Math.min(tExit, ty2);
			if (tExit < tEnter) {
				return Double.NaN;
			}
		}

		return tExit;
	}

	private static class Candidate implements Comparable<Candidate> {
		int index;
		double dist;
		double sortMetric;

		Candidate(int index, double dist, double sortMetric) {
			this.index = index;
			this.dist = dist;
			this.sortMetric = sortMetric;
		}

		@Override
		public int compareTo(Candidate o) {
			return Double.compare(this.sortMetric, o.sortMetric);
		}
	}

	private Envelope defaultClipEnvelope(Envelope dataEnv, double maxAbsW) {
		double diag = dataEnv.getDiameter();
		if (diag == 0) {
			diag = 1.0;
		}
		double half = opt.boundsScale * (diag + 2.0 * maxAbsW + 1.0);
		double cx = dataEnv.centre().x;
		double cy = dataEnv.centre().y;
		return new Envelope(cx - half, cx + half, cy - half, cy + half);
	}

	private static List<Coordinate> removeNearDuplicates(List<Coordinate> pts, double tol) {
		if (pts.isEmpty()) {
			return pts;
		}
		double tolSq = tol * tol;
		List<Coordinate> out = new ArrayList<>(pts.size());
		Coordinate prev = null;
		for (Coordinate c : pts) {
			if (prev == null || prev.distanceSq(c) > tolSq) {
				out.add(c);
				prev = c;
			}
		}
		if (out.size() >= 2 && out.get(0).distanceSq(out.get(out.size() - 1)) <= tolSq) {
			out.remove(out.size() - 1);
		}
		return out;
	}

	@SuppressWarnings("unchecked")
	private Precompute precomputePairs(double[] x, double[] y, double[] w) {
		int n = x.length;
		boolean[] empty = new boolean[n];
		List<Integer>[] cand = new ArrayList[n];
		for (int i = 0; i < n; i++) {
			cand[i] = new ArrayList<>();
		}

		double domTol = 10 * opt.eps;

		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double dx = x[j] - x[i];
				double dy = y[j] - y[i];
				double d = FastMath.hypot(dx, dy);

				double rij = w[j] - w[i];
				double rji = -rij;

				// dominance test: if wj - wi > dij then i is dominated everywhere
				if (rij > d + domTol) {
					empty[i] = true;
				}
				if (rji > d + domTol) {
					empty[j] = true;
				}

				// necessary condition for an in-front intersection along some ray: r + max(p·c)
				// > 0 => r + d > 0
				if (rij + d > opt.eps) {
					cand[i].add(j);
				}
				if (rji + d > opt.eps) {
					cand[j].add(i);
				}
			}
		}

		int[][] candidates = new int[n][];
		for (int i = 0; i < n; i++) {
			candidates[i] = cand[i].stream().mapToInt(Integer::intValue).toArray();
		}
		return new Precompute(empty, candidates);
	}

	static final class Precompute {
		final boolean[] empty;
		final int[][] candidates;

		Precompute(boolean[] empty, int[][] candidates) {
			this.empty = empty;
			this.candidates = candidates;
		}
	}
}