package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateSequence;
import org.locationtech.jts.geom.CoordinateSequences;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;

import net.jafama.FastMath;

/**
 * Line gaussian smoothing.
 * 
 * @author Julien Gaffuri
 * @author Michael Carleton
 *
 */
public class GaussianLineSmoothing {

	// based on https://github.com/locationtech/jts/pull/478

	// Tuning knobs
	private static final int SAMPLES_PER_SIGMA = 6; // 6..8 recommended
	private static final double CUTOFF_SIGMAS = 5.0; // truncate kernel at 5σ
	private static final int MAX_SAMPLES = 8192; // upper cap per geometry
	private static final int COUNT_QUANTUM = 32; // quantize grid size (stability)

	private GaussianLineSmoothing() {
	}

	/**
	 * Line gaussian smoothing. The position of each point is the average position
	 * of its neighbors, weighted by a gaussian kernel. For non-closed lines, the
	 * initial and final points are preserved.
	 *
	 * @param input The input line
	 * @param sigma The standard deviation of the gaussian kernel. The larger, the
	 *              more smoothed.
	 */
	public static LineString get(LineString line, double sigmaM) {
		if (line == null) {
			return null;
		}
		final int np = line.getNumPoints();
		if (np <= 2 || sigmaM <= 0) {
			return (LineString) line.copy();
		}

		boolean isClosed = line.isClosed();
		double length = line.getLength();
		if (length == 0) {
			return (LineString) line.copy();
		}

		// Extreme sigma fallback (same as your original idea)
		if (sigmaM > 0.333 * length) {
			GeometryFactory gf = line.getFactory();
			if (isClosed) {
				return gf.createLineString();
			} else {
				return line.getFactory().createLineString(new Coordinate[] { line.getCoordinateN(0), line.getCoordinateN(np - 1) });
			}
		}

		LineString src = line;
		if (isClosed) {
			/*
			 * LineString.norm() does not reorder coordinates so we norm with method below.
			 * Output of Gaussian smoothing is sensitive to vertex order!
			 */
			normalize((LinearRing) line, true);
		}

		// Build a STABLE uniform arc-length grid:
		// - density ~ SAMPLES_PER_SIGMA per sigma
		// - sample count quantized to COUNT_QUANTUM
		// - capped by MAX_SAMPLES
		Resampled rs = resampleUniformStable(src, sigmaM, isClosed);
		Coordinate[] samples = rs.coords;
		int n = samples.length;
		int M = isClosed ? (n - 1) : n; // number of unique samples to convolve
		double step = rs.step;

		// Pack into arrays
		double[] xs = new double[M];
		double[] ys = new double[M];
		for (int i = 0; i < M; i++) {
			xs[i] = samples[i].x;
			ys[i] = samples[i].y;
		}

		// Precompute Gaussian weights (truncated) using exact distances
		int halfWidth = Math.max(1, (int) Math.ceil(CUTOFF_SIGMAS * sigmaM / step));
		if (isClosed) {
			halfWidth = Math.min(halfWidth, (M - 1) / 2);
		}
		double[] w = new double[halfWidth + 1];
		final double twoSigma2 = 2.0 * sigmaM * sigmaM;
		w[0] = 1.0;
		for (int j = 1; j <= halfWidth; j++) {
			double d = j * step;
			w[j] = FastMath.exp(-(d * d) / twoSigma2);
		}

		// Convolve with per-point renormalization (exact discrete Gaussian on this
		// grid)
		double[] ox = new double[M];
		double[] oy = new double[M];

		if (isClosed) {
			for (int i = 0; i < M; i++) {
				double sx = xs[i] * w[0];
				double sy = ys[i] * w[0];
				double sw = w[0];
				for (int j = 1; j <= halfWidth; j++) {
					int il = wrap(i - j, M);
					int ir = wrap(i + j, M);
					double ww = w[j];
					sx += ww * (xs[il] + xs[ir]);
					sy += ww * (ys[il] + ys[ir]);
					sw += 2.0 * ww;
				}
				ox[i] = sx / sw;
				oy[i] = sy / sw;
			}
		} else {
			// Reflect at the ends; also preserve endpoints exactly
			for (int i = 0; i < M; i++) {
				if (i == 0 || i == M - 1) {
					ox[i] = xs[i];
					oy[i] = ys[i];
					continue;
				}
				double sx = xs[i] * w[0];
				double sy = ys[i] * w[0];
				double sw = w[0];
				for (int j = 1; j <= halfWidth; j++) {
					int il = reflect(i - j, M);
					int ir = reflect(i + j, M);
					double ww = w[j];
					sx += ww * (xs[il] + xs[ir]);
					sy += ww * (ys[il] + ys[ir]);
					sw += 2.0 * ww;
				}
				ox[i] = sx / sw;
				oy[i] = sy / sw;
			}
		}

		// Rebuild geometry
		Coordinate[] out;
		if (isClosed) {
			out = new Coordinate[M + 1];
			for (int i = 0; i < M; i++) {
				out[i] = new Coordinate(ox[i], oy[i]);
			}
			out[M] = new Coordinate(out[0]);
		} else {
			out = new Coordinate[M];
			for (int i = 0; i < M; i++) {
				out[i] = new Coordinate(ox[i], oy[i]);
			}
			out[0] = new Coordinate(samples[0]); // exact endpoints
			out[M - 1] = new Coordinate(samples[n - 1]);
		}
		return line.getFactory().createLineString(out);
	}

	// Stable resampling:
	// - samples-per-sigma fixed
	// - count quantized (COUNT_QUANTUM)
	// - independent of input vertex placement
	private static Resampled resampleUniformStable(LineString line, double sigmaM, boolean isClosed) {
		Coordinate[] in = line.getCoordinates();
		int n = in.length;
		if (n < 2) {
			return new Resampled(in, 1.0);
		}

		double total = 0.0;
		for (int i = 1; i < n; i++) {
			total += dist(in[i - 1], in[i]);
		}
		if (isClosed && !in[0].equals2D(in[n - 1])) {
			total += dist(in[n - 1], in[0]);
		}
		if (total == 0.0) {
			return new Resampled(in, 1.0);
		}

		// Target count based on samples-per-sigma, with caps and quantization for
		// stability
		int target = Math.max(8, (int) Math.ceil((SAMPLES_PER_SIGMA * total) / Math.max(sigmaM, 1e-12)));
		target = Math.min(target, MAX_SAMPLES);
		if (isClosed) {
			// ensure at least 4 samples
			target = Math.max(4, target);
		} else {
			target = Math.max(2, target);
		}
		// Quantize to a multiple of COUNT_QUANTUM to avoid flicker from tiny length
		// changes
		if (target > COUNT_QUANTUM) {
			int q = (target + COUNT_QUANTUM / 2) / COUNT_QUANTUM;
			target = Math.max(COUNT_QUANTUM, q * COUNT_QUANTUM);
			target = Math.min(target, MAX_SAMPLES);
		}

		double step = total / target;

		if (isClosed) {
			List<Coordinate> out = new ArrayList<>(target + 1);
			double targetS = 0.0;
			int segIdx = 0;
			double acc = 0.0;
			Coordinate a = in[0];
			Coordinate b = (n > 1 ? in[1] : in[0]);
			double segLen = dist(a, b);

			for (int k = 0; k < target; k++) {
				while (acc + segLen < targetS) {
					acc += segLen;
					segIdx++;
					int next = (segIdx + 1 < n) ? segIdx + 1 : 0;
					a = in[segIdx % n];
					b = in[next];
					segLen = dist(a, b);
				}
				double t = segLen > 0 ? (targetS - acc) / segLen : 0.0;
				out.add(new Coordinate(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)));
				targetS += step;
			}
			out.add(new Coordinate(out.get(0)));
			return new Resampled(out.toArray(new Coordinate[0]), step);
		} else {
			List<Coordinate> out = new ArrayList<>(target + 1);
			out.add(new Coordinate(in[0]));
			double targetS = step;
			int segIdx = 0;
			double acc = 0.0;
			Coordinate a = in[0];
			Coordinate b = in[1];
			double segLen = dist(a, b);

			while (targetS < total && segIdx < n - 1) {
				while (acc + segLen < targetS && segIdx < n - 1) {
					acc += segLen;
					segIdx++;
					a = in[segIdx];
					b = (segIdx + 1 < n) ? in[segIdx + 1] : in[segIdx];
					segLen = dist(a, b);
				}
				double t = segLen > 0 ? (targetS - acc) / segLen : 0.0;
				out.add(new Coordinate(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)));
				targetS += step;
			}
			out.add(new Coordinate(in[n - 1]));
			return new Resampled(out.toArray(new Coordinate[0]), step);
		}
	}

	private static int wrap(int i, int n) {
		int r = i % n;
		return r < 0 ? r + n : r;
	}

	private static int reflect(int i, int n) {
		if (n == 1) {
			return 0;
		}
		while (i < 0 || i >= n) {
			if (i < 0) {
				i = -i - 1;
			} else {
				i = 2 * n - i - 1;
			}
		}
		return i;
	}

	private static double dist(Coordinate a, Coordinate b) {
		return a.distance(b);
	}

	private static void normalize(LinearRing ring, boolean clockwise) {
		if (ring.isEmpty()) {
			return;
		}

		CoordinateSequence seq = ring.getCoordinateSequence();
		int minCoordinateIndex = CoordinateSequences.minCoordinateIndex(seq, 0, seq.size() - 2);
		CoordinateSequences.scroll(seq, minCoordinateIndex, true);
		if (Orientation.isCCW(seq) == clockwise) {
			CoordinateSequences.reverse(seq);
		}
	}

	private static final class Resampled {
		final Coordinate[] coords;
		final double step;

		Resampled(Coordinate[] coords, double step) {
			this.coords = coords;
			this.step = step;
		}
	}
}