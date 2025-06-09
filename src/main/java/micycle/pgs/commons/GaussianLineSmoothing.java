package micycle.pgs.commons;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;

/**
 * Line gaussian smoothing.
 *
 * @author Julien Gaffuri
 * @author Michael Carleton
 */
public class GaussianLineSmoothing {
	// from https://github.com/locationtech/jts/pull/478

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
	public static LineString get(LineString input, final double sigma) {
		if (input.getCoordinates().length <= 2) {
			return (LineString) input.copy();
		}

		final boolean isClosed = input.isClosed();
		input = (LineString) input.norm();
		final double Δx = sigma / 3;
//		double length = input.getLength();
//		if (Δx > length * 0.25) {
//		}

		// 1) Densify
		final LineString dens = (LineString) LittleThumblingDensifier.densify(input, Δx);
		final Coordinate[] coords = dens.getCoordinates();

		// 2) Build/normalize kernel
		int n = (int) Math.ceil(7 * sigma / Δx);
		n = Math.min(n, coords.length - 1);
		final double d = Δx * Δx, b = 2 * sigma * sigma;
		final double[] g = new double[n + 1];
		for (int i = 0; i <= n; i++) {
			g[i] = Math.exp(-i * (double) i * d / b);
		}
		// scale so ∑ g[i]*Δx = 1
		for (int i = 0; i <= n; i++) {
			g[i] *= Δx;
		}
		double sum = g[0];
		for (int i = 1; i <= n; i++) {
			sum += 2 * g[i];
		}
		for (int i = 0; i <= n; i++) {
			g[i] /= sum;
		}

		// Original densified coords always have coords[0]==coords[N] for a closed ring:
		final Coordinate[] densCoords = dens.getCoordinates();

		Coordinate[] work; // the unique points we will actually smooth
		int M; // number of unique points

		if (isClosed) {
			// toss the last duplicate in the densifier’s output
			M = densCoords.length - 1;
			work = new Coordinate[M];
			System.arraycopy(densCoords, 0, work, 0, M);
		} else {
			M = densCoords.length;
			work = densCoords;
		}

		// clamp radius so we never reflect/pad past our data
		n = Math.min(n, M - 1);

		final Coordinate[] pad = new Coordinate[M + 2 * n];
		if (isClosed) {
			// circular wrap
			for (int i = 0; i < n; i++) {
				pad[i] = work[(M - n + i) % M]; // last n points
				pad[n + M + i] = work[i]; // first n points
			}
			System.arraycopy(work, 0, pad, n, M);
		} else {
			for (int i = 0; i < n; i++) {
				pad[i] = work[n - i]; // reflect left
				pad[n + M + i] = work[M - 2 - i]; // reflect right
			}
			System.arraycopy(work, 0, pad, n, M);
		}

		Coordinate[] smooth;
		if (isClosed) {
			// we'll produce M+1 output coords: [0..M-1] + duplicate of 0
			smooth = new Coordinate[M + 1];
			for (int i = 0; i < M; i++) {
				double sx = 0, sy = 0;
				for (int j = -n; j <= n; j++) {
					final Coordinate c = pad[i + n + j];
					final double w = g[Math.abs(j)];
					sx += c.x * w;
					sy += c.y * w;
				}
				smooth[i] = new Coordinate(sx, sy);
			}
			// close the ring
			smooth[M] = smooth[0];

		} else {
			// open line: convolve interior and then clamp the two endpoints
			smooth = new Coordinate[M];
			// interior
			for (int i = 1; i < M - 1; i++) {
				double sx = 0, sy = 0;
				for (int j = -n; j <= n; j++) {
					final Coordinate c = pad[i + n + j];
					final double w = g[Math.abs(j)];
					sx += c.x * w;
					sy += c.y * w;
				}
				smooth[i] = new Coordinate(sx, sy);
			}
			// preserve original endpoints
			smooth[0] = work[0];
			smooth[M - 1] = work[M - 1];
		}

		final LineString result = input.getFactory().createLineString(smooth);
		return result;
	}

}