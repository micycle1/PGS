package micycle.pgs.commons;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;

import net.jafama.FastMath;

/**
 * Line gaussian smoothing.
 * 
 * @author Julien Gaffuri
 *
 */
public class GaussianLineSmoothing {

	// from https://github.com/locationtech/jts/pull/478

	private GaussianLineSmoothing() {
	}

	/**
	 * @param line
	 * @param sigmaM
	 */
	public static LineString get(LineString line, double sigmaM) {
		return get(line, sigmaM, -1);
	}

	/**
	 * Line gaussian smoothing. The position of each point is the average position
	 * of its neighbors, weighted by a gaussian kernel. For non-closed lines, the
	 * initial and final points are preserved.
	 * 
	 * @param input      The input line
	 * @param sigma     The standard deviation of the gaussian kernel. The larger,
	 *                   the more smoothed.
	 * @param resolution The target resolution of the geometry. This parameter is
	 *                   used to filter/simplify the final geometry.
	 */
	public static LineString get(LineString input, double sigma, double resolution) {
		if (input.getCoordinates().length <= 2) {
			return (LineString) input.copy();
		}

		boolean isClosed = input.isClosed();
		input = (LineString) input.norm();
		double length = input.getLength();
		double Δx = sigma / 3;
		if (Δx > length * 0.25) {
			// handle degenerate large‐sigma case…
		}

		// 1) Densify
		LineString dens = (LineString) LittleThumblingDensifier.densify(input, Δx);
		Coordinate[] coords = dens.getCoordinates();
		int N = coords.length - 1;

		// 2) Build/normalize kernel
		int n = (int) Math.ceil(7 * sigma / Δx);
		double d = Δx * Δx, b = 2 * sigma * sigma;
		double[] g = new double[n + 1];
		for (int i = 0; i <= n; i++) {
			g[i] = Math.exp(-i * (double) i * d / b);
		}
		// scale so ∑ g[i]*Δx = 1
		for (int i = 0; i <= n; i++)
			g[i] *= Δx;
		double sum = g[0];
		for (int i = 1; i <= n; i++)
			sum += 2 * g[i];
		for (int i = 0; i <= n; i++)
			g[i] /= sum;

		// 3) Pad coords on each side
		Coordinate[] pad = new Coordinate[N + 1 + 2 * n];
		// fill middle
		System.arraycopy(coords, 0, pad, n, N + 1);
		// fill ends
		if (isClosed) {
			for (int i = 0; i < n; i++) {
				pad[i] = coords[N - n + 1 + i]; // wrap tail
				pad[N + 1 + n + i] = coords[(i)]; // wrap head
			}
		} else {
			// reflect or clamp
			for (int i = 0; i < n; i++) {
				pad[i] = coords[n - i]; // reflect at start
				pad[N + 1 + n + i] = coords[N - 1 - i]; // reflect at end
			}
		}

		// 4) Convolve
		Coordinate[] out = new Coordinate[N + 1];
		for (int i = 0; i <= N; i++) {
			double x = 0, y = 0;
			for (int j = -n; j <= n; j++) {
				Coordinate c = pad[i + n + j];
				double w = g[Math.abs(j)];
				x += c.x * w;
				y += c.y * w;
			}
			out[i] = new Coordinate(x, y);
		}
		// preserve closure or endpoints
		if (isClosed) {
			out[N] = out[0];
		} else {
			out[0] = coords[0];
			out[N] = coords[N];
		}

		LineString smooth = input.getFactory().createLineString(out);

		return smooth;
	}

}