package micycle.pgs.commons;

import java.util.Arrays;

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
	 * @param line       The input line
	 * @param sigmaM     The standard deviation of the gaussian kernel. The larger,
	 *                   the more smoothed.
	 * @param resolution The target resolution of the geometry. This parameter is
	 *                   used to filter/simplify the final geometry.
	 */
	public static LineString get(LineString line, double sigmaM, double resolution) {
		if (line.getCoordinates().length <= 2) {
			return (LineString) line.copy();
		}

		// output is sensitive to vertex order, so normalise the line (note doesn't
		// normalise orientation)
		boolean isClosed = line.isClosed();
		if (isClosed) {
			line = normalise(line);
		}
		double length = line.getLength();
		double densifiedResolution = sigmaM / 3;

		// handle extreme cases: too large sigma resulting in too large densified
		// resolution.
		if (densifiedResolution > 0.25 * length) {
			if (isClosed) {
				// return tiny triangle nearby center point
				Coordinate c = line.getCentroid().getCoordinate();
				length *= 0.01;
				return line.getFactory()
						.createLineString(new Coordinate[] { new Coordinate(c.x - length, c.y - length), new Coordinate(c.x, c.y + length),
								new Coordinate(c.x + length, c.y - length), new Coordinate(c.x - length, c.y - length) });
			} else {
				// return segment
				return line.getFactory()
						.createLineString(new Coordinate[] { line.getCoordinateN(0), line.getCoordinateN(line.getNumPoints() - 1) });
			}
		}

		// compute densified line
		Coordinate[] densifiedCoords = LittleThumblingDensifier.densify(line, densifiedResolution).getCoordinates();

		// build output line structure
		int nb = (int) (length / densifiedResolution);
		Coordinate[] out = new Coordinate[nb + 1];

		// prepare gaussian coefficients
		int n = 7 * 3; // it should be: E(7*sigma/densifiedResolution), which is 7*3;
		double[] gcs = new double[n + 1];
		final double a = sigmaM * Math.sqrt(2 * Math.PI);
		final double b = sigmaM * sigmaM * 2;
		final double d = densifiedResolution * densifiedResolution;
		for (int i = 0; i < n + 1; i++) {
			gcs[i] = FastMath.exp(-i * i * d / b) / a;
		}

		final Coordinate c0 = densifiedCoords[0];
		final Coordinate cN = densifiedCoords[nb];
		for (int i = 0; i < nb; i++) {
			if (!isClosed && i == 0) {
				continue;
			}

			// compute coordinates of point i of the smoothed line (gauss mean)
			double x = 0.0, y = 0.0;
			for (int j = -n; j <= n; j++) {
				// index of the point to consider on the original densified line
				int q = i + j;
				// find coordinates (xq,yq) of point q
				double xq, yq;
				if (q < 0) {
					if (isClosed) {
						// make loop to get the right point
						q = q % nb;
						if (q < 0) {
							q += nb;
						}
						Coordinate c = densifiedCoords[q];
						xq = c.x;
						yq = c.y;
					} else {
						// get symetric point
						q = (-q) % nb;
						if (q == 0) {
							q = nb;
						}
						Coordinate c = densifiedCoords[q];
						xq = 2 * c0.x - c.x;
						yq = 2 * c0.y - c.y;
					}
				} else if (q > nb) {
					if (isClosed) {
						// make loop to get the right point
						q = q % nb;
						if (q == 0) {
							q = nb;
						}
						Coordinate c = densifiedCoords[q];
						xq = c.x;
						yq = c.y;
					} else {
						// get symetric point
						q = nb - q % nb;
						if (q == nb) {
							q = 0;
						}
						Coordinate c = densifiedCoords[q];
						xq = 2 * cN.x - c.x;
						yq = 2 * cN.y - c.y;
					}
				} else {
					// general case (most frequent)
					Coordinate c = densifiedCoords[q];
					xq = c.x;
					yq = c.y;
				}
				// get gaussian coefficient
				double gc = gcs[j >= 0 ? j : -j];
				// add contribution of point q to new position of point i
				x += xq * gc;
				y += yq * gc;
			}
			// assign smoothed position of point i
			out[i] = new Coordinate(x * densifiedResolution, y * densifiedResolution);
		}

		// handle start and end points
		if (isClosed) {
			// ensure start and end locations are the same
			out[nb] = out[0];
		} else {
			// ensure start and end points are at the same position as the initial geometry
			out[0] = densifiedCoords[0];
			out[nb] = densifiedCoords[densifiedCoords.length - 1];
		}

		LineString lsOut = line.getFactory().createLineString(out);
		return lsOut;
	}

	/**
	 * Normalises the LineString so that it starts from the coordinate with the
	 * smallest x value. In case of a tie on the x value, the smallest y value is
	 * used.
	 *
	 * @param line The open or closed LineString to be normalised.
	 * @return A new LineString with coordinates ordered starting from the smallest
	 *         x (and y, if tied).
	 */
	private static LineString normalise(LineString line) {
		boolean isClosed = line.isClosed();

		Coordinate[] originalCoords = line.getCoordinates();
		if (isClosed && originalCoords[0].equals2D(originalCoords[originalCoords.length - 1])) {
			// Remove last vertex if it is a duplicate of the first for closed lines
			originalCoords = Arrays.copyOf(originalCoords, originalCoords.length - 1);
		}

		// Find index of coordinate with smallest x value, tie by y value
		int minIndex = 0;
		for (int i = 1; i < originalCoords.length; i++) {
			if (originalCoords[i].x < originalCoords[minIndex].x
					|| (originalCoords[i].x == originalCoords[minIndex].x && originalCoords[i].y < originalCoords[minIndex].y)) {
				minIndex = i;
			}
		}

		// Rotate array to start from vertex with smallest x (and y, if tied)
		Coordinate[] rotatedCoords = new Coordinate[originalCoords.length + (isClosed ? 1 : 0)];
		System.arraycopy(originalCoords, minIndex, rotatedCoords, 0, originalCoords.length - minIndex);
		System.arraycopy(originalCoords, 0, rotatedCoords, originalCoords.length - minIndex, minIndex);

		if (isClosed) {
			rotatedCoords[rotatedCoords.length - 1] = rotatedCoords[0]; // Close the loop
		}

		return line.getFactory().createLineString(rotatedCoords);
	}

}