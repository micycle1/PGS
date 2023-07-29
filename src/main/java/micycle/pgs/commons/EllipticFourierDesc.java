package micycle.pgs.commons;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LinearRing;

import net.jafama.FastMath;

/**
 * Implements Elliptic Fourier Descriptors (EFD).
 * <p>
 * The EFD provides a normalized set of coefficients that are rotation,
 * translation and scale invariant. The very first coefficient relates to the
 * centroid of the input shape before the EFD is computed and can be ignored;
 * the second FD coefficient relates to a circle circumscribed about the
 * centroid before the EFD computation. After the EFD computation the second EFD
 * is always 2 and can be ignored. That leaves the remaining EFD coefficients
 * for use in comparing shapes.
 */
public class EllipticFourierDesc {

	// https://github.com/droiddeveloper1/android-wear-gestures-recognition/blob/master/mobile/src/main/java/com/dmt/gestureproject_1/

	private Coordinate[] xy; // the x and y coordinates
	/** The number of points on input contour */
	private int m;
	/** The number of FD coefficients */
	private int nFD;
	/** The Fourier Descriptors */
	private double[] efd_ax, efd_ay, efd_bx, efd_by;
	/** The normalized Elliptic Fourier Descriptors */
	private double[] efd;

	/**
	 * Constructs a descriptor object for a given polygon and computes the
	 * Elliptical Fourier Descriptors for it.
	 * 
	 * @param coords A two-dimensional array representing the vertices of the
	 *               polygon. Each row contains the x and y coordinates of a vertex.
	 * @param n      The number of Fourier descriptors to compute. This is also the
	 *               number of harmonics used in the Fourier series. >=2.
	 */
	public EllipticFourierDesc(LinearRing ring, int n) {
		this.xy = ring.getCoordinates();
		this.nFD = n;
		this.m = xy.length;
		computeEllipticFD();
	}

	/**
	 * Constructs a descriptor object for a given polygon and computes
	 * <code>#vertices/2</code> Elliptical Fourier Descriptors for it.
	 * 
	 * @param x the x coordinates of the contour
	 * @param y the y coordinates of the contour
	 */
	public EllipticFourierDesc(LinearRing ring) {
		this(ring, ring.getNumPoints() / 2);
	}

	/**
	 * Returns the elliptic fourier descriptors, which are computed upon
	 * initialisation.
	 * <p>
	 * The very first coefficient relates to the centroid of the input shape before
	 * the EFD is computed and can be ignored; the second FD coefficient relates to
	 * a circle circumscribed about the centroid before the EFD computation -- after
	 * the EFD computation this second EFD is always 2 and can be ignored.
	 * 
	 * @return list of coefficients/harmonics
	 */
	public double[] getEFD() {
		return efd;
	}

	/**
	 * Computes the Fourier and Elliptic Fourier Descriptors
	 */
	private void computeEllipticFD() {

		// the fourier descriptors
		efd_ax = new double[nFD];
		efd_ay = new double[nFD];
		efd_bx = new double[nFD];
		efd_by = new double[nFD];

		// preconfigure some values
		final double t = 2.0 * Math.PI / m;
		final double twoOverM = 2.0 / m;
		// step through each FD
		for (int k = 0; k < nFD; k++) {
			// and for each point
			for (int i = 0; i < m; i++) {
				final double p = k * t * i;
				final double sin = FastMath.sin(p);
				final double cos = FastMath.cos(p);
				efd_ax[k] += xy[i].x * cos;
				efd_bx[k] += xy[i].x * sin;
				efd_ay[k] += xy[i].y * cos;
				efd_by[k] += xy[i].y * sin;
			} // i-loop through the number of points

			efd_ax[k] *= twoOverM;
			efd_bx[k] *= twoOverM;
			efd_ay[k] *= twoOverM;
			efd_by[k] *= twoOverM;

		}

		// now compute the elliptic fourier descriptors as per Thomas Boudier's
		// implementation
		efd = new double[nFD];
		final int first = 1; // index of the normalization values
		// precompute the denominators
		final double denomA = (efd_ax[first] * efd_ax[first]) + (efd_ay[first] * efd_ay[first]);
		final double denomB = (efd_bx[first] * efd_bx[first]) + (efd_by[first] * efd_by[first]);
		for (int k = 0; k < nFD; k++) {
			efd[k] = Math.sqrt((efd_ax[k] * efd_ax[k] + efd_ay[k] * efd_ay[k]) / denomA)
					+ Math.sqrt((efd_bx[k] * efd_bx[k] + efd_by[k] * efd_by[k]) / denomB);
		}
	}

	/**
	 * Returns the polygon computed using all of this instance's EFD coefficients.
	 * 
	 * @return an array of (x,y) pairs that is the same length as the input polygon
	 *
	 */
	public Coordinate[] createPolygon() {
		return createPolygon(nFD);
	}

	/**
	 * Creates the polygon corresponding to the nth harmonic.
	 * 
	 * @param n >= 2.
	 * @return an array of coordinates that has the same length as the input polygon
	 *         vertices, and forms a closed ring
	 */
	public Coordinate[] createPolygon(int n) {
		if (m == 0) {
			return new Coordinate[0];
		}
		n = Math.min(n, nFD);
		Coordinate[] new_xy = new Coordinate[m + 1]; // NOTE +1 to close ring
		final double t = 2.0 * Math.PI / m;
		for (int i = 0; i < m; i++) {
			new_xy[i] = new Coordinate(efd_ax[0] / 2.0, efd_ay[0] / 2.0);

			for (int k = 1; k < n; k++) {
				final double p = t * k * i;
				final double sin = FastMath.sin(p);
				final double cos = FastMath.cos(p);
				new_xy[i].x += efd_ax[k] * cos + efd_bx[k] * sin;
				new_xy[i].y += efd_ay[k] * cos + efd_by[k] * sin;
			}
		}

		new_xy[m] = new_xy[0].copy(); // close ring
		return new_xy;
	}

	/**
	 * Computes pairwise euclidean distance between two descriptors (of equal
	 * length).
	 */
	public static double computeEFDDistance(double[] efd1, double[] efd2) {
		if (efd1.length != efd2.length) {
			throw new IllegalArgumentException(String.format("EFDs must be of the same length: %s vs %s", efd1.length, efd2.length));
		}
		double sum = 0.0;
		// ignore first 2 components
		for (int i = 2; i < efd1.length; i++) {
			double diff = efd1[i] - efd2[i];
			sum += diff * diff;
		}
		return Math.sqrt(sum);
	}

}