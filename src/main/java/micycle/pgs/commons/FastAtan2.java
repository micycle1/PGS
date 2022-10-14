package micycle.pgs.commons;

/**
 * Fast, high-quality polynomial-based atan2 approximation.
 * 
 * @author Michael Carleton
 *
 */
public final class FastAtan2 {

	private FastAtan2() {
	}

	private static final double PI = Math.PI;
	private static final float PI_F = (float) Math.PI;
	private static final double HALF_PI = (0.5 * Math.PI);
	private static final float HALF_PI_F = (float) (0.5 * Math.PI);
	private static final double THREE_QRTR_PI = (0.75 * Math.PI);
	private static final double QRTR_PI = (0.25 * Math.PI);
	private static final float QRTR_PI_F = (float) (0.25 * Math.PI);
	private static final double A = 0.0776509570923569;
	private static final double B = -0.287434475393028;
	private static final double C = (QRTR_PI - A - B);

	/**
	 * Maximum absolute error of ~0.00085 rad (~0.049º).
	 * 
	 * @return Angle from x axis positive side to (x,y) position, in radians, in
	 *         [-PI,PI].
	 */
	public static double atan2(final double y, final double x) {
		// Adapted from https://www.dsprelated.com/showarticle/1052.php
		final double ay = Math.abs(y), ax = Math.abs(x);
		final boolean invert = ay > ax;
		final double z = invert ? ax / ay : ay / ax; // [0,1]
		double th = atanFast2(z); // [0,π/4]
		if (invert) {
			th = HALF_PI - th; // [0,π/2]
		}
		if (x < 0) {
			th = PI - th; // [0,π]
		}
		return Math.copySign(th, y); // [-π,π]
	}

	/**
	 * Maximum absolute error of 0.0015 rad (0.086º).
	 * 
	 * @return Angle from x axis positive side to (x,y) position, in radians, in
	 *         [-PI,PI].
	 */
	public static float atan2(final float y, final float x) {
		// Adapted from https://www.dsprelated.com/showarticle/1052.php
		final float ay = Math.abs(y), ax = Math.abs(x);
		final boolean invert = ay > ax;
		final float z = invert ? ax / ay : ay / ax; // [0,1]
		float th = atanFast(z); // [0,π/4]
		if (invert) {
			th = HALF_PI_F - th; // [0,π/2]
		}
		if (x < 0) {
			th = PI_F - th; // [0,π]
		}
		return Math.copySign(th, y); // [-π,π]
	}

	/**
	 * Maximum absolute error of 0.01 rad (~0.57º). ~8% faster than
	 * {@link #atan2Fast(double, double)}.
	 */
	private static double atan2b(final double y, final double x) {
		// http://dspguru.com/dsp/tricks/fixed-point-atan2-with-self-normalization/
		double r, angle;
		final double abs_y = Math.abs(y) + 1e-10f; // kludge to prevent 0/0 condition

		if (x < 0.0f) {
			r = (x + abs_y) / (abs_y - x); // (3)
			angle = THREE_QRTR_PI; // (4)
		} else {
			r = (x - abs_y) / (x + abs_y); // (1)
			angle = QRTR_PI; // (2)
		}
		angle += (0.1963f * r * r - 0.9817f) * r; // (2 | 4)
		if (y < 0.0f) {
			return (-angle); // negate if in quad III or IV
		} else {
			return (angle);
		}
	}

	/**
	 * Maximum absolute error of 0.0015 rad (~0.086º).
	 */
	private static double atanFast(final double x) {
		// from 'Efficient Approximations for the Arctangent Function' by Rajan et al.
		return x * QRTR_PI - x * ((x - 1) * (0.2447 + 0.0663 * x));
	}

	/**
	 * Maximum absolute error of 0.0015 rad (~0.086º).
	 */
	private static float atanFast(final float x) {
		// from 'Efficient Approximations for the Arctangent Function' by Rajan et al.
		return x * QRTR_PI_F - x * ((x - 1) * (0.2447f + 0.0663f * x));
	}

	/**
	 * Maximum absolute error of ~0.00085 rad (~0.049º). Tested against
	 * {@link #atanFast(double)} here: https://stackoverflow.com/a/42542593/. Max
	 * relative error of ~0.48%.
	 */
	private static double atanFast2(final double x) {
		final double xx = x * x;
		return ((A * xx + B) * xx + C) * x;
	}

	/**
	 * Maximum absolute error of 0.0038 rad (0.22º). Provides better accuracy for
	 * the subintervals 1 > x > 0.5 and −1 < x < −0.5, where the maximum error is
	 * only about 0.001 rad (0.057º).
	 */
	private static double atanFast3(final double z) {
		// from 'Efficient Approximations for the Arctangent Function' by Rajan et al.
		return z * (QRTR_PI + 0.273 * (1 - z)); // NOTE removed Math.abs(); range now [0,1]
	}

}
