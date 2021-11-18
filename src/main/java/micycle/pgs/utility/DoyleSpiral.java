package micycle.pgs.utility;

import static java.lang.Math.pow;
import static java.lang.Math.sin;

import java.util.HashSet;
import java.util.Set;

import processing.core.PVector;

import static java.lang.Math.cos;
import static java.lang.Math.PI;

/**
 * Doyle spirals are patterns of non-crossing circles in the plane, each tangent
 * to six others.
 * 
 * @author Javascript implementation by Robin Houston
 * @author Java port by Michael Carleton
 *
 */
@java.lang.SuppressWarnings({"java:S2184", "java:S1905", "java:S100"})
public class DoyleSpiral {

	// https://bl.ocks.org/robinhouston/6096950
	// https://github.com/robinhouston/doyle-spirals

	/*
	 * NOTE In methods that calculate p / q, leaving this as integer division
	 * (versus double division), results in no missing circles around the center (as
	 * is the case with double division on some p & q) but some circles are
	 * duplicated (hence a set is returned).
	 */

	private double maxRadius = 500;
	private Set<PVector> circles;

	/**
	 * Creates a Doyle spiral having the given parameters. Parameters p and q
	 * control the number of arms the spiral has in each direction.
	 * 
	 * @param p         at least 2
	 * @param q         at least p + 1
	 * @param maxRadius the maximum radius of the packing arrangement (the maximum
	 *                  distance a circle centroid can be from the center of the
	 *                  arrangement)
	 */
	public DoyleSpiral(int p, int q, double maxRadius) {
		p = Math.max(p, 2);
		q = Math.max(q, p + 1);
		this.maxRadius = maxRadius;
		circles = new HashSet<>();

		Doyle root = new Doyle(p, q);
		double[] start = root.a;

		for (int i = 0; i < q; i++) {
			spiral(root.r, start, root.a);
			start = cmul(start, root.b);
		}
	}

	/**
	 * Centered on (0, 0).
	 * 
	 * @return
	 */
	public Set<PVector> getCircles() {
		return circles;
	}

	private void spiral(double r, double[] start_point, double[] delta) {
		final double mod_delta = modulus(delta);

		double[] q = start_point;
		for (double mod_q = modulus(q); mod_q < maxRadius; mod_q *= mod_delta) {
			circles.add(new PVector((float) q[0], (float) q[1], (float) (mod_q * r)));
			q = cmul(q, delta);
		}
	}

	private static double[] cmul(double[] w, double[] z) {
		return new double[] { w[0] * z[0] - w[1] * z[1], w[0] * z[1] + w[1] * z[0] };
	}

	private static double modulus(double[] p) {
		return Math.sqrt(p[0] * p[0] + p[1] * p[1]);
	}

	private static double _d(double z, double t, int p, int q) {
		// The square of the distance between z*e^(it) and z*e^(it)^(p/q).
		double w = pow(z, p / (int) q), s = (p * t + 2 * PI) / q;
		return (pow(z * cos(t) - w * cos(s), 2) + pow(z * sin(t) - w * sin(s), 2));
	}

	private static double ddz_d(double z, double t, int p, int q) {
		// The partial derivative of _d with respect to z.
		double w = pow(z, p / (int) q), s = (p * t + 2 * PI) / q, ddz_w = (p / (int) q) * pow(z, (p - q) / (int) q);
		return (2 * (w * cos(s) - z * cos(t)) * (ddz_w * cos(s) - cos(t)) + 2 * (w * sin(s) - z * sin(t)) * (ddz_w * sin(s) - sin(t)));
	}

	private static double ddt_d(double z, double t, int p, int q) {
		// The partial derivative of _d with respect to t.
		double w = pow(z, p / (int) q), s = (p * t + 2 * PI) / q, dds_t = (p / (int) q);
		return (2 * (z * cos(t) - w * cos(s)) * (-z * sin(t) + w * sin(s) * dds_t)
				+ 2 * (z * sin(t) - w * sin(s)) * (z * cos(t) - w * cos(s) * dds_t));
	}

	private static double _s(double z, double t, int p, int q) {
		// The square of the sum of the origin-distance of z*e^(it) and
		// the origin-distance of z*e^(it)^(p/q).
		return pow(z + pow(z, p / (int) q), 2);
	}

	private static double ddz_s(double z, double t, int p, int q) {
		// The partial derivative of _s with respect to z.
		double w = pow(z, p / (int) q), ddz_w = (p / (int) q) * pow(z, (p - q) / (int) q);
		return 2 * (w + z) * (ddz_w + 1);
	}

	private static double _r(double z, double t, int p, int q) {
		// The square of the radius-ratio implied by having touching circles
		// centred at z*e^(it) and z*e^(it)^(p/q).
		return _d(z, t, p, q) / _s(z, t, p, q);
	}

	private static double ddz_r(double z, double t, int p, int q) {
		// The partial derivative of _r with respect to z.
		return (ddz_d(z, t, p, q) * _s(z, t, p, q) - _d(z, t, p, q) * ddz_s(z, t, p, q)) / pow(_s(z, t, p, q), 2);
	}

	private static double ddt_r(double z, double t, int p, int q) {
		// The partial derivative of _r with respect to t.
		return (ddt_d(z, t, p, q) * _s(z, t, p, q)
		/* - _d(z,t,p,q) * ddt_s(z,t,p,q) */ // omitted because ddt_s is constant at zero
		) / pow(_s(z, t, p, q), 2);
	}

	private static class Doyle {

		private static double epsilon = 1e-10;

		int p, q;
		double[] a, b;

		private Doyle(int p, int q) {
			this.p = p;
			this.q = q;
			find_root(2, 0);
			if (!ok) {
				System.err.println("Failed to find root for p=" + p + ", q=" + q);
			}
			a = new double[] { z * cos(t), z * sin(t) };
			double corootZ = pow(z, p / (int) q);
			double corootT = (p * t + 2 * PI) / q;
			b = new double[] { corootZ * cos(corootT), corootZ * sin(corootT) };
		}

		boolean ok = false;
		double z, t, r;

		private void find_root(double z, double t) {
			while (true) {
				double v_f = _f(z, t), v_g = _g(z, t);
				if (-epsilon < v_f && v_f < epsilon && -epsilon < v_g && v_g < epsilon) {
					ok = true;
					this.z = z;
					this.t = t;
					this.r = Math.sqrt(_r(z, t, 0, 1));
					return;
				}

				double a = ddz_f(z, t), b = ddt_f(z, t), c = ddz_g(z, t), d = ddt_g(z, t);
				double det = a * d - b * c;
				if (-epsilon < det && det < epsilon) {
					ok = false;
					return;
				}

				z -= (d * v_f - b * v_g) / det;
				t -= (a * v_g - c * v_f) / det;

				if (z < epsilon) {
					ok = false;
					return;
				}
			}
		}

		private double _f(double z, double t) {
			return _r(z, t, 0, 1) - _r(z, t, p, q);
		}

		private double ddz_f(double z, double t) {
			return ddz_r(z, t, 0, 1) - ddz_r(z, t, p, q);
		}

		private double ddt_f(double z, double t) {
			return ddt_r(z, t, 0, 1) - ddt_r(z, t, p, q);
		}

		private double _g(double z, double t) {
			return _r(z, t, 0, 1) - _r(pow(z, p / (int) q), (p * t + 2 * PI) / q, 0, 1);
		}

		private double ddz_g(double z, double t) {
			return ddz_r(z, t, 0, 1) - ddz_r(pow(z, p / (int) q), (p * t + 2 * PI) / q, 0, 1) * (p / (int) q) * pow(z, (p - q) / (int) q);
		}

		private double ddt_g(double z, double t) {
			return ddt_r(z, t, 0, 1) - ddt_r(pow(z, p / (int) q), (p * t + 2 * PI) / q, 0, 1) * (p / (int) q);
		}
	}

}
