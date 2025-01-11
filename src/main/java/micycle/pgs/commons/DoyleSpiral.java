package micycle.pgs.commons;

import static java.lang.Math.PI;
import static net.jafama.FastMath.cos;
import static net.jafama.FastMath.pow;
import static net.jafama.FastMath.sin;

import java.util.ArrayList;
import java.util.List;

import processing.core.PVector;

/**
 * Doyle spirals are patterns of non-crossing circles in the plane, each tangent
 * to six others.
 * 
 * @author Javascript implementation by Robin Houston
 * @author Java port by Michael Carleton
 *
 */
@java.lang.SuppressWarnings({ "java:S100" })
public class DoyleSpiral {

	// https://bl.ocks.org/robinhouston/6096950
	// https://github.com/robinhouston/doyle-spirals

	private final double maxRadius;
	private List<PVector> circles;

	public DoyleSpiral(int p, int q, double maxRadius) {
		p = Math.max(p, 2);
		q = Math.max(q, p + 1);
		this.maxRadius = maxRadius;
		circles = new ArrayList<>();

		Doyle root = new Doyle(p, q);

		double[] Q;
		final double mod_b = hypot(root.b);

//		Q = new double[] { 1e-3, 0 }; // min radius of 1e-3
//		for (int i = 0; i < p; i++) {
//			double[] c = Q;
//			double d = hypot(c);
//			while (d < maxRadius) {
//				circles.add(makeCircle(c[0], c[1], d * root.r));
//				c = cmul(c, root.b);
//				d *= mod_b;
//			}
//			Q = cmul(Q, root.a);
//		}
		final float scale = (float) pow(root.z, -10);
		maxRadius = maxRadius / scale;
		Q = new double[] { 1, 0 };
		for (int i = 0; i < p; i++) {
			double[] c = Q;
			double d = hypot(c);
			while (d < maxRadius) {
				PVector circle = makeCircle(c[0], c[1], d * root.r).mult(scale);
				circles.add(circle);
				c = cmul(c, root.b);
				d *= mod_b;
			}
			Q = cmul(Q, root.a);
		}
	}

	/**
	 * Returns a list of circles comprising the spriral. The spiral is centered on
	 * (0, 0).
	 * 
	 * @return list of PVector circles where .z refers to radius of circle
	 */
	public List<PVector> getCircles() {
		return circles;
	}

	void spiral(double r, double[] start_point, double[] delta) {
		final double mod_delta = hypot(delta);

		double[] q = start_point;
		for (double mod_q = hypot(q); mod_q < maxRadius; mod_q *= mod_delta) {
			circles.add(makeCircle(q[0], q[1], mod_q * r));
			q = cmul(q, delta);
		}
	}

	private PVector makeCircle(double x, double y, double r) {
//		double[] mobius = transformCircle(x, y, r);
//		x = mobius[0];
//		y = mobius[1];
//		r = mobius[2];
		return new PVector((float) x, (float) y, (float) r);
	}

	// https://gist.github.com/robinhouston/6098310
	private static double[] transformCircle(double x, double y, double r) {
		// The image of a circle under the Möbius transformation
		double[] a = transformPoint(x - r, y);
		double[] b = transformPoint(x + r, y);
		double[] c = transformPoint(x, y + r);
		return circleThroughPoints(a, b, c);
	}

	private static double[] transformPoint(double x, double y) {
		// Möbius transformation that maps (0, 1, ∞) to (-1, 0, 1)
		double denom = (x + 1) * (x + 1) + y * y;
		return new double[] { (x * x - 1 + y * y) / denom, 2 * y / denom };
	}

	private static double[] circleThroughPoints(double[] a, double[] b, double[] c) {
		// The unique circle passing through three non-collinear points
		double na = a[0] * a[0] + a[1] * a[1];
		double nb = b[0] * b[0] + b[1] * b[1];
		double nc = c[0] * c[0] + c[1] * c[1];
		double y = ((a[0] - b[0]) * (nb - nc) - (b[0] - c[0]) * (na - nb))
				/ (2 * (b[1] - a[1]) * (b[0] - c[0]) - 2 * (a[0] - b[0]) * (c[1] - b[1])),
				x = (na - nb + 2 * (b[1] - a[1]) * y) / (2 * (a[0] - b[0])),
				r = Math.sqrt((x - a[0]) * (x - a[0]) + (y - a[1]) * (y - a[1]));
		return new double[] { x, y, r };
	}

	private static double[] cmul(double[] w, double[] z) {
		return new double[] { w[0] * z[0] - w[1] * z[1], w[0] * z[1] + w[1] * z[0] };
	}

	private static double[] crot(double[] z, double theta) {
		return cmul(z, new double[] { Math.cos(theta), Math.sin(theta) });
	}

	private static double hypot(double[] p) {
		return Math.sqrt(p[0] * p[0] + p[1] * p[1]);
	}

	private static double _d(double z, double t, double p, double q) {
		// The square of the distance between z*e^(it) and z*e^(it)^(p/q).
		double w = pow(z, p / q), s = (p * t + 2 * PI) / q;
		return (pow(z * cos(t) - w * cos(s), 2) + pow(z * sin(t) - w * sin(s), 2));
	}

	private static double ddz_d(double z, double t, double p, double q) {
		// The partial derivative of _d with respect to z.
		double w = pow(z, p / q), s = (p * t + 2 * PI) / q, ddz_w = (p / q) * pow(z, (p - q) / q);
		return (2 * (w * cos(s) - z * cos(t)) * (ddz_w * cos(s) - cos(t)) + 2 * (w * sin(s) - z * sin(t)) * (ddz_w * sin(s) - sin(t)));
	}

	private static double ddt_d(double z, double t, double p, double q) {
		// The partial derivative of _d with respect to t.
		double w = pow(z, p / q), s = (p * t + 2 * PI) / q, dds_t = (p / q);
		return (2 * (z * cos(t) - w * cos(s)) * (-z * sin(t) + w * sin(s) * dds_t)
				+ 2 * (z * sin(t) - w * sin(s)) * (z * cos(t) - w * cos(s) * dds_t));
	}

	private static double _s(double z, double t, double p, double q) {
		// The square of the sum of the origin-distance of z*e^(it) and
		// the origin-distance of z*e^(it)^(p/q).
		return pow(z + pow(z, p / q), 2);
	}

	private static double ddz_s(double z, double t, double p, double q) {
		// The partial derivative of _s with respect to z.
		double w = pow(z, p / q), ddz_w = (p / q) * pow(z, (p - q) / q);
		return 2 * (w + z) * (ddz_w + 1);
	}

	private static double _r(double z, double t, double p, double q) {
		// The square of the radius-ratio implied by having touching circles
		// centred at ze^(it) and ze^(it)^(p/q).
		return _d(z, t, p, q) / _s(z, t, p, q);
	}

	private static double ddz_r(double z, double t, double p, double q) {
		// The partial derivative of _r with respect to z.
		return (ddz_d(z, t, p, q) * _s(z, t, p, q) - _d(z, t, p, q) * ddz_s(z, t, p, q)) / pow(_s(z, t, p, q), 2);
	}

	private static double ddt_r(double z, double t, double p, double q) {
		// The partial derivative of _r with respect to t.
		return (ddt_d(z, t, p, q) * _s(z, t, p, q)
		/* - _d(z,t,p,q) * ddt_s(z,t,p,q) */ // This term is zero.
		) / pow(_s(z, t, p, q), 2);
	}

	private static class Doyle {

		private static final double EPSILON = 1e-10;

		final double p, q;

		boolean ok = false;
		double z;
		double t;
		double r;

		double[] a;
		double[] b;

		public Doyle(double p, double q) {
			this.p = p;
			this.q = q;

			findRoot(2, 0);

			if (!ok) {
				System.err.println("Failed to find root for p=" + p + ", q=" + q);
			}

			a = new double[] { z * cos(t), z * sin(t) };
			double corootZ = pow(z, p / q);
			double corootT = (p * t + 2 * PI) / q;
			b = new double[] { corootZ * cos(corootT), corootZ * sin(corootT) };
		}

		private void findRoot(double z, double t) {
			while (true) {
				double v_f = _f(z, t), v_g = _g(z, t);
				if (-EPSILON < v_f && v_f < EPSILON && -EPSILON < v_g && v_g < EPSILON) {
					this.z = z;
					this.t = t;
					this.r = Math.sqrt(_r(z, t, 0, 1));
					ok = true;
					return;
				}

				double a = ddz_f(z, t), b = ddt_f(z, t), c = ddz_g(z, t), d = ddt_g(z, t);
				double det = a * d - b * c;
				if (-EPSILON < det && det < EPSILON) {
					return;
				}

				z -= (d * v_f - b * v_g) / det;
				t -= (a * v_g - c * v_f) / det;

				if (z < EPSILON) {
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
			return _r(z, t, 0, 1) - _r(pow(z, p / q), (p * t + 2 * PI) / q, 0, 1);
		}

		private double ddz_g(double z, double t) {
			return ddz_r(z, t, 0, 1) - ddz_r(pow(z, p / q), (p * t + 2 * PI) / q, 0, 1) * (p / q) * pow(z, (p - q) / q);
		}

		private double ddt_g(double z, double t) {
			return ddt_r(z, t, 0, 1) - ddt_r(pow(z, p / q), (p * t + 2 * PI) / q, 0, 1) * (p / q);
		}
	}
}
