package micycle.pgs.commons;

import javax.vecmath.Point3d;
import javax.vecmath.Point4d;

/**
 * Computes 2D/3D weighted geometric median.
 *
 * @author Michael Carleton
 *
 */
public final class GeometricMedian {

	// java port of
	// https://github.com/postgis/postgis/blob/master/liblwgeom/lwgeom_median.c

	private static double DBL_EPSILON = 1E-11;

	/**
	 * Computes the median point of the input point set.
	 * 
	 * @param points   array of x,y,z,w where w is the weight; weights must be
	 *                 non-negative
	 * @param tol      tolerance
	 * @param max_iter max iterations
	 * @return median point of input
	 */
	public static Point3d median(Point4d[] points, double tol, int max_iter) {
		/*
		 * We need to count this ourselves so we can exclude empties and weightless
		 * points.
		 */
		int npoints = points.length;

		Point3d median = init_guess(points, npoints);

		iterate_4d(median, points, npoints, max_iter, tol);

		return median;
	}

	private static Point3d init_guess(Point4d[] points, int npoints) {
		Point3d guess = new Point3d();
		double mass = 0;
		int i;
		for (i = 0; i < npoints; i++) {
			guess.x += points[i].x * points[i].w;
			guess.y += points[i].y * points[i].w;
			guess.z += points[i].z * points[i].w;
			mass += points[i].w;
		}
		guess.x /= mass;
		guess.y /= mass;
		guess.z /= mass;
		return guess;
	}

	private static int iterate_4d(Point3d curr, final Point4d[] points, final int npoints, final int max_iter, final double tol) {
		int i, iter;
		double delta;
		double sum_curr = 0, sum_next = 0;
		boolean hit = false;
		double[] distances = new double[npoints];

		sum_curr = calc_weighted_distances_3d(curr, points, npoints, distances);

		for (iter = 0; iter < max_iter; iter++) {
			Point3d next = new Point3d();
			double denom = 0;

			/* Calculate denom to get the next point */
			for (i = 0; i < npoints; i++) {
				/*
				 * we need to use lower epsilon than in FP_IS_ZERO in the loop for calculation
				 * to converge
				 */
				if (distances[i] > DBL_EPSILON) {
					next.x += points[i].x / distances[i];
					next.y += points[i].y / distances[i];
					next.z += points[i].z / distances[i];
					denom += 1.0 / distances[i];
				} else {
					hit = true;
				}
			}

			if (denom < DBL_EPSILON) {
				/* No movement - Final point */
				break;
			}

			/* Calculate the new point */
			next.x /= denom;
			next.y /= denom;
			next.z /= denom;

			/*
			 * If any of the intermediate points in the calculation is found in the set of
			 * input points, the standard Weiszfeld method gets stuck with a divide-by-zero.
			 *
			 * To get ourselves out of the hole, we follow an alternate procedure to get the
			 * next iteration, as described in:
			 *
			 * Vardi, Y. and Zhang, C. (2011) "A modified Weiszfeld algorithm for the
			 * Fermat-Weber location problem." Math. Program., Ser. A 90: 559-566. DOI
			 * 10.1007/s101070100222
			 *
			 * Available online at the time of this writing at
			 * http://www.stat.rutgers.edu/home/cunhui/papers/43.pdf
			 */
			if (hit) {
				double dx = 0, dy = 0, dz = 0;
				double d_sqr;
				hit = false;

				for (i = 0; i < npoints; i++) {
					if (distances[i] > DBL_EPSILON) {
						dx += (points[i].x - curr.x) / distances[i];
						dy += (points[i].y - curr.y) / distances[i];
						dz += (points[i].z - curr.z) / distances[i];
					}
				}

				d_sqr = Math.sqrt(dx * dx + dy * dy + dz * dz);
				if (d_sqr > DBL_EPSILON) {
					double r_inv = Math.max(0, 1.0 / d_sqr); // note
					next.x = (1.0 - r_inv) * next.x + r_inv * curr.x;
					next.y = (1.0 - r_inv) * next.y + r_inv * curr.y;
					next.z = (1.0 - r_inv) * next.z + r_inv * curr.z;
				}
			}

			/* Check movement with next point */
			sum_next = calc_weighted_distances_3d(next, points, npoints, distances);
			delta = sum_curr - sum_next;
			if (delta < tol) {
				break;
			} else {
				curr.x = next.x;
				curr.y = next.y;
				curr.z = next.z;
				sum_curr = sum_next;
			}
		}

		return iter;
	}

	private static double calc_weighted_distances_3d(final Point3d curr, final Point4d[] points, int npoints,
			double[] distances) {
		int i;
		double weight = 0.0;
		for (i = 0; i < npoints; i++) {
			double dist = distance3d_pt_pt(curr, points[i]);
			distances[i] = dist / points[i].w;
			weight += dist * points[i].w;
		}

		return weight;
	}

	private static double distance3d_pt_pt(Point3d p0, Point4d p1) {
		double dx, dy, dz;

		dx = p0.x - p1.x;
		dy = p0.y - p1.y;
		dz = p0.z - p1.z;
		return Math.sqrt(dx * dx + dy * dy + dz * dz);
	}

}
