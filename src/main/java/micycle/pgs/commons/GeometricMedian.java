package micycle.pgs.commons;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateXYZM;

/**
 * Computes 2D/3D weighted geometric median.
 *
 * @author Michael Carleton
 *
 */
public final class GeometricMedian {

	// java port of
	// https://github.com/postgis/postgis/blob/master/liblwgeom/lwgeom_median.c

	private static final double DBL_EPSILON = 1E-11;

	private GeometricMedian() {
	}

	/**
	 * Computes the median point of the input point set.
	 * 
	 * @param points  array of x,y,z,w where w is the weight; weights must be
	 *                non-negative
	 * @param tol     tolerance
	 * @param maxIter max iterations
	 * @return median point of input
	 */
	public static Coordinate median(CoordinateXYZM[] points, double tol, int maxIter) {
		/*
		 * We need to count this ourselves so we can exclude empties and weightless
		 * points.
		 */
		int npoints = points.length;

		Coordinate median = initGuess(points, npoints);

		iterate4d(median, points, npoints, maxIter, tol);

		return median;
	}

	private static Coordinate initGuess(CoordinateXYZM[] points, int npoints) {
		Coordinate guess = new Coordinate();
		double mass = 0;
		int i;
		for (i = 0; i < npoints; i++) {
			final double weight = points[i].getM();
			guess.x += points[i].x * weight;
			guess.y += points[i].y * weight;
			guess.setZ(guess.getZ() + points[i].getZ() * weight);
			mass += weight;
		}
		guess.x /= mass;
		guess.y /= mass;
		guess.setZ(guess.getZ() / mass);
		return guess;
	}

	private static int iterate4d(Coordinate curr, final CoordinateXYZM[] points, final int npoints, final int maxIter, final double tol) {
		int i, iter;
		double delta;
		double sumCurr = 0, sumNext = 0;
		boolean hit = false;
		double[] distances = new double[npoints];

		sumCurr = calcWeightedDistances3d(curr, points, npoints, distances);

		for (iter = 0; iter < maxIter; iter++) {
			Coordinate next = new Coordinate();
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
					next.setZ(next.getZ() + points[i].getZ() / distances[i]);
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
			next.setZ(next.getZ() / denom);

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
				double dSqr;
				hit = false;

				for (i = 0; i < npoints; i++) {
					if (distances[i] > DBL_EPSILON) {
						dx += (points[i].x - curr.x) / distances[i];
						dy += (points[i].y - curr.y) / distances[i];
						dz += (points[i].getZ() - curr.getZ()) / distances[i];
					}
				}

				dSqr = Math.sqrt(dx * dx + dy * dy + dz * dz);
				if (dSqr > DBL_EPSILON) {
					double rInv = Math.max(0, 1.0 / dSqr); // note
					next.x = (1.0 - rInv) * next.x + rInv * curr.x;
					next.y = (1.0 - rInv) * next.y + rInv * curr.y;
					next.setZ((1.0 - rInv) * next.getZ() + rInv * curr.getZ());
				}
			}

			/* Check movement with next point */
			sumNext = calcWeightedDistances3d(next, points, npoints, distances);
			delta = sumCurr - sumNext;
			if (delta < tol) {
				break;
			} else {
				curr.x = next.x;
				curr.y = next.y;
				curr.setZ(next.getZ());
				sumCurr = sumNext;
			}
		}

		return iter;
	}

	private static double calcWeightedDistances3d(final Coordinate curr, final CoordinateXYZM[] points, int npoints, double[] distances) {
		int i;
		double weight = 0.0;
		for (i = 0; i < npoints; i++) {
			double dist = distance3dPtPt(curr, points[i]);
			distances[i] = dist / points[i].getM();
			weight += dist * points[i].getM();
		}

		return weight;
	}

	private static double distance3dPtPt(Coordinate p0, CoordinateXYZM p1) {
		double dx, dy, dz;

		dx = p0.x - p1.x;
		dy = p0.y - p1.y;
		dz = p0.getZ() - p1.getZ();
		return Math.sqrt(dx * dx + dy * dy + dz * dz);
	}

}
