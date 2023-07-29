package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import processing.core.PVector;

/**
 * An optimised implementation of Andrew's monotone chain algorithm for
 * constructing convex hulls.
 * <p>
 * Rougly 5x faster than a naive implementation of the algorithm.
 * 
 * @author Jernej Puc
 * @author Michael Carleton
 * @see <a href="https://github.com/JernejPuc/convex-hull/">Github Benchmark</a>
 *
 */
public class FastConvexHull {

	private FastConvexHull() {
	}

	/*
	 * A Java port of Jernej Puc's "optimised version" (Python) of A. M. Andrew's
	 * monotone chain algorithm for constructing convex hulls. Unlike the Python
	 * version, this version handles duplicate coordinates (that are otherwise a
	 * degenerate case).
	 */

	public static List<PVector> convexHull(List<PVector> P) {
		// Preprocess
		P.sort(new PVectorComparator());

		// Data structures
		List<PVector> upper = new ArrayList<>();
		List<PVector> lower = new ArrayList<>();

		// Endpoints
		PVector p0 = P.get(0);
		PVector pn = P.get(P.size() - 1);

		// Initial line(s) of separation
		double ku = (pn.y - p0.y) / (pn.x - p0.x + 1e-12);
		double nu = p0.y - ku * p0.x;
		double kl = ku;
		double nl = nu;

		// Add left endpoint
		upper.add(p0);
		lower.add(p0);

		// Construct the middle of the upper and lower hull sections
		for (int j = 1; j < P.size() - 1; j++) {
			PVector p = P.get(j);

			if (p.y > ku * p.x + nu) {
				while (upper.size() > 1 && isLeftTurn(upper.get(upper.size() - 2), upper.get(upper.size() - 1), p)) {
					upper.remove(upper.size() - 1);
				}
				upper.add(p);

				// Update upper line of separation
				ku = (pn.y - p.y) / (pn.x - p.x + 1e-12);
				nu = p.y - ku * p.x;

			} else if (p.y < kl * p.x + nl) {
				while (lower.size() > 1 && !isLeftTurn(lower.get(lower.size() - 2), lower.get(lower.size() - 1), p)) {
					lower.remove(lower.size() - 1);
				}

				lower.add(p);

				// Update lower line of separation
				kl = (pn.y - p.y) / (pn.x - p.x + 1e-12);
				nl = p.y - kl * p.x;
			}
		}

		// Add right endpoint (only once due to merging that follows)
		while (upper.size() > 1 && isLeftTurn(upper.get(upper.size() - 2), upper.get(upper.size() - 1), pn)) {
			upper.remove(upper.size() - 1);
		}

		while (lower.size() > 1 && !isLeftTurn(lower.get(lower.size() - 2), lower.get(lower.size() - 1), pn)) {
			lower.remove(lower.size() - 1);
		}

		upper.add(pn);

		// Reverse lower hull section
		Collections.reverse(lower);

		// Merge hull sections
		upper.addAll(lower);
		return upper;
	}

	private static class PVectorComparator implements Comparator<PVector> {
		@Override
		public int compare(final PVector p1, final PVector p2) {
			final int xComparison = Float.compare(p1.x, p2.x);
			if (xComparison == 0) {
				return Float.compare(p1.y, p2.y);
			}
			return xComparison;
		}
	}

	/**
	 * Returns True if the turn of pq -> qr is counter-clockwise and False
	 * otherwise.
	 */
	private static boolean isLeftTurn(final PVector p, final PVector q, final PVector r) {
		if (q.x == r.x && q.y == r.y) { // q==r
			return true; // handle duplicate point (degenerate case) -- return true to remove it
		}
		return (q.x - p.x) * (r.y - p.y) > (r.x - p.x) * (q.y - p.y);
	}

}
