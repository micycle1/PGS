package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.function.ToDoubleBiFunction;

/**
 * Gonzalez (Gon) heuristic for the k-center problem on an arbitrary vertex
 * collection with a distance function.
 * <p>
 * Chooses the first center uniformly at random, then repeatedly adds the
 * farthest vertex from the current center set (ties -> lowest index in the
 * input iteration order).
 * <p>
 * Time: O(k * |V|) distance evaluations Space: O(|V|)
 */
public final class GonHeuristic<V> {

	private final Random rng;

	public GonHeuristic(Random rng) {
		this.rng = Objects.requireNonNull(rng, "rng");
	}

	public List<V> getCenters(Collection<V> vertices, int k, ToDoubleBiFunction<V, V> distFunc) {
		Objects.requireNonNull(vertices, "vertices");
		Objects.requireNonNull(distFunc, "distFunc");

		final int n = vertices.size();
		if (k <= 0) {
			return new ArrayList<>();
		}
		if (n == 0) {
			throw new IllegalArgumentException("vertices must be non-empty");
		}
		if (n < k) {
			throw new IllegalArgumentException("number of vertices must be at least k");
		}
		if (n == k) {
			return new ArrayList<>();
		}

		// Indexable storage (fast scans, stable order for tie-breaking).
		@SuppressWarnings("unchecked")
		final V[] verts = (V[]) vertices.toArray();

		// Bookkeeping
		final boolean[] isCenter = new boolean[n];
		final double[] minDist = new double[n];
		Arrays.fill(minDist, Double.POSITIVE_INFINITY);

		// Centers by index, in selection order
		final int[] centers = new int[k];
		int centerCount = 0;

		// Pick first center randomly
		final int first = rng.nextInt(n);
		isCenter[first] = true;
		minDist[first] = 0.0;
		centers[centerCount++] = first;

		// Initialize minDist to first center
		final V firstV = verts[first];
		for (int i = 0; i < n; i++) {
			if (isCenter[i]) {
				continue;
			}
			final V vI = verts[i];
			final double d = distFunc.applyAsDouble(vI, firstV);
			if (Double.isNaN(d)) {
				throw new IllegalArgumentException("distFunc returned NaN for (" + vI + ", " + firstV + ")");
			}
			minDist[i] = d;
		}

		// Main loop
		while (centerCount < k) {
			// Find farthest non-center (ties -> lowest index)
			int farthest = -1;
			double maxD = -1.0;
			for (int i = 0; i < n; i++) {
				if (isCenter[i]) {
					continue;
				}
				final double d = minDist[i];
				if (d > maxD || (d == maxD && i < farthest)) {
					maxD = d;
					farthest = i;
				}
			}

			// Add farthest as new center
			isCenter[farthest] = true;
			minDist[farthest] = 0.0;
			centers[centerCount++] = farthest;

			// Update minDist using only the newly added center
			final V newC = verts[farthest];
			for (int i = 0; i < n; i++) {
				if (isCenter[i]) {
					continue;
				}
				final V vI = verts[i];
				final double d = distFunc.applyAsDouble(vI, newC);
				if (Double.isNaN(d)) {
					throw new IllegalArgumentException("distFunc returned NaN for (" + vI + ", " + newC + ")");
				}
				if (d < minDist[i]) {
					minDist[i] = d;
				}
			}
		}

		final List<V> result = new ArrayList<>(k);
		for (int i = 0; i < k; i++) {
			final V v = verts[centers[i]];
			result.add(v);
		}
		return result;
	}
}