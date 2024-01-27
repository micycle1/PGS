package micycle.pgs.commons;

import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.LineString;

import net.jafama.FastMath;

/**
 * This class simplifies curves while maintaining perceptual appearance. It does
 * this through the iterative removal of the least shape-relevant "kinks", which
 * are where consecutive line segments meet (at a vertex).
 * <p>
 * The shape relevance of a kink is calculated based on its turn angle and the
 * lengths of the line segments it comprises. Kinks with larger turn angles and
 * longer segments contribute more to the overall shape of the curve and are
 * thus more likely to be preserved during the simplification process.
 * <p>
 * This class uses an efficient implementation that reduces the time complexity
 * from a naive O(n<sup>2</sup>) to O(n log n), making it suitable for
 * processing large datasets. It handles both closed and open curves. The
 * algorithm is not guaranteed to be topology-preserving.
 *
 * @author Michael Carleton
 */
public class DiscreteCurveEvolution {

	/*-
	 * Relevance measure was lifted from github.com/DiegoCatalano/Catalano-Framework,
	 * which implements 'Convexity Rule for Shape Decomposition Based on Discrete
	 * Contour Evolution'.
	 * This implementation uses doubly-linked coordinates & ordered set for better 
	 * time complexity (nlogn vs n^2).
	 */

	/**
	 * The callback interface for determining the termination condition of the
	 * Discrete Curve Evolution (DCE) process.
	 * <p>
	 * This functional interface defines a single method that decides whether the
	 * DCE algorithm should terminate based on the current kink (having a candidate
	 * vertex), using its coordinates, relevance score, and the number of vertices
	 * remaining in the simplified geometry. Implementations can use this method to
	 * provide custom termination logic which may depend on various factors, such as
	 * a threshold relevance score, a specific number of vertices to preserve, or
	 * other criteria.
	 *
	 * @see #process(LineString, DCETerminationCallback)
	 */
	@FunctionalInterface
	public interface DCETerminationCallback {
		/**
		 * Determines whether the DCE process should terminate at the current step.
		 *
		 * @param currentVertex     The coordinates of the current kink being evaluated.
		 * @param relevance         The relevance score of the current kink. A value of
		 *                          30 or lower is generally imperceptible.
		 * @param verticesRemaining The number of vertices remaining in the simplified
		 *                          geometry. Note: this does not count the closing
		 *                          vertex of closed geometries.
		 * @return {@code true} if the DCE process should terminate; {@code false}
		 *         otherwise.
		 */
		boolean shouldTerminate(Coordinate currentVertex, double relevance, int verticesRemaining);
	}

	/**
	 * Applies the Discrete Curve Evolution (DCE) algorithm to a given polygonal
	 * curve defined by a {@link LineString}. This algorithm simplifies the curve by
	 * selectively removing vertices with the least shape-relevance while attempting
	 * to preserve the overall perceptual appearance of the shape. The algorithm
	 * uses a provided {@link DCETerminationCallback} to determine the termination
	 * condition for the simplification process.
	 * <p>
	 * If the input shape is closed (i.e., a polygon's LinearRing), the algorithm
	 * ensures that the resulting simplified shape remains closed. Vertices with the
	 * smallest turn angles and neighboring line segments are considered the least
	 * relevant and are candidates for removal.
	 * <p>
	 * The process continues iteratively, removing the least relevant vertex and
	 * recalculating the relevance of affected neighboring vertices, until the
	 * termination callback indicates that the process should stop.
	 *
	 * @param lineString          The {@link LineString} representing the original
	 *                            coordinates of the shape to be simplified.
	 * @param terminationCallback The callback used to determine when the
	 *                            simplification process should terminate.
	 * @return An array of coordinates representing the simplified shape, with a
	 *         potentially reduced number of vertices that maintains the perceptual
	 *         appearance of the original curve.
	 */
	public static Coordinate[] process(LineString lineString, DCETerminationCallback terminationCallback) {
		final boolean closed = lineString.isClosed();
		final Coordinate[] coords = lineString.getCoordinates();

		/*
		 * Initialise kinks (initially unlinked to each other), then link each one to
		 * its neighbours once all have been instantiated.
		 */
		final List<Kink> kinks = Arrays.asList(coords).stream().map(c -> new Kink(c)).collect(Collectors.toList());
		for (int i = 0; i < kinks.size(); i++) {
			Kink prev = kinks.get(i == 0 ? kinks.size() - 1 : i - 1);
			Kink next = kinks.get(i == kinks.size() - 1 ? 0 : i + 1);
			kinks.get(i).linkTo(prev, next);
		}

		if (!closed) { // unlink start and end if coords are unclosed
			kinks.get(0).markAsStart();
			kinks.get(kinks.size() - 1).markAsEnd();
		}

		/*
		 * The first element of the treeset is always the Kink having the least
		 * relevance. When it is removed, its neighbors must also be removed, re-linked
		 * then re-inserted with their new relevance values.
		 * 
		 * Topology preservation can be achieved by ignoring kinks whose triangle formed
		 * by its three vertices contains other vertices (unimplemented).
		 */
		final TreeSet<Kink> kinkRelevanceTree = new TreeSet<>(kinks);
		while (kinkRelevanceTree.size() > 2) {
			Kink candidate = kinkRelevanceTree.pollFirst();
			if (terminationCallback.shouldTerminate(candidate.c, candidate.relevance, kinkRelevanceTree.size() + 1)) {
				kinkRelevanceTree.add(candidate); // reinsert polled element
				break;
			}
			Kink previous = candidate.prev;
			Kink next = candidate.next;
			kinkRelevanceTree.remove(previous); // temporarily remove prev coordinate
			kinkRelevanceTree.remove(next); // temporarily remove next coordinate
			candidate.unlink(); // mutually link prev and next (and recompute their relevance)
			kinkRelevanceTree.add(previous); // re-add (with new relevance score) to reorder
			kinkRelevanceTree.add(next); // re-add (with new relevance score) to reorder
		}

		CoordinateList output = new CoordinateList();

		Kink first = closed ? kinkRelevanceTree.first() : kinks.get(0);
		Kink current = first;
		do {
			output.add(current.c);
		} while ((current = current.next) != first && current != null);
		if (closed) {
			output.closeRing();
		}

		return output.toCoordinateArray();
	}

	/**
	 * A doubly-linked coordinate. Models a "kink" in a curve where a pair of
	 * consecutive line segments meet.
	 */
	private static class Kink implements Comparable<Kink> {

		final Coordinate c;
		Kink prev, next;
		double relevance;

		Kink(Coordinate c) {
			this.c = c;
		}

		void linkTo(Kink prev, Kink next) {
			this.prev = prev;
			this.next = next;
			computeRelevance();
		}

		/**
		 * Mark this link as the starting one (relevant for non-closed inputs).
		 */
		void markAsStart() {
			this.prev = null;
			relevance = Double.MAX_VALUE - 1; // -1 because treeset values must be unique, can't be shared with end vertex
		}

		/**
		 * Mark this link as the terminating one (relevant for non-closed inputs).
		 */
		void markAsEnd() {
			this.next = null;
			relevance = Double.MAX_VALUE;
		}

		/**
		 * Removes this kink from the chain, links together its neighbours, and updates
		 * the relevancy measure of both.
		 */
		void unlink() {
			next.prev = prev;
			prev.next = next;
			prev.computeRelevance();
			next.computeRelevance();
		}

		void computeRelevance() {
			if (prev == null || next == null) {
				relevance = prev == null ? Double.MAX_VALUE - 1 : Double.MAX_VALUE;
				return;
			}
			double uv = prev.c.distance(c);
			double vw = c.distance(next.c);
			double uw = prev.c.distance(next.c);

			double alpha = FastMath.acos((vw * vw + uv * uv - uw * uw) / (2 * vw * uv));

			double a = 180 - alpha * 180 / Math.PI; // turning angle (0-180)
			relevance = a * vw * uv / (vw + uv);
		}

		@Override
		public int compareTo(Kink o) {
			return Double.compare(this.relevance, o.relevance);
		}

	}

}