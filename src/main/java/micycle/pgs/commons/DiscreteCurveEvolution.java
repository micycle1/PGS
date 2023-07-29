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
 * The DiscreteCurveEvolution class is used to simplify polygonal curves while
 * maintaining their perceptual appearance. It does this through the systematic
 * removal of the least shape-relevant vertices (kinks), where consecutive line
 * segments meet.
 * <p>
 * The shape relevance of a vertex is calculated based on its turn angle and the
 * lengths of its neighboring line segments. Vertices with larger turn angles
 * and neighboring segment lengths contribute more to the overall shape of the
 * curve and are thus more likely to be preserved during the simplification
 * process.
 * <p>
 * This class uses an efficient implementation that reduces the time complexity
 * from O(n^2) to O(n log n), making it suitable for processing large datasets.
 * It can handle both closed and open curves.
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

	private int vertexPreserveCount;

	/**
	 * Initialises an instance of the DiscreteCurveEvolution class. The provided
	 * parameter specifies the number of vertices to be preserved during the curve
	 * evolution process.
	 *
	 * @param vertices The number of vertices to be preserved. These vertices are
	 *                 chosen based on their shape relevance, which is calculated
	 *                 from their turn angle and the lengths of neighboring line
	 *                 segments. Larger relative lengths and turn angles make a
	 *                 vertex more likely to be preserved.
	 */
	public DiscreteCurveEvolution(int vertices) {
		this.vertexPreserveCount = vertices;
	}

	/**
	 * 
	 * @param lineString a closed or unclosed LineString
	 * @return
	 */
	public Coordinate[] process(LineString lineString) {
		return process(lineString.getCoordinates());
	}

	/**
	 * Processes an array of coordinates and applies the Discrete Curve Evolution
	 * algorithm. The algorithm simplifies the polygonal curve defined by the input
	 * coordinates while preserving the perceptual appearance by preserving a
	 * certain number of vertices with high shape relevance. If the shape is closed,
	 * the method handles it appropriately by preserving the closure.
	 * <p>
	 * Note: The vertices removed are always the least relevant vertices, those with
	 * the smallest turn angles and neighboring line segments.
	 *
	 * @param coords The original coordinates of the shape, an array of Coordinates
	 *               defining the polygonal curve.
	 *
	 * @return An array of Coordinates representing a simplified version of the
	 *         original shape. The number of vertices is reduced while preserving
	 *         the perceptual appearance of the shape.
	 */
	public Coordinate[] process(Coordinate[] coords) {
		boolean closed = coords[0].equals2D(coords[coords.length - 1]);
		int vertexRemoveCount = Math.min(coords.length - vertexPreserveCount, coords.length - (closed ? 4 : 2));

		/*
		 * Initialise LinkedCoordinates (initially unlinked to each other), then link
		 * each one to its neighbours once all have been instantiated.
		 */
		final List<LinkedCoordinate> linkedCoords = Arrays.asList(coords).stream().map(c -> new LinkedCoordinate(c))
				.collect(Collectors.toList());
		for (int i = 0; i < linkedCoords.size(); i++) {
			LinkedCoordinate prev = linkedCoords.get(i == 0 ? linkedCoords.size() - 1 : i - 1);
			LinkedCoordinate next = linkedCoords.get(i == linkedCoords.size() - 1 ? 0 : i + 1);
			linkedCoords.get(i).linkTo(prev, next);
		}

		if (!closed) { // unlink start and end if coords are unclosed
			linkedCoords.get(0).markAsStart();
			linkedCoords.get(linkedCoords.size() - 1).markAsEnd();
		}

		/*
		 * The first element of the treeset is always the LinkedCoordinate having the
		 * least relevance. When it is removed, its neighbors must also be removed,
		 * re-linked then re-inserted with their new relevance values.
		 */
		final TreeSet<LinkedCoordinate> t = new TreeSet<>(linkedCoords);
		for (int i = 0; i < vertexRemoveCount; i++) {
			LinkedCoordinate removed = t.pollFirst();
			LinkedCoordinate previous = removed.prev;
			LinkedCoordinate next = removed.next;
			t.remove(previous); // temporarily remove prev coordinate
			t.remove(next); // temporarily remove next coordinate
			previous.relinkToNext(next); // mutually link prev and next (and recompute relevance)
			t.add(previous); // re-add (with new relevance score)
			t.add(next); // re-add (with new relevance score)
		}

		CoordinateList output = new CoordinateList();

		LinkedCoordinate first = closed ? t.first() : linkedCoords.get(0);
		LinkedCoordinate current = first;
		do {
			output.add(current.c);
		} while ((current = current.next) != first && current != null);
		if (closed) {
			output.closeRing();
		}

		return output.toCoordinateArray();
	}

	/**
	 * A doubly-linked coordinate; models a "kink" in a linestring where consecutive
	 * line segments meet.
	 */
	private static class LinkedCoordinate implements Comparable<LinkedCoordinate> {

		final Coordinate c;
		LinkedCoordinate prev, next;
		double relevance;

		public LinkedCoordinate(Coordinate c) {
			this.c = c;
		}

		void linkTo(LinkedCoordinate prev, LinkedCoordinate next) {
			this.prev = prev;
			this.next = next;
			computeRelevance();
		}

		void markAsStart() {
			this.prev = null;
			relevance = Double.MAX_VALUE - 1; // -1 because treeset values must be unique
		}

		void markAsEnd() {
			this.next = null;
			relevance = Double.MAX_VALUE;
		}

		/**
		 * Mutually links this LC ('previous' in the chain) with the other LC ('next' in
		 * the chain'), and updates the relevancy measure of both.
		 * 
		 * @param newNext this LC's new 'next' LC
		 */
		void relinkToNext(LinkedCoordinate newNext) {
			this.next = newNext; // link me to other
			next.prev = this; // link other to me

			this.computeRelevance();
			newNext.computeRelevance();
		}

		private void computeRelevance() {
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
		public int compareTo(LinkedCoordinate o) {
			return Double.compare(this.relevance, o.relevance);
		}

	}

}