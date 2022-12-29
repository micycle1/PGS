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
 * Discrete Curve Evolution simplifies polygonal curves, neglecting minor
 * distortions while preserving perceptual appearance. The main idea of the
 * process is a stepwise elimination of kinks (vertex where consecutive line
 * segments meet) that are least relevant to the shape of the polygonal curve.
 * The relevance of kinks is intended to reflect their contribution to the
 * overall shape of the polygonal curve.
 * <p>
 * The shape relevance of every kink can be defined by the turn angle and the
 * lengths of the neighboring line segments; tthe larger both the relative
 * lengths and the turn angle of a kink, the greater is its contribution to the
 * shape of a curve.
 * 
 * @author Michael Carleton
 */
public class DiscreteCurveEvolution {

	/*-
	 * Relevance measure was lifted from github.com/DiegoCatalano/Catalano-Framework,
	 * which implements 'Convexity Rule for Shape Decomposition Based on Discrete
	 * Contour Evolution'.
	 * This implementation uses doubly-linked coordinates + ordered set for better 
	 * time complexity (nlogn vs n^2).
	 */

	private int vertexPreserveCount;

	/**
	 * 
	 * @param vertices the number of vertices to preserve
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

	public Coordinate[] process(Coordinate[] coords) {
		boolean closed = coords[0].equals2D(coords[coords.length - 1]);
		int vertexRemoveCount = Math.min(coords.length - vertexPreserveCount, coords.length - (closed ? 4 : 2));

		/*
		 * Create linked coordinates, then link each one to its neighbours once all have
		 * been instantiated.
		 */
		final List<LinkedCoordinate> linkedCoords = Arrays.asList(coords).stream().map(c -> new LinkedCoordinate(c))
				.collect(Collectors.toList());
		for (int i = 0; i < linkedCoords.size(); i++) {
			LinkedCoordinate prev = linkedCoords.get(i == 0 ? linkedCoords.size() - 1 : i - 1);
			LinkedCoordinate next = linkedCoords.get(i == linkedCoords.size() - 1 ? 0 : i + 1);
			linkedCoords.get(i).assignLinks(prev, next);
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
	 * Models a "kink".
	 */
	private static class LinkedCoordinate implements Comparable<LinkedCoordinate> {

		final Coordinate c;
		LinkedCoordinate prev, next;
		double relevance;

		public LinkedCoordinate(Coordinate c) {
			this.c = c;
		}

		void assignLinks(LinkedCoordinate prev, LinkedCoordinate next) {
			this.prev = prev;
			this.next = next;
			computeRelevance();
		}

		void markAsStart() {
			this.prev = null;
			relevance = Double.MAX_VALUE / 2;
		}

		void markAsEnd() {
			this.next = null;
			relevance = Double.MAX_VALUE;
		}

		void relinkToNext(LinkedCoordinate newNext) {
			this.next = newNext; // link me to other
			next.prev = this; // link other to me

			this.computeRelevance();
			newNext.computeRelevance();
		}

		void computeRelevance() {
			if (prev == null || next == null) {
				relevance = prev == null ? Double.MAX_VALUE / 2 : Double.MAX_VALUE;
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