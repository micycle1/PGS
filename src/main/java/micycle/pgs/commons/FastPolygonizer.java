package micycle.pgs.commons;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import micycle.pgs.color.RGB;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * This class solves the problem of finding polygons formed by a collection of
 * edges.
 * <p>
 * At least twice as fast as JTS'
 * {@link org.locationtech.jts.operation.polygonize.Polygonizer polygonizer}.
 * 
 * @author Michael Carleton
 *
 */
public class FastPolygonizer {

	// Implements John Hughes' answer at math.stackexchange.com/a/1809750/583690

	private FastPolygonizer() {
	}

	/**
	 * Polygonizes a set of edges which represent linework that forms some polygonal
	 * arrangement (a planar graph).
	 * 
	 * @param edges a collection of <b>NODED</b> (i.e. non intersecting / must only
	 *              meet at their endpoints) edges. The collection can contain
	 *              duplicates.
	 * @return a GROUP PShape, where each child shape represents a polygon face
	 *         formed by the edges input
	 */
	public static PShape polygonize(Collection<PEdge> edges) {
		// A “dart” is a directed edge.
		final Set<Dart> darts = new HashSet<>(edges.size() * 2);
		// A “star” for a given vertex is the set of darts that originate at that vertex
		final Map<PVector, HashSet<Dart>> stars = new HashMap<>(edges.size());

		edges.forEach(e -> {
			// For each input segment, compute the two darts associated to that segment and
			// place these two darts in a list
			final Dart e1 = new Dart(e.a, e.b);
			final Dart e2 = new Dart(e.b, e.a);
			e1.reverse = e2;
			e2.reverse = e1;
			darts.add(e1);
			darts.add(e2);

			// For each vertex, form a list of all darts in the dart-list that originate at
			// that vertex
			stars.computeIfAbsent(e.a, k -> new HashSet<>()).add(e1);
			stars.computeIfAbsent(e.b, k -> new HashSet<>()).add(e2);
		});

		// Remove dangles
		for (Iterator<Entry<PVector, HashSet<Dart>>> it = stars.entrySet().iterator(); it.hasNext();) {
			Entry<PVector, HashSet<Dart>> entry = it.next();
			if (entry.getValue().size() == 1) {
				it.remove();
				Dart dangle = entry.getValue().iterator().next();
				darts.remove(dangle);
				darts.remove(dangle.reverse);
				stars.get(dangle.b).remove(dangle.reverse);
			}
		}

		/*
		 * For each vertex, for each dart in the star of that vertex, compute the vector
		 * associated to that dart, and then find its angle from the x-axis; then sort
		 * the darts in the star by increasing angle.
		 */
		final Map<PVector, Map<Dart, Dart>> starsOrdered = new HashMap<>(stars.size());
		stars.forEach((vertex, star) -> {
			final List<Dart> starEdges = new ArrayList<>(star);
			starEdges.sort((e1, e2) -> Double.compare(e1.angle, e2.angle)); // star darts sorted by angle

			/*
			 * You can find a given dart within a given star in constant time by storing the
			 * darts of a star in a dictionary, in which the “value” for each dart is the
			 * “next” dart.
			 */
			final Map<Dart, Dart> nextEdgeMap = new HashMap<>(starEdges.size());
			for (int i = 0; i < starEdges.size(); i++) {
				nextEdgeMap.put(starEdges.get(i), starEdges.get((i + 1) % starEdges.size()));
			}

			starsOrdered.put(vertex, nextEdgeMap);
		});

		/*
		 * This algorithm natively produces a polygon which forms the union of all other
		 * faces (a "big outside polygon"), which is discovered using its vertex count
		 * (it is the face with the highest vertex count) and removed at the end.
		 */
		int largestVertexCount = -1;
		int bopIndex = -1; // "big outside polygon" index

		final PShape mesh = new PShape(PConstants.GROUP);

		int index = 0;
		for (Dart d : darts) { // Main algorithm (to find poly faces)
			if (d.seen) {
				continue; // do not process dart if seen before
			}

			Dart dart = d;
			final Deque<Dart> stack = new ArrayDeque<>();
			stack.push(dart);

			while (!dart.seen) {
				dart.seen = true;
				/*
				 * To ignore "dangling" edges, push each new dart onto the stack UNLESS the
				 * stack is non-empty and the top of the stack is the opposite dart, in which
				 * case you pop off the top of the stack.
				 */
				if (!stack.isEmpty() && dart != d) {
					if (stack.peek() == dart.reverse) {
						stack.pop();
					} else {
						stack.push(dart);
					}
				}

				/*
				 * In the star for vertex j, the dart [j,i] appears somewhere; set d to the next
				 * dart in the star after that one (we have already created a map for this, so
				 * the operation is O(1)).
				 */
				dart = starsOrdered.get(dart.b).get(dart.reverse);
			}

			if (stack.size() > 2) {

				final PShape polygon = new PShape(PShape.PATH);
				polygon.setFill(true);
				polygon.setStroke(true);
				polygon.setStrokeWeight(3);
				polygon.setStrokeCap(PConstants.ROUND);
				polygon.setFill(RGB.WHITE);
				polygon.setStroke(RGB.PINK);

				polygon.beginShape();
				while (!stack.isEmpty()) {
					final Dart e = stack.pop();
					polygon.vertex(e.b.x, e.b.y);
				}
				polygon.endShape(PConstants.CLOSE);

				mesh.addChild(polygon);

				if (polygon.getVertexCount() > largestVertexCount) {
					largestVertexCount = polygon.getVertexCount();
					bopIndex = index;
				}
				index++;
			}
		}

		if (bopIndex > -1) {
			mesh.removeChild(bopIndex);
		}

		return mesh;
	}

	/**
	 * Directional PEdge.
	 */
	static final class Dart {

		final PVector a, b;
		final double angle;
		Dart reverse;
		boolean seen = false;

		Dart(PVector a, PVector b) {
			this.a = a;
			this.b = b;
			PVector vector = PVector.sub(b, a);
			angle = atan2Quick(vector.y, vector.x);
		}

		@Override
		/**
		 * Directional hash
		 */
		public int hashCode() {
			return Float.floatToIntBits(a.x + a.y) ^ Float.floatToIntBits(b.x + b.y - 1);
		}

		@Override
		public boolean equals(Object obj) {
			if (obj instanceof Dart) {
				Dart other = (Dart) obj;
				return (other.a.equals(a) && other.b.equals(b));
			}
			return false;
		}

		@Override
		public String toString() {
			return a.toString() + " -> " + b.toString();
		}

		private static double atan2Quick(final double y, final double x) {
			final double THREE_QRTR_PI = Math.PI * 0.75;
			final double QRTR_PI = Math.PI * 0.25;

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

	}
}
