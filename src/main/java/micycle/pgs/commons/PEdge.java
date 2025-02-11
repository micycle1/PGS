package micycle.pgs.commons;

import processing.core.PVector;

/**
 * An undirected edge / line segment joining 2 PVectors.
 * <p>
 * Note: PEdges <code>PEdge(a, b)</code> and <code>PEdge(b, a)</code> are
 * considered equal.
 * 
 * @author Michael Carleton
 *
 */
public class PEdge implements Comparable<PEdge> {

	public final PVector a, b;

	public PEdge(PVector a, PVector b) {
		this.a = a;
		this.b = b;
	}

	public PEdge(double x1, double y1, double x2, double y2) {
		this(new PVector((float) x1, (float) y1), new PVector((float) x2, (float) y2));
	}

	/**
	 * Rounds (mutates) the vertex coordinates of this PEdge to their closest ints.
	 * 
	 * @return this PEdge
	 */
	public PEdge round() {
		a.x = Math.round(a.x);
		a.y = Math.round(a.y);
		b.x = Math.round(b.x);
		b.y = Math.round(b.y);
		return this;
	}

	public PVector midpoint() {
		return PVector.add(a, b).div(2);
	}

	/**
	 * Calculates the Euclidean distance of this PEdge.
	 * 
	 * @return
	 */
	public float length() {
		return a.dist(b);
	}

	/**
	 * Calculates the subsection of this PEdge as a new PEdge.
	 * 
	 * @param from the start of the subsection as a normalized value along the
	 *             length of this PEdge. 'from' should be less than or equal to
	 *             'to'.
	 * @param to   the end of the subsection as a normalized value along the length
	 *             of this PEdge. 'to' should be greater than or equal to 'from'.
	 * @return A new PEdge representing the subsection of this PEdge between 'from'
	 *         (from end <code>a</code>) and 'to' (towards end <code>b</code>).
	 */
	public PEdge slice(double from, double to) {
		final double TOLERANCE = 1e-7;
		if (Math.abs(from) < TOLERANCE) {
			from = 0;
		}
		if (Math.abs(to - 1) < TOLERANCE) {
			to = 1;
		}
		if (from < 0 || to > 1 || from > to) {
			throw new IllegalArgumentException("Parameters 'from' and 'to' must be between 0 and 1, and 'from' must be less than or equal to 'to'.");
		}

		final PVector pointFrom;
		if (from == 0) {
			pointFrom = a.copy();
		} else {
			pointFrom = PVector.lerp(a, b, (float) from);
		}

		final PVector pointTo;
		if (to == 1) {
			pointTo = b.copy();
		} else {
			pointTo = PVector.lerp(a, b, (float) to);
		}

		return new PEdge(pointFrom, pointTo);
	}

	@Override
	/**
	 * Direction-agnostic hash.
	 */
	public int hashCode() {
		return Float.floatToIntBits(b.y + a.y) ^ Float.floatToIntBits(b.x + a.x - 1);
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof PEdge) {
			PEdge other = (PEdge) obj;
			return (equals(a, other.a) && equals(b, other.b)) || (equals(b, other.a) && equals(a, other.b));
		}
		return false;
	}

	/**
	 * @return a deep copy of this PEdge
	 */
	public PEdge copy() {
		return new PEdge(a.copy(), b.copy());
	}

	@Override
	public String toString() {
		return a.toString() + " <-> " + b.toString();
	}

	private static boolean equals(PVector a, PVector b) {
		return a.x == b.x && a.y == b.y;
	}

	@Override
	public int compareTo(PEdge other) {
		PVector thisMidpoint = midpoint();
		PVector otherMidpoint = other.midpoint();
		return comparePVectors(thisMidpoint, otherMidpoint);
	}

	/**
	 * Helper method to compare two PVectors lexicographically.
	 */
	private int comparePVectors(PVector v1, PVector v2) {
		if (v1.x != v2.x) {
			return Float.compare(v1.x, v2.x);
		}
		if (v1.y != v2.y) {
			return Float.compare(v1.y, v2.y);
		}
		return Float.compare(v1.z, v2.z);
	}
}