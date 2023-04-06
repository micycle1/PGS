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
public class PEdge {

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
		return a.copy().add(b).div(2);
	}

	/**
	 * Calculates the Euclidean distance of this PEdge.
	 * 
	 * @return
	 */
	public float length() {
		return a.dist(b);
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
}