package micycle.pgs.commons;

import processing.core.PVector;

/**
 * Represents an adirectional edge between 2 PVectors.
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
	 * Rounds (mutates) the vertices of this PEdge.
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

	@Override
	/**
	 * Direction-agnostic hash
	 */
	public int hashCode() {
//			int x = Float.floatToIntBits(Math.min(a.x, b.x));
//			x = ((x >> 16) ^ Float.floatToIntBits(Math.max(a.x, b.x))) * 0x45d9f3b;
//			x = ((x >> 16) ^ Float.floatToIntBits(Math.min(a.y, b.y))) * 0x45d9f3b;
//			x = (x >> 16) ^ Float.floatToIntBits(Math.max(a.y, b.y));
//			return x;
		return Float.floatToIntBits(b.y + a.y) ^ Float.floatToIntBits(b.x + a.x - 1);
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof PEdge) {
			PEdge other = (PEdge) obj;
			return (other.a.equals(a) && other.b.equals(b)) || (other.a.equals(b) && other.b.equals(a));
		}
		return false;
	}

	@Override
	public PEdge clone() {
		return new PEdge(a.copy(), b.copy());
	}

	@Override
	public String toString() {
		return a.toString() + " <-> " + b.toString();
	}
}