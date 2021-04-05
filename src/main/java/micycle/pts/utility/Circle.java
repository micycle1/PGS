package micycle.pts.utility;

import processing.core.PVector;

public final class Circle {

	private static final float MULTIPLICATIVE_EPSILON = (float) (1 + 1e-14);

	/** center */
	public final PVector c;
	/** radius */
	public final float r;

	public Circle(PVector c, float r) {
		this.c = c;
		this.r = r;
	}

	public boolean contains(PVector p) {
		return c.dist(p) <= r * MULTIPLICATIVE_EPSILON;
	}
}
