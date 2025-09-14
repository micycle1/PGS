package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.algorithm.Distance;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

import processing.core.PVector;

/**
 * An undirected edge / line segment joining 2 PVectors.
 * <p>
 * Note: PEdges <code>PEdge(a, b)</code> and <code>PEdge(b, a)</code> are
 * considered equal (though the ordering for .a and .b is preserved).
 *
 * @author Michael Carleton
 *
 */
public class PEdge implements Comparable<PEdge> {

	public final PVector a, b;
	private final Coordinate aCoord, bCoord;

	/**
	 * A null PEdge.
	 */
	public PEdge() {
		this.a = null;
		this.b = null;
		aCoord = null;
		bCoord = null;
	}

	public PEdge(final PVector a, final PVector b) {
		this.a = a;
		this.b = b;
		aCoord = coordFromPVector(a);
		bCoord = coordFromPVector(b);
	}

	public PEdge(final double x1, final double y1, final double x2, final double y2) {
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
	 * Returns the point on the segment [a, b] at parameter t, where t=0 gives a,
	 * t=1 gives b, and values in between give the corresponding point on the line.
	 * If t is outside [0,1] it will be clamped.
	 */
	public PVector pointAt(double t) {
		if (t < 0) {
			t = 0;
		} else if (t > 1) {
			t = 1;
		}
		final float ft = (float) t;
		return new PVector(a.x + (b.x - a.x) * ft, a.y + (b.y - a.y) * ft);
	}

	/**
	 * Calculates the Euclidean distance of this PEdge.
	 */
	public float length() {
		return a.dist(b);
	}

	/**
	 * Computes the minimum distance between this and another edge.
	 */
	public double distance(final PEdge other) {
		return Distance.segmentToSegment(aCoord, bCoord, other.aCoord, other.bCoord);
	}

	/**
	 * Computes the distance from a point p to this edge.
	 */
	public double distance(final PVector point) {
		return Distance.pointToSegment(coordFromPVector(point), aCoord, bCoord);
	}

	public PVector closestPoint(final PVector point) {
		final LineSegment l = new LineSegment(aCoord, bCoord);
		return coordToPVector(l.closestPoint(coordFromPVector(point)));
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

	public List<PVector> sample(double d) {
		if (d <= 0) {
			throw new IllegalArgumentException("d must be > 0");
		}

		final double len = length(); // length() returns float, assign to double
		final List<PVector> samples = new ArrayList<>();

		// degenerate segment
		if (len == 0.0) {
			samples.add(new PVector(a.x, a.y));
			return samples;
		}

		// always include the first point
		samples.add(new PVector(a.x, a.y));

		final double step = d / len; // step in parametric [0..1] space
		// add intermediate samples at t = step, 2*step, ... while t < 1.0
		for (double t = step; t < 1.0; t += step) {
			samples.add(pointAt(t));
		}

		// always include the last point (avoid duplicate if last intermediate ~= b)
		final PVector last = samples.get(samples.size() - 1);
		final double eps = 1e-6;
		if (last.dist(b) > eps) {
			samples.add(new PVector(b.x, b.y));
		}

		return samples;
	}

	@Override
	/**
	 * Direction-agnostic hash.
	 */
	public int hashCode() {
		return Float.floatToIntBits(b.y + a.y) ^ Float.floatToIntBits(b.x + a.x - 1);
	}

	@Override
	public boolean equals(final Object obj) {
		if (obj instanceof final PEdge other) {
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

	private static boolean equals(final PVector a, final PVector b) {
		return a.x == b.x && a.y == b.y;
	}

	@Override
	public int compareTo(final PEdge other) {
		final PVector thisMidpoint = midpoint();
		final PVector otherMidpoint = other.midpoint();
		return comparePVectors(thisMidpoint, otherMidpoint);
	}

	/**
	 * Helper method to compare two PVectors lexicographically.
	 */
	private int comparePVectors(final PVector v1, final PVector v2) {
		if (v1.x != v2.x) {
			return Float.compare(v1.x, v2.x);
		}
		if (v1.y != v2.y) {
			return Float.compare(v1.y, v2.y);
		}
		return Float.compare(v1.z, v2.z);
	}

	private static final Coordinate coordFromPVector(final PVector p) {
		return new Coordinate(p.x, p.y);
	}

	private static final PVector coordToPVector(final Coordinate c) {
		return new PVector((float) c.x, (float) c.y);
	}
}