package micycle.pgs.utility;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.PrecisionModel;

import processing.core.PVector;

/**
 * Represents an adirectional edge between 2 PVectors.
 * 
 * @author Michael Carleton
 *
 */
public class PEdge {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	public final PVector a, b;

	public PEdge(PVector a, PVector b) {
		this.a = a;
		this.b = b;
	}

	public PEdge(double x1, double y1, double x2, double y2) {
		this(new PVector((float) x1, (float) y1), new PVector((float) x2, (float) y2));
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

	public LineString toLineString() {
		return GEOM_FACTORY.createLineString(new Coordinate[] { new Coordinate(a.x, a.y), new Coordinate(b.x, b.y) });
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