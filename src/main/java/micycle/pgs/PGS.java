package micycle.pgs;

import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;

import micycle.pgs.color.RGB;
import processing.core.PShape;
import processing.core.PVector;

/**
 * This class houses functions used by the library internally.
 * 
 * @author Michael Carleton
 */
final class PGS {

	/** Defines number of vertices that comprise constructed geometries. */
	static final int SHAPE_SAMPLES = 80;

	/**
	 * PGS global geometry factory (uses 32 bit float precision).
	 */
	public static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	private PGS() {
	}

	/**
	 * Create a LINES PShape, ready for vertices.
	 * 
	 * @param strokeColor  nullable
	 * @param strokeCap    nullable default = ROUND
	 * @param strokeWeight nullable. default = 2
	 * @return
	 */
	static final PShape prepareLinesPShape(Integer strokeColor, Integer strokeCap, Integer strokeWeight) {
		if (strokeColor == null) {
			strokeColor = RGB.PINK;
		}
		if (strokeCap == null) {
			strokeCap = ROUND;
		}
		if (strokeWeight == null) {
			strokeWeight = 2;
		}
		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(strokeCap);
		lines.setStroke(true);
		lines.setStrokeWeight(strokeWeight);
		lines.setStroke(strokeColor);
		lines.beginShape(LINES);
		return lines;
	}

	/**
	 * Euclidean distance between two coordinates
	 */
	static final double distance(Coordinate a, Coordinate b) {
		double deltaX = a.y - b.y;
		double deltaY = a.x - b.x;
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	/**
	 * Euclidean distance between two points
	 */
	static final double distance(Point a, Point b) {
		double deltaX = a.getY() - b.getY();
		double deltaY = a.getX() - b.getX();
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	static final LineString createLineString(PVector a, PVector b) {
		return GEOM_FACTORY.createLineString(new Coordinate[] { coordFromPVector(a), coordFromPVector(b) });
	}

	static final Point createPoint(double x, double y) {
		return GEOM_FACTORY.createPoint(new Coordinate(x, y));
	}

	static final Point pointFromPVector(PVector p) {
		return GEOM_FACTORY.createPoint(new Coordinate(p.x, p.y));
	}

	static final Coordinate coordFromPoint(Point p) {
		return new Coordinate(p.getX(), p.getY());
	}

	static final Coordinate coordFromPVector(PVector p) {
		return new Coordinate(p.x, p.y);
	}

	/**
	 * Reflection-based workaround to get the fill color of a PShape (this field is
	 * usually private).
	 */
	static final int getPShapeFillColor(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("fillColor");
			f.setAccessible(true);
			return f.getInt(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * Requires a closed hole
	 * 
	 * @param points
	 * @return
	 */
	static final boolean isClockwise(List<PVector> points) {
		boolean closed = true;
		if (points.get(0).equals(points.get(points.size() - 1))) {
			closed = false;
			points.add(points.get(0)); // mutate list
		}
		double area = 0;

		for (int i = 0; i < (points.size()); i++) {
			int j = (i + 1) % points.size();
			area += points.get(i).x * points.get(j).y;
			area -= points.get(j).x * points.get(i).y;
		}

		if (!closed) {
			points.remove(points.size() - 1); // revert mutation
		}

		return (area < 0);
	}

	/**
	 * Provides convenient iteration of exterior and linear rings (if any) of a
	 * polygonal JTS geometry. Supports MultiGeometries.
	 * 
	 * @author Michael Carleton
	 */
	static final class LinearRingIterator implements Iterable<LinearRing> {

		private LinearRing[] array;

		/**
		 * Constructs the iterator for the given geometry. The first ring returned by
		 * the iterator is the exterior ring; all other rings (if any) are interior
		 * rings.
		 * 
		 * @param g input geometry
		 */
		public LinearRingIterator(Geometry g) {
			ArrayList<LinearRing> rings = new ArrayList<>(g.getNumGeometries());
			for (int i = 0; i < g.getNumGeometries(); i++) {
				Polygon poly = (Polygon) g.getGeometryN(i);
				rings.add(poly.getExteriorRing());
				for (int j = 0; j < poly.getNumInteriorRing(); j++) {
					rings.add(poly.getInteriorRingN(j));
				}
			}
			array = rings.toArray(new LinearRing[rings.size()]);
		}

		@Override
		public Iterator<LinearRing> iterator() {
			return new Iterator<LinearRing>() {

				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < array.length;
				}

				@Override
				public LinearRing next() {
					if (!hasNext()) {
						throw new NoSuchElementException();
					}
					return array[currentIndex++];
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
		}
	}

	/**
	 * Represents an edge between 2 PVectors.
	 * 
	 * @author Michael Carleton
	 *
	 */
	static final class PEdge {

		final PVector a, b;

		public PEdge(PVector a, PVector b) {
			this.a = a;
			this.b = b;
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
	}

}
