package micycle.pgs;

import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

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
class PGS {

	/** Defines number of vertices that comprise constructed geometries. */
	protected static final int SHAPE_SAMPLES = 80;

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
	protected static PShape prepareLinesPShape(Integer strokeColor, Integer strokeCap, Integer strokeWeight) {
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
	protected static double distance(Coordinate a, Coordinate b) {
		double deltaX = a.y - b.y;
		double deltaY = a.x - b.x;
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	/**
	 * Euclidean distance between two points
	 */
	protected static double distance(Point a, Point b) {
		double deltaX = a.getY() - b.getY();
		double deltaY = a.getX() - b.getX();
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	protected static LineString createLineString(PVector a, PVector b) {
		return GEOM_FACTORY.createLineString(new Coordinate[] { coordFromPVector(a), coordFromPVector(b) });
	}

	protected static Point createPoint(double x, double y) {
		return GEOM_FACTORY.createPoint(new Coordinate(x, y));
	}

	protected static Point pointFromPVector(PVector p) {
		return GEOM_FACTORY.createPoint(new Coordinate(p.x, p.y));
	}

	protected static Coordinate coordFromPoint(Point p) {
		return new Coordinate(p.getX(), p.getY());
	}

	protected static Coordinate coordFromPVector(PVector p) {
		return new Coordinate(p.x, p.y);
	}

	/**
	 * Transforms a list of points into a POINTS PShape.
	 */
	protected static PShape toPointsPShape(Iterable<PVector> points) {
		PShape shape = new PShape();
		shape.setFamily(PShape.GEOMETRY);
		shape.setStrokeCap(ROUND);
		shape.beginShape(PShape.POINTS);
		points.forEach(p -> shape.vertex(p.x, p.y));
		shape.endShape();
		return shape;
	}

	/**
	 * Reflection-based workaround to get the fill color of a PShape (this field is
	 * usually private).
	 */
	protected static final int getPShapeFillColor(final PShape sh) {
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
	protected static boolean isClockwise(List<PVector> points) {
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
			points.remove(points.size() - 1); // undo mutation
		}

		return (area < 0);
	}

	/**
	 * Provides convenient iteration of exterior and linear rings (if any) of a JTS
	 * geometry.
	 * 
	 * @author Michael Carleton
	 */
	protected static class LinearRingIterator implements Iterable<LinearRing> {

		private LinearRing[] array;
		private int size;

		/**
		 * Constructs the iterator for the given geometry. The first ring returned by
		 * the iterator is the exterior ring; all other rings (if any) are interior
		 * rings.
		 * 
		 * @param g input geometry
		 */
		public LinearRingIterator(Geometry g) {
			Polygon poly = (Polygon) g;
			this.size = 1 + poly.getNumInteriorRing();
			this.array = new LinearRing[size];
			array[0] = poly.getExteriorRing();
			for (int i = 0; i < poly.getNumInteriorRing(); i++) {
				array[i + 1] = poly.getInteriorRingN(i);
			}
		}

		@Override
		public Iterator<LinearRing> iterator() {
			return new Iterator<LinearRing>() {

				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < size;
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

}
