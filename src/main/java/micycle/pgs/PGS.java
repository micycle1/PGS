package micycle.pgs;

import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.noding.BasicSegmentString;
import org.locationtech.jts.noding.Noder;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.snap.SnappingNoder;
import org.locationtech.jts.operation.polygonize.Polygonizer;

import micycle.pgs.color.RGB;
import micycle.pgs.utility.Nullable;
import micycle.pgs.utility.PEdge;
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
	 * Create a LINES PShape, ready for vertices (shape.vertex(x, y) calls).
	 * 
	 * @param strokeColor  nullable (default = {@link RGB#PINK})
	 * @param strokeCap    nullable (default = ROUND)
	 * @param strokeWeight nullable (default = 2)
	 * @return
	 */
	static final PShape prepareLinesPShape(@Nullable Integer strokeColor, @Nullable Integer strokeCap, @Nullable Integer strokeWeight) {
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

	static final SegmentString createSegmentString(PVector a, PVector b) {
		return new BasicSegmentString(new Coordinate[] { PGS.coordFromPVector(a), PGS.coordFromPVector(b) }, null);
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

	static final PVector toPVector(Coordinate c) {
		return new PVector((float) c.x, (float) c.y);
	}

	/**
	 * Reflection-based workaround to get the fill color of a PShape (this field is
	 * usually inaccessible).
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
	 * Reflection-based workaround to get the stroke color of a PShape (this field
	 * is usually inaccessible).
	 */
	static final int getPShapeStrokeColor(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("strokeColor");
			f.setAccessible(true);
			return f.getInt(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * Reflection-based workaround to get the stroke weight of a PShape (this field
	 * is usually inaccessible).
	 */
	static final float getPShapeStrokeWeight(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("strokeWeight");
			f.setAccessible(true);
			return f.getFloat(sh);
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
	 * Polygonizes a set of line segments via noding.
	 * 
	 * @param segments list of segments (noded or non-noded)
	 * @param node     whether to node the segments before polygonization. If the
	 *                 segments constitute a conforming mesh, then set this as
	 *                 false; otherwise true.
	 * @return
	 */
	@SuppressWarnings("unchecked")
	static final Collection<Geometry> polygonizeSegments(Collection<SegmentString> segments, boolean node) {
		final Polygonizer polygonizer = new Polygonizer(); // TODO use QuickPolygonizer?
		polygonizer.setCheckRingsValid(false);

		if (node) {
			segments = nodeSegmentStrings(segments);
		}

		final Set<PEdge> edges = new HashSet<>();
		segments.forEach(ss -> {
			/*
			 * If the same LineString is added more than once to the polygonizer, the string
			 * is "collapsed" and not counted as an edge. Therefore a set is used to ensure
			 * strings are added once only to the polygonizer. A PEdge is used to determine
			 * this (since LineString hashcode doesn't work).
			 */
			final PEdge e = new PEdge(toPVector(ss.getCoordinate(0)), toPVector(ss.getCoordinate(1)));
			if (edges.add(e)) {
				final LineString l = PGS.GEOM_FACTORY.createLineString(new Coordinate[] { ss.getCoordinate(0), ss.getCoordinate(1) });
				polygonizer.add(l);
			}
		});
		return polygonizer.getPolygons(); // NOTE rather slow method
	}

	/**
	 * Computes a robust noding for a collection of SegmentStrings.
	 * 
	 * @param segments
	 * @return
	 */
	@SuppressWarnings("unchecked")
	static final Collection<SegmentString> nodeSegmentStrings(Collection<SegmentString> segments) {
		/*
		 * Other noder implementations do not node correctly (fail to detect
		 * intersections) on many inputs; furthermore, using a very small tolerance
		 * (i.e. ~1e-10) on SnappingNoder noder on a small tolerance misses
		 * intersections too (hence 0.01 chosen as suitable). "Noding robustness issues
		 * are generally caused by nearly coincident line segments, or by very short
		 * line segments. Snapping mitigates both of these situations.".
		 */
		final Noder noder = new SnappingNoder(0.01);
		noder.computeNodes(segments);
		return noder.getNodedSubstrings();
	}

	/**
	 * Provides convenient iteration of the child geometries of a JTS MultiGeometry.
	 * This iterator does not recurse all geometries (as does
	 * {@link org.locationtech.jts.geom.GeometryCollectionIterator
	 * GeometryCollectionIterator}), but returns the first level geometries only.
	 * 
	 * @author Michael Carleton
	 */
	static final class GeometryIterator implements Iterable<Geometry> {

		private final Geometry g;

		public GeometryIterator(Geometry g) {
			this.g = g;
		}

		@Override
		public Iterator<Geometry> iterator() {
			return new Iterator<Geometry>() {
				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < g.getNumGeometries();
				}

				@Override
				public Geometry next() {
					if (!hasNext()) {
						throw new NoSuchElementException();
					}
					return g.getGeometryN(currentIndex++);
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
		}
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

}
