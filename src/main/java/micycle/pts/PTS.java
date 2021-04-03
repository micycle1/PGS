package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.geotools.geometry.jts.JTS;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.shape.random.RandomPointsBuilder;
import org.locationtech.jts.shape.random.RandomPointsInGridBuilder;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pts.color.RGB;
import micycle.pts.utility.PolygonDecomposition;
import micycle.pts.utility.RandomPolygon;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * PTS | Processing Topology Suite
 * <p>
 * <ul>
 * <li>https://github.com/IGNF/CartAGen
 * <li>https://ignf.github.io/CartAGen/docs/algorithms/others/spinalize.html
 * <li>https://discourse.processing.org/t/straight-skeleton-or-how-to-draw-a-center-line-in-a-polygon-or-shape/17208/9
 * </p>
 * </ul>
 * 
 * 
 * @author Michael Carleton
 */
public class PTS {

	// TODO check for getCoordinates() in loops (and replace) (if lots of child
	// geometries)
	// TODO replace line() in for-each with beginshape(LINES)...vertex...endShape()
	// TODO use LinearRingIterator when possible (refactor)
	// TODO https://ignf.github.io/CartAGen/docs/algorithms.html
	// TODO take into account strokweight (buffer by stroke amount?)
	// TODO maintain pshape fill, etc on output

	/**
	 * Calling Polygon#union repeatedly is one way to union several Polygons
	 * together. But here’s a trick that can be significantly faster (seconds rather
	 * than minutes) – add the Polygons to a GeometryCollection, then apply a buffer
	 * with zero distance
	 */

	/**
	 * Defines number of vertex samples per bezier line (during PShape->JTS)
	 */
	protected static final int CURVE_SAMPLES = 20;

	public static GeometryFactory GEOM_FACTORY = new GeometryFactory(
			new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	/**
	 * @return shape A - shape B
	 */
	public static PShape difference(PShape a, PShape b) {
		return toPShape(fromPShape(a).difference(fromPShape(b)));
	}

	/**
	 * Compute the parts that the shapes do not have in common.
	 * 
	 * @return A∪B - A∩B
	 */
	public static PShape symDifference(PShape a, PShape b) {
		return toPShape(fromPShape(a).symDifference(fromPShape(b)));
	}

	/**
	 * @return A∩B
	 */
	public static PShape intersection(PShape a, PShape b) {
		PShape out = toPShape(fromPShape(a).intersection(fromPShape(b)));
//		a.draw(p.getGraphics());
//		out.setFill(Blending.screen(getPShapeFillColor(a), getPShapeFillColor(b))); // TODO
		return out;
	}

	/**
	 * @return A∪B
	 * @see #union(PShape...)
	 */
	public static PShape union(PShape a, PShape b) {
//		OverlayNG.overlay(null, null, CURVE_SAMPLES) // TODO investigate
		return toPShape(fromPShape(a).union(fromPShape(b)));
	}

	public static PShape union(PShape... shapes) {
		// same as flatten?
		ArrayList<Geometry> geoms = new ArrayList<>();
		for (int i = 0; i < shapes.length; i++) {
			geoms.add(fromPShape(shapes[i]));
		}
		return toPShape(UnaryUnionOp.union(geoms));
	}

	/**
	 * Densifies a Geometry by inserting extra vertices along the line segments
	 * contained in the geometry. Specify maximum length of segments
	 */
	public static PShape densify(PShape shape, float distanceTolerance) {
		Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(distanceTolerance);
		d.setValidate(false);
		return toPShape(d.getResultGeometry());
	}

	/**
	 * Smoothes a geometry. The smoothing algorithm inserts new vertices which are
	 * positioned using Bezier splines.
	 * 
	 * @param shape
	 * @param fit   tightness of fit from 0 (loose) to 1 (tight)
	 * @return
	 */
	public static PShape smooth(PShape shape, float fit) {
		return toPShape(JTS.smooth(fromPShape(shape), fit));
	}

	/**
	 * Flatten, or merge/union, a shape's children
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape flatten(PShape shape) {
		// TODO iterate over geometries shapes, then create group PShape
		Polygon poly = (Polygon) fromPShape(shape).union().getGeometryN(0);
		return toPShape(poly.getExteriorRing());
	}

	/**
	 * Extracts a point from the perimeter of the given shape.
	 * 
	 * @param shape
	 * @param distance       0...1 around shape perimeter; or -1...0 (other
	 *                       direction)
	 * @param offsetDistance perpendicular offset distance, where 0 is exactly on
	 *                       the shape perimeter. The computed point is offset to
	 *                       the left of the line if the offset distance is
	 *                       positive, and to the right if negative
	 * @return
	 * @see #pointsOnPerimeter(PShape, int, float)
	 */
	public static PVector pointOnPerimeter(PShape shape, float distance, float offsetDistance) {
		distance %= 1;
		// TODO CHECK CAST (ITERATE OVER GROUP); apply to interior rings too?
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Extracts many points from the perimeter (faster than calling other method
	 * lots)
	 * 
	 * @param shape
	 * @param points         number of points to return; evenly distibuted around
	 *                       the perimeter of the shape
	 * @param offsetDistance offset distance along a line perpendicular to the
	 *                       perimeter
	 * @return
	 * @see #pointOnPerimeter(PShape, float, float)
	 * @see #pointsOnPerimeter(PShape, float, float)
	 */
	public static List<PVector> pointsOnPerimeter(PShape shape, int points, float offsetDistance) {
		// TODO another method that returns concave hull of returned points (when
		// offset)
		ArrayList<PVector> coords = new ArrayList<>(points);
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		final double increment = 1d / points;
		for (double distance = 0; distance < 1; distance += increment) {
			Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
			coords.add(new PVector((float) coord.x, (float) coord.y));
		}
		return coords;
	}

	/**
	 * 
	 * @param shape
	 * @param interPointDistance distance between each point on outline
	 * @param offsetDistance
	 * @return
	 */
	public static List<PVector> pointsOnPerimeter(PShape shape, float interPointDistance, float offsetDistance) {
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		if (interPointDistance > l.getEndIndex()) {
			System.err.println("Interpoint length greater than shape length");
			return new ArrayList<PVector>();
		}
		final int points = (int) Math.round(l.getEndIndex() / interPointDistance);

		ArrayList<PVector> coords = new ArrayList<>(points);

		final double increment = 1d / points;
		for (double distance = 0; distance < 1; distance += increment) {
			Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
			coords.add(new PVector((float) coord.x, (float) coord.y));
		}
		return coords;
	}

	/**
	 * Returns a hole-less version of the shape, or boundary of group of shapes
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape boundary(PShape shape) {
		return toPShape(fromPShape(shape).getBoundary());
	}

	/**
	 * Get N random points contained within the PShape region. Points are
	 * distributed randomly.
	 * 
	 * @param shape
	 * @param points number of points to generate
	 * @return
	 */
	public static ArrayList<PVector> generateRandomPoints(PShape shape, int points) {
		RandomPointsBuilder r = new RandomPointsBuilder();
		r.setExtent(fromPShape(shape));
		r.setNumPoints(points);

		ArrayList<PVector> vertices = new ArrayList<>();

		for (Coordinate coord : r.getGeometry().getCoordinates()) {
			vertices.add(new PVector((float) coord.x, (float) coord.y));
		}
		return vertices;
	}

	/**
	 * Random points generated in a grid of cells (one point randomly located in
	 * each cell) from the envelope of the shape
	 * 
	 * @param shape
	 * @param maxPoints           number of points, if this shape was its own
	 *                            envelope
	 * @param constrainedToCircle
	 * 
	 *                            Sets whether generated points are constrained to
	 *                            liewithin a circle contained within each grid
	 *                            cell.This provides greater separation between
	 *                            pointsin adjacent cells.
	 * @param gutterFraction      Sets the fraction of the grid cell side which will
	 *                            be treated asa gutter, in which no points will be
	 *                            created.The provided value is clamped to the range
	 *                            [0.0, 1.0].
	 * 
	 * @return
	 */
	public static ArrayList<PVector> generateRandomGridPoints(PShape shape, int maxPoints, boolean constrainedToCircle,
			double gutterFraction) {
		Geometry g = fromPShape(shape);
		IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);

		RandomPointsInGridBuilder r = new RandomPointsInGridBuilder();
		r.setConstrainedToCircle(constrainedToCircle);
		r.setExtent(g.getEnvelopeInternal());
		r.setNumPoints(maxPoints);
		r.setGutterFraction(gutterFraction);

		ArrayList<PVector> vertices = new ArrayList<>();

		for (Coordinate coord : r.getGeometry().getCoordinates()) {
			if (pointLocator.locate(coord) != Location.EXTERIOR) {
				vertices.add(new PVector((float) coord.x, (float) coord.y));
			}
		}
		return vertices;
	}

	/**
	 * Returns a copy of the shape with its small holes (i.e. inner rings with area
	 * < given threshold) removed.
	 * 
	 * @param polygon
	 * @return
	 */
	public static PShape removeSmallHoles(PShape shape, float area) {
		Polygon polygon = (Polygon) fromPShape(shape);
		Polygon noHolePol = GEOM_FACTORY.createPolygon(polygon.getExteriorRing());
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			LinearRing hole = polygon.getInteriorRingN(i);
			if (hole.getArea() < area)
				continue;
			noHolePol = (Polygon) noHolePol.difference(hole);

		}
		return toPShape(noHolePol);
	}

	/**
	 * 
	 * @param n    number of vertices
	 * @param xMax
	 * @param yMax
	 * @return
	 */
	public static PShape randomPolygon(int n, float xMax, float yMax) {
		return pshapeFromPVector(RandomPolygon.generateRandomConvexPolygon(n, xMax, yMax));
	}

	/**
	 * Same as getVertices for geometry PShapes; different (subdivides) for circles
	 * etc.
	 * 
	 * @param shape
	 * @return
	 */
	public static PVector[] vertices(PShape shape) {
		Coordinate[] coords = fromPShape(shape).getCoordinates();
		PVector[] vertices = new PVector[coords.length];
		for (int i = 0; i < coords.length; i++) {
			Coordinate coord = coords[i];
			vertices[i] = new PVector((float) coord.x, (float) coord.y);
		}
		return vertices;
	}

	/**
	 * aka envelope
	 * 
	 * @param shape
	 * @return float[] of [X,Y,W,H]
	 */
	public static float[] bounds(PShape shape) {
		// TODO move to ShapeMetrics?
		Envelope e = (Envelope) fromPShape(shape).getEnvelopeInternal();
		return new float[] { (float) e.getMinX(), (float) e.getMinY(), (float) e.getWidth(), (float) e.getHeight() };
	}

	/**
	 * aka envelope
	 * 
	 * @param shape
	 * @return float[] of [X1, Y1, X2, Y2]
	 */
	public static float[] boundCoords(PShape shape) {
		Envelope e = (Envelope) fromPShape(shape).getEnvelopeInternal();
		return new float[] { (float) e.getMinX(), (float) e.getMinY(), (float) e.getMaxX(), (float) e.getMaxY() };
	}

	/**
	 * 
	 * @param x      centre X
	 * @param y      center Y
	 * @param width  ?total? width
	 * @param height ?total? height
	 * @return
	 */
	public static PShape createSquircle(double x, double y, double width, double height) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 2);
		shapeFactory.setCentre(new Coordinate(x, y));
//		shapeFactory.setBase(new Coordinate(x, y));
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createSquircle());
	}

	/**
	 * Creates an elliptical arc polygon. The polygon is formed from the specified
	 * arc of an ellipse and the two radii connecting the endpoints to the centre of
	 * the ellipse.
	 * 
	 * @param x
	 * @param y
	 * @param width
	 * @param height
	 * @param angle  size of angle in radians
	 * @return
	 */
	public static PShape createArcPolygon(double x, double y, double width, double height, double angle) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 2);
		shapeFactory.setCentre(new Coordinate(x, y));
//		shapeFactory.setBase(new Coordinate(x, y));
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createArcPolygon(0, angle));
	}

	/**
	 * Splits a polygon into 4 equal quadrants
	 * 
	 * @param shape
	 * @return list containing the 4 split quadrants
	 */
	public static List<PShape> split(PShape shape) {
		// https://stackoverflow.com/questions/64252638/how-to-split-a-jts-polygon
		Geometry p = fromPShape(shape);
		ArrayList<PShape> ret = new ArrayList<>();

		final Envelope envelope = p.getEnvelopeInternal();
		double minX = envelope.getMinX();
		double maxX = envelope.getMaxX();
		double midX = minX + (maxX - minX) / 2.0;
		double minY = envelope.getMinY();
		double maxY = envelope.getMaxY();
		double midY = minY + (maxY - minY) / 2.0;

		Envelope llEnv = new Envelope(minX, midX, minY, midY);
		Envelope lrEnv = new Envelope(midX, maxX, minY, midY);
		Envelope ulEnv = new Envelope(minX, midX, midY, maxY);
		Envelope urEnv = new Envelope(midX, maxX, midY, maxY);
		Geometry ll = JTS.toGeometry(llEnv).intersection(p);
		Geometry lr = JTS.toGeometry(lrEnv).intersection(p);
		Geometry ul = JTS.toGeometry(ulEnv).intersection(p);
		Geometry ur = JTS.toGeometry(urEnv).intersection(p);
		ret.add(toPShape(ll));
		ret.add(toPShape(lr));
		ret.add(toPShape(ul));
		ret.add(toPShape(ur));

		return ret;
	}

	/**
	 * Partitions a shape into simple polygons using Mark Bayazit's algorithm.
	 * 
	 * @param shape
	 * @return list of convex (simple) polygons comprising the original shape
	 */
	public static ArrayList<PShape> partition(PShape shape) {
		// https://mpen.ca/406/bayazit
		// retry GreedyPolygonSplitter()?

		Geometry g = fromPShape(shape);

		ArrayList<PShape> out = new ArrayList<>();

		for (int i = 0; i < g.getNumGeometries(); i++) {
			Geometry child = g.getGeometryN(i);
			if (child.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) { // skip any linestrings etc
				List<Polygon> decomposed = PolygonDecomposition.decompose((Polygon) child);
				for (Polygon polygon : decomposed) {
					out.add(toPShape(polygon));
				}
			}
		}

		return out;
	}

	/**
	 * Slice a shape using a line given by its start and endpoints
	 * 
	 * @param shape
	 * @param p1    must be outside shape
	 * @param p2    must be outside shape
	 * @return Two PShapes
	 */
	public static ArrayList<PShape> slice(PShape shape, PVector p1, PVector p2) {
		// https://gis.stackexchange.com/questions/189976/
		Geometry poly = fromPShape(shape);
		LineString line = lineFromPVectors(p1, p2);
		Geometry nodedLinework = poly.getBoundary().union(line);
		Geometry polys = polygonize(nodedLinework);

		// Only keep polygons which are inside the input
		ArrayList<Polygon> slices = new ArrayList<>();
		for (int i = 0; i < polys.getNumGeometries(); i++) {
			Polygon candpoly = (Polygon) polys.getGeometryN(i);
			// TODO geometry cache for contains check?
			if (poly.contains(candpoly.getInteriorPoint())) {
				slices.add(candpoly);
			}
		}

		ArrayList<PShape> output = new ArrayList<>();

		for (Polygon polygon : slices) {
			output.add(toPShape(polygon));
		}
		return output;
	}

	public static List<PVector> toPVectorList(PShape shape) {
		final ArrayList<PVector> vertices = new ArrayList<>();
		for (int i = 0; i < shape.getVertexCount(); i++) {
			vertices.add(shape.getVertex(i));
		}
		return vertices;
	}

	/**
	 * Create a LINES PShape, ready for vertices.
	 * 
	 * @param strokeColor  nullable
	 * @param strokeCap    nullable default = ROUND
	 * @param strokeWeight nullable. default = 2
	 * @return
	 */
	static PShape prepareLinesPShape(Integer strokeColor, Integer strokeCap, Integer strokeWeight) {
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
	 * Used by slice()
	 */
	@SuppressWarnings("unchecked")
	private static Geometry polygonize(Geometry geometry) {
		List<LineString> lines = LineStringExtracter.getLines(geometry);
		Polygonizer polygonizer = new Polygonizer();
		polygonizer.add(lines);
		Collection<Polygon> polys = polygonizer.getPolygons();
		Polygon[] polyArray = GeometryFactory.toPolygonArray(polys);
		return geometry.getFactory().createGeometryCollection(polyArray);
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

	private static void removeCollinearVertices(Geometry g) {
		JTS.removeCollinearVertices(g);
	}

	static Point createPoint(float x, float y) {
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
	 * cirumcircle center must lie inside triangle
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	private static double smallestSide(Coordinate a, Coordinate b, Coordinate c) {
		double ab = Math.sqrt((b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x));
		double bc = Math.sqrt((c.y - b.y) * (c.y - b.y) + (c.x - b.x) * (c.x - b.x));
		double ca = Math.sqrt((a.y - c.y) * (a.y - c.y) + (a.x - c.x) * (a.x - c.x));
		return Math.min(Math.min(ab, bc), ca);
	}

	private static LineString lineFromPVectors(PVector a, PVector b) {
		return GEOM_FACTORY.createLineString(new Coordinate[] { new Coordinate(a.x, a.y), new Coordinate(b.x, b.y) });
	}

	/**
	 * Reflection-based workaround to get the fill color of a PShape (this field is
	 * usually private).
	 */
	private static final int getPShapeFillColor(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("fillColor");
			f.setAccessible(true);
			return f.getInt(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * Generate a simple polygon (no holes) from the given coordinate list. Used by
	 * randomPolygon().
	 */
	public static PShape pshapeFromPVector(List<PVector> coords) {
		PShape shape = new PShape();
		shape.setFamily(PShape.GEOMETRY);
		shape.setStroke(true);
		shape.setStrokeWeight(2);
		shape.setStroke(0);
		shape.setFill(-123712);
		shape.setFill(true);
		shape.beginShape();

		for (PVector v : coords) {
			shape.vertex(v.x, v.y);
		}

		shape.endShape(PConstants.CLOSE);
		return shape;
	}

	/**
	 * From com.badlogic.gdx
	 */
	private static boolean isClockwise(float[] polygon, int offset, int count) {
		if (count <= 2) {
			return false;
		}
		float area = 0;
		int last = offset + count - 2;
		float x1 = polygon[last], y1 = polygon[last + 1];
		for (int i = offset; i <= last; i += 2) {
			float x2 = polygon[i], y2 = polygon[i + 1];
			area += x1 * y2 - x2 * y1;
			x1 = x2;
			y1 = y2;
		}
		return area < 0;
	}

	/**
	 * Requires a closed hole
	 * 
	 * @param points
	 * @return
	 */
	static boolean isClockwise(List<PVector> points) {
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
			Iterator<LinearRing> it = new Iterator<LinearRing>() {

				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < size;
				}

				@Override
				public LinearRing next() {
					return array[currentIndex++];
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
			return it;
		}
	}

}
