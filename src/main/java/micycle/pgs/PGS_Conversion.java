package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.util.GeometricShapeFactory;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Conversion between Processing PShapes and JTS Geometries.
 * <p>
 * Methods in this class are used by the library but are kept accessible for
 * more advanced user use cases.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Conversion implements PConstants {

	/** Approximate distance between successive sample points on bezier curves */
	private static final float BEZIER_SAMPLE_DISTANCE = 2;

	private PGS_Conversion() {
	}

	/**
	 * Converts a JTS Geometry to an equivalent PShape. MultiGeometries (collections
	 * of geometries) become GROUP PShapes with children shapes.
	 * 
	 * @param g
	 * @param source PShape to copy fill/stroke details from
	 * @return
	 */
	public static PShape toPShape(final Geometry g) {
		if (g == null) {
			return new PShape(PShape.GEOMETRY);
		}

		PShape shape = new PShape();
		shape.setFill(true);
		shape.setFill(micycle.pgs.color.RGB.WHITE);
		shape.setStroke(true);
		shape.setStroke(micycle.pgs.color.RGB.PINK);
		shape.setStrokeWeight(4);

		switch (g.getGeometryType()) {
			case Geometry.TYPENAME_GEOMETRYCOLLECTION :
			case Geometry.TYPENAME_MULTIPOLYGON :
			case Geometry.TYPENAME_MULTILINESTRING :
				shape.setFamily(GROUP);
				for (int i = 0; i < g.getNumGeometries(); i++) {
					shape.addChild(toPShape(g.getGeometryN(i)));
				}
				break;
			case Geometry.TYPENAME_LINEARRING : // linearrings are closed by definition
			case Geometry.TYPENAME_LINESTRING :
				final LineString l = (LineString) g;
				final boolean closed = l.isClosed();
				shape.setFamily(PShape.PATH);
				shape.beginShape();
				Coordinate[] coords = l.getCoordinates();
				for (int i = 0; i < coords.length - (closed ? 1 : 0); i++) {
					shape.vertex((float) coords[i].x, (float) coords[i].y);
				}
				if (closed) { // closed vertex was skipped, so close the path
					shape.endShape(CLOSE);
				} else {
					shape.endShape();
				}
				break;
			case Geometry.TYPENAME_POLYGON :
				final Polygon polygon = (Polygon) g;
				shape.setFamily(PShape.PATH);
				shape.beginShape();

				/**
				 * Outer and inner loops are iterated up to length-1 to skip the point that
				 * closes the JTS shape (same as the first point).
				 */
				coords = polygon.getExteriorRing().getCoordinates();
				for (int i = 0; i < coords.length - 1; i++) {
					final Coordinate coord = coords[i];
					shape.vertex((float) coord.x, (float) coord.y);
				}

				for (int j = 0; j < polygon.getNumInteriorRing(); j++) { // holes
					shape.beginContour();
					coords = polygon.getInteriorRingN(j).getCoordinates();
					for (int i = 0; i < coords.length - 1; i++) {
						final Coordinate coord = coords[i];
						shape.vertex((float) coord.x, (float) coord.y);
					}
					shape.endContour();
				}
				shape.endShape(CLOSE);
				break;
			case Geometry.TYPENAME_POINT :
			case Geometry.TYPENAME_MULTIPOINT :
				coords = g.getCoordinates();
				shape.setFamily(PShape.GEOMETRY);
				shape.setStrokeCap(PConstants.ROUND);
				shape.beginShape(PShape.POINTS);
				for (int i = 0; i < coords.length; i++) {
					final Coordinate coord = coords[i];
					shape.vertex((float) coord.x, (float) coord.y);
				}
				shape.endShape();
				break;
			default :
				System.err.println("PGS_Conversion Error: " + g.getGeometryType() + " geometry types are unsupported.");
				break;
		}

		return shape;
	}

	/**
	 * Converts a collection of JTS Geometries to an equivalent GROUP PShape.
	 */
	public static PShape toPShape(Collection<Geometry> geometries) {
		System.out.println(geometries.size());
		PShape shape = new PShape(GROUP);
		shape.setFill(true);
		shape.setFill(micycle.pgs.color.RGB.WHITE);
		shape.setStroke(true);
		shape.setStroke(micycle.pgs.color.RGB.PINK);
		shape.setStrokeWeight(4);

		geometries.forEach(g -> shape.addChild(toPShape(g)));

		return shape;
	}

	/**
	 * Converts a PShape to an equivalent JTS Geometry.
	 * <p>
	 * PShapes with bezier curves are sampled at regular intervals (in which case
	 * the resulting geometry will have more vertices than the input)
	 * <p>
	 * For now, a PShape with multiple children is flattened/unioned since most
	 * library methods are not (yet) programmed to handle multi/disjoint geometries.
	 * 
	 * @param shape
	 * @return a JTS Polygon or MultiPolygon
	 */
	public static Geometry fromPShape(PShape shape) {

		Geometry g = null;

		switch (shape.getFamily()) {
			case PConstants.GROUP :
				final List<PShape> flatChildren = new ArrayList<>();
				getChildren(shape, flatChildren);
				flatChildren.removeIf(s -> s.getFamily() == PConstants.GROUP);
				if (flatChildren.isEmpty()) {
					return GEOM_FACTORY.createEmpty(2);
				}

				if (flatChildren.stream().allMatch(c -> c.getFamily() == PShape.PATH && c.isClosed() == false)) {
					// NOTE special case (for now): preserve paths as MultiLineString
					LineString[] children = new LineString[flatChildren.size()];
					for (int i = 0; i < children.length; i++) {
						children[i] = (LineString) fromPShape(flatChildren.get(i));
					}
					return GEOM_FACTORY.createMultiLineString(children);
				} else {
					List<Polygon> children = new ArrayList<>();
					for (int i = 0; i < flatChildren.size(); i++) {
						Geometry child = fromPShape(flatChildren.get(i));
						if (child.getGeometryType() == Geometry.TYPENAME_POLYGON
								|| child.getGeometryType() == Geometry.TYPENAME_LINEARRING) {
							children.add((Polygon) child);
						}
					}
					// NOTE for now, buffer/flatten multiple polygons into a single JTS polygon so
					// that methods handle them properly
					return (GEOM_FACTORY.createMultiPolygon(children.toArray(new Polygon[0])).buffer(0));
				}
			case PShape.GEOMETRY :
			case PShape.PATH :
				if (shape.getKind() == PConstants.POLYGON || shape.getKind() == PConstants.PATH || shape.getKind() == 0) {
					g = fromVertices(shape);
				} else {
					g = fromCreateShape(shape); // special paths (e.g. POINTS, LINES, etc.)
				}
				break;
			case PShape.PRIMITIVE :
				g = fromPrimitive(shape);
				break;
		}

		return g;
	}

	/**
	 * Converts a PShape made via beginShape(KIND), where KIND is not POLYGON.
	 * 
	 * @param shape
	 * @return
	 */
	private static Geometry fromCreateShape(PShape shape) {
		switch (shape.getKind()) {
			case PConstants.POINTS :
				final Coordinate[] coords = new Coordinate[shape.getVertexCount()];
				for (int i = 0; i < shape.getVertexCount(); i++) {
					coords[i] = PGS.coordFromPVector(shape.getVertex(i));
				}
				return GEOM_FACTORY.createMultiPointFromCoords(coords);
			case PConstants.LINES : // create multi line string consisting of each line
				final LineString[] lines = new LineString[shape.getVertexCount() / 2];
				for (int i = 0; i < lines.length; i++) {
					final Coordinate c1 = PGS.coordFromPVector(shape.getVertex(2 * i));
					final Coordinate c2 = PGS.coordFromPVector(shape.getVertex(2 * i + 1));
					lines[i] = GEOM_FACTORY.createLineString(new Coordinate[] { c1, c2 });
				}
				return GEOM_FACTORY.createMultiLineString(lines);
			default :
				System.err.println("PGS_Conversion Error: Unsupported PShape kind: " + shape.getKind());
				return GEOM_FACTORY.createEmpty(2);
		}
	}

	/**
	 * Creates a JTS Polygon from a geometry or path PShape, whose 'kind' is a
	 * polygon or path.
	 */
	private static Geometry fromVertices(PShape shape) {

		if (shape.getVertexCount() < 2) { // skip empty / point PShapes
			return GEOM_FACTORY.createPolygon();
		}
		
		final int[] contourGroups = getContourGroups(shape.getVertexCodes());
		final int[] vertexCodes = getVertexTypes(shape);

		final ArrayList<ArrayList<Coordinate>> coords = new ArrayList<>(); // list of coords representing rings

		int lastGroup = -1;
		for (int i = 0; i < shape.getVertexCount(); i++) {
			if (contourGroups[i] != lastGroup) {
				if (lastGroup == -1 && contourGroups[0] > 0) {
					lastGroup = 0;
				}
				lastGroup = contourGroups[i];
				coords.add(new ArrayList<>());
			}

			/**
			 * Sample bezier curves at regular intervals to produce smooth Geometry
			 */
			switch (vertexCodes[i]) { // VERTEX, BEZIER_VERTEX, CURVE_VERTEX, or BREAK
				case QUADRATIC_VERTEX :
					coords.get(lastGroup).addAll(getQuadraticBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
							shape.getVertex(i + 1), BEZIER_SAMPLE_DISTANCE));
					i += 1;
					continue;
				case BEZIER_VERTEX : // aka cubic bezier
					coords.get(lastGroup).addAll(getCubicBezierPoints(shape.getVertex(i - 1), shape.getVertex(i), shape.getVertex(i + 1),
							shape.getVertex(i + 2), BEZIER_SAMPLE_DISTANCE));
					i += 2;
					continue;
				default :
					coords.get(lastGroup).add(new Coordinate(shape.getVertexX(i), shape.getVertexY(i)));
					break;
			}
		}

		for (ArrayList<Coordinate> contour : coords) {
			final Iterator<Coordinate> iterator = contour.iterator();
			if (iterator.hasNext()) { // at least one vertex
				Coordinate previous = iterator.next();
				final List<Coordinate> duplicates = new ArrayList<>();

				while (iterator.hasNext()) { // find adjacent matching coordinates
					Coordinate current = iterator.next();
					if (current.equals2D(previous)) {
						duplicates.add(current);
					}
					previous = current;
				}
				contour.removeAll(duplicates); // remove adjacent matching coordinates

				if (!contour.get(0).equals2D(contour.get(contour.size() - 1)) && shape.isClosed()) {
					contour.add(contour.get(0)); // points of LinearRing must form a closed linestring
				}
			}
		}

		final Coordinate[] outerCoords = coords.get(0).toArray(new Coordinate[coords.get(0).size()]);
		
		if (outerCoords.length < 2) {
			return GEOM_FACTORY.createPolygon();
		}
		
		if (shape.isClosed() && outerCoords.length > 3) { // closed geometry or path
			LinearRing outer = GEOM_FACTORY.createLinearRing(outerCoords); // should always be valid

			LinearRing[] holes = new LinearRing[coords.size() - 1]; // Create linear ring for each hole in the shape
			for (int j = 1; j < coords.size(); j++) {
				final Coordinate[] innerCoords = coords.get(j).toArray(new Coordinate[coords.get(j).size()]);
				holes[j - 1] = GEOM_FACTORY.createLinearRing(innerCoords);
			}
			return GEOM_FACTORY.createPolygon(outer, holes);
		} else { // non-closed path (a string)
			return GEOM_FACTORY.createLineString(outerCoords);
		}
	}

	/**
	 * Creates a JTS Polygon from a primitive PShape. Primitive PShapes are those
	 * where createShape() is used to create them, and can take any of these types:
	 * POINT, LINE, TRIANGLE, QUAD, RECT, ELLIPSE, ARC, BOX, SPHERE. They do not
	 * have directly accessible vertex data.
	 */
	private static Polygon fromPrimitive(PShape shape) {
		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		switch (shape.getKind()) {
			case ELLIPSE :
				final double a = shape.getParam(2) / 2d;
				final double b = shape.getParam(3) / 2d;
				final double perimeter = Math.PI * (3 * (a + b) - Math.sqrt((3 * a + b) * (a + 3 * b)));
				shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
				shapeFactory.setWidth(a * 2);
				shapeFactory.setHeight(b * 2);
				shapeFactory.setNumPoints((int) (perimeter / BEZIER_SAMPLE_DISTANCE));
				return shapeFactory.createEllipse();
			case TRIANGLE :
				Coordinate[] coords = new Coordinate[3 + 1];
				Coordinate c1 = new Coordinate(shape.getParam(0), shape.getParam(1));
				coords[0] = c1;
				coords[1] = new Coordinate(shape.getParam(2), shape.getParam(3));
				coords[2] = new Coordinate(shape.getParam(4), shape.getParam(5));
				coords[3] = c1.copy(); // close loop
				return GEOM_FACTORY.createPolygon(coords);
			case RECT :
				shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
				shapeFactory.setWidth(shape.getParam(2));
				shapeFactory.setHeight(shape.getParam(3));
				return shapeFactory.createRectangle();
			case QUAD :
				coords = new Coordinate[4 + 1];
				c1 = new Coordinate(shape.getParam(0), shape.getParam(1));
				coords[0] = c1;
				coords[1] = new Coordinate(shape.getParam(2), shape.getParam(3));
				coords[2] = new Coordinate(shape.getParam(4), shape.getParam(5));
				coords[3] = new Coordinate(shape.getParam(6), shape.getParam(7));
				coords[4] = c1.copy(); // close loop
				return GEOM_FACTORY.createPolygon(coords);
			case ARC :
				shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
				shapeFactory.setWidth(shape.getParam(2));
				shapeFactory.setHeight(shape.getParam(3));
				// circumference (if it was full circle)
				final double circumference = Math.PI * Math.max(shape.getParam(2), shape.getParam(3));
				shapeFactory.setNumPoints((int) (circumference / BEZIER_SAMPLE_DISTANCE));
				return shapeFactory.createArcPolygon(-Math.PI / 2 + shape.getParam(4), shape.getParam(5));
			case LINE :
			case POINT :
				System.err.print("Non-polygon primitives are not supported.");
				break;
			case BOX :
			case SPHERE :
				System.err.print("3D primitives are not supported.");
				break;
			default :
				System.err.print(shape.getKind() + " primitives are not supported.");

		}
		return GEOM_FACTORY.createPolygon(); // empty polygon
	}

	/**
	 * Returns the vertices of a PShape as an unclosed list of PVector coordinates.
	 * 
	 * @param shape
	 * @return
	 */
	public static List<PVector> toPVector(PShape shape) {
		final ArrayList<PVector> vertices = new ArrayList<>();
		for (int i = 0; i < shape.getVertexCount(); i++) {
			vertices.add(shape.getVertex(i));
		}
		if (!vertices.isEmpty() && vertices.get(0).equals(vertices.get(vertices.size() - 1))) {
			vertices.remove(vertices.size() - 1);
		}
		return vertices;
	}

	/**
	 * Generate a simple polygon (no holes) from the given coordinate list. Used by
	 * randomPolygon().
	 */
	public static PShape fromPVector(List<PVector> coords) {
		PShape shape = new PShape();
		shape.setFamily(PShape.GEOMETRY);
		shape.setFill(micycle.pgs.color.RGB.WHITE);
		shape.setFill(true);
		shape.beginShape();

		for (PVector v : coords) {
			shape.vertex(v.x, v.y);
		}

		shape.endShape(PConstants.CLOSE);
		return shape;
	}

	/**
	 * Finds and returns all the children PShapes of a given PShape. All children
	 * (including the parent-most (input) shape) are put into the given list.
	 * <p>
	 * The output is flattened -- it does not respect a hierarchy of parent-child
	 * PShapes.
	 * 
	 * @param shape
	 * @param childrenOut user-provided list into which all child PShapes are placed
	 * @return
	 */
	public static PShape getChildren(PShape shape, List<PShape> childrenOut) {
		childrenOut.add(shape);

		if (shape.getChildCount() == 0 || shape.getKind() != GROUP) {
			return shape;
		}

		for (PShape child : shape.getChildren()) {
			getChildren(child, childrenOut);
		}
		return null;
	}

	/**
	 * Sets the fill color for the PShape and all of it's children recursively (and
	 * disables stroke).
	 * 
	 * @param shape
	 * @see #setAllStrokeColor(PShape, int, int)
	 */
	public static void setAllFillColor(PShape shape, int color) {
		List<PShape> all = new ArrayList<>();
		getChildren(shape, all);
		all.forEach(child -> {
			child.setStroke(false);
			child.setFill(true);
			child.setFill(color);
		});
	}

	/**
	 * Sets the stroke color for the PShape and all of it's children recursively.
	 * 
	 * @param shape
	 * @see {@link #setAllFillColor(PShape, int)}
	 */
	public static void setAllStrokeColor(PShape shape, int color, int strokeWeight) {
		List<PShape> all = new ArrayList<>();
		getChildren(shape, all);
		all.forEach(child -> {
			child.setStroke(true);
			child.setStroke(color);
			child.setStrokeWeight(strokeWeight);
		});
	}

	/**
	 * Calls setFill(false) on a PShape and all its children. This method mutates
	 * the input shape.
	 * 
	 * @param shape
	 */
	public static void disableAllFill(PShape shape) {
		ArrayList<PShape> all = new ArrayList<>();
		getChildren(shape, all);
		all.forEach(child -> child.setFill(false));
	}

	/**
	 * Calls setStrokefalse) on a PShape and all its children. This method mutates
	 * the input shape.
	 * 
	 * @param shape
	 */
	public static void disableAllStroke(PShape shape) {
		ArrayList<PShape> all = new ArrayList<>();
		getChildren(shape, all);
		all.forEach(child -> child.setStroke(false));
	}

	/**
	 * For every vertexcode, store the group (i.e. hole) it belongs to.
	 * 
	 * @param vertexCodes
	 * @return
	 */
	private static int[] getContourGroups(int[] vertexCodes) {

		int group = 0;

		ArrayList<Integer> groups = new ArrayList<>(vertexCodes.length * 2);

		for (int i = 0; i < vertexCodes.length; i++) {
			final int vertexCode = vertexCodes[i];
			switch (vertexCode) {
				case VERTEX :
					groups.add(group);
					break;

				case QUADRATIC_VERTEX :
					groups.add(group);
					groups.add(group);
					break;

				case BEZIER_VERTEX :
					groups.add(group);
					groups.add(group);
					groups.add(group);
					break;

				case CURVE_VERTEX :
					groups.add(group);
					break;

				case BREAK :
					// BREAK marks beginning/end of new contour, and should be proceeded by a VERTEX
					if (i > 0) {
						// In P2D, svg-loaded shapes begin with a break (so we don't want to increment)
						group++;
					}
					break;
				default :
					System.err.println("Unrecognised vertex code: " + vertexCode);
			}
		}

		final int[] vertexGroups = new int[groups.size()];
		Arrays.setAll(vertexGroups, groups::get);
		return vertexGroups;
	}

	/**
	 * Basically getVertexCodes, but returns the vertex type for every vertex
	 * 
	 * @param shape
	 * @return
	 */
	private static int[] getVertexTypes(PShape shape) {

		List<Integer> codes = new ArrayList<>(shape.getVertexCodeCount());

		for (int i = 0; i < shape.getVertexCodeCount(); i++) {
			int vertexCode = shape.getVertexCode(i);
			switch (vertexCode) {
				case VERTEX :
					codes.add(VERTEX);
					break;

				case QUADRATIC_VERTEX :
					codes.add(QUADRATIC_VERTEX);
					codes.add(QUADRATIC_VERTEX);
					break;

				case BEZIER_VERTEX :
					codes.add(BEZIER_VERTEX);
					codes.add(BEZIER_VERTEX);
					codes.add(BEZIER_VERTEX);
					break;

				case CURVE_VERTEX :
					codes.add(CURVE_VERTEX);
					break;

				case BREAK :
					break;

				default :
					System.err.println("Unrecognised vertex code: " + vertexCode);
			}
		}

		final int[] vertexGroups = new int[codes.size()];
		Arrays.setAll(vertexGroups, codes::get);
		return vertexGroups;
	}

	/**
	 * Subdivide/interpolate/discretise along a quadratic bezier curve, given by its
	 * start, end and control points
	 * 
	 * @return list of points along curve
	 */
	private static List<Coordinate> getQuadraticBezierPoints(PVector start, PVector controlPoint, PVector end, float sampleDistance) {
		List<Coordinate> coords = new ArrayList<>();

		if (start.dist(end) <= sampleDistance) {
			coords.add(new Coordinate(start.x, start.y));
			coords.add(new Coordinate(end.x, end.y));
			return coords;
		}

		final float length = bezierLengthQuadratic(start, controlPoint, end);
		final int samples = (int) Math.ceil(length / sampleDistance); // sample every x unit length (approximately)

		coords.add(new Coordinate(start.x, start.y));
		for (int j = 1; j < samples; j++) { // start at 1 -- don't sample at t=0
			final PVector bezierPoint = getQuadraticBezierCoordinate(start, controlPoint, end, j / (float) samples);
			coords.add(new Coordinate(bezierPoint.x, bezierPoint.y));
		}
		coords.add(new Coordinate(end.x, end.y));

		return coords;
	}

	/**
	 * 
	 * @param start
	 * @param controlPoint
	 * @param end
	 * @param t            0...1
	 * @return
	 */
	private static PVector getQuadraticBezierCoordinate(PVector start, PVector controlPoint, PVector end, float t) {
		float x = (1 - t) * (1 - t) * start.x + 2 * (1 - t) * t * controlPoint.x + t * t * end.x;
		float y = (1 - t) * (1 - t) * start.y + 2 * (1 - t) * t * controlPoint.y + t * t * end.y;
		return new PVector(x, y);
	}

	/**
	 * Approximate bezier length using Gravesen's approach. The insight is that the
	 * actual bezier length is always somewhere between the distance between the
	 * endpoints (the length of the chord) and the perimeter of the control polygon.
	 * For a quadratic Bézier, 2/3 the first + 1/3 the second is a reasonably good
	 * estimate.
	 * 
	 * @return
	 */
	private static float bezierLengthQuadratic(PVector start, PVector controlPoint, PVector end) {
		// https://raphlinus.github.io/curves/2018/12/28/bezier-arclength.html
		final float chord = PVector.sub(end, start).mag();
		final float cont_net = PVector.sub(start, controlPoint).mag() + PVector.sub(end, controlPoint).mag();
		return (2 * chord + cont_net) / 3f;

	}

	/**
	 * Generates a list of samples of a cubic bezier curve.
	 * 
	 * @param sampleDistance distance between successive samples on the curve
	 * @return
	 */
	private static List<Coordinate> getCubicBezierPoints(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end,
			float sampleDistance) {
		List<Coordinate> coords = new ArrayList<>();

		if (start.dist(end) <= sampleDistance) {
			coords.add(new Coordinate(start.x, start.y));
			coords.add(new Coordinate(end.x, end.y));
			return coords;
		}

		final float length = bezierLengthCubic(start, controlPoint1, controlPoint2, end);
		final int samples = (int) Math.ceil(length / sampleDistance); // sample every x unit length (approximately)

		coords.add(new Coordinate(start.x, start.y));
		for (int j = 1; j < samples; j++) { // start at 1 -- don't sample at t=0
			final PVector bezierPoint = getCubicBezierCoordinate(start, controlPoint1, controlPoint2, end, j / (float) samples);
			coords.add(new Coordinate(bezierPoint.x, bezierPoint.y));
		}
		coords.add(new Coordinate(end.x, end.y));
		return coords;
	}

	private static PVector getCubicBezierCoordinate(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end, float t) {
		final float t1 = 1.0f - t;
		float x = start.x * t1 * t1 * t1 + 3 * controlPoint1.x * t * t1 * t1 + 3 * controlPoint2.x * t * t * t1 + end.x * t * t * t;
		float y = start.y * t1 * t1 * t1 + 3 * controlPoint1.y * t * t1 * t1 + 3 * controlPoint2.y * t * t * t1 + end.y * t * t * t;
		return new PVector(x, y);
	}

	/**
	 * Approximate bezier length using Gravesen's approach. The insight is that the
	 * actual bezier length is always somewhere between the distance between the
	 * endpoints (the length of the chord) and the perimeter of the control polygon.
	 * 
	 * @return
	 */
	private static float bezierLengthCubic(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end) {
		// https://stackoverflow.com/a/37862545/9808792
		final float chord = PVector.sub(end, start).mag();
		final float cont_net = PVector.sub(start, controlPoint1).mag() + PVector.sub(controlPoint2, controlPoint1).mag()
				+ PVector.sub(end, controlPoint2).mag();
		return (cont_net + chord) / 2;

	}
}
