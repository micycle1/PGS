package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS.prepareLinesPShape;

import java.util.ArrayList;
import java.util.Arrays;

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
 * PShape<-->JTS Geometry conversion
 * 
 * @author MCarleton
 *
 */
public class PGS_Conversion implements PConstants {

	private PGS_Conversion() {

	}

	/**
	 * TODO GROUP, PRIMITIVE, PATH, or GEOMETRY TODO CACHE recent 5 calls? TODO
	 * split into voronoi, delaunay, bool algebra classes Morph class: smooth
	 * simplify. etc.
	 * 
	 * USE arraylist.toArray() where possible
	 * 
	 * @param shape
	 * @return
	 */
	public static Geometry fromPShape(PShape shape) {

		// TODO convert to switch statement

		if (shape.getFamily() == PShape.GROUP) {
			ArrayList<PShape> flatChildren = new ArrayList<PShape>();
			getChildren(shape, flatChildren);
			flatChildren.removeIf(s -> s.getFamily() == PShape.GROUP); // aka .remove(shape)
//			flatChildren.remove(shape);
			Polygon[] children = new Polygon[flatChildren.size()];
			for (int i = 0; i < children.length; i++) {
				children[i] = (Polygon) fromPShape(flatChildren.get(i));
			}
			// TODO return to multipoly instead to prevent some crashes
			return (GEOM_FACTORY.createMultiPolygon(children).buffer(0)); // TODO don't flatten?

		}

//		shape.getKind() // switch to get primitive, then == ELLIPSE
		if (shape.getFamily() == PShape.PRIMITIVE) {
			GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
			shapeFactory.setNumPoints(PGS.CURVE_SAMPLES * 4); // TODO magic constant
			switch (shape.getKind()) {
				case ELLIPSE:
					// TODO split into createCircleGeom method
					shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
					shapeFactory.setWidth(shape.getParam(2));
					shapeFactory.setHeight(shape.getParam(3));
					return shapeFactory.createEllipse();
				case TRIANGLE:
//					shapeFactor
					// TODO
					break;
				case QUAD:
					// TODO 4-sided polygon
					break;
				case RECT:
					shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
					shapeFactory.setWidth(shape.getParam(2));
					shapeFactory.setHeight(shape.getParam(3));
					return shapeFactory.createRectangle();
//				      * @param a x-coordinate of the ellipse
//				      * @param b y-coordinate of the ellipse
//				      * @param c width of the ellipse by default
//				      * @param d height of the ellipse by default
//					break;

				default:
					System.err.print("Primitive Shape" + shape.getKind() + " not implmented");
					return null;
			}
		}

		/// GEOMETRY PShape types ///

		final int[] contourGroups = getContourGroups(shape.getVertexCodes());
		final int[] vertexCodes = getVertexTypes(shape);

		final ArrayList<ArrayList<Coordinate>> coords = new ArrayList<>(); // list of coords representing rings

		int lastGroup = -1;

		for (int i = 0; i < shape.getVertexCount(); i++) {
			if (contourGroups[i] != lastGroup) {
				lastGroup = contourGroups[i];
				coords.add(new ArrayList<>());
			}

			/**
			 * Sample bezier curves at intervals to produce smooth Geometry
			 */
			switch (vertexCodes[i]) {

				case QUADRATIC_VERTEX:
					coords.get(lastGroup).addAll(getQuadraticBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
							shape.getVertex(i + 1), PGS.CURVE_SAMPLES));
					i += 1;
					continue;

				case BEZIER_VERTEX: // aka cubic bezier, untested
					coords.get(lastGroup).addAll(getCubicBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
							shape.getVertex(i + 1), shape.getVertex(i + 2), PGS.CURVE_SAMPLES));
					i += 2;
					continue;

				default:
					coords.get(lastGroup).add(new Coordinate(shape.getVertexX(i), shape.getVertexY(i)));
					break;
			}
		}

		for (ArrayList<Coordinate> contour : coords) {
			// TODO only add if first and last different
			if (!contour.get(0).equals2D(contour.get(contour.size() - 1))) {
				contour.add(contour.get(0)); // Points of LinearRing must form a closed linestring
			}
		}

		final Coordinate[] outerCoords = new Coordinate[coords.get(0).size()];
		Arrays.setAll(outerCoords, coords.get(0)::get);

		LinearRing outer = null;
		if (outerCoords.length >= 4 || outerCoords.length == 0) {
			outer = GEOM_FACTORY.createLinearRing(outerCoords);
		}
		else {
//			System.out.println("coords: " + outerCoords.length);
		}

		/**
		 * Create linear ring for each hole in the shape
		 */
		LinearRing[] holes = new LinearRing[coords.size() - 1];

		for (int j = 1; j < coords.size(); j++) {
			final Coordinate[] innerCoords = new Coordinate[coords.get(j).size()];
			Arrays.setAll(innerCoords, coords.get(j)::get);
			holes[j - 1] = GEOM_FACTORY.createLinearRing(innerCoords);
		}

		return GEOM_FACTORY.createPolygon(outer, holes);
	}

	/**
	 * broken for P2D (due to createshape)
	 * 
	 * @param polygon
	 * @return
	 */
	public static PShape toPShape(Polygon polygon) {
		PShape shape = new PShape();
		shape.setFamily(PShape.GEOMETRY);
		shape.setFill(true);
		shape.setFill(255);

		if (polygon == null) {
			return shape;
		}

		shape.beginShape();

		/**
		 * Draw one outer ring. Both inner and outer loops used to iterate upto length
		 * -1 to skip the point that closes the JTS shape (same as the first point).
		 * However calling buffer() produces broken results for some shapes.
		 */
		Coordinate[] coords = polygon.getExteriorRing().getCoordinates();
		for (int i = 0; i < coords.length; i++) {
			// -1: ignore last coord (a copy of the first) (TODO check)
			Coordinate coord = coords[i];
			shape.vertex((float) coord.x, (float) coord.y);
		}

		/**
		 * Draw any Contours/inner rings
		 */
		for (int j = 0; j < polygon.getNumInteriorRing(); j++) {
			shape.beginContour();
			for (int i = 0; i < polygon.getInteriorRingN(j).getCoordinates().length; i++) {
				// -1: ignore last coord (a copy of the first)
				Coordinate coord = polygon.getInteriorRingN(j).getCoordinates()[i];
				shape.vertex((float) coord.x, (float) coord.y);
			}
			shape.endContour();
		}
		shape.endShape(); // NOTE don't close (since points should be closed)

		return shape;
	}

	/**
	 * A geometry may include multiple geometries, so resulting PShape may contain
	 * multiple children.
	 * 
	 * @return
	 */
	public static PShape toPShape(Geometry geometry) {
		
		if (geometry == null) {
			return new PShape(GROUP);
		}
		
		if (geometry.getNumGeometries() == 1) {
			if (geometry.getNumPoints() == 1) { // single point
				// TODO
				PShape point = new PShape();
				point.setFamily(PShape.GEOMETRY);
				point.setStrokeCap(ROUND);
				point.setStroke(true);
				point.setStrokeWeight(5);
				point.setStroke(-1232222);
				point.beginShape(POINTS);
				point.vertex((float) geometry.getCoordinate().x, (float) geometry.getCoordinate().y);
				point.endShape();
				return point;

			} else {
				if (geometry.getNumPoints() == 2) { // line
					PShape line = prepareLinesPShape(null, null, 4);
					line.vertex((float) geometry.getCoordinates()[0].x, (float) geometry.getCoordinates()[0].y);
					line.vertex((float) geometry.getCoordinates()[1].x, (float) geometry.getCoordinates()[1].y);
					line.endShape();
					return line;
				} else {
					if (geometry.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) { // long linestring
						PShape line = new PShape();
						line.setFamily(PShape.PATH); // TODO check?
						line.setStroke(true);
						line.setStrokeWeight(4);
						line.setStroke(-1232222);
						line.beginShape();
						Coordinate[] coords = geometry.getCoordinates();
						for (int i = 0; i < coords.length - 1; i++) {
							line.vertex((float) coords[i].x, (float) coords[i].y);
						}
						line.endShape(CLOSE);
						return line;
					} else {
						return toPShape((Polygon) geometry);
					}
				}

			}
		} else { // not very primitive
			PShape parent = new PShape(GROUP);
			parent.setFill(-16711936); // TODO
			for (int i = 0; i < geometry.getNumGeometries(); i++) {

				PShape child = null;
				if (geometry.getGeometryN(i).getNumGeometries() > 1) { // is geom collection
					child = toPShape(geometry.getGeometryN(i));
				} else {

					// TODO switch case on: geometry.getGeometryType()

					if (geometry.getGeometryN(i) instanceof LineString) {
						LineString l = (LineString) geometry.getGeometryN(i);

						// TODO remove if-else here?
//						System.out.println("todo linestring: " + geometry.getCoordinates().length);
						child = new PShape();
						child.setFamily(PShape.GEOMETRY);
//						child.setStrokeCap(ROUND);
						child.setStroke(true);
						child.setStrokeWeight(2);
						child.setStroke(-1232222);
						child.beginShape();
						for (int j = 0; j < l.getCoordinates().length; j++) {
							float vx = (float) l.getCoordinates()[j].x;
							float vy = (float) l.getCoordinates()[j].y;
							child.vertex(vx, vy);
						}
						child.endShape();
						parent.addChild(child);
					} else {

						if (geometry.getGeometryN(i).getCoordinates().length > 2) { // not Point or linestring
							child = toPShape((Polygon) geometry.getGeometryN(i));
							child.setFill(-16711936); // TODO
							child.setStrokeCap(ROUND);
							child.setStroke(true);
							child.setStrokeWeight(2);
//					child.disableStyle(); // Inherit parent style; causes crash?

						} else {
//						System.out.println("debug line");
							child = toPShape(geometry.getGeometryN(i)); // Point or LineString
						}
						parent.addChild(child);
					}
				}
			}
			return parent;
		}
	}

//	/**
//	 * Converts PShapes to Gishur polygons. Only supports GEOMETRY-type PShapes so
//	 * far. Doesn't support holes.
//	 * 
//	 * @param shape
//	 * @return
//	 */
//	public static XPolygon toXPolygon(PShape shape) {
//
//		final int[] contourGroups = getContourGroups(shape.getVertexCodes());
//		final int[] vertexCodes = getVertexTypes(shape);
//
//		final ArrayList<ArrayList<XPoint>> coords = new ArrayList<>(); // list of coords representing rings
//
//		int lastGroup = -1;
//
//		for (int i = 0; i < shape.getVertexCount(); i++) {
//			if (contourGroups[i] != lastGroup) {
//				lastGroup = contourGroups[i];
//				coords.add(new ArrayList<>());
//			}
//
//			/**
//			 * Sample bezier curves at intervals to produce smooth Geometry
//			 */
//			switch (vertexCodes[i]) {
//
////				case QUADRATIC_VERTEX:
////					coords.get(lastGroup).addAll(getQuadraticBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
////							shape.getVertex(i + 1), PTS.CURVE_SAMPLES));
////					i += 1;
////					continue;
////
////				case BEZIER_VERTEX: // aka cubic bezier, untested
////					coords.get(lastGroup).addAll(getCubicBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
////							shape.getVertex(i + 1), shape.getVertex(i + 2), PTS.CURVE_SAMPLES));
////					i += 2;
////					continue;
//
//				default:
//					coords.get(lastGroup).add(new XPoint(shape.getVertexX(i), shape.getVertexY(i)));
////					System.out.println(new XPoint(shape.getVertexX(i), shape.getVertexY(i)).toString());
//			}
//		}
//
//		for (ArrayList<XPoint> contour : coords) {
//			// xpolys are not closed via vertices
////			contour.add(contour.get(0)); // Points of LinearRing must form a closed linestring
//		}
//
//		coords.get(0).remove(coords.size() - 1); // pshape is closed, so remove last vertex
//		final XPoint[] outerCoords = new XPoint[coords.get(0).size()];
//		Arrays.setAll(outerCoords, coords.get(0)::get);
//
//		if (outerCoords.length >= 3 || outerCoords.length == 0) {
//			return new XPolygon(outerCoords);
//		} else {
////			System.out.println("coords: " + outerCoords.length);
//			return new XPolygon();
//		}
//	}

	public static int[] getContourGroups(int[] vertexCodes) {

		int group = 0;

		ArrayList<Integer> groups = new ArrayList<>(vertexCodes.length * 2);

		for (int vertexCode : vertexCodes) {
			switch (vertexCode) {
				case VERTEX:
					groups.add(group);
					break;

				case QUADRATIC_VERTEX:
					groups.add(group);
					groups.add(group);
					break;

				case BEZIER_VERTEX:
					groups.add(group);
					groups.add(group);
					groups.add(group);
					break;

				case CURVE_VERTEX:
					groups.add(group);
					break;

				case BREAK:
					// Marks beginning/end of new contour, and should be proceeded by a VERTEX
					group++;
					break;
			}
		}

		final int[] vertexGroups = new int[groups.size()];
		Arrays.setAll(vertexGroups, groups::get);
		return vertexGroups;
	}

	/**
	 * Kinda recursive, caller must provide fresh arraylist. Output includes the
	 * parent-most (input) shape. Output is flattened, does not respect a hierarchy
	 * of parent-child PShapes.
	 * 
	 * @param shape
	 * @param visited
	 * @return
	 */
	public static PShape getChildren(PShape shape, ArrayList<PShape> visited) {
		visited.add(shape);

		if (shape.getChildCount() == 0 || shape.getKind() != GROUP) {
			return shape;
		}

		for (PShape child : shape.getChildren()) {
			getChildren(child, visited);
		}
		return null;
	}

	/**
	 * Sets the color of the PShape and all of it's children recursively (and
	 * disables stroke).
	 * 
	 * @param shape
	 */
	public static void setAllFillColor(PShape shape, int color) {
		ArrayList<PShape> all = new ArrayList<PShape>();
		getChildren(shape, all);
		all.forEach(child -> {
			child.setStroke(false);
			child.setFill(true);
			child.setFill(color);
		});
	}

	public static void setAllStrokeColor(PShape shape, int color, int strokeWeight) {
		ArrayList<PShape> all = new ArrayList<PShape>();
		getChildren(shape, all);
		all.forEach(child -> {
			child.setStroke(true);
			child.setStroke(color);
			child.setStrokeWeight(strokeWeight);
		});
	}

	public static void disableAllFill(PShape shape) {
		ArrayList<PShape> all = new ArrayList<PShape>();
		getChildren(shape, all);
		all.forEach(child -> {
			child.setFill(false);
		});
	}

	/**
	 * Basically getVertexCodes, but returns the vertex type for every vertex
	 * 
	 * @param shape
	 * @return
	 */
	private static int[] getVertexTypes(PShape shape) {

		ArrayList<Integer> codes = new ArrayList<>(shape.getVertexCodeCount());

		for (int i = 0; i < shape.getVertexCodeCount(); i++) {
			int vertexCode = shape.getVertexCode(i);
			switch (vertexCode) {
				case VERTEX:
					codes.add(VERTEX);
					break;

				case QUADRATIC_VERTEX:
					codes.add(QUADRATIC_VERTEX);
					codes.add(QUADRATIC_VERTEX);
					break;

				case BEZIER_VERTEX:
					codes.add(BEZIER_VERTEX);
					codes.add(BEZIER_VERTEX);
					codes.add(BEZIER_VERTEX);
					break;

				case CURVE_VERTEX:
					codes.add(CURVE_VERTEX);
					break;

				case BREAK:
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
	 * @param resolution points per spline
	 * @return list of points along curve
	 */
	private static ArrayList<Coordinate> getQuadraticBezierPoints(PVector start, PVector controlPoint, PVector end,
			int resolution) {

		ArrayList<Coordinate> coords = new ArrayList<>();

		for (int j = 0; j < resolution; j++) {
			PVector bezierPoint = getQuadraticBezierCoordinate(start, controlPoint, end, j / (float) resolution);
			coords.add(new Coordinate(bezierPoint.x, bezierPoint.y));
		}
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

	private static ArrayList<Coordinate> getCubicBezierPoints(PVector start, PVector controlPoint1,
			PVector controlPoint2, PVector end, int resolution) {

		ArrayList<Coordinate> coords = new ArrayList<>();

		for (int j = 0; j < resolution; j++) {
			PVector bezierPoint = getCubicBezierCoordinate(start, controlPoint1, controlPoint2, end,
					j / (float) resolution);
			coords.add(new Coordinate(bezierPoint.x, bezierPoint.y));
		}
		return coords;
	}

	private static PVector getCubicBezierCoordinate(PVector start, PVector controlPoint1, PVector controlPoint2,
			PVector end, float t) {
		final float t1 = 1.0f - t;
		float x = start.x * t1 * t1 * t1 + 3 * controlPoint1.x * t * t1 * t1 + 3 * controlPoint2.x * t * t * t1
				+ end.x * t * t * t;
		float y = start.y * t1 * t1 * t1 + 3 * controlPoint1.y * t * t1 * t1 + 3 * controlPoint2.y * t * t * t1
				+ end.y * t * t * t;
		return new PVector(x, y);
	}
}