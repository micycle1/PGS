package micycle.pts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.log4j.Level;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.geotools.geometry.jts.JTS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.geomgraph.Edge;
import org.locationtech.jts.geomgraph.GeometryGraph;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.locationtech.jts.shape.random.RandomPointsBuilder;
import org.locationtech.jts.shape.random.RandomPointsInGridBuilder;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.triangulate.ConformingDelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;
import org.locationtech.jts.triangulate.quadedge.EdgeConnectedTriangleTraversal;
import org.locationtech.jts.util.GeometricShapeFactory;

import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineSegment;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineString;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IPolygon;
import fr.ign.cogit.geoxygene.contrib.delaunay.TriangleDelaunay;
import fr.ign.cogit.geoxygene.spatial.coordgeom.GM_Polygon;
import fr.ign.cogit.geoxygene.util.conversion.AdapterFactory;
import micycle.pts.color.Blending;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.concavehull.ConcaveHull;
import uk.osgb.algorithm.concavehull.TriCheckerChi;
import uk.osgb.algorithm.minkowski_sum.Minkowski_Sum;

/**
 * PTS | Processing Topology Suite
 * 
 * http://thecloudlab.org/processing/library.html
 * https://github.com/Rogach/jopenvoronoi https://github.com/IGNF/CartAGen
 * https://ignf.github.io/CartAGen/docs/algorithms/others/spinalize.html
 * https://github.com/twak/campskeleton
 * https://discourse.processing.org/t/straight-skeleton-or-how-to-draw-a-center-line-in-a-polygon-or-shape/17208/9
 * http://lastresortsoftware.blogspot.com/2010/12/smooth-as.html#note1
 * 
 * TODO https://ignf.github.io/CartAGen/docs/algorithms.html TODO take into
 * account strokweight (buffer by stroke amount?)
 * 
 * @author MCarleton
 *
 */
public class PTS implements PConstants {

	private static final int CURVE_SAMPLES = 20;

	protected static GeometryFactory geometryFactory = new GeometryFactory();

	static {
		// stop spina/skeleton console logging
		List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
		loggers.add(LogManager.getRootLogger());
		for (Logger logger : loggers) {
			logger.setLevel(Level.OFF);
		}
	}

	/**
	 * TODO GROUP, PRIMITIVE, PATH, or GEOMETRY TODO CACHE recent 5 calls? TODO
	 * split into voronoi, delaunay, bool algebra classes Morph class: smooth
	 * simplify. etc.
	 * 
	 * @param shape
	 * @return
	 */
	public static Polygon fromPShape(PShape shape) {

		// TODO convert to switch statement

		if (shape.getFamily() == PShape.GROUP) {
			ArrayList<PShape> flatChildren = new ArrayList<PShape>();
			getChildren(shape, flatChildren);
			flatChildren.removeIf(s -> s.getFamily() == PShape.GROUP); // aka .remove(shape)
//			flatChildren.remove(shape);
			Polygon[] children = new Polygon[flatChildren.size()];
			for (int i = 0; i < children.length; i++) {
				children[i] = fromPShape(flatChildren.get(i));
			}
			return (Polygon) (geometryFactory.createMultiPolygon(children).union()); // TODO don't flatten?
		}

//		shape.getKind() // switch to get primitive, then == ELLIPSE
		if (shape.getFamily() == PShape.PRIMITIVE) {
			GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
			shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
			switch (shape.getKind()) {
				case ELLIPSE:
					shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
					shapeFactory.setWidth(shape.getParam(2));
					shapeFactory.setHeight(shape.getParam(3));
					return shapeFactory.createEllipse();
				case TRIANGLE:
//					shapeFactor
					// TODO
					break;
				case QUAD:
					// TODO
					break;
				case RECT:
					// TODO
					break;
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
							shape.getVertex(i + 1), CURVE_SAMPLES));
					i += 1;
					continue;

				case BEZIER_VERTEX: // aka cubic bezier, untested
					coords.get(lastGroup).addAll(getCubicBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
							shape.getVertex(i + 1), shape.getVertex(i + 2), CURVE_SAMPLES));
					i += 2;
					continue;

				default:
					coords.get(lastGroup).add(new Coordinate(shape.getVertexX(i), shape.getVertexY(i)));
					break;
			}
		}

		for (ArrayList<Coordinate> contour : coords) {
			// TODO only add if first and last different
			contour.add(contour.get(0)); // Points of LinearRing must form a closed linestring
		}

		final Coordinate[] outerCoords = new Coordinate[coords.get(0).size()];
		Arrays.setAll(outerCoords, coords.get(0)::get);

		GeometryFactory geometryFactory = new GeometryFactory();
		LinearRing outer = geometryFactory.createLinearRing(outerCoords);

		/**
		 * Create linear ring for each hole in the shape
		 */
		LinearRing[] holes = new LinearRing[coords.size() - 1];

		for (int j = 1; j < coords.size(); j++) {
			final Coordinate[] innerCoords = new Coordinate[coords.get(j).size()];
			Arrays.setAll(innerCoords, coords.get(j)::get);
			holes[j - 1] = geometryFactory.createLinearRing(innerCoords);
		}

		return geometryFactory.createPolygon(outer, holes);
	}

	/**
	 * Uses other package version (for GeOxygene compatibility)
	 * 
	 * @param shape
	 * @return
	 */
	public static com.vividsolutions.jts.geom.Polygon fromPShapeVivid(PShape shape) {

//		shape.getKind() // switch to get primitive, then == ELLIPSE
		if (shape.getFamily() == PShape.PRIMITIVE) {
			com.vividsolutions.jts.util.GeometricShapeFactory shapeFactory = new com.vividsolutions.jts.util.GeometricShapeFactory();
			shapeFactory.setNumPoints(40); // TODO magic constant
			switch (shape.getKind()) {
				case ELLIPSE:
					shapeFactory.setCentre(
							new com.vividsolutions.jts.geom.Coordinate(shape.getParam(0), shape.getParam(1)));
					shapeFactory.setWidth(shape.getParam(2));
					shapeFactory.setHeight(shape.getParam(3));
					return shapeFactory.createEllipse();
				case TRIANGLE:
//					shapeFactor
					// TODO
					break;
				case QUAD:
					// TODO
					break;
				case RECT:
					// TODO
					break;
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

		// GEOMETRY PShape types:

		final int[] contourGroups = getContourGroups(shape.getVertexCodes());
		final int[] vertexCodes = getVertexTypes(shape);

		final ArrayList<ArrayList<com.vividsolutions.jts.geom.Coordinate>> coords = new ArrayList<>(); // list of coords
																										// representing
																										// rings

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
// TODO
//				case QUADRATIC_VERTEX:
//					coords.get(lastGroup).addAll(getQuadraticBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
//							shape.getVertex(i + 1), 20));
//					i += 1;
//					continue;
//
//				case BEZIER_VERTEX: // aka cubic bezier, untested
//					coords.get(lastGroup).addAll(getCubicBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
//							shape.getVertex(i + 1), shape.getVertex(i + 2), 20));
//					i += 2;
//					continue;

				default:
					coords.get(lastGroup)
							.add(new com.vividsolutions.jts.geom.Coordinate(shape.getVertexX(i), shape.getVertexY(i)));
					break;
			}
		}

		for (ArrayList<com.vividsolutions.jts.geom.Coordinate> contour : coords) {
			contour.add(contour.get(0)); // Points of LinearRing must form a closed linestring
		}

		final com.vividsolutions.jts.geom.Coordinate[] outerCoords = new com.vividsolutions.jts.geom.Coordinate[coords
				.get(0).size()];
		Arrays.setAll(outerCoords, coords.get(0)::get);

		com.vividsolutions.jts.geom.GeometryFactory geometryFactory = new com.vividsolutions.jts.geom.GeometryFactory();
		com.vividsolutions.jts.geom.LinearRing outer = geometryFactory.createLinearRing(outerCoords);

		/**
		 * Create linear ring for each hole in the shape
		 */
		com.vividsolutions.jts.geom.LinearRing[] holes = new com.vividsolutions.jts.geom.LinearRing[coords.size() - 1];

		for (int j = 1; j < coords.size(); j++) {
			final com.vividsolutions.jts.geom.Coordinate[] innerCoords = new com.vividsolutions.jts.geom.Coordinate[coords
					.get(j).size()];
			Arrays.setAll(innerCoords, coords.get(j)::get);
			holes[j - 1] = geometryFactory.createLinearRing(innerCoords);
		}

		return geometryFactory.createPolygon(outer, holes);
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

		shape.beginShape();

		/**
		 * Draw one outer ring. Both inner and outer loops used to iterate upto length
		 * -1 to skip the point that closes the JTS shape (same as the first point).
		 * However calling buffer() produces broken results for some shapes.
		 */
		for (int i = 0; i < polygon.getExteriorRing().getCoordinates().length; i++) {
			// -1: ignore last coord (a copy of the first)
			Coordinate coord = polygon.getExteriorRing().getCoordinates()[i];
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
		shape.endShape(CLOSE);

		return shape;
	}

	/**
	 * A geometry may include multiple geometries, so resulting PShape may contain
	 * multiple children.
	 * 
	 * @return
	 */
	public static PShape toPShape(Geometry geometry) {
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
					PShape line = new PShape();
					line.setFamily(PShape.GEOMETRY);
					line.setStrokeCap(ROUND);
					line.setStroke(true);
					line.setStrokeWeight(5);
					line.setStroke(-1232222);
					line.beginShape(LINES);
					line.vertex((float) geometry.getCoordinates()[0].x, (float) geometry.getCoordinates()[0].y);
					line.vertex((float) geometry.getCoordinates()[1].x, (float) geometry.getCoordinates()[1].y);
					line.endShape();
					return line;
				} else {
					return toPShape((Polygon) geometry);
				}

			}
		} else {
			PShape parent = new PShape(GROUP);
			parent.setFill(-16711936); // TODO
			for (int i = 0; i < geometry.getNumGeometries(); i++) {

				PShape child = null;
				if (geometry.getGeometryN(i).getNumGeometries() > 1) { // geom collection
					child = toPShape(geometry.getGeometryN(i));
				} else {

					if (geometry.getGeometryN(i).getCoordinates().length > 2) { // not Point or linestring
						child = toPShape((Polygon) geometry.getGeometryN(i));
						child.setFill(-16711936); // TODO
						child.setStrokeCap(ROUND);
						child.setStroke(true);
						child.setStrokeWeight(2);
//					child.disableStyle(); // Inherit parent style; causes crash?
						parent.addChild(child);
					}
				}
			}
			return parent;
		}
	}

	public static void triangulate(PShape shape) {
		// TODO
//		fromPShape(shape).getExteriorRing().get
//		org.locationtech.jts.triangulate.DelaunayTriangulationBuilder.
		// see delaunay method
	}

	/**
	 * 
	 * @param shape
	 * @param fit   0...1
	 * @return
	 */
	public static Polygon smooth(Polygon shape, float fit) {
		return (Polygon) JTS.smooth(shape, fit);
	}

	public static Geometry smooth(Geometry shape, float fit) {
		return JTS.smooth(shape, fit);
	}

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
		out.setFill(Blending.screen(getPShapeFillColor(a), getPShapeFillColor(b)));
		return out;
	}

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
	 * @return A∪B
	 */
	public static PShape union(PShape a, PShape b) {
		return toPShape(fromPShape(a).union(fromPShape(b)));
	}

	public static PShape buffer(PShape shape, float buffer) {
		// TODO read
		// https://locationtech.github.io/jts/javadoc/org/locationtech/jts/operation/buffer/BufferOp.html
		return toPShape(fromPShape(shape).buffer(buffer, 4));
	}

	/**
	 * Simplifies a PShape using the Douglas-Peucker algorithm.
	 * 
	 * @param shape
	 * @param distanceTolerance
	 * @return
	 */
	public static PShape simplify(PShape shape, float distanceTolerance) {
		return toPShape(DouglasPeuckerSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Smoothes a geometry
	 * 
	 * @param shape
	 * @param buffer 0...1
	 * @return
	 */
	public static PShape smooth(PShape shape, float fit) {
		return toPShape(smooth(fromPShape(shape), fit));
	}

	public static PShape flatten(PShape shape) {
		return toPShape(fromPShape(shape).getExteriorRing());
	}

	public static PShape convexHull(PShape... shapes) {
		Geometry g = fromPShape(shapes[0]);
		for (int i = 1; i < shapes.length; i++) {
			g = g.union(fromPShape(shapes[i]));
		}
		return toPShape(g.convexHull());
	}

	public static PShape concaveHull(ArrayList<PVector> points, float threshold) {

		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
			points.add(points.get(0)); // close geometry
		} else { // already closed
			coords = new Coordinate[points.size()];
		}

		for (int i = 0; i < coords.length; i++) {
			coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
		}

		Geometry g = geometryFactory.createPolygon(coords);
		ConcaveHull hull = new ConcaveHull(g);
		return toPShape(hull.getConcaveHullBFS(new TriCheckerChi(threshold), false, false).get(0));
	}

	/**
	 * Has a more "organic" structure compared to other concave method.
	 * 
	 * @param points
	 * @param threshold 0...1
	 * @return
	 */
	public static PShape concaveHull2(ArrayList<PVector> points, float threshold) {

		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
			points.add(points.get(0)); // close geometry
		} else { // already closed
			coords = new Coordinate[points.size()];
		}

		for (int i = 0; i < coords.length; i++) {
			coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
		}

		Geometry g = geometryFactory.createPolygon(coords);
		org.geodelivery.jap.concavehull.ConcaveHull hull = new org.geodelivery.jap.concavehull.ConcaveHull(threshold);

		return toPShape(hull.transform(g));
	}

	/**
	 * Minkowski addition a.k.a dilation
	 * 
	 * @return
	 */
	public static PShape minkSum(PShape source, PShape addition) {
//		Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
		Geometry sum = Minkowski_Sum.compMinkSum(fromPShape(source), fromPShape(addition), false, false);
		return toPShape(sum);
	}

	/**
	 * TODO check
	 * 
	 * @param source
	 * @param addition
	 * @return
	 */
	public static PShape minkDifference(PShape source, PShape addition) {
//			Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
		Geometry sum = Minkowski_Sum.compMinkDiff(fromPShape(source), fromPShape(addition), false, false);
		return toPShape(sum);
	}

	/**
	 * tolerance = 0
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape delaunayTriangulation(PShape shape) {
		Geometry g = fromPShape(shape);
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setSites(g);
		Geometry out = d.getTriangles(geometryFactory); // triangulates convex hull of points
		out = out.intersection(g); // get concave hull
		return toPShape(out);
	}

	private static Coordinate coordFromPoint(Point p) {
		return new Coordinate(p.getX(), p.getY());
	}

	/**
	 * cirumcircle center must lie inside triangle
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	public static double smallestSide(Coordinate a, Coordinate b, Coordinate c) {
		double ab = Math.sqrt((b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x));
		double bc = Math.sqrt((c.y - b.y) * (c.y - b.y) + (c.x - b.x) * (c.x - b.x));
		double ca = Math.sqrt((a.y - c.y) * (a.y - c.y) + (a.x - c.x) * (a.x - c.x));
		return Math.min(Math.min(ab, bc), ca);
	}

	/**
	 * Calc from vertices of PShape
	 * 
	 * @param shape
	 * @param tolerance
	 * @return
	 */
	public static PShape delaunayTriangulation(PShape shape, float tolerance) {

		/**
		 * Assume each pair of points on exterior ring to be segments,
		 */
		// http://lin-ear-th-inking.blogspot.com/2011/04/polygon-triangulation-via-ear-clipping.html
		Geometry g = fromPShape(shape);
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setTolerance(tolerance);
		d.setSites(g);
		Geometry out = d.getTriangles(geometryFactory); // triangulates concave hull of points

//		Coordinate

		ArrayList<Coordinate> coords = new ArrayList<>();
		// add from g
		for (int i = 0; i < g.getCoordinates().length; i++) {
			coords.add(g.getCoordinates()[i]);
		}

		for (int refinement = 0; refinement < 1; refinement++) {
			// refinement: add new points (centroids)
			for (int i = 0; i < out.getNumGeometries(); i++) {
				if (out.getGeometryN(i).getArea() > tolerance) {
					coords.add(coordFromPoint(out.getGeometryN(i).getCentroid()));
				}
			}
			d = new DelaunayTriangulationBuilder();
			d.setSites(coords);
			out = d.getTriangles(geometryFactory); // triangulates concave hull of points
		}

		// use d.getSubdivision() to get connected to a given segment to test
		// encroachment

		out = out.intersection(g); // get convex hull
		return toPShape(out);

//		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
//		b.setTolerance(tolerance);
//		b.setSites(g);
//		b.setConstraints(fromPShape(createSquircle(ThreadLocalRandom.current().nextInt(400, 450), 400, 200, 200)));

//		d.getSubdivision().getVoronoiCellPolygons(geometryFactory)
//		out = b.getTriangles(geometryFactory);
	}

	/**
	 * Create a triangulation with that forces certain required segments into the
	 * triangulation from a 'constrain' shape.
	 * 
	 * @param shape
	 * @param tolerance
	 * @return
	 */
	public static PShape constrainedDelaunayTriangulation(PShape shape, PShape constraints, float tolerance) {
		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
		Geometry g = fromPShape(shape);
		b.setSites(g);
		b.setTolerance(tolerance);
		b.setConstraints(fromPShape(constraints));
		Geometry out = b.getTriangles(geometryFactory); // triangulates concave hull of points
		out = out.intersection(g); // get convex hull
		return toPShape(out);
	}

	public static PShape constrainedDelaunayTriangulation(ArrayList<PVector> points, PShape constraints,
			float tolerance) {
		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
		Coordinate[] coords = new Coordinate[points.size()];
		for (int j = 0; j < points.size(); j++) {
			coords[j] = new Coordinate(points.get(j).x, points.get(j).y);
		}
		b.setSites(geometryFactory.createPolygon(coords));
		b.setTolerance(tolerance);
//		b.setConstraints(fromPShape(constraints));
		Geometry out = b.getTriangles(geometryFactory); // triangulates concave hull of points
		out = out.intersection(fromPShape(constraints)); // get convex hull
		return toPShape(out);
	}

	/**
	 * 
	 * @param points
	 * @param tolerance default is 10, snaps pixels to the nearest N. higher values
	 *                  result in more coarse triangulation
	 * @return A PShape whose children are triangles making up the triangulation
	 */
	public static PShape delaunayTriangulation(ArrayList<PVector> points, float tolerance) {
		ArrayList<Coordinate> coords = new ArrayList<>(points.size());
		for (PVector p : points) {
			coords.add(new Coordinate(p.x, p.y));
		}
		/**
		 * Use ConformingDelau... for better result? (get constrained segments from edge
		 * ring) Code example: https://github.com/locationtech/jts/issues/320
		 */
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setTolerance(tolerance);
		d.setSites(coords);
		Geometry out = d.getTriangles(geometryFactory);
		return toPShape(out);
	}

	/**
	 * TODO set clip envelope?
	 * 
	 * @param shape
	 * @param tolerance
	 * @return
	 */
	public static PShape voronoiDiagram(PShape shape, float tolerance) {
		Geometry g = fromPShape(shape);
		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setTolerance(tolerance);
		v.setSites(g);
		Geometry out = v.getDiagram(geometryFactory);
		return toPShape(out);
	}

	public static PShape voronoiDiagram(ArrayList<PVector> points, float tolerance) {
		ArrayList<Coordinate> coords = new ArrayList<>(points.size());
		for (PVector p : points) {
			coords.add(new Coordinate(p.x, p.y));
		}
		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setTolerance(tolerance);
//		v.setClipEnvelope(new Envelope(0, 0, 500, 500));
		v.setSites(coords);
		Geometry out = v.getDiagram(geometryFactory);
		return toPShape(out);
	}

	public static void prims(PShape shape) { // TODO
//		PrimMinimumSpanningTree<Vertex, Edge> prims = new PrimMinimumSpanningTree<>(null);
		GeometryGraph graph = new GeometryGraph(0, fromPShape(shape));
		for (Iterator<Edge> iterator = graph.getEdgeIterator(); iterator.hasNext();) {
			Edge type = iterator.next();
//		 System.out.println(type.getCoordinate().getX());
//		 fromPShape(shape).get
		}
//	Graph<Vertex, Edge> graph = new DefaultUndirectedGraph<Vertex, Edge>(Edge.class);
//	graph.
	}

	/**
	 * Not working: "Skeleton screw up"
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape spinalize(PShape shape) {
		// also see
		// https://github.com/IGNF/geoxygene/blob/master/geoxygene-spatial/src/main/java/fr/ign/cogit/geoxygene/util/algo/geometricAlgorithms/morphomaths/MorphologyTransform.java
//		GM_Polygon p = new GM_Polygon(fromPShape(null));
		try {
			GM_Polygon polygon = (GM_Polygon) AdapterFactory.toGM_Object(fromPShapeVivid(shape));

//			AdapterFactory.log
			List<IPolygon> l = new ArrayList<>();
			l.add(polygon);
			List<ILineString> lines = Spinalize.spinalize(l, 5, 100, false);
			System.out.println(lines.size());
			PShape skeleton = new PShape();
			skeleton.setFamily(PShape.GEOMETRY);
			skeleton.setStroke(true);
			skeleton.setStrokeWeight(3);
//			skeleton.setstrok
			skeleton.beginShape(PShape.LINES);
//			for (ILineString segment : lines) {
////				segment.get
////				System.out.println(segment.getStartPoint().getX());
//				skeleton.vertex((float) segment.getStartPoint().getX(), (float) segment.getStartPoint().getY());
//				skeleton.vertex((float) segment.getEndPoint().getX(), (float) segment.getEndPoint().getY());
//			}
			skeleton.endShape();
			return skeleton;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
//			return new ArrayList<ILineString>();
		}
//		JtsGeOxygene.makeGeOxygeneGeom(fromPShapeVivid(shape));
//		JTSGeomFactory f = new JTSGeomFactory();
//		f.createIPolygon((IRing) fromPShape(shape).getExteriorRing());
	}

	/**
	 * Works, but sloooow
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape skeletonize(PShape shape) {
		try {
			GM_Polygon polygon = (GM_Polygon) AdapterFactory.toGM_Object(fromPShapeVivid(shape));
			Set<ILineSegment> lines = Skeletonize.skeletonizeStraightSkeleton(polygon);
//			System.out.println(lines.size());
			PShape skeleton = new PShape();
			skeleton.setFamily(PShape.GEOMETRY);
			skeleton.setStroke(true);
			skeleton.setStrokeWeight(3);
//			skeleton.setstrok
			skeleton.beginShape(PShape.LINES);
			for (ILineSegment segment : lines) {
//				System.out.println(segment.getStartPoint().getX());
				skeleton.vertex((float) segment.getStartPoint().getX(), (float) segment.getStartPoint().getY());
				skeleton.vertex((float) segment.getEndPoint().getX(), (float) segment.getEndPoint().getY());
			}
			skeleton.endShape();
			return skeleton;
//			
//			for (ILineString line : lines) {
//				System.out.println(line.startPoint().getX());
//			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * TODO: outline of exterior ring only
	 * 
	 * @param shape
	 * @param points
	 * @param offsetDistance the distance the point is offset from the
	 *                       segment(positive is to the left, negative is to the
	 *                       right)
	 * @return
	 * @see #equidistantOutlineByDistance(PShape, float, float)
	 */
	public static PVector[] equidistantOutline(PShape shape, int points, float offsetDistance) {
		Polygon p = fromPShape(shape);

		LengthIndexedLine l = new LengthIndexedLine(p.getExteriorRing());
		PVector[] outlinePoints = new PVector[points * (p.getNumInteriorRing() + 1)];
		// exterior ring
		for (int i = 0; i < points; i++) {
			Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
			outlinePoints[i] = new PVector((float) q.x, (float) q.y);
		}

		for (int j = 0; j < p.getNumInteriorRing(); j++) {
			l = new LengthIndexedLine(p.getInteriorRingN(j));
			for (int i = 0; i < points; i++) {
				Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
				outlinePoints[(j + 1) * points + i] = new PVector((float) q.x, (float) q.y);
			}
		}

		return outlinePoints;
	}

	/**
	 * 
	 * @param shape
	 * @param interPointDistance Distance between each point on outline
	 * @param offsetDistance
	 * @return nearest distance such that every distance is equal (maybe be
	 *         different due to rounding)
	 */
	public static PVector[] equidistantOutlineByDistance(PShape shape, float interPointDistance, float offsetDistance) {
		Polygon p = fromPShape(shape);

		LengthIndexedLine l = new LengthIndexedLine(p.getExteriorRing());
		if (interPointDistance > l.getEndIndex()) {
			System.err.println("Interpoint greater than shape length");
			return null;
		}
		int points = (int) Math.round(l.getEndIndex() / interPointDistance);
		ArrayList<PVector> outlinePoints = new ArrayList<>();
		// exterior ring
		for (int i = 0; i < points; i++) {
			Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
			outlinePoints.add(new PVector((float) q.x, (float) q.y));
		}

		for (int j = 0; j < p.getNumInteriorRing(); j++) {
			l = new LengthIndexedLine(p.getInteriorRingN(j));
			points = (int) Math.round(l.getEndIndex() / interPointDistance);
			for (int i = 0; i < points; i++) {
				Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
				outlinePoints.add(new PVector((float) q.x, (float) q.y));
			}
		}

		final PVector[] out = new PVector[outlinePoints.size()];
		Arrays.setAll(out, outlinePoints::get);
		return out;
	}

	/**
	 * 
	 * @param shape
	 * @param distance 0...1 around shape perimeter; or -1...0 (other direction)
	 * @return
	 */
	public static PVector pointOnOutline(PShape shape, float distance, float offsetDistance) {
		distance %= 1;
		LengthIndexedLine l = new LengthIndexedLine(fromPShape(shape).getExteriorRing());
		Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
		return new PVector((float) coord.x, (float) coord.y);
	}

	public static PVector getCentroid(PShape shape) {
		Point point = fromPShape(shape).getCentroid();
		return new PVector((float) point.getX(), (float) point.getY());
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
	 * @param points
	 * @return
	 */
	public static ArrayList<PVector> generateRandomPoints(PShape shape, int points) {
		// or g.getInteriorPoint()?
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
	 * random points generated in a grid of cells (one point randomly located in
	 * each cell) from the envelope of the shape
	 * 
	 * @param shape
	 * @param points number of points, if this shape was its own envelope
	 * @return
	 */
	public static ArrayList<PVector> generateRandomGridPoints(PShape shape, int maxPoints) {
		Geometry g = fromPShape(shape);

		RandomPointsInGridBuilder r = new RandomPointsInGridBuilder();
		r.setConstrainedToCircle(true);
		r.setExtent(g.getEnvelopeInternal());
		r.setNumPoints(maxPoints);

		ArrayList<PVector> vertices = new ArrayList<>();

		for (Coordinate coord : r.getGeometry().getCoordinates()) {
			// manually prune envelope
			if (g.contains(geometryFactory.createPoint(coord))) {
				vertices.add(new PVector((float) coord.x, (float) coord.y));
			}
		}
		return vertices;
	}

	/**
	 * TODO return list of points when shape is group
	 * 
	 * @param shape
	 * @param point
	 * @return
	 */
	public static PVector closestVertexToPoint(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		Coordinate coord = DistanceOp.nearestPoints(g, pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

	public static boolean contains(PShape outer, PShape inner) {
		return fromPShape(outer).contains(fromPShape(inner));
	}

	public static boolean containsPoint(PShape shape, PVector point) {
		return fromPShape(shape).contains(pointFromPVector(point));
	}

	public static float distance(PShape a, PShape b) {
		return (float) fromPShape(a).distance(fromPShape(b));
	}

	public static Point getPoint(float x, float y) {
		return geometryFactory.createPoint(new Coordinate(x, y));
	}

	public static float area(PShape shape) {
		return (float) fromPShape(shape).getArea();
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
	 * @return X,Y,W,H
	 */
	public static float[] bound(PShape shape) {
		Envelope e = (Envelope) fromPShape(shape).getEnvelopeInternal();
		return new float[] { (float) e.getMinX(), (float) e.getMinY(), (float) e.getWidth(), (float) e.getHeight() };
	}

	/**
	 * aka envelope
	 * 
	 * @param shape
	 * @return X1, Y1, X2, Y2
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
	public static PShape createSquircle(float x, float y, float width, float height) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(40);
		shapeFactory.setCentre(new Coordinate(x, y));
//		shapeFactory.setBase(new Coordinate(x, y));
//		shapeFactory.setRotation(1);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createSquircle());
	}

	/**
	 * https://stackoverflow.com/questions/64252638/how-to-split-a-jts-polygon Split
	 * a polygon into 4 equal quadrants
	 * 
	 * @param p
	 * @return
	 */
	public static ArrayList<PShape> split(PShape shape) {
		Polygon p = fromPShape(shape);
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
	 * Sets the color of the PShape and all of it's children recursively
	 * 
	 * @param shape
	 */
	public static void setAllFillColor(PShape shape, int color) {
		ArrayList<PShape> all = new ArrayList<PShape>();
		getChildren(shape, all);
		all.forEach(child -> {
			child.setFill(true);
			child.setFill(color);
		});
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
	private static PShape getChildren(PShape shape, ArrayList<PShape> visited) {
		visited.add(shape);

		if (shape.getChildCount() == 0) {
			return shape;
		}

		for (PShape child : shape.getChildren()) {
			getChildren(child, visited);
		}
		return null;
	}

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

	private static Point pointFromPVector(PVector p) {
		return geometryFactory.createPoint(new Coordinate(p.x, p.y));
	}

}
