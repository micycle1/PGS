package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS.coordFromPVector;
import static micycle.pgs.color.RGB.decomposeclrRGB;
import static processing.core.PConstants.GROUP;
import static processing.core.PConstants.QUADRATIC_VERTEX;
import static processing.core.PConstants.BEZIER_VERTEX;
import static processing.core.PConstants.CURVE_VERTEX;

import java.awt.Shape;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import org.jgrapht.alg.drawing.IndexedFRLayoutAlgorithm2D;
import org.jgrapht.alg.drawing.LayoutAlgorithm2D;
import org.jgrapht.alg.drawing.model.Box2D;
import org.jgrapht.alg.drawing.model.LayoutModel2D;
import org.jgrapht.alg.drawing.model.MapLayoutModel2D;
import org.jgrapht.alg.drawing.model.Point2D;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.awt.ShapeReader;
import org.locationtech.jts.awt.ShapeWriter;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKBReader;
import org.locationtech.jts.io.WKBWriter;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.io.WKTWriter;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pgs.commons.PEdge;
import processing.core.PConstants;
import processing.core.PMatrix;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Conversion between <i>Processing</i> <code>PShapes</code> and <i>JTS</i>
 * <code>Geometries</code> (amongst other formats). Also includes helper/utility
 * methods for PShapes.
 * <p>
 * Methods in this class are used by the library internally but are kept
 * accessible for more advanced user use cases.
 * <p>
 * Notably, JTS geometries do not support bezier curves so any bezier curves are
 * finely subdivided into straight linestrings during <code>PShape</code> -> JTS
 * <code>Geometry</code> conversion.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Conversion {

	/** Approximate distance between successive sample points on bezier curves */
	private static final float BEZIER_SAMPLE_DISTANCE = 2;
	private static Field MATRIX_FIELD;
	/**
	 * Boolean flag that affects whether a PShape's style (fillColor, strokeColor,
	 * strokeWidth) is preserved during PShape->Geometry->PShape conversion (i.e.
	 * when <code>toPShape(fromPShape(myPShape))</code> is called). Default = true.
	 */
	public static boolean PRESERVE_STYLE = true;

	static {
		try {
			MATRIX_FIELD = PShape.class.getDeclaredField("matrix");
			MATRIX_FIELD.setAccessible(true);
		} catch (NoSuchFieldException e) {
			System.err.println(e.getLocalizedMessage());
		}
	}

	private PGS_Conversion() {
	}

	/**
	 * Converts a JTS Geometry to an equivalent PShape. MultiGeometries (collections
	 * of geometries) become GROUP PShapes containing the appropriate children
	 * PShapes.
	 * 
	 * @param g JTS geometry to convert
	 * @return
	 */
	public static PShape toPShape(final Geometry g) {
		if (g == null) {
			return new PShape();
		}
		PShape shape = new PShape();
		// apply PGS style by default
		shape.setFill(true);
		shape.setFill(micycle.pgs.color.RGB.WHITE);
		shape.setStroke(true);
		shape.setStroke(micycle.pgs.color.RGB.PINK);
		shape.setStrokeWeight(4);

		switch (g.getGeometryType()) {
			case Geometry.TYPENAME_GEOMETRYCOLLECTION :
			case Geometry.TYPENAME_MULTIPOLYGON :
			case Geometry.TYPENAME_MULTILINESTRING :
				if (g.getNumGeometries() == 1) {
					shape = toPShape(g.getGeometryN(0));
				} else {
					shape.setFamily(GROUP);
					for (int i = 0; i < g.getNumGeometries(); i++) {
						shape.addChild(toPShape(g.getGeometryN(i)));
					}
				}
				break;
			case Geometry.TYPENAME_LINEARRING : // LinearRings are closed by definition
			case Geometry.TYPENAME_LINESTRING : // LineStrings may be open
				final LineString l = (LineString) g;
				final boolean closed = l.isClosed();
				shape.setFamily(PShape.PATH);
				shape.beginShape();
				Coordinate[] coords = l.getCoordinates();
				for (int i = 0; i < coords.length - (closed ? 1 : 0); i++) {
					shape.vertex((float) coords[i].x, (float) coords[i].y);
				}
				if (closed) { // closed vertex was skipped, so close the path
					shape.endShape(PConstants.CLOSE);
				} else {
					// shape is more akin to an unconnected line: keep as PATH shape, but don't fill
					// visually
					shape.endShape();
					shape.setFill(false);
				}
				break;
			case Geometry.TYPENAME_POLYGON :
				final Polygon polygon = (Polygon) g;
				shape.setFamily(PShape.PATH);
				shape.beginShape();

				/*
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
				shape.endShape(PConstants.CLOSE);
				break;
			case Geometry.TYPENAME_POINT :
			case Geometry.TYPENAME_MULTIPOINT :
				coords = g.getCoordinates();
				shape.setFamily(PShape.GEOMETRY);
				shape.setFill(false);
				shape.setStrokeCap(PConstants.ROUND);
				shape.beginShape(PConstants.POINTS);
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

		if (PRESERVE_STYLE && g.getUserData() != null && g.getUserData() instanceof PShapeData) {
			PShapeData style = (PShapeData) g.getUserData();
			style.applyTo(shape);
		}

		return shape;
	}

	/**
	 * Converts a collection of JTS Geometries to an equivalent GROUP PShape. If the
	 * collection contains only one geometry, an equivalent PShape will be output
	 * directly (not a GROUP shape).
	 */
	public static PShape toPShape(Collection<? extends Geometry> geometries) {
		PShape shape = new PShape(GROUP);
		shape.setFill(true);
		shape.setFill(micycle.pgs.color.RGB.WHITE);
		shape.setStroke(true);
		shape.setStroke(micycle.pgs.color.RGB.PINK);
		shape.setStrokeWeight(4);

		geometries.forEach(g -> shape.addChild(toPShape(g)));
		if (shape.getChildCount() == 1) {
			return shape.getChild(0);
		}

		return shape;
	}

	/**
	 * Converts a PShape to an equivalent JTS Geometry.
	 * <p>
	 * PShapes with bezier curves are sampled at regular intervals (in which case
	 * the resulting geometry will have more vertices than the input PShape).
	 * 
	 * @param shape
	 * @return a JTS Geometry equivalent to the input PShape
	 */
	public static Geometry fromPShape(PShape shape) {
		Geometry g = GEOM_FACTORY.createEmpty(2);

		switch (shape.getFamily()) {
			case PConstants.GROUP :
				final List<PShape> flatChildren = getChildren(shape);
				List<Geometry> geoChildren = flatChildren.stream().map(PGS_Conversion::fromPShape).collect(Collectors.toList());
				g = GEOM_FACTORY.buildGeometry(geoChildren);
				break;
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

		if (PRESERVE_STYLE && g != null) {
			g.setUserData(new PShapeData(shape));
		}

		/*
		 * Finally, apply PShape's affine transformations (which are not applied to its
		 * vertices directly).
		 */
		try {
			final PMatrix matrix = (PMatrix) MATRIX_FIELD.get(shape);
			if (matrix != null) { // is null if no affine transformations have been applied to shape
				final float[] affine = matrix.get(null);
				if (affine.length == 6) { // process 2D shape matrix only
					AffineTransformation t = new AffineTransformation(affine[0], affine[1], affine[2], affine[3], affine[4], affine[5]);
					return t.transform(g);
				}
			}
		} catch (Exception e) {
		}

		return g;
	}

	/**
	 * Converts a PShape made via beginShape(KIND), where KIND is not POLYGON, to
	 * its equivalent JTS Geometry.
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
			case PConstants.TRIANGLES :
				final Polygon[] triangles = new Polygon[shape.getVertexCount() / 3];
				for (int i = 0; i < triangles.length; i++) {
					final Coordinate c1 = PGS.coordFromPVector(shape.getVertex(3 * i));
					final Coordinate c2 = PGS.coordFromPVector(shape.getVertex(3 * i + 1));
					final Coordinate c3 = PGS.coordFromPVector(shape.getVertex(3 * i + 2));
					triangles[i] = GEOM_FACTORY.createPolygon(new Coordinate[] { c1, c2, c3, c1 });
				}
				return GEOM_FACTORY.createMultiPolygon(triangles);
			case PConstants.QUADS :
				final Polygon[] quads = new Polygon[shape.getVertexCount() / 4];
				for (int i = 0; i < quads.length; i++) {
					final Coordinate c1 = PGS.coordFromPVector(shape.getVertex(4 * i));
					final Coordinate c2 = PGS.coordFromPVector(shape.getVertex(4 * i + 1));
					final Coordinate c3 = PGS.coordFromPVector(shape.getVertex(4 * i + 2));
					final Coordinate c4 = PGS.coordFromPVector(shape.getVertex(4 * i + 3));
					quads[i] = GEOM_FACTORY.createPolygon(new Coordinate[] { c1, c2, c3, c4, c1 });
				}
				return GEOM_FACTORY.createMultiPolygon(quads);
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

		int[] rawVertexCodes = shape.getVertexCodes();
		/*
		 * getVertexCodes() is null for P2D shapes with no contours created via
		 * createShape(), so need to instantiate array here.
		 */
		if (rawVertexCodes == null) {
			rawVertexCodes = new int[shape.getVertexCount()];
			Arrays.fill(rawVertexCodes, PConstants.VERTEX);
		}

		final int[] contourGroups = getContourGroups(rawVertexCodes);
		final int[] vertexCodes = getVertexTypes(rawVertexCodes);

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
				default : // VERTEX
					coords.get(lastGroup).add(coordFromPVector(shape.getVertex(i)));
					break;
			}
		}

		for (ArrayList<Coordinate> contour : coords) {
			final Iterator<Coordinate> iterator = contour.iterator();
			if (iterator.hasNext()) { // has at least one vertex
				final List<Coordinate> contourNoDupes = new ArrayList<>(contour.size());
				Coordinate previous = iterator.next();
				contourNoDupes.add(previous);

				/*
				 * Remove consecutive duplicate coordinates
				 */
				while (iterator.hasNext()) {
					Coordinate current = iterator.next();
					if (!current.equals2D(previous)) {
						contourNoDupes.add(current);
					}
					previous = current;
				}

				// mutate contour list
				contour.clear();
				contour.addAll(contourNoDupes);

				if (!contour.get(0).equals2D(contour.get(contour.size() - 1)) && shape.isClosed()) {
					contour.add(contour.get(0)); // close LinearRing: "points of LinearRing must form a closed linestring"
				}
			}
		}

		final Coordinate[] outerCoords = coords.get(0).toArray(new Coordinate[coords.get(0).size()]);

		if (outerCoords.length == 0) {
			return GEOM_FACTORY.createPolygon(); // empty polygon
		} else if (outerCoords.length == 1) {
			return GEOM_FACTORY.createPoint(outerCoords[0]);
		} else if (outerCoords.length == 2) {
			return GEOM_FACTORY.createLineString(outerCoords);
		} else if (shape.isClosed()) { // closed geometry or path
			LinearRing outer = GEOM_FACTORY.createLinearRing(outerCoords); // should always be valid
			LinearRing[] holes = new LinearRing[coords.size() - 1]; // Create linear ring for each hole in the shape
			for (int j = 1; j < coords.size(); j++) {
				final Coordinate[] innerCoords = coords.get(j).toArray(new Coordinate[coords.get(j).size()]);
				holes[j - 1] = GEOM_FACTORY.createLinearRing(innerCoords);
			}
			return GEOM_FACTORY.createPolygon(outer, holes);
		} else { // not closed
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
			case PConstants.ELLIPSE :
				final double a = shape.getParam(2) / 2d;
				final double b = shape.getParam(3) / 2d;
				final double perimeter = Math.PI * (3 * (a + b) - Math.sqrt((3 * a + b) * (a + 3 * b)));
				if ((int) Math.ceil(perimeter / BEZIER_SAMPLE_DISTANCE) < 4) {
					return GEOM_FACTORY.createPolygon();
				}
				shapeFactory.setNumPoints((int) Math.ceil(perimeter / BEZIER_SAMPLE_DISTANCE));
				shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
				shapeFactory.setWidth(a * 2);
				shapeFactory.setHeight(b * 2);
				return shapeFactory.createEllipse();
			case PConstants.TRIANGLE :
				Coordinate[] coords = new Coordinate[3 + 1];
				Coordinate c1 = new Coordinate(shape.getParam(0), shape.getParam(1));
				coords[0] = c1;
				coords[1] = new Coordinate(shape.getParam(2), shape.getParam(3));
				coords[2] = new Coordinate(shape.getParam(4), shape.getParam(5));
				coords[3] = c1.copy(); // close loop
				return GEOM_FACTORY.createPolygon(coords);
			case PConstants.RECT :
				final float w = shape.getParam(2);
				final float h = shape.getParam(3);
				shapeFactory.setCentre(new Coordinate(shape.getParam(0) + w / 2, shape.getParam(1) + h / 2));
				shapeFactory.setNumPoints(4);
				shapeFactory.setWidth(w);
				shapeFactory.setHeight(h);
				return shapeFactory.createRectangle();
			case PConstants.QUAD :
				coords = new Coordinate[4 + 1];
				c1 = new Coordinate(shape.getParam(0), shape.getParam(1));
				coords[0] = c1;
				coords[1] = new Coordinate(shape.getParam(2), shape.getParam(3));
				coords[2] = new Coordinate(shape.getParam(4), shape.getParam(5));
				coords[3] = new Coordinate(shape.getParam(6), shape.getParam(7));
				coords[4] = c1.copy(); // close loop
				return GEOM_FACTORY.createPolygon(coords);
			case PConstants.ARC :
				shapeFactory.setCentre(new Coordinate(shape.getParam(0), shape.getParam(1)));
				shapeFactory.setWidth(shape.getParam(2));
				shapeFactory.setHeight(shape.getParam(3));
				// circumference (if it was full circle)
				final double circumference = Math.PI * Math.max(shape.getParam(2), shape.getParam(3));
				shapeFactory.setNumPoints((int) Math.ceil(circumference / BEZIER_SAMPLE_DISTANCE));
				return shapeFactory.createArcPolygon(-Math.PI / 2 + shape.getParam(4), shape.getParam(5));
			case PConstants.LINE :
			case PConstants.POINT :
				System.err.print("Non-polygon primitives are not supported.");
				break;
			case PConstants.BOX :
			case PConstants.SPHERE :
				System.err.print("3D primitives are not supported.");
				break;
			default :
				System.err.print(shape.getKind() + " primitives are not supported.");

		}
		return GEOM_FACTORY.createPolygon(); // empty polygon
	}

	/**
	 * Transforms a list of points into a POINTS PShape.
	 * 
	 * @since 1.2.0
	 */
	public static final PShape toPointsPShape(Collection<PVector> points) {
		PShape shape = new PShape();
		shape.setFamily(PShape.GEOMETRY);
		shape.setStrokeCap(PConstants.ROUND);
		shape.setStroke(true);
		shape.setStroke(micycle.pgs.color.RGB.WHITE);
		shape.setStrokeWeight(5);
		shape.beginShape(PShape.POINTS);
		points.forEach(p -> shape.vertex(p.x, p.y));
		shape.endShape();
		return shape;
	}

	/**
	 * Returns the vertices of a PShape as an <b>unclosed</b> list of PVector
	 * coordinates.
	 * 
	 * @param shape a non-GROUP PShape
	 * @return
	 */
	public static List<PVector> toPVector(PShape shape) {
		if (shape.getFamily() == PShape.PRIMITIVE) {
			// getVertex() doesn't work on PShape primitives
			shape = toPShape(fromPrimitive(shape));
		}
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
	 * Converts a shape into a simple graph; graph vertices represent shape
	 * vertices, and graph edges represent shape edges (formed from adjacent
	 * vertices in polygonal shapes).
	 * 
	 * @param shape the shape to convert
	 * @return graph representation of the input shape
	 * @since 1.2.1
	 * @see #toDualGraph(PShape)
	 */
	public static SimpleGraph<PVector, PEdge> toGraph(PShape shape) {
		final SimpleGraph<PVector, PEdge> graph = new SimpleWeightedGraph<>(PEdge.class);
		for (PShape face : getChildren(shape)) {
			for (int i = 0; i < face.getVertexCount(); i++) {
				final PVector a = face.getVertex(i);
				final PVector b = face.getVertex((i + 1) % face.getVertexCount());
				if (a.equals(b)) {
					continue;
				}
				final PEdge e = new PEdge(a, b);

				graph.addVertex(a);
				graph.addVertex(b);
				graph.addEdge(a, b, e);
				graph.setEdgeWeight(e, e.length());
			}
		}

		return graph;
	}

	/**
	 * Takes as input a graph and computes a layout for the graph vertices using a
	 * Force-Directed placement algorithm (not vertex coordinates, if any exist).
	 * Vertices are joined by their edges.
	 * <p>
	 * The output is a rather abstract representation of the input graph, and not a
	 * geometric equivalent (unlike most other conversion methods in the class).
	 * 
	 * @param <V>                 any vertex type
	 * @param <E>                 any edge type
	 * @param graph               the graph whose edges and vertices to lay out
	 * @param normalizationFactor normalization factor for the optimal distance,
	 *                            between 0 and 1.
	 * @param boundsX             horizontal vertex bounds
	 * @param boundsY             vertical vertex bounds
	 * @return a GROUP PShape consisting of 2 children; child 0 is the linework
	 *         (LINES) depicting edges and child 1 is the points (POINTS) depicting
	 *         vertices. The bounds of the layout are anchored at (0, 0);
	 * @since 1.2.1
	 */
	public static <V, E> PShape fromGraph(SimpleGraph<V, E> graph, double normalizationFactor, double boundsX, double boundsY) {
		normalizationFactor = Math.min(Math.max(normalizationFactor, 0.001), 1);
		LayoutAlgorithm2D<V, E> layout;
		layout = new IndexedFRLayoutAlgorithm2D<>(50, 0.7, normalizationFactor, new Random(1337));
		LayoutModel2D<V> model = new MapLayoutModel2D<>(new Box2D(boundsX, boundsY));
		layout.layout(graph, model);

		NeighborCache<V, E> cache = new NeighborCache<>(graph);
		Set<PEdge> edges = new HashSet<>(graph.edgeSet().size());
		Map<V, PVector> pointMap = new HashMap<>();
		model.forEach(a -> {
			Point2D point = a.getValue();
			pointMap.put(a.getKey(), new PVector((float) point.getX(), (float) point.getY()));
		});
		pointMap.keySet().forEach(v -> cache.neighborsOf(v).forEach(n -> edges.add(new PEdge(pointMap.get(v), pointMap.get(n)))));

		PShape lines = PGS.prepareLinesPShape(null, null, null);
		edges.forEach(e -> {
			lines.vertex(e.a.x, e.a.y);
			lines.vertex(e.b.x, e.b.y);
		});
		lines.endShape();

		PShape pointsS = toPointsPShape(pointMap.values());
		pointsS.setStrokeWeight(10);

		return flatten(lines, pointsS);
	}

	/**
	 * Converts a mesh-like PShape into its undirected, unweighted dual-graph.
	 * <p>
	 * The output is a <i>dual graph</i> of the input; it has a vertex for each face
	 * (PShape) of the input, and an edge for each pair of faces that are adjacent.
	 * 
	 * @param mesh a GROUP PShape, whose children constitute the polygonal faces of
	 *             a <b>conforming mesh</b>. A conforming mesh consists of adjacent
	 *             cells that not only share edges, but every pair of shared edges
	 *             are identical (having the same coordinates) (such as a
	 *             triangulation).
	 * 
	 * @return the dual graph of the input mesh; an undirected graph containing no
	 *         graph loops or multiple edges.
	 * @since 1.2.1
	 * @see #toGraph(PShape)
	 */
	public static SimpleGraph<PShape, DefaultEdge> toDualGraph(PShape mesh) {
		return toDualGraph(getChildren(mesh));
	}

	static SimpleGraph<PShape, DefaultEdge> toDualGraph(Collection<PShape> meshFaces) {
		final SimpleGraph<PShape, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
		// map of which edge belong to each face; used to detect half-edges
		final HashMap<PEdge, PShape> edgesMap = new HashMap<>(meshFaces.size() * 4);

		for (PShape face : meshFaces) {
			graph.addVertex(face); // always add child so disconnected shapes are colored
			for (int i = 0; i < face.getVertexCount(); i++) {
				final PVector a = face.getVertex(i);
				final PVector b = face.getVertex((i + 1) % face.getVertexCount());
				if (a.equals(b)) {
					continue;
				}
				final PEdge e = new PEdge(a, b);
				final PShape neighbour = edgesMap.get(e);

				if (neighbour != null) {
					// edge seen before, so faces must be adjacent; create edge between faces
					if (neighbour.equals(face)) { // probably bad input (3 edges the same)
						System.err.println("toGraph: Bad input — saw the same edge 3 times.");
						continue; // continue to prevent self-loop in graph
					}
					graph.addEdge(neighbour, face);
				} else {
					edgesMap.put(e, face); // edge is new
				}
			}
		}
		return graph;
	}

	/**
	 * Writes the <i>Well-Known Text</i> representation of a shape. The
	 * <i>Well-Known Text</i> format is defined in the OGC Simple Features
	 * Specification for SQL.
	 * 
	 * @param shape shape to process
	 * @return a Geometry Tagged Text string
	 * @since 1.2.1
	 * @see #fromWKT(String)
	 */
	public static String toWKT(PShape shape) {
		WKTWriter writer = new WKTWriter(2);
		writer.setPrecisionModel(new PrecisionModel(PrecisionModel.FIXED)); // 1 d.p.
//		writer.setMaxCoordinatesPerLine(1);
		return writer.writeFormatted(fromPShape(shape));
	}

	/**
	 * Converts a geometry in <i>Well-Known Text</i> format into a PShape.
	 * 
	 * @param textRepresentation one or more Geometry Tagged Text strings, separated
	 *                           by whitespace
	 * @return a PShape specified by the text
	 * @since 1.2.1
	 * @see #toWKT(PShape)
	 */
	public static PShape fromWKT(String textRepresentation) {
		WKTReader reader = new WKTReader(GEOM_FACTORY);
		reader.setFixStructure(true);
		try {
			return toPShape(reader.read(textRepresentation));
		} catch (ParseException e) {
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Writes a shape into <i>Well-Known Binary</i> format. The WKB format is
	 * specified in the OGC <i>Simple Features for SQL specification</i></a>.
	 * 
	 * @param shape shape to process
	 * @return WKB byte representation of shape
	 * @since 1.2.1
	 * @see #fromWKB(byte[])
	 * @see #toHexWKB(PShape)
	 */
	public static byte[] toWKB(PShape shape) {
		WKBWriter writer = new WKBWriter();
		return writer.write(fromPShape(shape));
	}

	/**
	 * Converts a geometry in <i>Well-Known Binary</i> format into a PShape.
	 * 
	 * @param shapeWKB byte representation of shape to process
	 * @return a PShape specified by the WKB
	 * @since 1.2.1
	 * @see #toWKB(PShape)
	 */
	public static PShape fromWKB(byte[] shapeWKB) {
		WKBReader reader = new WKBReader();
		try {
			return toPShape(reader.read(shapeWKB));
		} catch (ParseException e) {
			return new PShape();
		}
	}

	/**
	 * Writes a shape into the hexadecimal string representation of its
	 * <i>Well-Known Binary</i> format.
	 * 
	 * @param shape shape to process
	 * @return hexadecimal string representation of shape WKB
	 * @since 1.2.1
	 * @see #toWKB(PShape)
	 */
	public static String toHexWKB(PShape shape) {
		return WKBWriter.toHex(toWKB(shape));
	}

	/**
	 * Converts a geometry in <i>Well-Known Binary</i> hex format into a PShape.
	 * 
	 * @param shapeWKB hex string WKB representation of shape to process
	 * @return a PShape specified by the WKB
	 * @since 1.2.1
	 * @see #toWKB(PShape)
	 */
	public static PShape fromHexWKB(String shapeWKB) {
		return fromWKB(WKBReader.hexToBytes(shapeWKB));
	}

	/**
	 * Creates a Java2D/java.awt Shape representing a PShape.
	 * 
	 * @param shape the PShape to convert
	 * @return a Java2D shape representing the PShape
	 * @since 1.2.1
	 */
	public static Shape toJava2D(PShape shape) {
		return new ShapeWriter().toShape(fromPShape(shape));
	}

	/**
	 * Converts a Java2D/java.awt Shape to a Processing PShape.
	 * <p>
	 * If the shape contains bezier components (such as <code>CubicCurve2D</code>,
	 * these are decomposed into straight-line segments in the output.
	 * 
	 * @param shape the Java2D shape to convert
	 * @return a PShape representing the Java2D shape
	 * @since 1.2.1
	 */
	public static PShape fromJava2D(Shape shape) {
		return toPShape(ShapeReader.read(shape, BEZIER_SAMPLE_DISTANCE, GEOM_FACTORY));
	}

	/**
	 * Generates a shape from a list of vertices. If the list of vertices is closed
	 * (first and last vertices are the same), the vertices are interpreted as a
	 * closed polygon (having no holes); if the list is unclosed, they are treated
	 * as a linestring.
	 * 
	 * @param vertices list of (un)closed shape vertices
	 * @return a PATH PShape (either open linestring or closed polygon)
	 * @see #fromPVector(PVector...)
	 */
	public static PShape fromPVector(Collection<PVector> vertices) {
		List<PVector> verticesList = new ArrayList<>(vertices);
		boolean closed = false;
		if (!vertices.isEmpty() && verticesList.get(0).equals(verticesList.get(vertices.size() - 1))) {
			closed = true;
		}

		PShape shape = new PShape();
		shape.setFamily(PShape.PATH);
		shape.setFill(micycle.pgs.color.RGB.WHITE);
		shape.setFill(closed);
		shape.setStroke(!closed);
		shape.setStroke(micycle.pgs.color.RGB.WHITE);
		shape.setStrokeWeight(2);

		shape.beginShape();
		for (int i = 0; i < verticesList.size() - (closed ? 1 : 0); i++) {
			PVector v = verticesList.get(i);
			shape.vertex(v.x, v.y);
		}
		shape.endShape(closed ? PConstants.CLOSE : PConstants.OPEN);

		return shape;
	}

	/**
	 * Generates a simple closed polygon (assumes no holes) from the list of
	 * vertices (varargs).
	 * 
	 * @param vertices list of (un)closed shape vertices
	 * @see #fromPVector(List)
	 */
	public static PShape fromPVector(PVector... vertices) {
		return fromPVector(Arrays.asList(vertices));
	}

	/**
	 * Flattens a collection of PShapes into a single GROUP PShape which has the
	 * input shapes as its children.
	 * 
	 * @since 1.2.0
	 * @see #flatten(PShape...)
	 */
	public static PShape flatten(Collection<PShape> shapes) {
		PShape group = new PShape(GROUP);
		shapes.forEach(group::addChild);
		return group;
	}

	/**
	 * Flattens a collection of PShapes into a single GROUP PShape which has the
	 * input shapes as its children.
	 * 
	 * @since 1.2.1
	 * @see #flatten(Collection)
	 */
	public static PShape flatten(PShape... shapes) {
		return flatten(Arrays.asList(shapes));
	}

	/**
	 * Recurses a GROUP PShape, finding all of its non-GROUP child PShapes.
	 * <p>
	 * This method differs from PShape.getChildren(): that method will return GROUP
	 * child shapes, whereas this method will recurse such shapes, returing their
	 * non-group children (in other words, this method explores the whole tree of
	 * shapes, returning non-group shapes only).
	 * 
	 * @param shape
	 * @return a list of non-GROUP PShapes
	 * @since 1.2.0
	 */
	public static List<PShape> getChildren(PShape shape) {
		final List<PShape> children = new ArrayList<>();
		final ArrayDeque<PShape> parents = new ArrayDeque<>();

		if (shape.getFamily() == GROUP) {
			parents.add(shape);
		} else {
			children.add(shape);
			return children;
		}

		while (!parents.isEmpty()) {
			final PShape parent = parents.pop(); // will always be a GROUP PShape
			if (parent.getChildCount() > 0) { // avoid NPE on .getChildren()
				for (PShape child : parent.getChildren()) {
					if (child.getFamily() == GROUP) {
						parents.add(child);
					} else {
						children.add(child);
					}
				}
			}
		}

		return children;
	}

	/**
	 * Creates a single GROUP shape whose children shapes are the list given.
	 * 
	 * @param children
	 * @return a GROUP PShape consisting of the given children
	 */
	public static PShape fromChildren(List<PShape> children) {
		final PShape parent = new PShape(GROUP);
		children.forEach(parent::addChild);
		return parent;
	}

	/**
	 * Sets the fill color for the PShape and all of its children recursively (and
	 * disables stroke).
	 * 
	 * @param shape
	 * @return the input object (having now been mutated)
	 * @see #setAllStrokeColor(PShape, int, float)
	 */
	public static PShape setAllFillColor(PShape shape, int color) {
		getChildren(shape).forEach(child -> {
			child.setStroke(false);
			child.setFill(true);
			child.setFill(color);
		});
		return shape;
	}

	/**
	 * Sets the stroke color for the PShape and all of its children recursively.
	 * 
	 * @param shape
	 * @return the input object (having now been mutated)
	 * @see {@link #setAllFillColor(PShape, int)}
	 */
	public static PShape setAllStrokeColor(PShape shape, int color, float strokeWeight) {
		getChildren(shape).forEach(child -> {
			child.setStroke(true);
			child.setStroke(color);
			child.setStrokeWeight(strokeWeight);
		});
		return shape;
	}

	/**
	 * Sets the stroke color equal to the fill color for the PShape and all of its
	 * descendent shapes individually (that is, each child shape belonging to the
	 * shape (if any) will have its stroke color set to <b>its own fill color</b>,
	 * and not the parent-most shape's fill color).
	 * 
	 * @param shape
	 * @return the input object (having now been mutated)
	 * @since 1.2.0
	 */
	public static PShape setAllStrokeToFillColor(PShape shape) {
		getChildren(shape).forEach(child -> {
			child.setStroke(true);
			child.setStroke(PGS.getPShapeFillColor(child));
		});
		return shape;
	}

	/**
	 * Calls setFill(false) on a PShape and all its children. This method mutates
	 * the input shape.
	 * 
	 * @param shape
	 * @return the input object (having now been mutated)
	 */
	public static PShape disableAllFill(PShape shape) {
		getChildren(shape).forEach(child -> child.setFill(false));
		return shape;
	}

	/**
	 * Calls setStroke(false) on a PShape and all its children. This method mutates
	 * the input shape.
	 * 
	 * @param shape
	 * @return the input object (having now been mutated)
	 */
	public static PShape disableAllStroke(PShape shape) {
		getChildren(shape).forEach(child -> child.setStroke(false));
		return shape;
	}

	/**
	 * Rounds the x and y coordinates (to the closest int) of all vertices belonging
	 * to the shape, <b>mutating</b> the shape. This can sometimes fix a visual
	 * problem in Processing where narrow gaps can appear between otherwise flush
	 * shapes.
	 * 
	 * @return the input object (having now been mutated)
	 * @since 1.1.3
	 */
	public static PShape roundVertexCoords(PShape shape) {
		getChildren(shape).forEach(c -> {
			for (int i = 0; i < c.getVertexCount(); i++) {
				final PVector v = c.getVertex(i);
				c.setVertex(i, Math.round(v.x), Math.round(v.y));
			}
		});
		return shape;
	}

	/**
	 * Produces a deep copy / clone of the input shape. Handles GROUP, PRIMITIVE,
	 * GEOMETRY and PATH PShapes.
	 * 
	 * @param shape the PShape to copy
	 * @return a deep copy of the given shape
	 * @since 1.2.0
	 */
	public static PShape copy(PShape shape) {
		final PShape copy = new PShape();
		copy.setName(shape.getName());

		try {
			Method method;
			switch (shape.getFamily()) {
				case GROUP :
					copy.setFamily(GROUP);
					getChildren(shape).forEach(child -> copy.addChild(copy(child)));
					return copy;
				case PShape.PRIMITIVE :
					copy.setFamily(PShape.PRIMITIVE);
					method = PShape.class.getDeclaredMethod("copyPrimitive", PShape.class, PShape.class);
					break;
				case PShape.GEOMETRY :
					copy.setFamily(PShape.GEOMETRY);
					method = PShape.class.getDeclaredMethod("copyGeometry", PShape.class, PShape.class);
					break;
				case PShape.PATH :
					copy.setFamily(PShape.PATH);
					method = PShape.class.getDeclaredMethod("copyPath", PShape.class, PShape.class);
					break;
				default :
					return copy;
			}
			method.setAccessible(true);
			method.invoke(null, shape, copy);
		} catch (NoSuchMethodException | SecurityException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException e) {
			e.printStackTrace();
		}

		return copy;
	}

	/**
	 * For every vertexcode, store the group (i.e. hole) it belongs to.
	 * 
	 * @param vertexCodes
	 * @return
	 */
	private static int[] getContourGroups(int[] vertexCodes) {

		int group = 0;
		List<Integer> groups = new ArrayList<>(vertexCodes.length * 2);

		for (int i = 0; i < vertexCodes.length; i++) {
			final int vertexCode = vertexCodes[i];
			switch (vertexCode) {
				case PConstants.VERTEX :
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

				case PConstants.BREAK :
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
	private static int[] getVertexTypes(int[] rawVertexCodes) {

		List<Integer> codes = new ArrayList<>(rawVertexCodes.length);

		for (int i = 0; i < rawVertexCodes.length; i++) {
			int vertexCode = rawVertexCodes[i];
			switch (vertexCode) {
				case PConstants.VERTEX :
					codes.add(PConstants.VERTEX);
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

				case PConstants.BREAK :
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
			coords.add(coordFromPVector(start));
			coords.add(coordFromPVector(end));
			return coords;
		}

		final float length = bezierLengthQuadratic(start, controlPoint, end);
		final int samples = (int) Math.ceil(length / sampleDistance); // sample every x unit length (approximately)

		coords.add(coordFromPVector(start));
		for (int j = 1; j < samples; j++) { // start at 1 -- don't sample at t=0
			final PVector bezierPoint = getQuadraticBezierCoordinate(start, controlPoint, end, j / (float) samples);
			coords.add(coordFromPVector(bezierPoint));
		}
		coords.add(coordFromPVector(end));

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
			coords.add(coordFromPVector(start));
			coords.add(coordFromPVector(end));
			return coords;
		}

		final float length = bezierLengthCubic(start, controlPoint1, controlPoint2, end);
		final int samples = (int) Math.ceil(length / sampleDistance); // sample every x unit length (approximately)

		coords.add(coordFromPVector(start));
		for (int j = 1; j < samples; j++) { // start at 1 -- don't sample at t=0
			final PVector bezierPoint = getCubicBezierCoordinate(start, controlPoint1, controlPoint2, end, j / (float) samples);
			coords.add(coordFromPVector(bezierPoint));
		}
		coords.add(coordFromPVector(end));
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

	static class PShapeData {

		private static Field fillColorF, fillF, strokeColorF, strokeWeightF, strokeF;

		static {
			try {
				fillColorF = PShape.class.getDeclaredField("fillColor");
				fillColorF.setAccessible(true);
				fillF = PShape.class.getDeclaredField("fill");
				fillF.setAccessible(true);
				strokeColorF = PShape.class.getDeclaredField("strokeColor");
				strokeColorF.setAccessible(true);
				strokeWeightF = PShape.class.getDeclaredField("strokeWeight");
				strokeWeightF.setAccessible(true);
				strokeF = PShape.class.getDeclaredField("stroke");
				strokeF.setAccessible(true);
			} catch (NoSuchFieldException | SecurityException e) {
				e.printStackTrace();
			}
		}

		int fillColor, strokeColor;
		float strokeWeight;
		boolean fill, stroke;

		private PShapeData(PShape shape) {
			try {
				fillColor = fillColorF.getInt(shape);
				fill = fillF.getBoolean(shape);
				stroke = strokeF.getBoolean(shape);
				strokeColor = strokeColorF.getInt(shape);
				strokeWeight = strokeWeightF.getFloat(shape);
			} catch (IllegalArgumentException | IllegalAccessException e) {
				e.printStackTrace();
			}
		}

		/**
		 * Apply this shapedata to a given PShape.
		 * 
		 * @param other
		 */
		void applyTo(PShape other) {
			other.setFill(fill);
			other.setFill(fillColor);
			other.setStroke(stroke);
			other.setStroke(strokeColor);
			other.setStrokeWeight(strokeWeight);
		}

		@Override
		public String toString() {
			return String.format("fillColor: %s; strokeColor: %s; strokeWeight: %.1f", Arrays.toString(decomposeclrRGB(fillColor)),
					Arrays.toString(decomposeclrRGB(strokeColor)), strokeWeight);
		}
	}
}
