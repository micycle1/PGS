package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.color.ColorUtils.decomposeclrRGB;
import static processing.core.PConstants.BEZIER_VERTEX;
import static processing.core.PConstants.CURVE_VERTEX;
import static processing.core.PConstants.GROUP;
import static processing.core.PConstants.QUADRATIC_VERTEX;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.PathIterator;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
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
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.awt.ShapeReader;
import org.locationtech.jts.awt.ShapeWriter;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
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
import org.locationtech.jts.io.geojson.GeoJsonReader;
import org.locationtech.jts.io.geojson.GeoJsonWriter;
import org.locationtech.jts.util.GeometricShapeFactory;
import org.scoutant.polyline.PolylineDecoder;

import com.github.micycle1.betterbeziers.CubicBezier;

import it.rambow.master.javautils.PolylineEncoder;
import it.rambow.master.javautils.Track;
import it.rambow.master.javautils.Trackpoint;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import micycle.pgs.color.Colors;
import micycle.pgs.commons.Nullable;
import micycle.pgs.commons.PEdge;
import net.jafama.FastMath;
import processing.core.PConstants;
import processing.core.PMatrix;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Facilitates conversion between <i>Processing's</i> {@code PShapes} and
 * <i>JTS's</i> {@code Geometries}, along with various other formats. It also
 * offers additional utility methods to assist with handling {@code PShapes}.
 * <p>
 * Though certain conversion methods are utilised internally by the library,
 * they have been kept public to cater to more complex user requirements.
 * <p>
 * Note: JTS {@code Geometries} do not provide support for bezier curves. As
 * such, bezier curves are linearised/divided into straight line segments during
 * the conversion process from {@code PShape} to JTS {@code Geometry}.
 * Furthermore, any consecutive duplicate vertices are not preserved.
 * <p>
 * Two configurable boolean flags influence the conversion process:
 * {@link #PRESERVE_STYLE} (set to true by default), and
 * {@link #HANDLE_MULTICONTOUR} (set to false by default). Users are encouraged
 * to review these flags as part of more complicated workflows with this class.
 * 
 * @author Michael Carleton
 */
public final class PGS_Conversion {

	/** Approximate distance between successive sample points on bezier curves */
	static final float BEZIER_SAMPLE_DISTANCE = 2;
	private static Field MATRIX_FIELD, PSHAPE_FILL_FIELD;
	/**
	 * A boolean flag that affects whether a PShape's style (fillColor, strokeColor,
	 * strokeWidth) is preserved during <code>PShape->Geometry->PShape</code>
	 * conversion (i.e. when <code>toPShape(fromPShape(myPShape))</code> is called).
	 * Default = <code>true</code>.
	 */
	public static boolean PRESERVE_STYLE = true;
	/**
	 * A boolean flag that, when true, enables a specialised subroutine during the
	 * {@link #fromPShape(PShape) fromPShape()} conversion to correctly convert
	 * <b>single</b> PShapes comprised of multiple contours, each representing a
	 * separate shape. If set to <code>false</code>, the {@link #fromPShape(PShape)
	 * fromPShape()} method assumes that in multi-contour shapes, every contour
	 * beyond the first represents a hole, which is generally an adequate
	 * assumption.
	 * <p>
	 * This feature is disabled by default because it necessitates additional
	 * computation, such as determining polygon ring orientations, and is seldom
	 * required (unless one is dealing with fonts, I have observed). Default =
	 * <code>false</code>.
	 * <p>
	 * For more information, refer to the discussion on this topic at
	 * <a href="https://github.com/micycle1/PGS/issues/67">GitHub</a>.
	 */
	public static boolean HANDLE_MULTICONTOUR = false;

	static {
		try {
			MATRIX_FIELD = PShape.class.getDeclaredField("matrix");
			MATRIX_FIELD.setAccessible(true);
			PSHAPE_FILL_FIELD = PShape.class.getDeclaredField("fillColor");
			PSHAPE_FILL_FIELD.setAccessible(true);
		} catch (NoSuchFieldException e) {
			System.err.println(e.getLocalizedMessage());
		}

		try {
			// use JTS "next gen" overlay. configurable via system property but will often
			// load before, hence use reflection
			Class<?> geometryOverlayClass = Class.forName("org.locationtech.jts.geom.GeometryOverlay");
			Field isOverlayNGField = geometryOverlayClass.getDeclaredField("isOverlayNG");
			isOverlayNGField.setAccessible(true);
			isOverlayNGField.setBoolean(null, true);
		} catch (ClassNotFoundException | NoSuchFieldException | IllegalAccessException e) {
			e.printStackTrace();
		}
	}

	private PGS_Conversion() {
	}

	/**
	 * Converts a JTS Geometry into a corresponding PShape. In the case of
	 * MultiGeometries (which include collections of geometries), the result is a
	 * GROUP PShape containing the appropriate child PShapes.
	 * <p>
	 * The conversion process follows the geometry types supported by JTS, namely:
	 * <ul>
	 * <li>{@link Geometry#TYPENAME_GEOMETRYCOLLECTION GEOMETRYCOLLECTION}:
	 * Converted to a GROUP PShape if it contains multiple geometries. For single
	 * geometry collections, it extracts and converts the single geometry.</li>
	 * <li>{@link Geometry#TYPENAME_MULTIPOLYGON MULTIPOLYGON}: Similar to
	 * GeometryCollection, MultiPolygons are converted to a GROUP PShape, with each
	 * polygon converted to a child PShape.</li>
	 * <li>{@link Geometry#TYPENAME_MULTILINESTRING MULTILINESTRING}:
	 * MultiLineStrings are handled in the same way as GeometryCollections and
	 * MultiPolygons, converted to a GROUP PShape containing child PShapes.</li>
	 * <li>{@link Geometry#TYPENAME_LINEARRING LINEARRING} and
	 * {@link Geometry#TYPENAME_LINESTRING LINESTRING}: These are converted to a
	 * PATH PShape, preserving the closed or open nature of the original
	 * LineString.</li>
	 * <li>{@link Geometry#TYPENAME_POLYGON POLYGON}: Converted to a PATH PShape,
	 * with each contour of the polygon represented as a series of vertices in the
	 * PShape. Inner contours, or 'holes' in the polygon, are handled
	 * separately.</li>
	 * <li>{@link Geometry#TYPENAME_POINT POINT} and
	 * {@link Geometry#TYPENAME_MULTIPOINT MULTIPOINT}: These are converted to a
	 * GEOMETRY PShape with each point represented as a vertex.</li>
	 * </ul>
	 * <p>
	 * Please note that any unsupported geometry types will result in an error
	 * message.
	 * <p>
	 * If {@link #PRESERVE_STYLE} is enabled and the geometry includes user data in
	 * the form of PShapeData, the style from the data is applied to the resulting
	 * PShape.
	 * 
	 * @param g The JTS geometry to convert.
	 * @return A PShape that represents the input geometry, or a new, empty PShape
	 *         if the input is null.
	 */
	public static PShape toPShape(final Geometry g) {
		if (g == null) {
			return new PShape();
		}
		PShape shape = new PShape();
		// apply PGS style by default
		shape.setFill(true);
		shape.setFill(micycle.pgs.color.Colors.WHITE);
		shape.setStroke(true);
		shape.setStroke(micycle.pgs.color.Colors.PINK);
		shape.setStrokeWeight(4);
		shape.setStrokeJoin(PConstants.ROUND);
		shape.setStrokeCap(PConstants.ROUND);

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
			// TODO treat closed linestrings as unfilled & unclosed paths?
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
				for (final Coordinate coord : coords) {
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
	 * 
	 * Converts a collection of JTS Geometries into a corresponding GROUP PShape.
	 * This method loops through the provided geometries, converting each individual
	 * geometry into a PShape, and then adds it as a child to the GROUP PShape.
	 * <p>
	 * In case the collection only contains a single geometry, this method will
	 * instead return a PShape that directly corresponds to that single geometry. It
	 * will not be wrapped in a GROUP shape in this case.
	 * 
	 * @param geometries A collection of JTS Geometries to convert into a PShape.
	 * @return A PShape that represents the collection of input geometries. If the
	 *         collection contains only a single geometry, the return is a PShape
	 *         directly equivalent to that geometry. Otherwise, the return is a
	 *         GROUP PShape containing child PShapes for each geometry in the
	 *         collection.
	 */
	public static PShape toPShape(Collection<? extends Geometry> geometries) {
		PShape shape = new PShape(GROUP);
		shape.setFill(true);
		shape.setFill(micycle.pgs.color.Colors.WHITE);
		shape.setStroke(true);
		shape.setStroke(micycle.pgs.color.Colors.PINK);
		shape.setStrokeWeight(4);

		geometries.forEach(g -> shape.addChild(toPShape(g)));
		if (shape.getChildCount() == 1) {
			return shape.getChild(0);
		}

		return shape;
	}

	/**
	 * 
	 * Transforms a PShape into a corresponding JTS Geometry.
	 * <p>
	 * During this transformation, any bezier curve elements within the input PShape
	 * are linearised, meaning they are sampled at regular, equidistant intervals.
	 * This process results in the created geometry having a greater number of
	 * vertices than the original PShape.
	 * <p>
	 * Additionally, please be aware that the method does not preserve multi-level
	 * child shape hierarchies present in the input PShape. All child shapes are
	 * flattened to the same level in the output geometry.
	 * <p>
	 * The conversion process depends on the PShape's type and can be broadly
	 * categorized as follows:
	 * <ul>
	 * <li>{@link PConstants#GROUP}: The method recursively converts each child of
	 * the PShape into a corresponding Geometry and groups these into a
	 * GeometryCollection.</li>
	 * <li>{@link PShape#GEOMETRY} and {@link PShape#PATH}: Here, the method further
	 * distinguishes between the kinds of the shape. For POLYGON, PATH and
	 * unspecified kinds, it creates a Geometry from the vertices of the PShape. For
	 * special paths (e.g., POINTS, LINES), it uses a separate conversion
	 * routine.</li>
	 * <li>{@link PShape#PRIMITIVE}: It converts the PShape using a separate routine
	 * dedicated to primitive shapes.</li>
	 * </ul>
	 * <p>
	 * If {@link #PRESERVE_STYLE} is enabled, the method preserves the style of the
	 * PShape in the output Geometry as user data.
	 * <p>
	 * Lastly, any affine transformations applied to the PShape (which do not
	 * directly affect its vertices) are also applied to the resulting Geometry
	 * (baked into its coordinates).
	 * 
	 * @param shape The PShape to convert into a JTS Geometry.
	 * @return A JTS Geometry that corresponds to the input PShape, or an empty 2D
	 *         Geometry if the PShape is null or the type of the PShape is
	 *         unsupported.
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
					g = fromVertices(shape);// holes possible
				} else {
					g = fromCreateShape(shape); // special paths (e.g. POINTS, LINES, etc.) (no holes)
				}
				break;
			case PShape.PRIMITIVE :
				g = fromPrimitive(shape); // (no holes)
				break;
		}

		if (PRESERVE_STYLE && g != null) {
			g.setUserData(new PShapeData(shape));
		}

		/*
		 * Finally, apply PShape's affine transformations, if present, to geometry
		 * vertices (which are not applied to its vertices directly).
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
			case PConstants.TRIANGLE :
				final Coordinate[] triangle = new Coordinate[3 + 1];
				final Coordinate a = PGS.coordFromPVector(shape.getVertex(0));
				triangle[0] = a;
				triangle[1] = PGS.coordFromPVector(shape.getVertex(1));
				triangle[2] = PGS.coordFromPVector(shape.getVertex(2));
				triangle[3] = a.copy();
				return GEOM_FACTORY.createPolygon(triangle);
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
	 * Extracts the contours from a <code>POLYGON</code> or <code>PATH</code>
	 * PShape, represented as lists of PVector points. It extracts both the exterior
	 * contour (perimeter) and interior contours (holes). For such PShape types, all
	 * contours after the first are guaranteed to be holes.
	 * <p>
	 * Background: The PShape data structure stores all vertices in a single array,
	 * with contour breaks designated in a separate array of vertex codes. This
	 * design makes it challenging to directly access holes and their vertices. This
	 * method overcomes that limitation by analyzing the vertex codes to identify
	 * and extract both the exterior contour and interior contours (holes).
	 *
	 * @param shape The PShape object (of kind POLYGON or PATH) from which to
	 *              extract contours.
	 * @return A List of Lists of PVector coordinates, where each inner List
	 *         represents a contour. <b>The first contour is always the exterior,
	 *         while all subsequent contours represent holes</b>.
	 * @since 2.0
	 */
	public static List<List<PVector>> toContours(PShape shape) {
		List<List<PVector>> rings = new ArrayList<>();

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

		int currentGroup = -1;
		PVector lastVertex = null;

		for (int i = 0; i < shape.getVertexCount(); i++) {
			if (contourGroups[i] != currentGroup) {
				currentGroup = contourGroups[i];
				rings.add(new ArrayList<>());
				lastVertex = null;
			}

			/*
			 * Sample bezier curves at regular intervals to produce smooth Geometry
			 */
			switch (vertexCodes[i]) { // VERTEX, BEZIER_VERTEX, CURVE_VERTEX, or BREAK
				case QUADRATIC_VERTEX :
					rings.get(currentGroup)
							.addAll(getQuadraticBezierPoints(shape.getVertex(i - 1), shape.getVertex(i), shape.getVertex(i + 1), BEZIER_SAMPLE_DISTANCE));
					i += 1;
					continue;
				case BEZIER_VERTEX : // aka cubic bezier
					rings.get(currentGroup).addAll(getCubicBezierPoints(shape.getVertex(i - 1), shape.getVertex(i), shape.getVertex(i + 1),
							shape.getVertex(i + 2), BEZIER_SAMPLE_DISTANCE));
					i += 2;
					continue;
				default : // VERTEX
					PVector v = shape.getVertex(i).copy();
					// skip consecutive duplicate vertices
					if (lastVertex == null || !(v.x == lastVertex.x && v.y == lastVertex.y)) {
						rings.get(currentGroup).add(v);
						lastVertex = v;
					}
					break;
			}
		}

		// close rings
		if (shape.isClosed()) {
			rings.forEach(ring -> ring.add(ring.get(0).copy()));
		}

		return rings;
	}

	/**
	 * Creates a JTS Polygon from a lineal/polygonal PShape, whose 'kind' is a
	 * <code>POLYGON</code> or <code>PATH</code>.
	 * <p>
	 * Note that repeated vertices are not preserved during conversion (to maximise
	 * compatibility with geometric algorithms).
	 */
	private static Geometry fromVertices(PShape shape) {
		final List<List<PVector>> contours = toContours(shape);

		if (contours.isEmpty()) { // skip empty / pointal PShapes
			return GEOM_FACTORY.createPolygon();
		}

		final List<Coordinate[]> rings = contours.stream().map(PGS::toCoords).toList();
		final Coordinate[] outerRing = rings.get(0);

		if (outerRing.length == 0) {
			return GEOM_FACTORY.createPolygon(); // empty polygon
		} else if (outerRing.length == 1) {
			return GEOM_FACTORY.createPoint(outerRing[0]);
		} else if (outerRing.length == 2) {
			return GEOM_FACTORY.createLineString(outerRing);
		} else if (shape.isClosed()) { // closed geometry or path
			if (HANDLE_MULTICONTOUR) { // handle single shapes that *may* represent multiple shapes over many contours
				return fromMultiContourShape(rings, false, false);
			} else { // assume all contours beyond the first represent holes
				LinearRing outer = GEOM_FACTORY.createLinearRing(outerRing); // should always be valid
				LinearRing[] holes = new LinearRing[rings.size() - 1]; // Create linear ring for each hole in the shape
				for (int j = 1; j < rings.size(); j++) {
					final Coordinate[] innerCoords = rings.get(j);
					holes[j - 1] = GEOM_FACTORY.createLinearRing(innerCoords);
				}
				return GEOM_FACTORY.createPolygon(outer, holes);
			}
		} else { // not closed
			return GEOM_FACTORY.createLineString(outerRing);
		}
	}

	/**
	 * <p>
	 * Transforms a {@code PShape} object, which might contain multiple contours
	 * representing several polygons (including nested polygons with holes), into a
	 * valid JTS representation.
	 * <p>
	 * The input shape may not be directly producible using JTS because JTS
	 * represents multiple polygons as separate objects within a
	 * {@code MultiGeometry}, resulting in a {@code GROUP} shape. However,
	 * Processing allows <b>single</b> {@code PATH} PShapes to include multiple
	 * contour groups, representing multiple polygons despite being considered a
	 * "single" shape.
	 * <p>
	 * The {@code PShape} does not carry information about which contours represent
	 * holes or how contours should be grouped to represent the same shapes. This
	 * method determines this internally.
	 * 
	 * @param contours A {@code List} of contours/rings. The contour at
	 *                 {@code index==0} is assumed to be polygonal and considered an
	 *                 exterior ring.
	 * @param sort     A boolean value determining whether to sort the contours by
	 *                 orientation (clockwise contours first), which can solve some
	 *                 complex cases.
	 * @param reverse  A boolean value indicating whether to reverse the contour
	 *                 collection, which can also solve complex cases.
	 * @return Returns either a {@code Polygon} or {@code MultiPolygon} depending on
	 *         the input shape.
	 * 
	 * @see Geometry
	 * @see CoordinateList
	 */
	private static Geometry fromMultiContourShape(List<Coordinate[]> contours, boolean sort, boolean reverse) {
		if (reverse) {
			Collections.reverse(contours);
		}
		if (sort) {
			contours.sort((a, b) -> {
				boolean aCW = !Orientation.isCCWArea(a);
				boolean bCW = !Orientation.isCCWArea(b);
				if (aCW == bCW) {
					return 0;
				} else if (aCW && !bCW) {
					return -1;
				} else {
					return 1;
				}
			});
		}
		if (contours.isEmpty()) {
			return GEOM_FACTORY.createPolygon();
		} else if (contours.size() == 1) {
			return GEOM_FACTORY.createPolygon(contours.get(0));
		} else {
			boolean previousRingIsCCW = false;
			boolean previousRingIsHole = false;

			/*
			 * Each list entry (ring group) holds the exterior polygon ring at index 0,
			 * followed by any number of hole rings. Created rings follow JTS expected
			 * orientation: clockwise as polygon exterior ring orientation; anti-clockwise
			 * for holes.
			 */
			final List<List<LinearRing>> polygonRingGroups = new ArrayList<>();

			/*
			 * Ideally, the hole that succeeds a polygon exterior would belong to that
			 * exterior. However, sometimes contours in valid PShapes aren't ordered in this
			 * fashion: it is possible to have four contours: a,b,1,2; where the exteriors
			 * (a & b) are followed by two holes (1 & 2), where hole 1 belongs to a and hole
			 * 2 belongs to b. Hence (to be really thorough), we must check that a hole is
			 * actually contained by the last contour, otherwise we associate it with the
			 * first exterior (working backwards) that contains it.
			 */
			for (int j = 0; j < contours.size(); j++) {
				Coordinate[] contourCoords = contours.get(j);
				LinearRing ring;
				/*
				 * Measure contour orientation to determine whether the contour represents a
				 * hole (orientated opposite to previous polygon exterior) or a polygon exterior
				 * (orientated in same direction to previous contour).
				 */
				final boolean ringIsCCW = Orientation.isCCWArea(contourCoords);
				final boolean switched = previousRingIsCCW != ringIsCCW;
				previousRingIsCCW = ringIsCCW;

				if (!ringIsCCW) {
					contourCoords = reversedCopy(contourCoords); // alt to .toCoordinateArray(ringIsCCW)
				}
				if (((switched && !previousRingIsHole) || (!switched && previousRingIsHole)) && j > 0) { // this ring is hole
					if (switched) {
						previousRingIsHole = true;
					}

					/*
					 * Find exterior that contains the hole (usually is most recent one).
					 */
					ring = GEOM_FACTORY.createLinearRing(contourCoords);
					int checkContainsPolygonIndex = polygonRingGroups.size() - 1;
					while (checkContainsPolygonIndex >= 0
							&& !GEOM_FACTORY.createPolygon(polygonRingGroups.get(checkContainsPolygonIndex).get(0)).contains(ring)) {
						checkContainsPolygonIndex--;
					}
					if (checkContainsPolygonIndex >= 0) {
						polygonRingGroups.get(checkContainsPolygonIndex).add(ring);
					} else {
						if (!sort) {
							if (!reverse) {
								// retry and reverse contours, which can fix some cases
								return fromMultiContourShape(contours, false, true);
							}
							// retry and sort by orientation, which can fix some cases (assuming CW exterior
							// rings)
							return fromMultiContourShape(contours, true, false);
						} else {
							System.err.println(String.format(
									"PGS_Conversion Error: Shape contour #%s was identified as a hole but no existing exterior rings contained it.", j));
						}
					}
				} else { // this ring is new polygon (or explictly contour #1)
					ring = GEOM_FACTORY.createLinearRing(contourCoords);
					if (previousRingIsHole) {
						previousRingIsHole = false;
					}
					polygonRingGroups.add(new ArrayList<>());
					polygonRingGroups.get(polygonRingGroups.size() - 1).add(ring);
				}
			}

			/*
			 * Convert each ring group to their representative polygon.
			 */
			final List<Polygon> polygons = new ArrayList<>();
			polygonRingGroups.forEach(perPolygonRings -> {
				LinearRing[] holes = null;
				if (perPolygonRings.size() > 1) { // has holes
					holes = perPolygonRings.subList(1, perPolygonRings.size()).toArray(new LinearRing[0]);
				}
				polygons.add(GEOM_FACTORY.createPolygon(perPolygonRings.get(0), holes));
			});

			if (polygons.size() > 1) {
				return GEOM_FACTORY.createMultiPolygon(polygons.toArray(new Polygon[0]));
			} else {
				return polygons.get(0);
			}
		}
	}

	/**
	 * Creates a JTS Polygon from a primitive PShape. Primitive PShapes are
	 * generated by invoking the <code>createShape()</code> method from Processing.
	 * <p>
	 * Supported primitives include POINT, LINE, TRIANGLE, QUAD, RECT, ELLIPSE, and
	 * ARC. Other types such as BOX and SPHERE (3D primitives), or any non-polygon
	 * primitives are not supported and will print an error message.
	 * <p>
	 * Primitive PShapes do not have directly accessible vertex data. This method
	 * generates equivalent JTS Polygons by accessing the shape's parameters (e.g.,
	 * width, height, coordinates) via the <code>getParam()</code> method.
	 * <p>
	 * In the event where a non-supported primitive is supplied, an empty JTS
	 * Polygon is returned.
	 *
	 * @param shape The primitive PShape to be converted to a JTS Polygon
	 * @return A JTS Geometry (Polygon) representing the input primitive PShape. For
	 *         non-supported or non-polygon primitives, an empty JTS Polygon is
	 *         returned.
	 */
	private static Geometry fromPrimitive(PShape shape) {
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
			case PConstants.POINT :
				return GEOM_FACTORY.createPoint(new Coordinate(shape.getParam(0), shape.getParam(1)));
			case PConstants.LINE :
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
	 * Transforms a variable arg list of points into a POINTS PShape.
	 * 
	 * @param vertices
	 * @return a POINTS PShape
	 * @since 1.4.0
	 */
	public static final PShape toPointsPShape(PVector... vertices) {
		return toPointsPShape(Arrays.asList(vertices));
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
		shape.setStroke(micycle.pgs.color.Colors.PINK);
		shape.setStrokeWeight(6);
		shape.beginShape(PConstants.POINTS);
		points.forEach(p -> shape.vertex(p.x, p.y));
		shape.endShape();
		return shape;
	}

	/**
	 * Creates a PShape having circle geometries representing a collection of
	 * circles.
	 * 
	 * @param circles The collection of PVector objects representing the circles.
	 *                The x and y components represent the center of the circle, and
	 *                the z component represents the <b>radius</b>.
	 * @return The GROUP PShape object representing the collection of circles.
	 * @since 1.4.0
	 */
	public static final PShape toCircles(Collection<PVector> circles) {
		return toPShape(circles.stream().map(c -> PGS_Construction.createCirclePoly(c.x, c.y, c.z)).collect(Collectors.toList()));
	}

	/**
	 * 
	 * Extracts the vertices of a PShape into a list of PVectors.
	 * <p>
	 * The function navigates through all children of the given shape if it is of
	 * the GROUP type, recursively flattening their vertices and adding them to the
	 * list. In the case of PShape primitives, where the <code>getVertex()</code>
	 * method fails, the shape is converted to its equivalent path representation
	 * before vertex extraction.
	 * <p>
	 * If the input shape represents a closed polygon, the method returns an
	 * "unclosed" version of the shape. This means that the duplicate vertex that
	 * closes the shape (which is identical to the first vertex) is omitted from the
	 * output.
	 * <p>
	 * The resulting list contains all vertices from the input PShape in the order
	 * they appear in the shape.
	 * 
	 * @param shape the PShape from which vertices are to be extracted
	 * @return a list of PVector objects representing the vertices of the input
	 *         shape
	 */
	public static List<PVector> toPVector(PShape shape) {
		// use getChildren() incase shape is GROUP
		final List<PVector> vertices = new ArrayList<>();
		getChildren(shape).forEach(s -> {
			if (s.getFamily() == PShape.PRIMITIVE) {
				// getVertex() doesn't work on PShape primitives
				s = toPShape(fromPrimitive(s));
			}
			for (int i = 0; i < s.getVertexCount(); i++) {
				vertices.add(s.getVertex(i));
			}
		});
		if (!vertices.isEmpty() && shape.getChildCount() > 0 && vertices.get(0).equals(vertices.get(vertices.size() - 1))) {
			vertices.remove(vertices.size() - 1);
		}
		return vertices;
	}

	/**
	 * Transforms a given PShape into a simple graph representation. In this
	 * representation, the vertices of the graph correspond to the vertices of the
	 * shape, and the edges of the graph correspond to the edges of the shape.
	 * <p>
	 * The edge weights in the graph are set to the length (euclidean distance) of
	 * the corresponding geometric edge in the shape.
	 * 
	 * @param shape the PShape to convert into a graph. LINES and polygonal shapes
	 *              are accepted (and GROUP shapes thereof).
	 * @return A SimpleGraph object that represents the structure of the input shape
	 * @since 1.3.0
	 * @see #toDualGraph(PShape)
	 */
	public static SimpleGraph<PVector, PEdge> toGraph(PShape shape) {
		final SimpleGraph<PVector, PEdge> graph = new SimpleWeightedGraph<>(PEdge.class);
		for (PShape child : getChildren(shape)) {
			final int stride = child.getKind() == PShape.LINES ? 2 : 1;
			// Handle other child shapes (e.g., faces)
			for (int i = 0; i < child.getVertexCount() - (child.isClosed() ? 0 : 1); i += stride) {
				final PVector a = child.getVertex(i);
				final PVector b = child.getVertex((i + 1) % child.getVertexCount());
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
	 * Converts a given SimpleGraph consisting of PVectors and PEdges into a PShape
	 * by polygonizing its edges. If the graph represented a shape with holes, these
	 * will not be preserved during the conversion.
	 * 
	 * @param graph the graph to be converted into a PShape.
	 * @return a PShape representing the polygonized edges of the graph.
	 * @since 1.4.0
	 */
	public static PShape fromGraph(SimpleGraph<PVector, PEdge> graph) {
		return PGS.polygonizeNodedEdges(graph.edgeSet());
	}

	/**
	 * Computes a layout for the vertices of a graph using a Force-Directed
	 * placement algorithm. The algorithm generates vertex coordinates based on the
	 * graph's topology, preserving its structure (i.e., connectivity and
	 * relationships between vertices and edges). Existing vertex coordinates, if
	 * any, are ignored.
	 * <p>
	 * The output is an abstract representation of the input graph, not a geometric
	 * equivalent (unlike most other conversion methods in this class). The layout
	 * is bounded by the specified dimensions and anchored at (0, 0).
	 *
	 * @param <V>                 the type of vertices in the graph
	 * @param <E>                 the type of edges in the graph
	 * @param graph               the graph whose vertices and edges are to be laid
	 *                            out
	 * @param normalizationFactor the normalization factor for the optimal distance
	 *                            between vertices, clamped between 0.001 and 1
	 * @param boundsX             the horizontal bounds for the layout
	 * @param boundsY             the vertical bounds for the layout
	 * @return a GROUP PShape containing two children: child 0 represents the edges
	 *         as linework (LINES), and child 1 represents the vertices as points
	 *         (POINTS)
	 * @since 1.3.0
	 */
	public static <V, E> PShape fromGraph(SimpleGraph<V, E> graph, double normalizationFactor, double boundsX, double boundsY) {
		normalizationFactor = Math.min(Math.max(normalizationFactor, 0.001), 1);
		LayoutAlgorithm2D<V, E> layout;
		layout = new IndexedFRLayoutAlgorithm2D<>(50, 0.7, normalizationFactor, new XoRoShiRo128PlusRandom(1337));
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
	 * @since 1.3.0
	 * @see #toGraph(PShape)
	 */
	public static SimpleGraph<PShape, DefaultEdge> toDualGraph(PShape mesh) {
		return toDualGraph(getChildren(mesh));
	}

	/**
	 * Converts a mesh-like PShape into its centroid-based undirected dual-graph.
	 * <p>
	 * The output is a <i>dual graph</i> of the input; it has a vertex for each
	 * centroid of the face of the input, and an edge (connecting the centroids) for
	 * each pair of faces that are adjacent. Each vertex represents the geometric
	 * center or centroid of the respective face in the mesh.
	 *
	 * @param mesh a GROUP PShape, whose children constitute the polygonal faces of
	 *             a <b>conforming mesh</b>. A conforming mesh consists of adjacent
	 *             cells that not only share edges, but every pair of shared edges
	 *             are identical (having the same coordinates) (such as a
	 *             triangulation).
	 * @return the centroid-based dual graph of the input mesh; an undirected graph
	 *         containing no graph loops or multiple edges. Each vertex in the graph
	 *         represents the centroid of a face in the input mesh, and each edge
	 *         represents adjacency between two faces.
	 * @since 1.4.0
	 * @see #toDualGraph(PShape)
	 * @see PGS_ShapePredicates#centroid(PShape)
	 */
	public static SimpleGraph<PVector, PEdge> toCentroidDualGraph(PShape mesh) {
		final SimpleGraph<PShape, DefaultEdge> toplogy = toDualGraph(getChildren(mesh));
		Map<PShape, PVector> centroids = toplogy.vertexSet().stream().collect(Collectors.toMap(x -> x, PGS_ShapePredicates::centroid));

		SimpleGraph<PVector, PEdge> graph = new SimpleGraph<>(PEdge.class);
		centroids.values().forEach(v -> graph.addVertex(v));
		toplogy.edgeSet().forEach(e -> {
			PVector c1 = centroids.get(toplogy.getEdgeSource(e));
			PVector c2 = centroids.get(toplogy.getEdgeTarget(e));
			PEdge edge = new PEdge(c1, c2);
			graph.addEdge(c1, c2, edge);
		});

		return graph;
	}

	/**
	 * @param meshFaces collection of faces comprising a conforming mesh.
	 * @return
	 */
	static SimpleGraph<PShape, DefaultEdge> toDualGraph(Collection<PShape> meshFaces) {
		final SimpleGraph<PShape, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
		final Map<PEdge, List<PShape>> edgeFacesMap = new HashMap<>();

		// Phase 1: Collect edges and their associated faces
		for (PShape face : meshFaces) {
			graph.addVertex(face);
			for (int i = 0; i < face.getVertexCount(); i++) {
				PVector a = face.getVertex(i);
				PVector b = face.getVertex((i + 1) % face.getVertexCount());
				if (a.equals(b)) {
					continue;
				}

				PEdge edge = new PEdge(a, b);
				edgeFacesMap.computeIfAbsent(edge, k -> new ArrayList<>()).add(face);
			}
		}

		// Phase 2: Process edges in sorted order for graph iteration consistency
		edgeFacesMap.entrySet().stream().sorted(Comparator.comparing(e -> e.getKey())) // Sort edges to ensure deterministic processing
				.forEach(entry -> {
					List<PShape> faces = entry.getValue();
					if (faces.size() == 2) {
						// If exactly two faces share this edge, connect them in the dual graph
						PShape f1 = faces.get(0);
						PShape f2 = faces.get(1);
						if (!f1.equals(f2)) {
							graph.addEdge(f1, f2); // Avoid self-loops
						} else {
							// Handle case where the same face is associated with the edge twice
							System.err.println("toDualGraph(): Bad input — saw the same edge 3+ times for face: " + f1);
						}
					} else if (faces.size() > 2) {
						// Handle edges shared by more than two faces
						System.err.println("toDualGraph(): Bad input — edge shared by more than two faces: " + entry.getKey().toString());
					}
				});

		return graph;
	}

	/**
	 * Writes the <i>Well-Known Text</i> representation of a shape. The
	 * <i>Well-Known Text</i> format is defined in the OGC Simple Features
	 * Specification for SQL.
	 *
	 * @param shape shape to process
	 * @return a Geometry Tagged Text string
	 * @since 1.3.0
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
	 * @since 1.3.0
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
	 * specified in the OGC <i>Simple Features for SQL specification</i>.
	 *
	 * @param shape shape to process
	 * @return WKB byte representation of shape
	 * @since 1.3.0
	 * @see #fromWKB(byte[])
	 * @see #toHexWKB(PShape)
	 */
	public static byte[] toWKB(PShape shape) {
		WKBWriter writer = new WKBWriter();
		return writer.write(fromPShape(shape));
	}

	/**
	 * Converts a shape into <i>Well-Known Binary</i> format and writes the bytes to
	 * a file.
	 * 
	 * @param shape    shape to process
	 * @param filename Absolute file path (with filename and extension). Prefix with
	 *                 "./" for a relative path.
	 * @since 1.4.0
	 */
	public static void toWKB(PShape shape, String filename) {
		WKBWriter writer = new WKBWriter();
		byte[] bytes = writer.write(fromPShape(shape));
		try {
			FileUtils.writeByteArrayToFile(new File(filename), bytes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Converts a geometry in <i>Well-Known Binary</i> format into a PShape.
	 *
	 * @param shapeWKB byte representation of shape to process
	 * @return a PShape specified by the WKB
	 * @since 1.3.0
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
	 * Reads a shape from a (binary) file containing the <i>Well-Known Binary</i>
	 * representation of it.
	 * 
	 * @param filename Absolute file path (with filename and extension). Prefix with
	 *                 "./" for a relative path.
	 * @return a PShape specified by the WKB in the file
	 */
	public static PShape fromWKB(String filename) {
		byte[] shapeWKB;
		try {
			shapeWKB = FileUtils.readFileToByteArray(new File(filename));
			WKBReader reader = new WKBReader();
			return toPShape(reader.read(shapeWKB));
		} catch (IOException | ParseException e) {
			e.printStackTrace();
			return new PShape();
		}

	}

	/**
	 * Writes a shape into the hexadecimal string representation of its
	 * <i>Well-Known Binary</i> format.
	 *
	 * @param shape shape to process
	 * @return hexadecimal string representation of shape WKB
	 * @since 1.3.0
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
	 * @since 1.3.0
	 * @see #toWKB(PShape)
	 */
	public static PShape fromHexWKB(String shapeWKB) {
		return fromWKB(WKBReader.hexToBytes(shapeWKB));
	}

	/**
	 * Writes a <b>single holeless</b> shape into the string representation of its
	 * Google <i>Encoded Polyline</i> format.
	 * 
	 * @param shape single (holeless) polygon or line
	 * @return String with the encoded polyline representing the shape
	 * @since 1.3.0
	 */
	public static String toEncodedPolyline(PShape shape) {
		Track line = new Track();
		toPVector(shape).forEach(p -> line.addTrackpoint(new Trackpoint(p.x, p.y)));
		line.addTrackpoint(line.getTrackpoints().get(0)); // close
		// PolylineEncoder.createEncodings() writes to console, so suppress that...
		PrintStream old = System.out;
		System.setOut(new PrintStream(new OutputStream() {
			@Override
			public void write(int b) throws IOException {
			}
		}));
		String encoding = (String) PolylineEncoder.createEncodings(line, 0, 1).get("encodedPoints");
		System.setOut(old);
		return encoding;
	}

	/**
	 * Converts a geometry in <i>Encoded Polyline</i> format into a PShape.
	 * 
	 * @param encodedPolyline an encoded polyline string representing a shape
	 * @return a PShape represented by the encoded polyline string
	 * @since 1.3.0
	 */
	public static PShape fromEncodedPolyline(String encodedPolyline) {
		PolylineDecoder decoder = new PolylineDecoder();
		CoordinateList coords = new CoordinateList();

		decoder.decode(encodedPolyline).forEach(p -> {
			double x = p.getLat();
			double y = p.getLng();
			coords.add(new Coordinate(x, y));
		});

		return toPShape(GEOM_FACTORY.createLineString(coords.toCoordinateArray()));
	}

	/**
	 * Writes a shape into the string representation of its <i>GeoJSON</i> format.
	 * 
	 * @param shape
	 * @return json JSON string
	 * @since 1.3.0
	 */
	public static String toGeoJSON(PShape shape) {
		final GeoJsonWriter writer = new GeoJsonWriter(1);
		writer.setForceCCW(true);
		return writer.write(fromPShape(shape));
	}

	/**
	 * Converts a GeoJSON representation of a shape into its PShape counterpart.
	 * 
	 * @param json GeoJSON string
	 * @return PShape represented by the GeoJSON
	 * @since 1.3.0
	 */
	public static PShape fromGeoJSON(String json) {
		final GeoJsonReader reader = new GeoJsonReader(GEOM_FACTORY);
		try {
			return toPShape(reader.read(json));
		} catch (ParseException e) {
			System.err.println("Error occurred when converting json to shape.");
			return new PShape();
		}
	}

	/**
	 * Creates a Java2D/java.awt Shape representing a PShape.
	 *
	 * @param shape the PShape to convert
	 * @return a Java2D shape representing the PShape
	 * @since 1.3.0
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
	 * @since 1.3.0
	 */

	public static PShape fromJava2D(Shape shape) {
		if (shape != null) {
			PathIterator pathIt = shape.getPathIterator(AffineTransform.getScaleInstance(1, 1), BEZIER_SAMPLE_DISTANCE);
			return toPShape(ShapeReader.read(pathIt, GEOM_FACTORY));
		} else {
			return new PShape();
		}
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
	 * @see #fromContours(List, List)
	 * @see #toPointsPShape(Collection)
	 */
	public static PShape fromPVector(Collection<PVector> vertices) {
		List<PVector> verticesList = new ArrayList<>(vertices);
		boolean closed = false;
		if (!vertices.isEmpty() && verticesList.get(0).equals(verticesList.get(vertices.size() - 1))) {
			closed = true;
		}

		PShape shape = new PShape();
		shape.setFamily(PShape.PATH);
		shape.setFill(Colors.WHITE);
		shape.setFill(closed);
		shape.setStroke(true);
		shape.setStroke(closed ? Colors.PINK : Colors.WHITE);
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
	 * Generates a shape from a list of vertices. If the list of vertices is closed
	 * (first and last vertices are the same), the vertices are interpreted as a
	 * closed polygon (having no holes); if the list is unclosed, they are treated
	 * as a linestring.
	 *
	 * @param vertices list of (un)closed shape vertices
	 * @return a PATH PShape (either open linestring or closed polygon)
	 * @see #fromPVector(PVector...)
	 */
	public static PShape fromPVector(PVector... vertices) {
		return fromPVector(Arrays.asList(vertices));
	}

	/**
	 * Generates a polygonal shape from lists of vertices representing its shell and
	 * holes.
	 * <p>
	 * According to Processing shape definitions, the shell should be orientated CW
	 * and the holes CCW, but this method will detect orientation and handle it
	 * accordingly, so the orientation of input does not matter.
	 * 
	 * @param shell vertices of the shell of the polygon
	 * @param holes (optional) list of holes
	 * @since 1.4.0
	 */
	public static PShape fromContours(List<PVector> shell, @Nullable List<List<PVector>> holes) {
		boolean closed = false;
		if (!shell.isEmpty() && shell.get(0).equals(shell.get(shell.size() - 1))) {
			closed = true;
		}

		PShape shape = new PShape();
		shape.setFamily(PShape.PATH);
		shape.setFill(true);
		shape.setFill(Colors.WHITE);
		shape.setStroke(true);
		shape.setStroke(Colors.PINK);
		shape.setStrokeWeight(4);

		shape.beginShape();
		if (!PGS.isClockwise(shell)) {
			Collections.reverse(shell);
		}
		for (int i = 0; i < shell.size() - (closed ? 1 : 0); i++) {
			PVector v = shell.get(i);
			shape.vertex(v.x, v.y);
		}

		if (holes != null) {
			holes.forEach(hole -> {
				if (hole.size() < 3) {
					return;
				}
				final boolean holeClosed = hole.get(0).equals(hole.get(hole.size() - 1));
				if (PGS.isClockwise(hole)) {
					Collections.reverse(hole);
				}
				shape.beginContour();
				for (int i = 0; i < hole.size() - (holeClosed ? 1 : 0); i++) {
					PVector v = hole.get(i);
					shape.vertex(v.x, v.y);
				}
				shape.endContour();
			});
		}

		shape.endShape(PConstants.CLOSE);

		return shape;
	}

	/**
	 * Converts a simple PShape into an array of its coordinates.
	 * 
	 * @param shape      a simple shape (closed polygon or open line) represented by
	 *                   a coordinate array [[x1, y1], [x2, y2]...]
	 * @param keepClosed flag to determine whether to keep the (last) closing vertex
	 *                   in the output if the input forms a closed polygon
	 * @return coordinate array in the form [[x1, y1], [x2, y2]]
	 * @since 1.4.0 an array of coordinates representing the PShape
	 * @see #fromArray(double[][], boolean)
	 */
	public static double[][] toArray(PShape shape, boolean keepClosed) {
		List<PVector> points = toPVector(shape); // unclosed
		if (shape.isClosed() && keepClosed) {
			points.add(points.get(0).copy()); // since toPVector returns unclosed view
		}
		double[][] out = new double[points.size()][2];
		for (int i = 0; i < points.size(); i++) {
			PVector point = points.get(i);
			out[i][0] = point.x;
			out[i][1] = point.y;
		}
		return out;
	}

	/**
	 * Creates a PShape from an array of doubles representing coordinates.
	 * 
	 * @param shape coordinate array representing a simple shape (closed polygon or
	 *              open line) [[x1, y1], [x2, y2]...]
	 * @param close close the coordinates (if unclosed)
	 * @return a PShape represented by the coordinates
	 * @since 1.4.0
	 * @see #toArray(PShape, boolean)
	 */
	public static PShape fromArray(double[][] shape, boolean close) {
		List<PVector> points = new ArrayList<>(shape.length);
		for (double[] p : shape) {
			points.add(new PVector((float) p[0], (float) p[1]));
		}
		// add closing vertex if close==true and data isn't already closed
		if (close && !points.get(0).equals(points.get(points.size() - 1))) {
			points.add(new PVector((float) shape[0][0], (float) shape[0][1]));
		}
		return fromPVector(points);
	}

	/**
	 * Flattens a collection of PShapes into a single GROUP PShape which has the
	 * input shapes as its children. If the collection contains only one shape, it
	 * is returned directly as a non-GROUP shape.
	 *
	 * @since 1.2.0
	 * @see #flatten(PShape...)
	 */
	public static PShape flatten(Collection<PShape> shapes) {
		PShape group = new PShape(GROUP);
		shapes.stream().filter(Objects::nonNull).forEach(group::addChild);
		if (group.getChildCount() == 1) {
			return group.getChild(0);
		}
		return group;
	}

	/**
	 * Flattens a collection of PShapes into a single GROUP PShape which has the
	 * input shapes as its children.
	 *
	 * @since 1.3.0
	 * @see #flatten(Collection)
	 */
	public static PShape flatten(PShape... shapes) {
		return flatten(Arrays.asList(shapes));
	}

	/**
	 * Recurses a GROUP PShape, finding <b>all</b> of its non-GROUP child PShapes.
	 * If the shape type is not GROUP, only it is returned.
	 * <p>
	 * Note: this method differs from {@link processing.core.PShape#getChildren()
	 * PShape.getChildren()}. That method will return GROUP child shapes, whereas
	 * this method will recurse such shapes, returing their non-group children (in
	 * other words, this method explores the whole tree of shapes, returning
	 * non-group shapes only).
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
				for (PShape candidate : parent.getChildren()) {
					if (candidate.getFamily() == GROUP) {
						parents.add(candidate);
					} else {
						children.add(candidate);
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
	public static PShape fromChildren(Collection<PShape> children) {
		final PShape parent = new PShape(GROUP);
		children.forEach(parent::addChild);
		return parent;
	}

	/**
	 * Retrieves the fill color of a PShape.
	 * 
	 * @param shape The PShape object for which to retrieve the fill color.
	 * @return The integer representation of the fill color in ARGB format (32-bit).
	 * @since 1.4.0
	 */
	public static int getFillColor(PShape shape) {
		try {
			return PSHAPE_FILL_FIELD.getInt(shape);
		} catch (IllegalArgumentException | IllegalAccessException e) {
			e.printStackTrace();
			return 0;
		}
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
	 * @see #setAllStrokeColor(PShape, int, double, int)
	 * @see #setAllFillColor(PShape, int)
	 */
	public static PShape setAllStrokeColor(PShape shape, int color, double strokeWeight) {
		getChildren(shape).forEach(child -> {
			child.setStroke(true);
			child.setStroke(color);
			child.setStrokeWeight((float) strokeWeight);
		});
		return shape;
	}

	/**
	 * Sets the stroke color and cap style for the PShape and all of its children
	 * recursively.
	 *
	 * @param strokeCap either <code>SQUARE</code>, <code>PROJECT</code>, or
	 *                  <code>ROUND</code>
	 * @return the input object (having now been mutated)
	 * @see #setAllStrokeColor(PShape, int, double)
	 * @see #setAllFillColor(PShape, int)
	 * @since 2.1
	 */
	public static PShape setAllStrokeColor(PShape shape, int color, double strokeWeight, int strokeCap) {
		getChildren(shape).forEach(child -> {
			child.setStroke(true);
			child.setStroke(color);
			child.setStrokeWeight((float) strokeWeight);
			child.setStrokeCap(strokeCap);
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
	 * Sets the stroke color equal to the fill color, and the strokeWeight to the
	 * specified value, for the PShape and all of its descendent shapes individually
	 * (that is, each child shape belonging to the shape (if any) will have its
	 * stroke color set to <b>its own fill color</b>, and not the parent-most
	 * shape's fill color).
	 *
	 * @param shape
	 * @return the input object (having now been mutated)
	 * @since 1.3.0
	 */
	public static PShape setAllStrokeToFillColor(PShape shape, double strokeWeight) {
		getChildren(shape).forEach(child -> {
			child.setStroke(true);
			child.setStrokeWeight((float) strokeWeight);
			child.setStroke(PGS.getPShapeFillColor(child));
		});
		return shape;
	}

	/**
	 * Retrieves the styling data associated with the specified PShape object. The
	 * method creates an instance of PShapeData containing the fill, stroke, fill
	 * color, stroke color, and stroke weight extracted from the PShape.
	 *
	 * @param shape the PShape instance from which to extract styling information.
	 * @return a PShapeData object with the extracted styling properties of the
	 *         shape.
	 * @since 2.0
	 */
	public static PShapeData getShapeStylingData(PShape shape) {
		return new PShapeData(shape);
	}

	/**
	 * Reorders the child shapes of a given shape.
	 * <p>
	 * Creates a new GROUP shape, having the same children as the input, but in a
	 * different order; child shapes of the new shape are ordered according to the
	 * given comparator.
	 * 
	 * @param shape      a GROUP shape
	 * @param comparator PShape comparison function
	 * @return a new GROUP PShape object having its children in a different order.
	 *         Child shapes reference the same objects as the input.
	 * @since 1.4.0
	 */
	public static PShape reorderChildren(PShape shape, Comparator<PShape> comparator) {
		List<PShape> children = getChildren(shape);
		children.sort(comparator);
		return flatten(children);
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
	 * to the shape. This can sometimes fix a visual problem in Processing where
	 * narrow gaps can appear between otherwise flush shapes. If the shape is a
	 * GROUP, the rounding is applied to all child shapes.
	 *
	 * @param shape the PShape to round vertex coordinates for.
	 * @return a rounded copy of the input shape.
	 * @see #roundVertexCoords(PShape, int)
	 * @since 1.1.3
	 */
	public static PShape roundVertexCoords(PShape shape) {
		return roundVertexCoords(shape, 0);
	}

	/**
	 * Rounds the x and y coordinates (to <code>n</code> decimal places) of all
	 * vertices belonging to the shape. This can sometimes fix a visual problem in
	 * Processing where narrow gaps can appear between otherwise flush shapes. If
	 * the shape is a GROUP, the rounding is applied to all child shapes.
	 *
	 * @param shape the PShape to round vertex coordinates for.
	 * @param n     The number of decimal places to which the coordinates should be
	 *              rounded.
	 * @return a rounded copy of the input shape.
	 * @since 2.1
	 */
	public static PShape roundVertexCoords(PShape shape, int n) {
		return PGS_Processing.transform(shape, s -> {
			var c = copy(s);
			for (int i = 0; i < c.getVertexCount(); i++) {
				final PVector v = c.getVertex(i);
				c.setVertex(i, round(v.x, n), round(v.y, n));
			}
			return c;
		});
	}

	/**
	 * Produces a deep copy / clone of the input shape. Handles GROUP, PRIMITIVE,
	 * GEOMETRY and PATH PShapes. Clones both geometry and styling.
	 *
	 * @param shape the PShape to copy
	 * @return a deep copy of the given shape
	 * @since 1.2.0
	 */
	public static PShape copy(PShape shape) {
		final PShape copy = new PShape();
		copy.setName(shape.getName());
		final PShapeData style = new PShapeData(shape);

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
					if (shape.getKind() == PConstants.POLYGON) { // kind = POLYGON by default
						copy.setFamily(PShape.PATH);
						method = PShape.class.getDeclaredMethod("copyPath", PShape.class, PShape.class);
					} else {
						copy.setFamily(PShape.GEOMETRY);
						method = PShape.class.getDeclaredMethod("copyGeometry", PShape.class, PShape.class);
					}
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
		} catch (NoSuchMethodException | SecurityException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
			e.printStackTrace();
		}

		return style.applyTo(copy);
	}

	/**
	 * Creates a PATH PShape representing a quadratic bezier curve, given by its
	 * parameters.
	 *
	 * @param start        starting point of the bezier curve
	 * @param controlPoint control point of the curve
	 * @param end          end point of the bezier curve
	 * @return a PShape representing the quadratic bezier curve as a PATH, sampled
	 *         every 2 units along its length
	 * @since 1.4.0
	 */
	public static PShape fromQuadraticBezier(PVector start, PVector controlPoint, PVector end) {
		// convert to cubic bezier form
		PVector cp1 = start.copy().add(controlPoint.copy().sub(start).mult(2 / 3f));
		PVector cp2 = end.copy().add(controlPoint.copy().sub(end).mult(2 / 3f));
		return fromCubicBezier(start, cp1, cp2, end);
	}

	/**
	 * Creates a PATH PShape representing a cubic bezier curve, given by its
	 * parameters.
	 *
	 * @param start         starting point of the bezier curve
	 * @param controlPoint1 first control point of the curve
	 * @param controlPoint2 second control point of the curve
	 * @param end           end point of the bezier curve
	 * @return a PShape representing the cubic bezier curve as a PATH, sampled every
	 *         2 units along its length
	 * @since 1.4.0
	 */
	public static PShape fromCubicBezier(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end) {
		CubicBezier bezier = new CubicBezier(start.x, start.y, controlPoint1.x, controlPoint1.y, controlPoint2.x, controlPoint2.y, end.x, end.y);
		double[][] samples = bezier.sampleEquidistantPoints(BEZIER_SAMPLE_DISTANCE);
		final List<PVector> coords = new ArrayList<>(samples.length);
		for (double[] sample : samples) {
			coords.add(new PVector((float) sample[0], (float) sample[1]));
		}
		return fromPVector(coords);
	}

	/**
	 * Subdivide/interpolate/discretise along a quadratic bezier curve, given by its
	 * start, end and control points
	 *
	 * @return list of points along curve
	 */
	private static List<PVector> getQuadraticBezierPoints(PVector start, PVector controlPoint, PVector end, float sampleDistance) {
		// convert to cubic form
		PVector cp1 = start.copy().add(controlPoint.copy().sub(start).mult(2 / 3f));
		PVector cp2 = end.copy().add(controlPoint.copy().sub(end).mult(2 / 3f));
		return getCubicBezierPoints(start, cp1, cp2, end, sampleDistance);
	}

	/**
	 * Generates a list of equidistant samples along a cubic bezier curve.
	 *
	 * @param sampleDistance distance between successive samples on the curve
	 * @return
	 */
	private static List<PVector> getCubicBezierPoints(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end, float sampleDistance) {
		CubicBezier bezier = new CubicBezier(start.x, start.y, controlPoint1.x, controlPoint1.y, controlPoint2.x, controlPoint2.y, end.x, end.y);
		double[][] samples = bezier.sampleEquidistantPoints(sampleDistance);
		final List<PVector> coords = new ArrayList<>(samples.length);
		for (double[] sample : samples) {
			coords.add(new PVector((float) sample[0], (float) sample[1]));
		}

		return coords;
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

		for (int vertexCode : rawVertexCodes) {
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

	private static <T> T[] reversedCopy(T[] original) {
		@SuppressWarnings("unchecked")
		T[] reversed = (T[]) Array.newInstance(original.getClass().getComponentType(), original.length);
		for (int i = 0; i < original.length; i++) {
			reversed[i] = original[original.length - 1 - i];
		}
		return reversed;
	}

	private static float round(float x, float n) {
		float m = (float) FastMath.pow(10, n);

		return FastMath.floor(m * x + 0.5f) / m;
	}

	/**
	 * A utility class for storing and manipulating the visual properties of PShapes
	 * from the Processing library. It encapsulates the stroke, fill, stroke color,
	 * stroke weight, and fill color attributes by directly accessing and modifying
	 * the corresponding fields of a given PShape using reflection.
	 */
	public static class PShapeData {

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

		public int fillColor, strokeColor;
		public float strokeWeight;
		public boolean fill, stroke;
		private final PShape source;

		PShapeData(PShape shape) {
			source = shape;
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
		 * @return other (fluent interface)
		 */
		public PShape applyTo(PShape other) {
			if (source.getFamily() == GROUP && other.getFamily() != GROUP) {
				// opinionated -- don't apply group styling to non-group
				return other;
			}

			if (other.getFamily() == GROUP) {
				// ONLY IF child fill/stroke aren't defaults!?
				getChildren(other).forEach(c -> applyTo(c));
			}
			other.setFill(fill);
			other.setFill(fillColor);
			other.setStroke(stroke);
			other.setStroke(strokeColor);
			other.setStrokeWeight(strokeWeight);

			return other;
		}

		@Override
		public String toString() {
			return String.format("fillColor: %s; strokeColor: %s; strokeWeight: %.1f", Arrays.toString(decomposeclrRGB(fillColor)),
					Arrays.toString(decomposeclrRGB(strokeColor)), strokeWeight);
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (fill ? 1231 : 1237);
			result = prime * result + fillColor;
			result = prime * result + (stroke ? 1231 : 1237);
			result = prime * result + strokeColor;
			result = prime * result + Float.floatToIntBits(strokeWeight);
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null || getClass() != obj.getClass()) {
				return false;
			}
			PShapeData other = (PShapeData) obj;
			return fillColor == other.fillColor && strokeColor == other.strokeColor && strokeWeight == other.strokeWeight && fill == other.fill
					&& stroke == other.stroke;
		}
	}
}
