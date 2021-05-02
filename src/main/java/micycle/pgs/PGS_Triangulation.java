package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static processing.core.PConstants.TRIANGLES;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Location;
import org.tinfour.common.IConstraint;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;

import earcut4j.Earcut;
import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.color.RGB;
import micycle.pgs.utility.PoissonDistribution;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Delaunay and earcut triangulation of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Triangulation {

	private PGS_Triangulation() {
	}

	/**
	 * Performs constrained + refined Delaunay Triangulation of a shape.
	 * 
	 * @param shape
	 * @param steinerPoints A list of additional points to triangulate in addition
	 *                      to the input shape. Can be null.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to both maintain the Delaunay nature when
	 *                      constraining and check circumcircle location during
	 *                      refinement. When true, triangles in the triangulation
	 *                      may be slightly more regular in shape/size. There is a
	 *                      small performance overhead which becomes more
	 *                      considerable at higher refinement levels.
	 * @return
	 * @see #delaunayTriangulationPoints(PShape, List, boolean, int, boolean)
	 * @see #delaunayTriangulationTin(PShape, List, boolean, int, boolean)
	 */
	public static PShape delaunayTriangulation(PShape shape, List<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		IncrementalTin tin = delaunayTriangulationTin(shape, steinerPoints, constrain, refinements, pretty);

		PShape triangulation = new PShape(PShape.GEOMETRY);

		PGS_Conversion.setAllFillColor(triangulation, RGB.composeColor(50, 50, 50));
		PGS_Conversion.setAllStrokeColor(triangulation, RGB.PINK, 2);
		triangulation.beginShape(TRIANGLES);

		Consumer<Vertex[]> triangleVertexConsumer = t -> {
			triangulation.vertex((float) t[0].x, (float) t[0].y);
			triangulation.vertex((float) t[1].x, (float) t[1].y);
			triangulation.vertex((float) t[2].x, (float) t[2].y);
		};
		if (constrain) {
			TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		} else {
			TriangleCollector.visitTriangles(tin, triangleVertexConsumer);
		}

		triangulation.endShape();
		return triangulation;
	}

	/**
	 * Performs Constrained + refined Delaunay Triangulation of a shape. This method returns
	 * the triangulation as a list of points.
	 * 
	 * @param shape
	 * @param steinerPoints A list of additional points to triangulate in addition
	 *                      to the input shape. Can be null.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to both maintain the Delaunay nature when
	 *                      constraining and check circumcircle location during
	 *                      refinement. When true, triangles in the triangulation
	 *                      may be slightly more regular in shape/size. There is a
	 *                      small performance overhead which becomes more
	 *                      considerable at higher refinement levels.
	 * @return List of coordinates, where each consecutive triplet of coordinates
	 *         are the vertices belonging to one triangle
	 * @see #delaunayTriangulationPoints(PShape, List, boolean, int, boolean)
	 * @see #delaunayTriangulationTin(PShape, List, boolean, int, boolean)
	 */
	public static List<PVector> delaunayTriangulationPoints(PShape shape, List<PVector> steinerPoints,
			boolean constrain, int refinements, boolean pretty) {
		IncrementalTin tin = delaunayTriangulationTin(shape, steinerPoints, constrain, refinements, pretty);

		ArrayList<PVector> triangles = new ArrayList<>();
		Consumer<Vertex[]> triangleVertexConsumer = t -> {
			triangles.add(toPVector(t[0]));
			triangles.add(toPVector(t[1]));
			triangles.add(toPVector(t[2]));
		};
		if (constrain) {
			TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		} else {
			TriangleCollector.visitTriangles(tin, triangleVertexConsumer);
		}
		return triangles;
	}

	/**
	 * Performs Delaunay Triangulation of a shape. This method returns the raw Tinfour
	 * triangulated network object.
	 * 
	 * @param shape
	 * @param steinerPoints A list of additional points to triangulate in addition
	 *                      to the input shape. Can be null.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to both maintain the Delaunay nature when
	 *                      constraining and check circumcircle location during
	 *                      refinement. When true, triangles in the triangulation
	 *                      may be slightly more regular in shape/size. There is a
	 *                      small performance overhead which becomes more
	 *                      considerable at higher refinement levels.
	 * @return Triangulated Irregular Network object
	 * @see #delaunayTriangulation(PShape, List, boolean, int, boolean)
	 * @see #delaunayTriangulationPoints(PShape, List, boolean, int, boolean)
	 */
	public static IncrementalTin delaunayTriangulationTin(PShape shape, List<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		final Geometry g = fromPShape(shape);
		final IncrementalTin tin = new IncrementalTin(10);

		final ArrayList<Vertex> vertices = new ArrayList<>();
		Coordinate[] coords = g.getCoordinates();
		for (int i = 0; i < coords.length; i++) {
//			PVector v = shape.getVertex(i);
			vertices.add(new Vertex(coords[i].x, coords[i].y, 0));
		}

		tin.add(vertices, null); // initial triangulation

		if (steinerPoints != null) {
			steinerPoints.forEach(v -> tin.add(new Vertex(v.x, v.y, 0))); // add steiner points
		}

		if (refinements > 0) {

			final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
			final ArrayList<Vertex> refinementVertices = new ArrayList<>();
			final Consumer<SimpleTriangle> triangleConsumer = t -> {
				if (t.getArea() > 85) { // don't refine small triangles
					final Coordinate center = centroid(t); // use centroid rather than circumcircle center
					if (pretty || pointLocator.locate(center) != Location.EXTERIOR) {
						refinementVertices.add(new Vertex(center.x, center.y, 0));
					}
				}
			};

			/**
			 * Possible optimisation is to recursely split within each triangle upto the
			 * refinement depth (in one pass), so perform many less location checks. Another
			 * is to rasterise the PShape and check pixel[] array.
			 */
			for (int i = 0; i < refinements; i++) {
				refinementVertices.clear();
				TriangleCollector.visitSimpleTriangles(tin, triangleConsumer);
				tin.add(refinementVertices, null); // add refinement (steiner) points
			}
		}

		if (constrain) {
			List<IConstraint> constraints = new ArrayList<>();
			for (int n = 0; n < g.getNumGeometries(); n++) {
				LinearRingIterator lri = new LinearRingIterator(g.getGeometryN(n));
				lri.forEach(ring -> {
					ArrayList<Vertex> points = new ArrayList<>();
					Coordinate[] c = ring.getCoordinates();
					if (Orientation.isCCW(c)) {
						for (int i = 0; i < c.length; i++) {
							points.add(new Vertex(c[i].x, c[i].y, 0));
						}
					} else {
						for (int i = c.length - 1; i >= 0; i--) { // iterate backwards if CW
							points.add(new Vertex(c[i].x, c[i].y, 0));
						}
					}
					constraints.add(new PolygonConstraint(points));
				});
			}
			tin.addConstraints(constraints, pretty); // true/false is negligible?
		}

		return tin;
	}

	/**
	 * Creates a Delaunay triangulation with steiner points within the shape
	 * populated by poisson sampling.
	 * 
	 * @param shape
	 * @param spacing (Minimum) spacing between poisson points
	 * @return
	 */
	public static List<PVector> poissonTriangulation(PShape shape, double spacing) {

		final Geometry g = fromPShape(shape).buffer(-spacing);
		Envelope e = (Envelope) g.getEnvelopeInternal();

		PoissonDistribution pd = new PoissonDistribution(0);
		List<PVector> poissonPoints = pd.generate(e.getMinX(), e.getMinY(), e.getMinX() + e.getWidth(),
				e.getMinY() + e.getHeight(), spacing, 8);
//		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
//		List<PVector> pp = poissonPoints.parallelStream()
//				.filter(p -> pointLocator.locate(PGS.coordFromPVector(p)) != Location.EXTERIOR)
//				.collect(Collectors.toList());

		IncrementalTin tin = delaunayTriangulationTin(shape, poissonPoints, true, 0, true);

		ArrayList<PVector> triangles = new ArrayList<>();
		Consumer<Vertex[]> triangleVertexConsumer = t -> {
			triangles.add(toPVector(t[0]));
			triangles.add(toPVector(t[1]));
			triangles.add(toPVector(t[2]));
		};
		TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		return triangles;
	}

	/**
	 * Computes a trianglation of the points according to the ear clipping ("earcut") method.
	 * @param points
	 * @return
	 */
	public static PShape earCutTriangulation(List<PVector> points) {
		double[] arrCoords = new double[points.size() * 2];

		for (int i = 0; i < points.size(); i++) {
			arrCoords[2 * i] = points.get(i).x;
			arrCoords[2 * i + 1] = points.get(i).y;
		}

		List<Integer> triangles = Earcut.earcut(arrCoords, null, 2);

		PShape triangulation = new PShape();
		triangulation.setFamily(PShape.GEOMETRY);
		// triangulation.setStrokeCap(ROUND);
		triangulation.setStroke(true);
		triangulation.setStrokeWeight(2);
		triangulation.setStroke(RGB.PINK);
		triangulation.setFill(true);
		triangulation.setFill(micycle.pgs.color.RGB.composeColor(255, 255, 255, 255));
		triangulation.beginShape(TRIANGLES);
		for (int i = 0; i < triangles.size(); i += 3) {
			final int v1 = 2 * triangles.get(i);
			final int v2 = 2 * triangles.get(i + 1);
			final int v3 = 2 * triangles.get(i + 2);
			triangulation.vertex((float) arrCoords[v1], (float) arrCoords[v1 + 1]);
			triangulation.vertex((float) arrCoords[v2], (float) arrCoords[v2 + 1]);
			triangulation.vertex((float) arrCoords[v3], (float) arrCoords[v3 + 1]);
		}
		triangulation.endShape();
		return triangulation;
	}

	private static PVector toPVector(Vertex v) {
		return new PVector((float) v.getX(), (float) v.getY());
	}

	private static Coordinate centroid(SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new Coordinate(x, y);
	}
}
