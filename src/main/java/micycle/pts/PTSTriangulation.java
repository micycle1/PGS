package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
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
import micycle.pts.PTS.LinearRingIterator;
import micycle.pts.utility.PoissonDistribution;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Delaunay and earcut triangulation of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public class PTSTriangulation {

	/**
	 * Constrained + refined Delaunay Triangulation of a shape.
	 * 
	 * @param shape
	 * @param steinerPoints can be null
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set)
	 * @param refinements   number of triangulation refinements to perform. 0 or
	 *                      greater
	 * @param pretty        false = enforce delaunay triangulation during
	 *                      constraining and don't check circumcircle during
	 *                      refinement. when true output triangles may be slightly
	 *                      more regular in shape/size. small performance overhead
	 *                      -- becomes more considerable at higher refinement levels
	 * @return array of coordinates, each triplet of coordinates are the vertices of
	 *         one triangle
	 */
	public static List<PVector> delaunayTriangulation(PShape shape, List<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		IncrementalTin tin = delaunayTriangulationTin(shape, steinerPoints, constrain, refinements, pretty);

		ArrayList<PVector> triangles = new ArrayList<>();
		if (constrain) {
			TriangleCollector.visitTrianglesConstrained(tin, new Consumer<Vertex[]>() {
				public void accept(Vertex[] t) {
					triangles.add(pVectorFromVertex(t[0]));
					triangles.add(pVectorFromVertex(t[1]));
					triangles.add(pVectorFromVertex(t[2]));
				}
			});
		} else {
			TriangleCollector.visitTriangles(tin, new Consumer<Vertex[]>() {
				public void accept(Vertex[] t) {
					triangles.add(pVectorFromVertex(t[0]));
					triangles.add(pVectorFromVertex(t[1]));
					triangles.add(pVectorFromVertex(t[2]));
				}
			});
		}
		return triangles;
	}

	public static IncrementalTin delaunayTriangulationTin(PShape shape, List<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		final Geometry g = fromPShape(shape);
		final IncrementalTin tin = new IncrementalTin(10);

		final ArrayList<Vertex> vertices = new ArrayList<>();
		for (int i = 0; i < shape.getVertexCount(); i++) {
			PVector v = shape.getVertex(i);
			vertices.add(new Vertex(v.x, v.y, 0));
		}
		tin.add(vertices, null); // initial triangulation

		if (steinerPoints != null) {
			steinerPoints.forEach(v -> tin.add(new Vertex(v.x, v.y, 0))); // add steiner points
		}
		/**
		 * Possible optimisation is to recursely split within each triangle upto the
		 * refinement depth (in one pass), so perform many less location checks. Another
		 * is to rasterise the PShape and check pixel[] array.
		 */
		for (int i = 0; i < refinements; i++) {
			IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
			ArrayList<Vertex> refinementVertices = new ArrayList<>();
			TriangleCollector.visitSimpleTriangles(tin, new Consumer<SimpleTriangle>() {
				public void accept(SimpleTriangle t) {
					if (t.getArea() > 85) { // don't refine small triangles
						// use centroid rather than circumcircle center
						final Coordinate center = centroid(t);
						if (pretty || pointLocator.locate(center) != Location.EXTERIOR) {
							refinementVertices.add(new Vertex(center.x, center.y, 0));
						}
					}
				}
			});
			tin.add(refinementVertices, null); // add refinement (steiner) points
		}

		if (constrain) {
			List<IConstraint> constraints = new ArrayList<>();
			LinearRingIterator lri = new LinearRingIterator(g);
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
			tin.addConstraints(constraints, pretty); // true/false is negligible?
		}
		return tin;
	}

	/**
	 * Creates a triangulation with steiner points within the shape populated by
	 * poisson sampling.
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
				e.getMinY() + e.getHeight(),
				spacing, 10);
//		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
//		List<PVector> pp = poissonPoints.parallelStream()
//				.filter(p -> pointLocator.locate(PTS.coordFromPVector(p)) != Location.EXTERIOR)
//				.collect(Collectors.toList());

		IncrementalTin tin = delaunayTriangulationTin(shape, poissonPoints, true, 0, true);

		ArrayList<PVector> triangles = new ArrayList<>();
		TriangleCollector.visitTrianglesConstrained(tin, new Consumer<Vertex[]>() {
			public void accept(Vertex[] t) {
				triangles.add(pVectorFromVertex(t[0]));
				triangles.add(pVectorFromVertex(t[1]));
				triangles.add(pVectorFromVertex(t[2]));
			}
		});
		return triangles;
	}

	public static PShape earCutTriangulation(ArrayList<PVector> points) {
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
		triangulation.setStroke(-123222);
		triangulation.setFill(true);
		triangulation.setFill(micycle.pts.color.RGB.composeclr(255, 255, 255, 255));
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

	private static PVector pVectorFromVertex(Vertex v) {
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
