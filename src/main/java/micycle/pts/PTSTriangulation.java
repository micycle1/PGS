package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import static micycle.pts.PTS.GEOM_FACTORY;
import static processing.core.PConstants.TRIANGLES;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.function.Consumer;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.triangulate.ConformingDelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.IncrementalDelaunayTriangulator;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;
import org.tinfour.common.IConstraint;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;

import earcut4j.Earcut;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Rename to Mesh?
 * 
 * @author Michael Carleton
 *
 */
public class PTSTriangulation {

	public static PShape earCutTriangulation(ArrayList<PVector> points) {
		double[] arrCoords = new double[points.size() * 2];

		for (int i = 0; i < points.size(); i++) {
			arrCoords[2 * i] = points.get(i).x;
			arrCoords[2 * i + 1] = points.get(i).y;
		}

		List<Integer> triangles = Earcut.earcut(arrCoords, null, 2);

		PShape triangulation = new PShape();
		triangulation.setFamily(PShape.GEOMETRY);
//		triangulation.setStrokeCap(ROUND);
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
		Geometry out = d.getTriangles(GEOM_FACTORY); // triangulates convex hull of points
		out = out.intersection(g); // get concave hull
		return toPShape(out);
	}

	public static List<PVector> delaunayTinFour(PShape shape, PShape constraint) {
		Geometry g = fromPShape(shape);
//		IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
		IncrementalTin tin = new IncrementalTin(10);
		Coordinate[] coords = g.getCoordinates();
		final ArrayList<Vertex> vertices = new ArrayList<>();
		for (int i = 0; i < coords.length; i++) {
			vertices.add(new Vertex(coords[i].x, coords[i].y, 0));
		}
		tin.add(vertices, null); // insert point set; points are triangulated upon insertion
		for (int i = 0; i < constraint.getVertexCount(); i++) {
			PVector v = constraint.getVertex(i);
			tin.add(new Vertex(v.x, v.y, 0));
		}
		Collections.reverse(vertices);
		PolygonConstraint pc = new PolygonConstraint(vertices);
		List<IConstraint> c = new ArrayList<IConstraint>();
		c.add(pc);
		tin.addConstraints(c, true);
//		tin.addConstraints(, false);
		ArrayList<PVector> triangles = new ArrayList<>(coords.length / 2);
		TriangleCollector.visitSimpleTriangles(tin, new Consumer<SimpleTriangle>() {
			@Override
			public void accept(SimpleTriangle t) {
				IConstraint constraint = t.getContainingRegion();
				if (constraint != null && constraint.definesConstrainedRegion()) {
//				if (pointLocator.locate(circumcircle(t)) != Location.EXTERIOR) {
					triangles.add(pVectorFromVertex(t.getVertexA()));
					triangles.add(pVectorFromVertex(t.getVertexB()));
					triangles.add(pVectorFromVertex(t.getVertexC()));
				}
//				}
			}
		});
		return triangles;
	}

	/**
	 * 
	 * @param points
	 * @return array of PVectors. each triplet of vectors are the vertices of one
	 *         triangle
	 */
	public static IncrementalTin delaunayTinFour(List<PVector> points, List<PVector> steinerPoints) {
		// TODO points orientation test (for constraints)
		final IncrementalTin tin = new IncrementalTin(10);
//		Collections.reverse(points);
//		Collections.reverse(steinerPoints);
		final ArrayList<Vertex> vertices = new ArrayList<>();
		for (PVector v : points) {
			vertices.add(new Vertex(v.x, v.y, 0));
		}
		tin.add(vertices, null); // insert point set; points are triangulated upon insertion
		if (steinerPoints != null) {
			for (PVector v : steinerPoints) {
				tin.add(new Vertex(v.x, v.y, 0));
			}
		}

		if (isClockwise(points)) {
			Collections.reverse(vertices); // constraint should be CCW
		}
//		
		PolygonConstraint pc = new PolygonConstraint(vertices);
		List<IConstraint> c = new ArrayList<IConstraint>();
		c.add(pc);
		tin.addConstraints(c, true); // true/false is negligible?
		return tin;
//		ArrayList<PVector> triangles = new ArrayList<>(vertices.size() / 2);
//		TriangleCollector.visitSimpleTriangles(tin, new Consumer<SimpleTriangle>() {
//			@Override
//			public void accept(SimpleTriangle t) {
//				IConstraint constraint = t.getContainingRegion();
//				if (constraint != null && constraint.definesConstrainedRegion()) {
//					triangles.add(pVectorFromVertex(t.getVertexA()));
//					triangles.add(pVectorFromVertex(t.getVertexB()));
//					triangles.add(pVectorFromVertex(t.getVertexC()));
//				}
//			}
//		});
//		return triangles;
	}

	public static IncrementalTin delaunayTinFourNoConstraints(List<PVector> points, List<PVector> steinerPoints) {
		// TODO points orientation test (for constraints)
		final IncrementalTin tin = new IncrementalTin(10);
//		Collections.reverse(points);
//		Collections.reverse(steinerPoints);
		final ArrayList<Vertex> vertices = new ArrayList<>();
		for (PVector v : points) {
			vertices.add(new Vertex(v.x, v.y, 0));
		}
		tin.add(vertices, null); // insert point set; points are triangulated upon insertion
		if (steinerPoints != null) {
			for (PVector v : steinerPoints) {
				tin.add(new Vertex(v.x, v.y, 0));
			}
		}
		return tin;
	}

	static Geometry refinedTriangulation(Geometry g, int nRefinements, double tolerance) {

		DelaunayTriangulationBuilder builder = new DelaunayTriangulationBuilder();
		builder.setSites(g); // set vertex sites
		builder.setTolerance(tolerance); // set tolerance for initial triangulation only

		Geometry triangulation = builder.getTriangles(GEOM_FACTORY); // initial triangulation

		HashSet<Coordinate> sites = new HashSet<>();
		for (int i = 0; i < triangulation.getCoordinates().length; i++) {
			sites.add(triangulation.getCoordinates()[i]);
		}

		for (int refinement = 0; refinement < nRefinements; refinement++) {
			for (int i = 0; i < triangulation.getNumGeometries(); i++) {
				Polygon triangle = (Polygon) triangulation.getGeometryN(i);

				if (triangle.getArea() > 100) { // skip small triangles
					sites.add(new Coordinate(triangle.getCentroid().getX(), triangle.getCentroid().getY()));
				}
			}
			builder = new DelaunayTriangulationBuilder();
			builder.setSites(sites);
			triangulation = builder.getTriangles(GEOM_FACTORY); // re-triangulate using new centroid sites
		}

		triangulation = triangulation.intersection(g); // restore concave hull and any holes
		return triangulation;
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
//		Geometry g = fromPShape(shape);
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
//		d.setTolerance(tolerance);
//		d.setSites(g);
//		Geometry out = d.getTriangles(geometryFactory); // triangulates concave hull of points
//
////		Coordinate
//
////		var x = d.getSubdivision();
//		HashSet<Coordinate> coords = new HashSet<>();
//
//		// add from g
//		for (int i = 0; i < out.getCoordinates().length; i++) {
//			coords.add(out.getCoordinates()[i]);
//		}
//
//		for (int refinement = 0; refinement < 3; refinement++) {
//			// refinement: add new points (centroids)
//			for (int i = 0; i < out.getNumGeometries(); i++) {
//				Polygon triangle = (Polygon) out.getGeometryN(i);
//
//				if (triangle.getArea() > 50) { // skip small triangles
//					coords.add(coordFromPoint(triangle.getCentroid()));
//				}
//			}
//			d = new DelaunayTriangulationBuilder();
//			d.setSites(coords);
//			out = d.getTriangles(geometryFactory); // triangulates concave hull of points
//		}
//
//		// use d.getSubdivision() to get connected to a given segment to test
//		// encroachment
//
//		out = out.intersection(g); // get convex hull
		return toPShape(refinedTriangulation(fromPShape(shape), 3, 0));

//		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
//		b.setTolerance(tolerance);
//		b.setSites(g);
//		b.setConstraints(fromPShape(createSquircle(ThreadLocalRandom.current().nextInt(400, 450), 400, 200, 200)));

//		d.getSubdivision().getVoronoiCellPolygons(geometryFactory)
//		out = b.getTriangles(geometryFactory);
	}

	public static void incrementalDelaunay() {

		// TODO it's the same?
		IncrementalDelaunayTriangulator i = new IncrementalDelaunayTriangulator(new QuadEdgeSubdivision(null, 5));
		i.insertSite(null);
	}

	/**
	 * Create a triangulation with that forces certain required segments into the
	 * triangulation from a 'constrain' shape.
	 * 
	 * Steiner point is a point that is not part of the input to a geometric
	 * optimization problem but is added during the solution of the problem
	 * 
	 * @param shape
	 * @param constraints points from this shape are used as 'steiner points',
	 *                    applied to the input
	 * @param tolerance
	 * @return
	 */
	public static PShape constrainedDelaunayTriangulation(PShape shape, PShape constraints, float tolerance) {
		// TODO point-set array constraints
		// TODO AS list of triangles (Vertices?)
		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
		Geometry g = fromPShape(shape);
		b.setSites(g);
		b.setTolerance(tolerance);
		b.setConstraints(fromPShape(constraints));
		Geometry out = b.getTriangles(GEOM_FACTORY); // triangulates concave hull of points
		out = out.intersection(g); // get convex hull TODO use other method SLOW
		return toPShape(out);
	}

	public static PShape constrainedDelaunayTriangulation(ArrayList<PVector> points, PShape constraints,
			float tolerance) {
		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
		Coordinate[] coords = new Coordinate[points.size()];
		for (int j = 0; j < points.size(); j++) {
			coords[j] = new Coordinate(points.get(j).x, points.get(j).y);
		}
		b.setSites(GEOM_FACTORY.createPolygon(coords));
		b.setTolerance(tolerance);
//		b.setConstraints(fromPShape(constraints)); // TODO?
		Geometry out = b.getTriangles(GEOM_FACTORY); // triangulates concave hull of points
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
		Geometry out = d.getTriangles(GEOM_FACTORY);
		return toPShape(out);
	}

	private static PVector pVectorFromVertex(Vertex v) {
		return new PVector((float) v.getX(), (float) v.getY());
	}

	private static Coordinate circumcircle(SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();

		double D = (a.getX() - c.getX()) * (b.getY() - c.getY()) - (b.getX() - c.getX()) * (a.getY() - c.getY());
		double px = (((a.getX() - c.getX()) * (a.getX() + c.getX()) + (a.getY() - c.getY()) * (a.getY() + c.getY())) / 2
				* (b.getY() - c.getY())
				- ((b.getX() - c.getX()) * (b.getX() + c.getX()) + (b.getY() - c.getY()) * (b.getY() + c.getY())) / 2
						* (a.getY() - c.getY()))
				/ D;

		double py = (((b.getX() - c.getX()) * (b.getX() + c.getX()) + (b.getY() - c.getY()) * (b.getY() + c.getY())) / 2
				* (a.getX() - c.getX())
				- ((a.getX() - c.getX()) * (a.getX() + c.getX()) + (a.getY() - c.getY()) * (a.getY() + c.getY())) / 2
						* (b.getX() - c.getX()))
				/ D;

		return new Coordinate(px, py);
	}

	/**
	 * Requires a closed hole
	 * 
	 * @param points
	 * @return
	 */
	private static boolean isClockwise(List<PVector> points) {
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

}
