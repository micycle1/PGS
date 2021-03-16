package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import static micycle.pts.PTS.GEOM_FACTORY;
import static processing.core.PConstants.TRIANGLES;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.triangulate.ConformingDelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.IncrementalDelaunayTriangulator;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;

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
//		arrCoords = new double[] { 0, 0, 100, 0, 100, 100, 0, 100, 20, 20, 80, 20, 80, 80, 20, 80 };

		List<Integer> triangles = Earcut.earcut(arrCoords, null, 2);

		PShape triangulation = new PShape();
		triangulation.setFamily(PShape.GEOMETRY);
//		triangulation.setStrokeCap(ROUND);
		triangulation.setStroke(true);
		triangulation.setStrokeWeight(2);
		triangulation.setStroke(-123222);
		triangulation.setFill(true);
		triangulation.setFill(micycle.pts.color.RGB.composeclr(0, 0, 0, 25));
		triangulation.beginShape(TRIANGLES);
		for (int i = 0; i < triangles.size(); i += 3) {
			int v1 = triangles.get(i);
			int v2 = triangles.get(i + 1);
			int v3 = triangles.get(i + 2);
//			triangulation.vertex((float) arrCoords[v1], (float) arrCoords[v1 + 1]);
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
		return toPShape(refinedTriangulation(fromPShape(shape), 3, 10));

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
		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
		Geometry g = fromPShape(shape);
		b.setSites(g);
		b.setTolerance(tolerance);
		b.setConstraints(fromPShape(constraints));
		Geometry out = b.getTriangles(GEOM_FACTORY); // triangulates concave hull of points
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

}
