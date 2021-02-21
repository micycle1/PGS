package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.quadedge.QuadEdge;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;
import org.locationtech.jts.triangulate.quadedge.Vertex;

import processing.core.PShape;

/**
 * Delaunay triangulation refinement for JTS geometry using Ruppert's method.
 * <p>
 * The central operation of Chew’s and Ruppert’s Delaunay refinement algorithms
 * is the insertion of a vertex at the circumcenter of an element of poor
 * quality. y. The Delaunay property is maintained, using Lawson’s algorithm or
 * the Bowyer– Watson algorithm for the incremental update of Delaunay
 * triangulations.
 * 
 * <p>
 * "poor quality" includes only triangles that have a circumradius-to-shortest
 * edge ratio larger than some appropriate bound B. Ruppert’s algorithm employs
 * a bound of B = √ 2, and Chew’s second refinement algorithm employs a bound of
 * B = 1.
 * 
 * @author MCarleton
 *
 */
public class RuppertTriangulation {

	// consider Triangle.circum() etc...
	// also see https://people.eecs.berkeley.edu/~jrs/meshpapers/delnotes.pdf
	// https://people.eecs.berkeley.edu/~jrs/papers/2dj.pdf

	public static PShape delaunayTriangulation(PShape shape) {
		Geometry g = fromPShape(shape);
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setSites(g);
		Geometry out = d.getTriangles(PTS.GEOM_FACTORY); // triangulates convex hull of points
//		d.getSubdivision().fr
//		var z = ((Polygon) g).getExteriorRing(); // segments
		out = out.intersection(g); // get concave hull
		return toPShape(out);
	}

	// @formatter:off
	/**
	 * 1. Compute delaunay triangulation first of all
	 * 2. measure triangles for "poorness"
	 * 3. if poor T: if circle of T encroaches, then split S
	 * 		else: 
	 * 4. insert points into triangulation (either retriangulate or incremental update using Lawson’s algorithm [19] or the Bowyer–Watson algorithm .
	 */
	// @formatter:on
	public static PShape ruppert(PShape shape, float x, float y) {
		Geometry g = fromPShape(shape);
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setSites(g);

		ArrayList<Coordinate> coords = new ArrayList<>();
		// add from g
		for (int i = 0; i < g.getCoordinates().length; i++) {
			coords.add(g.getCoordinates()[i]);
		}

		Geometry out = d.getTriangles(PTS.GEOM_FACTORY);

		int r = 0;

		QuadEdgeSubdivision qed = d.getSubdivision();
//		qed.

		for (int i = 0; i < out.getNumGeometries(); i++) {
			Coordinate[] z = out.getGeometryN(i).getCoordinates();
			double[] ssa = smallestSideAndAngle(z[0], z[1], z[2]);
			if ((ssa[1] * (180 / Math.PI)) < 40 && ssa[0] > 10) {
				double[] circumcircle = circumcircle(z[0], z[1], z[2]);
				qed.insertSite(new Vertex(circumcircle[0], circumcircle[1]));
				coords.add(new Coordinate(circumcircle[0], circumcircle[1]));
				d = new DelaunayTriangulationBuilder();
				d.setSites(coords);
				out = d.getTriangles(PTS.GEOM_FACTORY); // retriangulates concave hull of points (After one change)
//						coords.clear();
//						for (int j = 0; j < out.getCoordinates().length; j++) {
//							coords.add(out.getCoordinates()[j]);
//						}
				r++;
//						break loop;
			}

		}

//		for (int refinement = 0; refinement < 1; refinement++) {
//			// refinement: add new points (centroids)
//			for (int i = 0; i < out.getNumGeometries(); i++) {
////				if (out.getGeometryN(i).getArea() > tolerance) {
//				var z = out.getGeometryN(i).getCoordinates();
//				double[] circumcircle = circumcircle(z[0], z[1], z[2]);
//				double[] ssa = smallestSideAndAngle(z[0], z[1], z[2]);
//				if ((ssa[1] * (180 / Math.PI)) < 60) {
//					coords.add(new Coordinate(circumcircle[0], circumcircle[1]));
//				}
////				coords.add(PTS.coordFromPoint(out.getGeometryN(i).getCentroid()));
////				}
//			}
//			d = new DelaunayTriangulationBuilder();
//			d.setSites(coords);
//			out = d.getTriangles(PTS.geometryFactory); // triangulates concave hull of points
//			coords.clear();
//			for (int i = 0; i < g.getCoordinates().length; i++) {
//				coords.add(g.getCoordinates()[i]);
//			}
//		}

//		qed.delete(
//				(QuadEdge) qed.getPrimaryEdges(false)
//						.get(ThreadLocalRandom.current().nextInt(qed.getEdges().size() - 1)));
//		for (int i = 0; i < 2; i++) {
//			QuadEdgeSubdivision qed = d.getSubdivision();
//		ArrayList<Vertex[]> verts = (ArrayList<Vertex[]>) qed.getTriangleVertices(true);
//		for (Vertex[] c : verts) {
//			double[] circumcircle = circumcircle(c[0], c[1], c[2]);
//			double[] ssa = smallestSideAndAngle(c[0], c[1], c[2]);
////			System.out.println(ssa[1] * (180 / Math.PI));
////			if ((ssa[1] * (180 / Math.PI)) < 40) {
//			coords.add(new Coordinate(circumcircle[0], circumcircle[1]));
////			}
////			try {
////				qed.insertSite(new Vertex(circumcircle[0], circumcircle[1]));
////			} catch (Exception e) {
////				// TODO: handle exception
////			}
//		}
//	}
//		qed.insertSite(new Vertex(x, y));
//		Geometry out = d.getTriangles(PTS.geometryFactory); // triangulates convex hull of points
//		d = new DelaunayTriangulationBuilder();
//		d.setSites(coords);
//		out = d.getTriangles(PTS.geometryFactory); // triangulates concave hull of points
//		d.getSubdivision().fr
//		var z = ((Polygon) g).getExteriorRing(); // segments
		out = out.intersection(g); // get concave hull
		return toPShape(out);
	}

	/**
	 * smallest angle is opposite the smallest side
	 * 
	 * calc sides then use law of cosines to calc angle
	 * 
	 * @return [smallestSideLength, smallestAngle(rads)]
	 */
	public static double[] smallestSideAndAngle(Coordinate a, Coordinate b, Coordinate c) {
		// (lengths squared)
		double ab = (b.getY() - a.getY()) * (b.getY() - a.getY()) + (b.getX() - a.getX()) * (b.getX() - a.getX());
		double bc = (c.getY() - b.getY()) * (c.getY() - b.getY()) + (c.getX() - b.getX()) * (c.getX() - b.getX());
		double ca = (a.getY() - c.getY()) * (a.getY() - c.getY()) + (a.getX() - c.getX()) * (a.getX() - c.getX());

		// real lengths
		double abS = Math.sqrt(ab);
		double bcS = Math.sqrt(bc);
		double caS = Math.sqrt(ca);
		double min = Math.min(Math.min(abS, bcS), caS);

		double theta;
		if (min == abS) {
			theta = Math.acos((bc + ca - ab) / (2 * bcS * caS));

		} else {
			if (min == bcS) {
				theta = Math.acos((ab + ca - bc) / (2 * abS * caS));
			} else {
				theta = Math.acos((ab + bc - ca) / (2 * abS * bcS));
			}
		}
		return new double[] { min, theta };

	}

	/**
	 * smallest angle is opposite the smallest side
	 * 
	 * calc sides then use law of cosines to calc angle
	 * 
	 * @return [smallestSideLength, smallestAngle(rads)]
	 */
	public static double[] smallestSideAndAngle(Vertex a, Vertex b, Vertex c) {
		// (lengths squared)
		double ab = (b.getY() - a.getY()) * (b.getY() - a.getY()) + (b.getX() - a.getX()) * (b.getX() - a.getX());
		double bc = (c.getY() - b.getY()) * (c.getY() - b.getY()) + (c.getX() - b.getX()) * (c.getX() - b.getX());
		double ca = (a.getY() - c.getY()) * (a.getY() - c.getY()) + (a.getX() - c.getX()) * (a.getX() - c.getX());

		// real lengths
		double abS = Math.sqrt(ab);
		double bcS = Math.sqrt(bc);
		double caS = Math.sqrt(ca);
		double min = Math.min(Math.min(abS, bcS), caS);

		double theta;
		if (min == abS) {
			theta = Math.acos((bc + ca - ab) / (2 * bcS * caS));

		} else {
			if (min == bcS) {
				theta = Math.acos((ab + ca - bc) / (2 * abS * caS));
			} else {
				theta = Math.acos((ab + bc - ca) / (2 * abS * bcS));
			}
		}
		return new double[] { min, theta };
	}

	/**
	 * https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
	 * 
	 * @return [centerX, centerY, radius]
	 */
	public static double[] circumcircle(Coordinate a, Coordinate b, Coordinate c) {

		double D = (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y);
		double px = (((a.x - c.x) * (a.x + c.x) + (a.y - c.y) * (a.y + c.y)) / 2 * (b.y - c.y)
				- ((b.x - c.x) * (b.x + c.x) + (b.y - c.y) * (b.y + c.y)) / 2 * (a.y - c.y)) / D;

		double py = (((b.x - c.x) * (b.x + c.x) + (b.y - c.y) * (b.y + c.y)) / 2 * (a.x - c.x)
				- ((a.x - c.x) * (a.x + c.x) + (a.y - c.y) * (a.y + c.y)) / 2 * (b.x - c.x)) / D;

		double rs = (c.x - px) * (c.x - px) + (c.y - py) * (c.y - py);

		return new double[] { px, py, Math.sqrt(rs) };

	}

	/**
	 * @return [centerX, centerY, radius]
	 */
	public static double[] circumcircle(Vertex a, Vertex b, Vertex c) {

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

		double rs = (c.getX() - px) * (c.getX() - px) + (c.getY() - py) * (c.getY() - py);

		return new double[] { px, py, Math.sqrt(rs) };

	}

	/**
	 * The diametral circle of a segment is the (unique) smallest circle that
	 * contains the segment.
	 * 
	 * <p>
	 * A segment is said to be encroached if a point lies within its diametral
	 * circle. Any encroached segment that arises is immediately split by inserting
	 * a vertex at its midpoint. The two resulting subsegments have smaller
	 * diametral circles, and may or may not be encroached themselves.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static double[] diametralCircle(Coordinate a, Coordinate b) {
		// https://www.cs.cmu.edu/~quake/tripaper/triangle3.html
		double d = Math.sqrt((b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x));
		return new double[] { (a.x + b.x) / 2, (a.y + b.y) / 2, d / 2 };
	}

}
