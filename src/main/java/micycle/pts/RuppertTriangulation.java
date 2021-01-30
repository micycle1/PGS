package micycle.pts;

import static micycle.pts.PTS.fromPShape;
import static micycle.pts.PTS.toPShape;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;

import processing.core.PShape;

public class RuppertTriangulation {

	public static PShape delaunayTriangulation(PShape shape) {
		Geometry g = fromPShape(shape);
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setSites(g);
		Geometry out = d.getTriangles(PTS.geometryFactory); // triangulates convex hull of points
//		d.getSubdivision().fr
		var z = ((Polygon) g).getExteriorRing();
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
	void ruppert() {

	}

	/**
	 * smallest angle is opposite the smallest side
	 * 
	 * calc sides then use law of cosines to calc angle
	 * 
	 * @return [smallestSideLength, smallestAngle(rads)]
	 */
	public static double[] smallestSideAndAngle(Coordinate a, Coordinate b, Coordinate c) {
		// don't sqrt lengths yet
		double ab = (b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x);
		double bc = (c.y - b.y) * (c.y - b.y) + (c.x - b.x) * (c.x - b.x);
		double ca = (a.y - c.y) * (a.y - c.y) + (a.x - c.x) * (a.x - c.x);
		double min = Math.min(Math.min(ab, bc), ca);

		double theta;
		if (min == ab) {
			theta = Math.acos((bc + ca - ab) / (2 * bc * ca));

		} else {
			if (min == bc) {
				theta = Math.acos((ab + ca - bc) / (2 * ab * ca));
			} else {
				theta = Math.acos((ab + bc - ca) / (2 * ab * bc));
			}
		}
		return new double[] { Math.sqrt(min), theta };

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
