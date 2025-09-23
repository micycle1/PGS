package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.SegmentStringUtil;
import org.locationtech.jts.operation.distance.DistanceOp;

/**
 * Compute a uniform scale about the shape centroid so the shape first touches a
 * boundary geometry (shrink or grow to first contact).
 * <p>
 * Uses analytic vertex–edge events (shape-vertex → boundary-edge and
 * boundary-vertex → shape-edge) to compute the minimal positive scale factor.
 * </p>
 * 
 * @author Michael Carleton
 */
public class TouchScale {

	/**
	 * <p>
	 * Scales a shape about its centroid so it first contacts boundary.
	 * </p>
	 * 
	 * <ul>
	 * <li><b>Inputs:</b> shape and boundary (may be polygons, lines,
	 * multiparts).</li>
	 * <li><b>Behavior:</b> returns a new Geometry = scaled shape. The centroid of
	 * gShape is used as the scale center</li>
	 * </ul>
	 * 
	 * @param shape    the geometry to scale
	 * @param boundary the boundary geometry to touch
	 * @return scaled geometry (new Geometry instance)
	 */
	@SuppressWarnings("unchecked")
	public static Geometry scale(Geometry shape, Geometry boundary) {
		/*
		 * Uniform scaling about the centroid c sends every point x to c + s·(x − c).
		 * First contact happens at the smallest s > 0 where either: a shape-vertex
		 * (scaled along its ray from c) hits a boundary edge, or a boundary-vertex hits
		 * a shape edge (scaled). Solve these two cases in closed form for all
		 * vertex–edge pairs; take the minimal positive s.
		 */
		if (shape.isEmpty() || boundary.isEmpty())
			return shape;

		final Coordinate c = shape.getCentroid().getCoordinate();
		final double EPS = 1e-12;

		// Vertices (duplicates are fine)
		final Coordinate[] shapeVerts = shape.getCoordinates();
		final Coordinate[] boundVerts = boundary.getCoordinates();

		// Edges from boundaries (exterior + holes)
		final List<SegmentString> ssB = SegmentStringUtil.extractSegmentStrings(boundary.getBoundary());
		final List<LineSegment> edgesB = toLineSegments(ssB);

		final List<SegmentString> ssS = SegmentStringUtil.extractSegmentStrings(shape.getBoundary());
		final List<LineSegment> edgesS = toLineSegments(ssS);

		double sBest = Double.POSITIVE_INFINITY;

		// 1) shape-vertex -> boundary-edge
		for (Coordinate v : shapeVerts) {
			double dvx = v.x - c.x, dvy = v.y - c.y;
			if (Math.abs(dvx) + Math.abs(dvy) < EPS)
				continue;

			for (LineSegment seg : edgesB) {
				double ex = seg.p1.x - seg.p0.x, ey = seg.p1.y - seg.p0.y;
				double denom = cross(dvx, dvy, ex, ey);
				if (Math.abs(denom) < EPS)
					continue;

				double s = -cross(c.x - seg.p0.x, c.y - seg.p0.y, ex, ey) / denom;
				if (!(s > 0))
					continue;

				double yx = c.x + s * dvx, yy = c.y + s * dvy;
				double txNum = ((yx - seg.p0.x) * ex + (yy - seg.p0.y) * ey);
				double txDen = ex * ex + ey * ey;
				if (txDen <= EPS)
					continue;
				double t = txNum / txDen;
				if (t >= -1e-9 && t <= 1 + 1e-9) {
					if (s < sBest)
						sBest = s;
				}
			}
		}

		// 2) boundary-vertex -> shape-edge
		for (Coordinate w : boundVerts) {
			double wcx = w.x - c.x, wcy = w.y - c.y;

			for (LineSegment seg : edgesS) {
				double rax = seg.p0.x - c.x, ray = seg.p0.y - c.y;
				double rbx = seg.p1.x - c.x, rby = seg.p1.y - c.y;
				double mx = rbx - rax, my = rby - ray;

				double denom = cross(rax, ray, mx, my);
				if (Math.abs(denom) < EPS)
					continue;

				double s1 = cross(wcx, wcy, mx, my) / denom;
				if (!(s1 > 0))
					continue;

				double s2 = cross(rax, ray, wcx, wcy) / denom;
				double u = s2 / s1;
				if (u >= -1e-9 && u <= 1 + 1e-9) {
					if (s1 < sBest)
						sBest = s1;
				}
			}
		}

		// Fallback if degenerate
		if (!Double.isFinite(sBest)) {
			Coordinate nb = nearestPointOnBoundary(boundary, c);
			double ux = nb.x - c.x, uy = nb.y - c.y;
			double ulen = Math.hypot(ux, uy);
			if (ulen < EPS)
				return shape;
			ux /= ulen;
			uy /= ulen;

			double rB = firstRayHitDistance(c.x, c.y, ux, uy, edgesB);
			double rS = maxDotAlong(shapeVerts, c, ux, uy);
			if (rB > 0 && rS > EPS)
				sBest = rB / rS;
			else
				return shape;
		}

		// Tiny relative backoff
		double sApply = Math.max(0.0, sBest * (1 - 1e-9));

		AffineTransformation T = AffineTransformation.scaleInstance(sApply, sApply, c.x, c.y);
		Geometry out = T.transform(shape);
		return out;
	}

	/* -------- helpers -------- */

	private static List<LineSegment> toLineSegments(List<org.locationtech.jts.noding.SegmentString> ssList) {
		ArrayList<LineSegment> out = new ArrayList<>();
		for (SegmentString ss : ssList) {
			Coordinate[] cs = ss.getCoordinates();
			for (int i = 0; i < cs.length - 1; i++) {
				if (!cs[i].equals2D(cs[i + 1])) {
					out.add(new LineSegment(cs[i], cs[i + 1]));
				}
			}
		}
		return out;
	}

	private static double cross(double ax, double ay, double bx, double by) {
		return ax * by - ay * bx;
	}

	private static double firstRayHitDistance(double cx, double cy, double rx, double ry, List<LineSegment> edges) {
		final double EPS = 1e-12;
		double best = Double.POSITIVE_INFINITY;
		for (LineSegment seg : edges) {
			double ex = seg.p1.x - seg.p0.x, ey = seg.p1.y - seg.p0.y;
			double denom = cross(rx, ry, ex, ey);
			if (Math.abs(denom) < EPS)
				continue;
			double dx = seg.p0.x - cx, dy = seg.p0.y - cy;
			double t = cross(dx, dy, ex, ey) / denom; // along ray
			double u = cross(dx, dy, rx, ry) / denom; // along segment
			if (t > 0 && u >= -1e-9 && u <= 1 + 1e-9) {
				if (t < best)
					best = t;
			}
		}
		return best;
	}

	private static double maxDotAlong(Coordinate[] verts, Coordinate c, double ux, double uy) {
		double best = -Double.MAX_VALUE;
		for (Coordinate v : verts) {
			double dx = v.x - c.x, dy = v.y - c.y;
			double d = dx * ux + dy * uy;
			if (d > best)
				best = d;
		}
		return best;
	}

	private static Coordinate nearestPointOnBoundary(Geometry boundary, Coordinate c) {
		Coordinate[] pair = DistanceOp.nearestPoints(boundary, boundary.getFactory().createPoint(c));
		return pair[0];
	}

}
