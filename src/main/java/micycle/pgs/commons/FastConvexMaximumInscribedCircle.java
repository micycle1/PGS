package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;

/**
 * Computes the exact maximum inscribed circle inside a convex polygon by
 * solving a linear program.
 * 
 * @author Michael Carleton
 */
public final class FastConvexMaximumInscribedCircle {

	/*
	 * From 'An Efficient Algorithm to Calculate the Center of the Biggest Inscribed
	 * Circle in an Irregular Polygon' by OSCAR MARTINEZ.
	 */

	/**
	 * Solve for the maximum inscribed circle from a Geometry.
	 *
	 * This method requires a single Polygon. It throws IllegalArgumentException if
	 * the input is null, not polygonal, or not a single Polygon (e.g. a
	 * MultiPolygon is rejected). For an empty polygon the delegated solve(Polygon)
	 * may return null.
	 *
	 * @param geom input geometry (must be a single Polygon)
	 * @return Coordinate (x,y,r) with r in z, or null on failure / empty polygon
	 * @throws IllegalArgumentException if geom is null or not a single polygonal
	 *                                  Polygon
	 */
	public static Coordinate getCircle(Geometry geom) {
		if (geom == null) {
			throw new IllegalArgumentException("geometry must not be null");
		}
		if (!(geom instanceof Polygonal)) {
			throw new IllegalArgumentException("geometry must be polygonal");
		}
		if (!(geom instanceof Polygon)) {
			throw new IllegalArgumentException("expected a single Polygon (not a MultiPolygon)");
		}
		return getCircle((Polygon) geom);
	}

	/**
	 * Solve for the maximum inscribed circle inside the given polygon.
	 * <p>
	 * The method formulates a linear program that maximizes r subject to half-plane
	 * constraints derived from the polygon edges. For each edge the constraint is
	 * n.x * x + n.y * y - r >= n · p0 where n is the inward unit normal and p0 is
	 * an edge endpoint.
	 * <p>
	 * Important: - The polygon should be convex for the solution to be valid. - If
	 * a solution cannot be found or an error occurs, null is returned.
	 * <p>
	 * The returned Coordinate encodes (x,y,r) with r in the z component.
	 *
	 * @param polygon polygon to inscribe the maximum inscribed circle into (may
	 *                contain holes)
	 * @return Coordinate (x,y,r) with r in z, or null on failure or invalid input
	 */
	public static Coordinate getCircle(Polygon polygon) {
		if (polygon == null || polygon.isEmpty()) {
			return null;
		}

		// maximise r
		final LinearObjectiveFunction obj = new LinearObjectiveFunction(new double[] { 0, 0, 1 }, 0);

		final List<LinearConstraint> constraints = new ArrayList<>();

//		 Envelope bounds to help solver stability (keeps x,y near polygon)
//		final Envelope env = polygon.getEnvelopeInternal();
//		constraints.add(new LinearConstraint(new double[] { 1, 0, 0 }, Relationship.GEQ, env.getMinX()));
//		constraints.add(new LinearConstraint(new double[] { 1, 0, 0 }, Relationship.LEQ, env.getMaxX()));
//		constraints.add(new LinearConstraint(new double[] { 0, 1, 0 }, Relationship.GEQ, env.getMinY()));
//		constraints.add(new LinearConstraint(new double[] { 0, 1, 0 }, Relationship.LEQ, env.getMaxY()));
//		// r >= 0
//		constraints.add(new LinearConstraint(new double[] { 0, 0, 1 }, Relationship.GEQ, 0.0));

		// Add constraints for outer ring and any holes
		addRingConstraints(polygon.getExteriorRing(), polygon, constraints);
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			addRingConstraints(polygon.getInteriorRingN(i), polygon, constraints);
		}

		try {
			SimplexSolver solver = new SimplexSolver(1e-3);
			PointValuePair sol = solver.optimize(obj, new LinearConstraintSet(constraints), GoalType.MAXIMIZE, new NonNegativeConstraint(false),
					new MaxIter(1000));
			if (sol == null) {
				return null;
			}

			double[] p = sol.getPoint();
			double r = p[2];
			Coordinate c = new Coordinate(p[0], p[1], r);
			return c;
		} catch (Exception e) {
			return null;
		}
	}

	/**
	 * Convert a ring (edge list) into half-plane linear constraints and append them
	 * to the provided list.
	 *
	 * @param ring   edge ring (exterior or interior)
	 * @param parent parent polygon (used to detect exterior ring)
	 * @param out    list to append LinearConstraint objects to
	 */
	private static void addRingConstraints(LineString ring, Polygon parent, List<LinearConstraint> out) {
		final Coordinate[] coords = ring.getCoordinates();
		if (coords.length < 2) {
			return;
		}

		final boolean isShell = ring == parent.getExteriorRing();
		final boolean ccw = Orientation.isCCW(coords);
		// If isShell == ccw, left normal points inward; otherwise flip it
		final double inwardSign = (isShell == ccw) ? +1.0 : -1.0;

		for (int i = 0; i < coords.length - 1; i++) {
			final Coordinate p0 = coords[i];
			final Coordinate p1 = coords[i + 1];

			final double len = p0.distance(p1);
			if (len == 0) {
				continue;
			}

			final double dx = p1.x - p0.x;
			final double dy = p1.y - p0.y;
			// Left unit normal
			double nx = (-dy / len) * inwardSign;
			double ny = (dx / len) * inwardSign;

			// Half-plane: nx*x + ny*y - r >= nx*p0.x + ny*p0.y
			final double rhs = nx * p0.x + ny * p0.y;
			out.add(new LinearConstraint(new double[] { nx, ny, -1.0 }, Relationship.GEQ, rhs));
		}
	}

}