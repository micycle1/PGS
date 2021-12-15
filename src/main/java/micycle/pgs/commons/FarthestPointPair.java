package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collection;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;

import micycle.pgs.PGS_Conversion;
import processing.core.PVector;

/**
 * The {@code FarthestPair} data type computes the farthest pair of points in a
 * set of <em>n</em> points in the plane and provides accessor methods for
 * getting the farthest pair of points and the distance between them. The
 * distance between two points is their Euclidean distance.
 * <p>
 * This implementation computes the convex hull of the set of points and uses
 * the rotating calipers method to find all antipodal point pairs and the
 * farthest pair. It runs in O(<em>n</em> log <em>n</em>) time in the worst case
 * and uses O(<em>N</em>) extra space.
 *
 * @author Robert Sedgewick
 * @author Kevin Wayne
 * @author Adapeted by Michael Carleton
 */
public class FarthestPointPair {

	// https://algs4.cs.princeton.edu/99hull/FarthestPair.java.html

	// farthest pair of points and distance
	private PVector best1, best2;
	private double bestDistanceSquared = Double.NEGATIVE_INFINITY;

	/**
	 * Computes the farthest pair of points in the specified array of points.
	 *
	 * @param points an array of points
	 */
	public FarthestPointPair(Collection<PVector> points) {

		final Geometry convexHull = PGS_Conversion.fromPShape(PGS_Conversion.fromPVector(new ArrayList<>(points))).convexHull();
		Coordinate[] coords = convexHull.getCoordinates();
		if (!Orientation.isCCW(coords)) {
			coords = convexHull.reverse().getCoordinates();
		}

		// number of points on the hull
		int m = coords.length;

		// single point
		if (m <= 1) {
			return;
		}

		// the hull, in counterclockwise order hull[1] to hull[m]
		PVector[] hull = new PVector[m + 1];

		m = 1;
		for (int i = m; i < coords.length; i++) {
			hull[m++] = new PVector((float) coords[i].x, (float) coords[i].y);
		}
		m--;

		// all points are equal
		if (m == 1) {
			return;
		}

		// points are collinear
		if (m == 2) {
			best1 = hull[1];
			best2 = hull[2];
			bestDistanceSquared = best1.dist(best2);
			return;
		}

		// k = farthest vertex from edge from hull[1] to hull[m]
		int k = 2;
		while (area2(hull[m], hull[1], hull[k + 1]) > area2(hull[m], hull[1], hull[k])) {
			k++;
		}

		int j = k;
		for (int i = 1; i <= k && j <= m; i++) {
			if (hull[i].dist(hull[j]) > bestDistanceSquared) {
				best1 = hull[i];
				best2 = hull[j];
				bestDistanceSquared = hull[i].dist(hull[j]);
			}
			while ((j < m) && area2(hull[i], hull[i + 1], hull[j + 1]) > area2(hull[i], hull[i + 1], hull[j])) {
				j++;
				double distanceSquared = hull[i].dist(hull[j]);
				if (distanceSquared > bestDistanceSquared) {
					best1 = hull[i];
					best2 = hull[j];
					bestDistanceSquared = hull[i].dist(hull[j]);
				}
			}
		}
	}

	/**
	 * Returns one of the points in the farthest pair of points.
	 *
	 * @return one of the two points in the farthest pair of points; {@code null} if
	 *         no such point (because there are fewer than 2 points)
	 */
	public PVector either() {
		return best1;
	}

	/**
	 * Returns the other point in the farthest pair of points.
	 *
	 * @return the other point in the farthest pair of points {@code null} if no
	 *         such point (because there are fewer than 2 points)
	 */
	public PVector other() {
		return best2;
	}

	/**
	 * Returns twice the signed area of the triangle a-b-c.
	 * 
	 * @param a first point
	 * @param b second point
	 * @param c third point
	 * @return twice the signed area of the triangle a-b-c
	 */
	private static double area2(PVector a, PVector b, PVector c) {
		return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	}

}
