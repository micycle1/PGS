package micycle.pgs.commons;

import org.locationtech.jts.geom.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Geometry smoothing using Lane-Riesenfeld (LR) curve subdivision with 4-point
 * refinement to reduce contraction.
 * <p>
 * The LR algorithm is a generalization of the Chaikin subdivision which
 * generates splines with variable continuity. The 4-point (Dyn-Levin-Gregory)
 * refinement is a variant that interpolates the control points. A combination
 * of LR and 4-point refinement can be used to reduce contraction when smoothing
 * a geometry.
 * <p>
 * This class provides a utility method to subdivide geometries using the LR
 * algorithm with 4-point refinement. The algorithm can be applied to both open
 * and closed geometries (e.g., LineString and LinearRing).
 * 
 * @author Michael Carleton
 */
public class LaneRiesenfeldSmoothing {

	// https://observablehq.com/@esperanc/lane-riesenfeld-subdivision
	// https://tiborstanko.sk/teaching/geo-num-2017/tp5.html

	/**
	 * Subdivides the input geometry using the Lane-Riesenfeld algorithm with
	 * 4-point refinement.
	 *
	 * @param geometry              The lineal input geometry to subdivide.
	 * @param degree                The degree of the LR algorithm. Higher degrees
	 *                              influence the placement of vertices and the
	 *                              overall shape of the curve, but only slightly
	 *                              increase the number of vertices generated.
	 *                              Increasing the degree also increases the
	 *                              contraction of the curve toward its control
	 *                              points. The degree does not directly control the
	 *                              smoothness of the curve. A value of 3 or 4 is
	 *                              usually sufficient for most applications.
	 * @param subdivisions          The number of times the subdivision process is
	 *                              applied. More subdivisions result in finer
	 *                              refinement and visually smoother curves between
	 *                              vertices. A value of 3 or 4 is usually
	 *                              sufficient for most applications.
	 * @param antiContractionFactor The weight parameter for the 4-point refinement.
	 *                              Controls the interpolation strength. A value of
	 *                              0 effectively disables the contraction
	 *                              reduction. Generally suitable values are in
	 *                              [0...0.1]. Larger values may create
	 *                              self-intersecting geometry.
	 * @return A new subdivided geometry (LineString or LinearRing).
	 */
	public static LineString subdivide(LineString geometry, int degree, int subdivisions, double antiContractionFactor) {
		Coordinate[] coords = geometry.getCoordinates();
		boolean closed = geometry.isClosed();
		if (closed && coords.length > 0) {
			coords = Arrays.copyOf(coords, coords.length - 1); // Remove the last coordinate if closed
		}
		Coordinate[] subdivided = lr4(coords, degree, closed, antiContractionFactor, subdivisions);
		GeometryFactory factory = geometry.getFactory();
		return createGeometry(factory, subdivided, closed);
	}

	private static LineString createGeometry(GeometryFactory factory, Coordinate[] coords, boolean closed) {
		if (closed && coords.length > 0) {
			List<Coordinate> coordList = new ArrayList<>(Arrays.asList(coords));
			coordList.add(new Coordinate(coordList.get(0)));
			return factory.createLinearRing(coordList.toArray(new Coordinate[0]));
		} else {
			return factory.createLineString(coords);
		}
	}

	/**
	 * Applies the Lane-Riesenfeld algorithm with 4-point refinement to an array of
	 * coordinates. This replaces the vanilla midpoint averaging with four-point
	 * averaging.
	 *
	 * @param points       The input array of coordinates.
	 * @param degree       The degree of the Lane-Riesenfeld algorithm.
	 * @param closed       Whether the geometry is closed.
	 * @param w            The weight parameter for the 4-point refinement.
	 * @param subdivisions The number of times the subdivision process is applied.
	 * @return A new array of subdivided coordinates.
	 */
	private static Coordinate[] lr4(Coordinate[] points, int degree, boolean closed, double w, int subdivisions) {
		if (degree < 1 || subdivisions < 1) {
			return Arrays.copyOf(points, points.length);
		}

		Coordinate[] v = points;

		for (int s = 0; s < subdivisions; s++) {
			v = fourPoint(v, closed, w);

			for (int d = 1; d < degree; d++) {
				int n = v.length;
				List<Coordinate> u = new ArrayList<>();

				for (int i = 0; i < n; i++) {
					int prevIndex = getPreviousIndex(i, n, closed);
					int nextIndex = getNextIndex(i, n, closed);
					int nextNextIndex = getNextNextIndex(i, n, closed);

					Coordinate p0 = v[prevIndex];
					Coordinate p1 = v[i];
					Coordinate p2 = v[nextIndex];
					Coordinate p3 = v[nextNextIndex];

					double qx = computeQ(p0.x, p1.x, p2.x, p3.x, w);
					double qy = computeQ(p0.y, p1.y, p2.y, p3.y, w);
					u.add(new Coordinate(qx, qy));
				}

				if (closed) {
					v = u.toArray(new Coordinate[0]);
				} else {
					List<Coordinate> newV = new ArrayList<>();
					if (v.length > 0) {
						newV.add(v[0]);
					}
					if (!u.isEmpty()) {
						newV.addAll(u.subList(0, u.size() - 1));
					}
					if (v.length > 0) {
						newV.add(v[v.length - 1]);
					}
					v = newV.toArray(new Coordinate[0]);
				}
			}
		}

		return v;
	}

	/**
	 * Applies the 4-point (Dyn-Levin-Gregory) refinement to an array of
	 * coordinates.
	 *
	 * @param points The input array of coordinates.
	 * @param closed Whether the geometry is closed.
	 * @param w      The weight parameter for the 4-point refinement.
	 * @return A new array of refined coordinates.
	 */
	private static Coordinate[] fourPoint(Coordinate[] points, boolean closed, double w) {
		int n = points.length;
		if (n == 0) {
			return new Coordinate[0];
		}

		List<Coordinate> result = new ArrayList<>(points.length);

		for (int i = 0; i < n; i++) {
			int prevIndex = getPreviousIndex(i, n, closed);
			int nextIndex = getNextIndex(i, n, closed);
			int nextNextIndex = getNextNextIndex(i, n, closed);

			Coordinate p0 = points[prevIndex];
			Coordinate p1 = points[i];
			Coordinate p2 = points[nextIndex];
			Coordinate p3 = points[nextNextIndex];

			double qx = computeQ(p0.x, p1.x, p2.x, p3.x, w);
			double qy = computeQ(p0.y, p1.y, p2.y, p3.y, w);
			Coordinate q = new Coordinate(qx, qy);

			result.add(p1);
			result.add(q);
		}

		if (!closed && !result.isEmpty()) {
			result.remove(result.size() - 1);
		}

		return result.toArray(new Coordinate[0]);
	}

	/**
	 * Computes a new coordinate value using the 4-point refinement formula.
	 *
	 * @param p0 The first coordinate value.
	 * @param p1 The second coordinate value.
	 * @param p2 The third coordinate value.
	 * @param p3 The fourth coordinate value.
	 * @param w  The weight parameter.
	 * @return The computed coordinate value.
	 */
	private static double computeQ(double p0, double p1, double p2, double p3, double w) {
		return -w * p0 + (0.5 + w) * p1 + (0.5 + w) * p2 - w * p3;
	}

	private static int getPreviousIndex(int i, int n, boolean closed) {
		if (closed) {
			return (i - 1 + n) % n;
		} else {
			return Math.max(0, i - 1);
		}
	}

	private static int getNextIndex(int i, int n, boolean closed) {
		if (closed) {
			return (i + 1) % n;
		} else {
			return Math.min(n - 1, i + 1);
		}
	}

	private static int getNextNextIndex(int i, int n, boolean closed) {
		if (closed) {
			return (i + 2) % n;
		} else {
			return Math.min(n - 1, i + 2);
		}
	}
}