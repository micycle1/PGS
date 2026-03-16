package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.betterbeziers.CubicBezier;

/**
 * Piecewise cubic Bezier curve fitting + resampling for JTS geometries.
 *
 * <p>
 * Implements Philip J. Schneider’s “An Algorithm for Automatically Fitting
 * Digitized Curves” (Graphics Gems, 1990): given a polyline (digitized
 * vertices) and an error tolerance, it fits one or more <b>cubic Bezier
 * segments</b> that approximate the input, then returns a JTS
 * {@link org.locationtech.jts.geom.LineString LineString} by sampling those
 * Beziers at a fixed arc-length spacing.
 */
public final class SchneiderBezierFitter {

	private static final int MAX_FIT_ITERS = 4;
	private static final double DET_EPS = 1e-12;
	private static final double U_EPS = 1e-12;
	private static final double POINT_EPS2 = 1e-24;
	/**
	 * Upper bound on the distance of each inner Bezier control point from its
	 * corresponding endpoint, expressed as a multiple of the endpoint chord length.
	 *
	 * <p>
	 * Used as a stability clamp in {@code generateBezier(...)} to prevent
	 * near-singular least-squares solves from producing extremely long control
	 * handles, which can create large overshoots, self-intersections, or wild
	 * looping segments even when pointwise fitting error appears small.
	 *
	 * <p>
	 * For a segment with endpoints {@code v1} and {@code v4}, the solved handle
	 * lengths {@code alphaL} and {@code alphaR} are limited to at most
	 * {@code MAX_HANDLE_SCALE * |v4 - v1|}. If either handle exceeds this bound,
	 * the solver falls back to the standard heuristic placing control points at
	 * one-third of the chord length along the endpoint tangents.
	 *
	 * <p>
	 * Smaller values are more conservative and reduce overshoot at the cost of
	 * potentially splitting the input into more Bezier segments. Larger values
	 * allow looser fits but increase the risk of pathological curves.
	 */
	private static final double MAX_HANDLE_SCALE = 2.0;

	private SchneiderBezierFitter() {
	}

	/**
	 * Fits a piecewise cubic Bezier approximation to a JTS {@link LineString} and
	 * returns a new {@link LineString} created by sampling the fitted curve at a
	 * fixed arc-length spacing.
	 *
	 * <p>
	 * Uses the input line's {@link GeometryFactory}. If the input is closed, the
	 * output will be closed as well.
	 *
	 * @param line                input polyline/ring (must have at least 2
	 *                            vertices)
	 * @param error               maximum allowed distance error used by Schneider's
	 *                            fitter
	 * @param interSampleDistance target distance between successive output vertices
	 *                            along the fitted curve (must be &gt; 0)
	 * @return fitted-and-sampled LineString (closed if input was closed)
	 */
	public static LineString fitAndSample(LineString line, double error, double interSampleDistance) {
		Objects.requireNonNull(line, "line");
		return fitAndSample(line, error, interSampleDistance, line.getFactory());
	}

	/**
	 * Fits a piecewise cubic Bezier approximation to a JTS {@link LineString} and
	 * returns a new {@link LineString} created by sampling the fitted curve at a
	 * fixed arc-length spacing.
	 *
	 * <p>
	 * If the input is closed, the output will be closed as well.
	 *
	 * @param line                input polyline/ring (must have at least 2
	 *                            vertices)
	 * @param error               maximum allowed distance error used by Schneider's
	 *                            fitter
	 * @param interSampleDistance target distance between successive output vertices
	 *                            along the fitted curve (must be &gt; 0)
	 * @param gf                  geometry factory used to build the output
	 *                            LineString
	 * @return fitted-and-sampled LineString (closed if input was closed)
	 */
	public static LineString fitAndSample(LineString line, double error, double interSampleDistance, GeometryFactory gf) {
		Objects.requireNonNull(line, "line");
		Objects.requireNonNull(gf, "gf");

		line = (LineString) line.norm();

		Coordinate[] coords = line.getCoordinates();
		if (coords.length < 2)
			throw new IllegalArgumentException("LineString must have at least 2 coordinates");

		// if closed, keep the duplicate last coordinate so the closing edge is fitted
		List<Coordinate> pts = new ArrayList<>(coords.length);
		for (Coordinate coord : coords)
			pts.add(coord);

		// Fit+sample directly; if input is closed, output will already end where it
		// starts
		return fitAndSample(pts, error, interSampleDistance, gf);
	}

	/**
	 * Fits a piecewise cubic Bezier approximation to an ordered list of coordinates
	 * and returns a new {@link LineString} created by sampling the fitted curve at
	 * a fixed arc-length spacing.
	 *
	 * <p>
	 * This overload does not infer “closed-ness” from the input list; if you need
	 * ring-closure preservation, prefer the
	 * {@link #fitAndSample(LineString, double, double)} overload or ensure closure
	 * yourself after sampling.
	 *
	 * @param points              input vertices (must contain at least 2 points)
	 * @param error               maximum allowed distance error used by Schneider's
	 *                            fitter
	 * @param interSampleDistance target distance between successive output vertices
	 *                            along the fitted curve (must be &gt; 0)
	 * @param gf                  geometry factory used to build the output
	 *                            LineString
	 * @return fitted-and-sampled LineString
	 */
	private static LineString fitAndSample(List<Coordinate> points, double error, double interSampleDistance, GeometryFactory gf) {
		Objects.requireNonNull(gf, "gf");
		if (points == null || points.size() < 2) {
			throw new IllegalArgumentException("The number of points must be greater than 1");
		}
		if (!(interSampleDistance > 0)) {
			throw new IllegalArgumentException("interSampleDistance must be > 0");
		}
		if (!(error > 0)) {
			throw new IllegalArgumentException("error must be > 0");
		}

		// Detect closed ring by duplicate last==first
		final boolean closed = points.size() >= 4 && points.get(0) != null && points.get(points.size() - 1) != null
				&& points.get(0).equals2D(points.get(points.size() - 1));

		// Remove consecutive duplicate points; they cause zero-length tangents and bad
		// u values.
		List<Vector2D> pts = sanitizePoints(points, closed);

		if (!closed && pts.size() < 2) {
			throw new IllegalArgumentException("Need at least 2 distinct points");
		}
		if (closed && pts.size() < 4) {
			// Degenerate closed input after sanitization; just return the original
			// geometry.
			return gf.createLineString(points.toArray(new Coordinate[0]));
		}

		final Vector2D leftTangent;
		final Vector2D rightTangent;

		if (closed) {
			int n = pts.size();
			Vector2D pPrev = pts.get(n - 2); // last unique vertex
			Vector2D pNext = pts.get(1);

			Vector2D t = safeNormalize(pNext.subtract(pPrev));
			leftTangent = t;
			rightTangent = t.negate();
		} else {
			leftTangent = computeLeftTangent(pts);
			rightTangent = computeRightTangent(pts);
		}

		MultiBezierCurve fitted = fitCurves(0, pts.size() - 1, error, new MultiBezierCurve(), leftTangent, rightTangent, pts);

		List<Coordinate> out = new ArrayList<>();
		boolean firstSeg = true;

		for (BezierSeg seg : fitted.segments) {
			CubicBezier cb = seg.toBetterBezier();
			double[][] samples = cb.sampleEquidistantPoints(interSampleDistance);

			for (int i = 0; i < samples.length; i++) {
				if (!firstSeg && i == 0) {
					continue;
				}
				out.add(new Coordinate(samples[i][0], samples[i][1]));
			}
			firstSeg = false;
		}

		if (out.isEmpty()) {
			Coordinate a = points.get(0);
			Coordinate b = points.get(points.size() - 1);
			return gf.createLineString(new Coordinate[] { new Coordinate(a), new Coordinate(b) });
		}

		return gf.createLineString(out.toArray(new Coordinate[0]));
	}

	private static MultiBezierCurve fitCurves(int first, int last, double error, MultiBezierCurve curve, Vector2D leftTangent, Vector2D rightTangent,
			List<Vector2D> points) {

		int nPts = last - first + 1;

		if (nPts == 2) {
			double distance = dist(points.get(first), points.get(last)) / 3.0;
			BezierSeg bez = new BezierSeg();
			bez.v1 = points.get(first);
			bez.v4 = points.get(last);
			bez.v2 = bez.v1.add(leftTangent.multiply(distance));
			bez.v3 = bez.v4.add(rightTangent.multiply(distance));
			curve.segments.add(bez);
			return curve;
		}

		// Important: don't try to fit a single cubic across coincident endpoints.
		// This is common at the top of closed rings and can generate pathological
		// loops.
		if (samePoint(points.get(first), points.get(last))) {
			int split = findSplitForCoincidentEnds(points, first, last);
			Vector2D centerTangent = computeCenterTangent(points, split);
			fitCurves(first, split, error, curve, leftTangent, centerTangent, points);
			fitCurves(split, last, error, curve, centerTangent.negate(), rightTangent, points);
			return curve;
		}

		double[] u = chordLengthParameterize(points, first, last);
		BezierSeg bezier = generateBezier(points, first, last, u, leftTangent, rightTangent);

		int[] splitIndex = new int[] { (first + last) >>> 1 };
		double maxError = computeMaxError(points, first, last, bezier, u, splitIndex);

		if (maxError < error) {
			curve.segments.add(bezier);
			return curve;
		}

		// Try reparameterization only when we're already somewhat close.
		double iterationError = error * 4.0;
		if (maxError < iterationError) {
			double prevError = maxError;

			for (int i = 0; i < MAX_FIT_ITERS; i++) {
				double[] uPrime = reparameterize(points, first, last, u, bezier);

				// Critical stability check: if parameters are not strictly increasing,
				// the fit can double back and form loops.
				if (!isStrictlyIncreasing(uPrime)) {
					break;
				}

				BezierSeg candidate = generateBezier(points, first, last, uPrime, leftTangent, rightTangent);

				int[] candidateSplit = new int[] { splitIndex[0] };
				double candidateError = computeMaxError(points, first, last, candidate, uPrime, candidateSplit);

				if (candidateError < error) {
					curve.segments.add(candidate);
					return curve;
				}

				// If the iteration is no longer improving, stop and split instead.
				if (!(candidateError + 1e-15 < prevError)) {
					break;
				}

				bezier = candidate;
				u = uPrime;
				prevError = candidateError;
				splitIndex[0] = candidateSplit[0];
			}
		}

		Vector2D centerTangent = computeCenterTangent(points, splitIndex[0]);
		fitCurves(first, splitIndex[0], error, curve, leftTangent, centerTangent, points);
		fitCurves(splitIndex[0], last, error, curve, centerTangent.negate(), rightTangent, points);
		return curve;
	}

	private static BezierSeg generateBezier(List<Vector2D> points, int first, int last, double[] u, Vector2D leftVector, Vector2D rightVector) {
		int size = last - first + 1;

		Vector2D[][] A = new Vector2D[size][2];
		for (int i = 0; i < size; i++) {
			A[i][0] = leftVector.multiply(B1(u[i]));
			A[i][1] = rightVector.multiply(B2(u[i]));
		}

		Vector2D firstPoint = points.get(first);
		Vector2D lastPoint = points.get(last);

		double[][] C = new double[][] { { 0, 0 }, { 0, 0 } };
		double[] X = new double[] { 0, 0 };

		for (int i = 0; i < size; i++) {
			C[0][0] += A[i][0].dot(A[i][0]);
			C[0][1] += A[i][0].dot(A[i][1]);
			C[1][0] = C[0][1];
			C[1][1] += A[i][1].dot(A[i][1]);

			double ui = u[i];
			Vector2D tmp = points.get(first + i).subtract(firstPoint.multiply(B0(ui) + B1(ui)).add(lastPoint.multiply(B2(ui) + B3(ui))));

			X[0] += A[i][0].dot(tmp);
			X[1] += A[i][1].dot(tmp);
		}

		double detC0C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
		double detC0X = C[0][0] * X[1] - C[1][0] * X[0];
		double detXC1 = X[0] * C[1][1] - X[1] * C[0][1];

		double segLength = dist(firstPoint, lastPoint);
		double epsilon = 1.0e-6 * Math.max(segLength, 1.0);
		double maxHandle = MAX_HANDLE_SCALE * segLength;

		BezierSeg bez = new BezierSeg();
		bez.v1 = firstPoint;
		bez.v4 = lastPoint;

		boolean fallback = (segLength == 0.0);

		// Treat near-singular systems as singular; exact == 0.0 is not enough.
		double detTol = DET_EPS * (Math.abs(C[0][0] * C[1][1]) + Math.abs(C[0][1] * C[1][0]) + 1.0);

		double alphaL = 0.0;
		double alphaR = 0.0;

		if (!fallback) {
			if (Math.abs(detC0C1) <= detTol) {
				fallback = true;
			} else {
				alphaL = detXC1 / detC0C1;
				alphaR = detC0X / detC0C1;

				// Reject bad handle lengths; giant handles are what create the "crazy loops".
				if (!Double.isFinite(alphaL) || !Double.isFinite(alphaR) || alphaL < epsilon || alphaR < epsilon || alphaL > maxHandle || alphaR > maxHandle) {
					fallback = true;
				}
			}
		}

		if (fallback) {
			double d = segLength / 3.0;
			bez.v2 = bez.v1.add(leftVector.multiply(d));
			bez.v3 = bez.v4.add(rightVector.multiply(d));
		} else {
			bez.v2 = bez.v1.add(leftVector.multiply(alphaL));
			bez.v3 = bez.v4.add(rightVector.multiply(alphaR));
		}

		return bez;
	}

	private static double computeMaxError(List<Vector2D> points, int first, int last, BezierSeg bezier, double[] u, int[] splitIndex) {
		splitIndex[0] = (first + last) >>> 1;
		double maxDist2 = -1.0;

		CubicBezier cb = bezier.toBetterBezier();

		for (int i = first + 1; i < last; i++) {
			double t = u[i - first];
			double[] pt = cb.getPointAtParameter(t);

			double dx = pt[0] - points.get(i).getX();
			double dy = pt[1] - points.get(i).getY();
			double dist2 = dx * dx + dy * dy;

			if (dist2 >= maxDist2) {
				maxDist2 = dist2;
				splitIndex[0] = i;
			}
		}
		return Math.sqrt(maxDist2);
	}

	private static double[] reparameterize(List<Vector2D> points, int first, int last, double[] u, BezierSeg bezier) {
		int n = last - first + 1;
		double[] out = new double[n];
		CubicBezier cb = bezier.toBetterBezier();

		// Keep endpoints fixed. Refining them is unnecessary and can introduce drift.
		out[0] = 0.0;
		out[n - 1] = 1.0;

		for (int i = first + 1; i < last; i++) {
			out[i - first] = newtonRaphsonRootFind(cb, points.get(i), u[i - first]);
		}
		return out;
	}

	/**
	 * Newton-Raphson root refinement. Since your
	 * CubicBezier#getGradientAtParameter(u) returns a single double (not a
	 * derivative vector), Q'(u) and Q''(u) are computed by finite differences of
	 * getPointAtParameter(u).
	 */
	private static double newtonRaphsonRootFind(CubicBezier curve, Vector2D p, double u) {
		final double h0 = 1e-4;

		double um = Math.max(0.0, u - h0);
		double up = Math.min(1.0, u + h0);
		if (up == um) {
			return u;
		}

		double h = (up - um) / 2.0;

		double[] qm = curve.getPointAtParameter(um);
		double[] q = curve.getPointAtParameter(u);
		double[] qp = curve.getPointAtParameter(up);

		double q1x = (qp[0] - qm[0]) / (2.0 * h);
		double q1y = (qp[1] - qm[1]) / (2.0 * h);

		double h2 = h * h;
		double q2x = (qp[0] - 2.0 * q[0] + qm[0]) / h2;
		double q2y = (qp[1] - 2.0 * q[1] + qm[1]) / h2;

		double dx = q[0] - p.getX();
		double dy = q[1] - p.getY();

		double numerator = dx * q1x + dy * q1y;
		double denominator = (q1x * q1x + q1y * q1y) + (dx * q2x + dy * q2y);

		if (denominator == 0.0) {
			return u;
		}

		double uPrime = u - (numerator / denominator);
		return Math.max(0.0, Math.min(1.0, uPrime));
	}

	private static double[] chordLengthParameterize(List<Vector2D> points, int first, int last) {
		int n = last - first + 1;
		double[] u = new double[n];

		for (int i = 1; i < n; i++) {
			u[i] = u[i - 1] + dist(points.get(first + i), points.get(first + i - 1));
		}

		double total = u[n - 1];
		if (total == 0.0) {
			for (int i = 1; i < n; i++) {
				u[i] = i / (double) (n - 1);
			}
			return u;
		}

		for (int i = 1; i < n; i++) {
			u[i] /= total;
		}
		return u;
	}

	private static Vector2D computeLeftTangent(List<Vector2D> points) {
		return safeNormalize(points.get(1).subtract(points.get(0)));
	}

	private static Vector2D computeRightTangent(List<Vector2D> points) {
		int n = points.size();
		return safeNormalize(points.get(n - 2).subtract(points.get(n - 1)));
	}

	private static Vector2D computeCenterTangent(List<Vector2D> points, int centerIndex) {
		Vector2D v1 = points.get(centerIndex - 1).subtract(points.get(centerIndex));
		Vector2D v2 = points.get(centerIndex).subtract(points.get(centerIndex + 1));
		return safeNormalize(v1.add(v2).multiply(0.5));
	}

	private static Vector2D safeNormalize(Vector2D v) {
		double len = v.length();
		return len == 0.0 ? new Vector2D(0, 0) : v.multiply(1.0 / len);
	}

	private static double dist(Vector2D a, Vector2D b) {
		return a.distance(b);
	}

	private static List<Vector2D> sanitizePoints(List<Coordinate> points, boolean closed) {
		List<Vector2D> pts = new ArrayList<>(points.size());

		int limit = closed ? points.size() - 1 : points.size();
		Coordinate prev = null;

		for (int i = 0; i < limit; i++) {
			Coordinate c = points.get(i);
			if (c == null) {
				throw new IllegalArgumentException("points contains null Coordinate");
			}
			if (prev == null || !prev.equals2D(c)) {
				pts.add(new Vector2D(c.x, c.y));
				prev = c;
			}
		}

		if (closed && !pts.isEmpty()) {
			Vector2D p0 = pts.get(0);
			pts.add(new Vector2D(p0.getX(), p0.getY()));
		}

		return pts;
	}

	private static boolean isStrictlyIncreasing(double[] u) {
		for (int i = 1; i < u.length; i++) {
			if (!(u[i] > u[i - 1] + U_EPS)) {
				return false;
			}
		}
		return true;
	}

	private static boolean samePoint(Vector2D a, Vector2D b) {
		double dx = a.getX() - b.getX();
		double dy = a.getY() - b.getY();
		return dx * dx + dy * dy <= POINT_EPS2;
	}

	private static int findSplitForCoincidentEnds(List<Vector2D> points, int first, int last) {
		Vector2D a = points.get(first);

		int split = (first + last) >>> 1;
		double maxD2 = -1.0;

		for (int i = first + 1; i < last; i++) {
			double dx = points.get(i).getX() - a.getX();
			double dy = points.get(i).getY() - a.getY();
			double d2 = dx * dx + dy * dy;
			if (d2 > maxD2) {
				maxD2 = d2;
				split = i;
			}
		}

		if (split <= first)
			split = first + 1;
		if (split >= last)
			split = last - 1;
		return split;
	}

	private static double B0(double t) {
		double tmp = 1 - t;
		return tmp * tmp * tmp;
	}

	private static double B1(double t) {
		double tmp = 1 - t;
		return 3 * t * tmp * tmp;
	}

	private static double B2(double t) {
		double tmp = 1 - t;
		return 3 * t * t * tmp;
	}

	private static double B3(double t) {
		return t * t * t;
	}

	private static final class MultiBezierCurve {
		final List<BezierSeg> segments = new ArrayList<>();
	}

	private static final class BezierSeg {
		Vector2D v1, v2, v3, v4;

		CubicBezier toBetterBezier() {
			return new CubicBezier(v1.getX(), v1.getY(), v2.getX(), v2.getY(), v3.getX(), v3.getY(), v4.getX(), v4.getY());
		}
	}
}