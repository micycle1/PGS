package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SplittableRandom;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;

import micycle.betterbeziers.CubicBezier;
import net.jafama.FastMath;

/**
 * 2D random shape generator using Bezier curves.
 * 
 * @author Michael Carleton
 * @author Jonathan Viquerat
 *
 */
public class BezierShapeGenerator {

	/*
	 * Ports https://github.com/jviquerat/shapes which is in turn inspired by
	 * https://stackoverflow.com/a/50751932/9808792.
	 */

	private double[][] controlPts;
	private final int nControlPts;
	private final double bezierSampleDistance;
	private double[] radius;
	private double[] spikiness;
	private double[][] curvePts;

	/**
	 * Constructs a BezierShapeGenerator object with the specified number of control
	 * points, bezier sample distance, radius and spikiness.
	 *
	 * @param nControlPts          The number of control points for the bezier
	 *                             curve.
	 * @param bezierSampleDistance The sample distance along the bezier curve.
	 * @param radius               The radius relative to the distance between
	 *                             adjacent points. The radius is used to position
	 *                             the control points of the bezier curve, should be
	 *                             between 0 and 1. Larger values result in sharper
	 *                             features on the curve.
	 * @param spikiness            A measure of the curve's smoothness. If 0, the
	 *                             curve's angle through each point will be the
	 *                             average between the direction to adjacent points.
	 *                             As it increases, the angle will be determined
	 *                             mostly by one adjacent point, making the curve
	 *                             more "spiky".
	 */
	public BezierShapeGenerator(int nControlPts, double bezierSampleDistance, double radius, double spikiness) {
		this.nControlPts = nControlPts;
		this.bezierSampleDistance = bezierSampleDistance;

		this.radius = new double[nControlPts];
		for (int i = 0; i < nControlPts; i++) {
			this.radius[i] = radius;
		}

		this.spikiness = new double[nControlPts];
		for (int i = 0; i < nControlPts; i++) {
			this.spikiness[i] = spikiness;
		}

	}

	/**
	 * Generates a shape's curve based on the control points. The shape's curve can
	 * be scaled, centered and sorted in counter-clockwise order.
	 *
	 * @param centering If true, the generated shape will be centered at the origin
	 *                  (0, 0).
	 * @param cylinder  If true, the control points will be initialized in a
	 *                  cylindrical formation. If false, the control points are
	 *                  generated randomly.
	 * @param scale     The factor by which the initial control points are
	 *                  magnified.
	 * @param seed      The seed used for random control point generation
	 *                  (applicable when cylinder is false).
	 * @return An array of Coordinates representing the generated shape's vertices.
	 */
	public Coordinate[] generate(boolean centering, boolean cylinder, double scale, long seed) {
		
		if (nControlPts < 3) {
			return new Coordinate[] {};
		}

		// Generate random control points if empty
		if (controlPts == null || controlPts.length == 0) {
			if (cylinder) {
				controlPts = generateCylinderPts(nControlPts);
			} else {
				controlPts = generateRandomPts(nControlPts, seed);
			}
		}

		// Magnify
		for (int i = 0; i < controlPts.length; i++) {
			controlPts[i][0] *= scale;
			controlPts[i][1] *= scale;
		}

		// Center set of points
		if (centering) {
			double centerX = 0;
			double centerY = 0;
			for (int i = 0; i < controlPts.length; i++) {
				centerX += controlPts[i][0];
				centerY += controlPts[i][1];
			}
			centerX /= controlPts.length;
			centerY /= controlPts.length;

			for (int i = 0; i < controlPts.length; i++) {
				controlPts[i][0] -= centerX;
				controlPts[i][1] -= centerY;
			}
		}

		// Sort points counter-clockwise to avoid shape self-intersections
		double[][] sortedPts = new double[controlPts.length][2];
		double[] sortedRadius = new double[radius.length];
		double[] sortedEdgy = new double[spikiness.length];
		ccwSort(controlPts, radius, spikiness, sortedPts, sortedRadius, sortedEdgy);
		controlPts = sortedPts;
		radius = sortedRadius;
		spikiness = sortedEdgy;

		// Other local variables
		double[][] delta = new double[nControlPts][2];
		double[][] radii = new double[nControlPts][2];
		double[][] deltaB = new double[nControlPts][2];

		// Compute curve information
		for (int i = 0; i < nControlPts; i++) {

			int prev = (i - 1 + nControlPts) % nControlPts;
			int next = (i + 1) % nControlPts;

			double[] ptPrev = controlPts[prev];
			double[] ptCurr = controlPts[i];
			double[] ptNext = controlPts[next];

			// Compute delta vector
			double[] diff = subtract(ptNext, ptPrev);
			double norm = norm(diff);
			delta[i][0] = diff[0] / norm;
			delta[i][1] = diff[1] / norm;

			// Compute spikiness vector
			double[] deltaBpt = scale(0.5, add(ptPrev, ptNext));
			deltaB[i] = subtract(deltaBpt, ptCurr);

			// Compute radii
			double dist = distance(ptPrev, ptCurr);
			radii[i][0] = 0.5 * dist * radius[i];

			dist = distance(ptCurr, ptNext);
			radii[i][1] = 0.5 * dist * radius[i];

		}

		// Generate curves
		List<double[][]> localCurves = new ArrayList<>();
		for (int i = 0; i < nControlPts; i++) {
			final int next = (i + 1) % nControlPts;

			final double[] ptCurr = controlPts[i];
			final double[] ptNext = controlPts[next];

			double[][] localCurve = generateBezierCurve(ptCurr, ptNext, delta[i], delta[next], deltaB[i], deltaB[next], radii[i][1],
					radii[next][0], spikiness[i], spikiness[next]);

			localCurves.add(localCurve);

		}
		curvePts = localCurves.stream().flatMap(Arrays::stream).toArray(double[][]::new);

		// Recenter if needed
		if (centering) {
			double[] center = mean(curvePts);
			for (int i = 0; i < curvePts.length; i++) {
				curvePts[i][0] -= center[0];
				curvePts[i][1] -= center[1];
			}
			for (int i = 0; i < controlPts.length; i++) {
				controlPts[i][0] -= center[0];
				controlPts[i][1] -= center[1];
			}
		}

		final Coordinate[] coords = Arrays.stream(curvePts).map(p -> new Coordinate(p[0], p[1])).toArray(Coordinate[]::new);

		return new CoordinateList(coords, false).toCoordinateArray();

	}

	// Compute distance between two points
	private static double distance(double[] p1, double[] p2) {
		double dx = p1[0] - p2[0];
		double dy = p1[1] - p2[1];
		return Math.sqrt(dx * dx + dy * dy);
	}

	// Generate random points
	private static double[][] generateRandomPts(int nPts, long seed) {
		SplittableRandom r = new SplittableRandom(seed);
		double[][] pts = new double[nPts][2];
		for (int i = 0; i < nPts; i++) {
			pts[i][0] = r.nextDouble();
			pts[i][1] = r.nextDouble();
		}
		return pts;
	}

	// Generate cylinder points
	private static double[][] generateCylinderPts(int n_pts) {
		if (n_pts < 4) {
			System.out.println("Not enough points to generate cylinder");
		}

		double[][] pts = new double[n_pts][2];
		double ang = 2.0 * Math.PI / n_pts;
		for (int i = 0; i < n_pts; i++) {
			pts[i][0] = FastMath.cos(i * ang);
			pts[i][1] = FastMath.sin(i * ang);
		}

		return pts;
	}

	// Compute mean of array
	private static double[] mean(double[][] x) {
		double[] mean = new double[x[0].length];
		for (double[] xi : x) {
			for (int j = 0; j < xi.length; j++) {
				mean[j] += xi[j];
			}
		}
		for (int j = 0; j < mean.length; j++) {
			mean[j] /= x.length;
		}
		return mean;
	}

	// Compute ccwise ordering
	private static void ccwSort(double[][] pts, double[] radius, double[] edgy, double[][] sortedPts, double[] sortedRadius,
			double[] sortedEdgy) {

		// Compute geometric center
		double[] center = mean(pts);

		// Translate points
		double[][] translatedPts = new double[pts.length][2];
		for (int i = 0; i < pts.length; i++) {
			translatedPts[i][0] = pts[i][0] - center[0];
			translatedPts[i][1] = pts[i][1] - center[1];
		}

		// Compute angles
		double[] angles = new double[pts.length];
		for (int i = 0; i < pts.length; i++) {
			angles[i] = FastMath.atan2(translatedPts[i][1], translatedPts[i][0]);
		}

		// Sort by ascending angles
		Integer[] index = new Integer[pts.length];
		for (int i = 0; i < pts.length; i++) {
			index[i] = i;
		}
		Arrays.sort(index, (i1, i2) -> Double.compare(angles[i1], angles[i2]));

		// Copy to output arrays
		for (int i = 0; i < pts.length; i++) {
			sortedPts[i] = pts[index[i]];
			sortedRadius[i] = radius[index[i]];
			sortedEdgy[i] = edgy[index[i]];
		}

	}

	// Generate bezier curve
	private double[][] generateBezierCurve(double[] p1, double[] p2, double[] delta1, double[] delta2, double[] deltaB1, double[] deltaB2,
			double radius1, double radius2, double edgy1, double edgy2) {

		// Control points
		double[][] controlPts = new double[4][2];
		controlPts[0] = p1;
		controlPts[3] = p2;

		// Baseline intermediate control points
		double[] ctrlPt1Base = scale(radius1, delta1);
		double[] ctrlPt2Base = scale(-radius2, delta2);

		// Edge intermediate control points
		double[] ctrlPt1Edgy = deltaB1;
		double[] ctrlPt2Edgy = deltaB2;

		// Compute intermediate control points
		controlPts[1] = add(p1, add(ctrlPt1Base, scale(edgy1, ctrlPt1Edgy)));
		controlPts[2] = add(p2, add(ctrlPt2Base, scale(edgy2, ctrlPt2Edgy)));

		CubicBezier curve = new CubicBezier(controlPts);
		return curve.sampleEquidistantPoints(bezierSampleDistance);

	}

	// Utility functions
	private static double[] add(double[] x, double[] y) {
		double[] out = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			out[i] = x[i] + y[i];
		}
		return out;
	}

	private static double[] subtract(double[] x, double[] y) {
		double[] out = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			out[i] = x[i] - y[i];
		}
		return out;
	}

	private static double[] scale(double a, double[] x) {
		double[] out = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			out[i] = a * x[i];
		}
		return out;
	}

	private static double norm(double[] x) {
		double sumSq = 0;
		for (double v : x) {
			sumSq += v * v;
		}
		return Math.sqrt(sumSq);
	}

}