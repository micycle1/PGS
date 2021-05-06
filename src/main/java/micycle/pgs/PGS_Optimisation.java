package micycle.pgs;

import static micycle.pgs.PGS.SHAPE_SAMPLES;
import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pgs.color.RGB;
import micycle.pgs.utility.MaximumInscribedRectangle;
import micycle.pgs.utility.MinimumBoundingEllipse;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Solve geometric optimisation problems, such as bounding volumes, inscribed
 * areas, optimal distances, etc.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Optimisation {

	private PGS_Optimisation() {
	}

	/**
	 * Computes the shape's envelope (bounding box).
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape envelope(PShape shape) {
		return toPShape(fromPShape(shape).getEnvelope());
	}

	/**
	 * The Maximum Inscribed Circle is determined by a point in the interior of the
	 * area which has the farthest distance from the area boundary, along with a
	 * boundary point at that distance.
	 * 
	 * @param shape
	 * @param tolerance the distance tolerance for computing the centre point
	 *                  (around 1)
	 */
	public static PShape maximumInscribedCircle(PShape shape, double tolerance) {
		MaximumInscribedCircle mic = new MaximumInscribedCircle(fromPShape(shape), tolerance);

		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(SHAPE_SAMPLES);
		shapeFactory.setCentre(new Coordinate(mic.getCenter().getX(), mic.getCenter().getY()));
		shapeFactory.setWidth(mic.getRadiusLine().getLength() * 2); // r*2 for total width & height
		shapeFactory.setHeight(mic.getRadiusLine().getLength() * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());

	}

	/**
	 * Return the maximum circle (at a given centerpoint inside/outside the circle)
	 * 
	 * @param shape
	 * @param centerPoint
	 * @return A circular PShape
	 */
	public static PShape maximumInscribedCircle(PShape shape, PVector centerPoint) {
		Geometry g = fromPShape(shape);
		Point p = PGS.pointFromPVector(centerPoint);
		Coordinate closestEdgePoint = DistanceOp.nearestPoints(g.getBoundary(), p)[0];

		double radius = PGS.distance(GEOM_FACTORY.createPoint(closestEdgePoint), p);
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(SHAPE_SAMPLES);
		shapeFactory.setCentre(p.getCoordinate());
		shapeFactory.setWidth(radius * 2); // r*2 for total width & height
		shapeFactory.setHeight(radius * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * Compute the minimum inscribed rectangle for a shape. The method computes the
	 * MIR for convex shapes only; if a concave shape is passed in, the resulting
	 * rectangle will be computed based on its convex hull.
	 * 
	 * <p>
	 * This method uses a brute force algorithm to perform an exhaustive search for
	 * a solution (therefore it is slow relative to other optimisation methods).
	 * 
	 * @param shape
	 * @param fast  whether to compute MIR based on a lower resolution input. When
	 *              true processing is ~6 times faster but potentially a little
	 *              inaccurate
	 */
	public static PShape maximumInscribedRectangle(PShape shape, boolean fast) {
		double f = fast ? 5 : 2;
		final MaximumInscribedRectangle mir = new MaximumInscribedRectangle(fromPShape(shape), f);
		int[] r = mir.getInscribedRectangle();

		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setCentre(new Coordinate((r[0] + r[2] / 2d) * f, (r[1] + r[3] / 2d) * f));
		shapeFactory.setWidth(r[2] * f);
		shapeFactory.setHeight(r[3] * f);
		return toPShape(shapeFactory.createRectangle());
	}

	/**
	 * Computes the Minimum Bounding Circle (MBC) for the points in a Geometry. The
	 * MBC is the smallest circle which covers all the vertices of the input shape
	 * (this is also known as the Smallest Enclosing Circle). This is equivalent to
	 * computing the Maximum Diameter of the input vertex set.
	 */
	public static PShape minimumBoundingCircle(PShape shape) {
		MinimumBoundingCircle mbc = new MinimumBoundingCircle(fromPShape(shape));
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(SHAPE_SAMPLES);
		shapeFactory.setCentre(new Coordinate(mbc.getCentre().getX(), mbc.getCentre().getY()));
		shapeFactory.setWidth(mbc.getRadius() * 2); // r*2 for total width & height
		shapeFactory.setHeight(mbc.getRadius() * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * Computes the minimum bounding rectangle that encloses a shape. Unlike the
	 * envelope for a shape, the rectangle returned by this method can have any
	 * orientation (it's not axis-aligned).
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumBoundingRectangle(PShape shape) {
		Polygon md = (Polygon) MinimumDiameter.getMinimumRectangle(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Computes the minimum bounding ellipse that encloses a shape.
	 * 
	 * @param shape
	 * @param errorTolerance Mean-squared error tolerance (this value does not
	 *                       correspond to a pixel distance). 0.001 to 0.01
	 *                       recommended. Higher values are a looser (yet quicker)
	 *                       fit.
	 * @return
	 */
	public static PShape minimumBoundingEllipse(PShape shape, double errorTolerance) {
		final Geometry hull = fromPShape(shape).convexHull();
		final Coordinate[] coords = hull.getCoordinates();

		double[][] points = new double[coords.length][2];
		for (int i = 0; i < points.length; i++) {
			points[i][0] = coords[i].x;
			points[i][1] = coords[i].y;
		}

		final MinimumBoundingEllipse e = new MinimumBoundingEllipse(points, Math.max(errorTolerance, 0.001));
		double[][] eEoords = e.getBoundingCoordinates(100);

		final PShape ellipse = new PShape(PShape.PATH);
		ellipse.setFill(true);
		ellipse.setFill(RGB.WHITE);
		ellipse.beginShape();
		for (int i = 0; i < eEoords.length; i++) {
			ellipse.vertex((float) eEoords[i][0], (float) eEoords[i][1]);
		}
		ellipse.endShape();

		return ellipse;
	}

	/**
	 * Computes the minimum diameter of a shape.
	 * <p>
	 * The minimum diameter is defined to be the width of the smallest band that
	 * contains the shape, where a band is a strip of the plane defined by two
	 * parallel lines. This can be thought of as the smallest hole that the geometry
	 * can bemoved through, with a single rotation.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumDiameter(PShape shape) {
		LineString md = (LineString) MinimumDiameter.getMinimumDiameter(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Returns the nearest point of the shape to the given point. If the shape is
	 * has multiple children/geometries, the single closest point is returned.
	 * 
	 * @param shape
	 * @param point
	 * @return
	 * @see #closestPoints(PShape, PVector)
	 */
	public static PVector closestPoint(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		Coordinate coord = DistanceOp.nearestPoints(g, PGS.pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Returns the nearest point for each "island" in the input shape.
	 * 
	 * @param shape
	 * @param point
	 * @return list of closest points for each child shape. Output is identical to
	 *         {@link #closestPoint(PShape, PVector)} if the input shape is a single
	 *         polygon
	 * @see #closestPoint(PShape, PVector)
	 */
	public static List<PVector> closestPoints(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		ArrayList<PVector> points = new ArrayList<>();
		for (int i = 0; i < g.getNumGeometries(); i++) {
			Coordinate coord = DistanceOp.nearestPoints(g.getGeometryN(i), PGS.pointFromPVector(point))[0];
			points.add(new PVector((float) coord.x, (float) coord.y));
		}
		return points;
	}

	/**
	 * Solves the Problem of Apollonius (finding a circle tangent to three other
	 * circles in the plane). Circles are represented by PVectors, where the z
	 * coordinate is interpreted as radius.
	 * 
	 * @param c1 One of the circles in the problem
	 * @param c2 One of the circles in the problem
	 * @param c3 One of the circles in the problem
	 * @param s1 An indication if the solution should be externally or internally
	 *           tangent (+1/-1) to c1
	 * @param s2 An indication if the solution should be externally or internally
	 *           tangent (+1/-1) to c2
	 * @param s3 An indication if the solution should be externally or internally
	 *           tangent (+1/-1) to c3
	 * @return The circle (as a PVector) that is tangent to c1, c2 and c3.
	 */
	public static PVector solveApollonius(PVector c1, PVector c2, PVector c3, int s1, int s2, int s3) {

		// https://github.com/DIKU-Steiner/ProGAL/blob/master/src/ProGAL/geom2d/ApolloniusSolver.java

		double x1 = c1.x;
		double y1 = c1.y;
		double r1 = c1.z;
		double x2 = c2.x;
		double y2 = c2.y;
		double r2 = c2.z;
		double x3 = c3.x;
		double y3 = c3.y;
		double r3 = c3.z;

		// Currently optimized for fewest multiplications. Should be optimized for
		// readability
		double v11 = 2 * x2 - 2 * x1;
		double v12 = 2 * y2 - 2 * y1;
		double v13 = x1 * x1 - x2 * x2 + y1 * y1 - y2 * y2 - r1 * r1 + r2 * r2;
		double v14 = 2 * s2 * r2 - 2 * s1 * r1;

		double v21 = 2 * x3 - 2 * x2;
		double v22 = 2 * y3 - 2 * y2;
		double v23 = x2 * x2 - x3 * x3 + y2 * y2 - y3 * y3 - r2 * r2 + r3 * r3;
		double v24 = 2 * s3 * r3 - 2 * s2 * r2;

		double w12 = v12 / v11;
		double w13 = v13 / v11;
		double w14 = v14 / v11;

		double w22 = v22 / v21 - w12;
		double w23 = v23 / v21 - w13;
		double w24 = v24 / v21 - w14;

		double P = -w23 / w22;
		double Q = w24 / w22;
		double M = -w12 * P - w13;
		double N = w14 - w12 * Q;

		double a = N * N + Q * Q - 1;
		double b = 2 * M * N - 2 * N * x1 + 2 * P * Q - 2 * Q * y1 + 2 * s1 * r1;
		double c = x1 * x1 + M * M - 2 * M * x1 + P * P + y1 * y1 - 2 * P * y1 - r1 * r1;

		// Find a root of a quadratic equation. This requires the circle centers not
		// to be e.g. colinear
		double D = b * b - 4 * a * c;
		double rs = (-b - Math.sqrt(D)) / (2 * a);
		double xs = M + N * rs;
		double ys = P + Q * rs;
		return new PVector((float) xs, (float) ys, (float) rs);
	}

}
