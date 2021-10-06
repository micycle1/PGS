package micycle.pgs;

import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pgs.color.RGB;
import micycle.pgs.utility.RandomPolygon;
import micycle.pgs.utility.Star;
import processing.core.PConstants;
import processing.core.PShape;

/**
 * Construct uncommon/interesting 2D geometries.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Construction {

	private PGS_Construction() {
	}

	/**
	 * Generates a random simple convex polygon (n-gon).
	 * 
	 * @param n         number of polygon vertices/sides
	 * @param maxWidth  maximum width of generated random polygon
	 * @param maxHeight maximum height of generated random polygon
	 * @return
	 * @see {@link #createRandomPolygonExact(int, double, double)} to specify exact
	 *      dimensions
	 */
	public static PShape createRandomPolygon(int n, double maxWidth, double maxHeight) {
		return PGS_Transformation.translateTo(PGS_Conversion.fromPVector(RandomPolygon.generateRandomConvexPolygon(n, maxWidth, maxHeight)),
				maxWidth / 2, maxHeight / 2);
	}

	/**
	 * Generates a random simple convex polygon (n-gon), where the output's bounding
	 * box has the dimensions of those specified.
	 * 
	 * @param n      number of polygon vertices/sides
	 * @param width  width of generated random polygon
	 * @param height height of generated random polygon
	 * @return
	 */
	public static PShape createRandomPolygonExact(int n, double width, double height) {
		return PGS_Transformation.resize(PGS_Conversion.fromPVector(RandomPolygon.generateRandomConvexPolygon(n, width, height)), width,
				height);
	}

	/**
	 * Creates a supercircle PShape.
	 * 
	 * @param x      centre point X
	 * @param y      centre point Y
	 * @param width
	 * @param height
	 * @param power  circularity of super circle. Values less than 1 create
	 *               star-like shapes; power=1 is a square;
	 * @return
	 */
	public static PShape createSupercircle(double x, double y, double width, double height, double power) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(PGS.SHAPE_SAMPLES);
		shapeFactory.setCentre(new Coordinate(x, y));
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createSupercircle(power));
	}

	/**
	 * Creates a supershape PShape. The parameters feed into the superformula, which
	 * is a simple 2D analytical expression allowing to draw a wide variety of
	 * geometric and natural shapes (starfish, petals, snowflakes) by choosing
	 * suitable values relevant to few parameters.
	 * 
	 * @param x      centre point X
	 * @param y      centre point Y
	 * @param radius maximum radius
	 * @param m      specifies the rotational symmetry of the shape (3 = 3 sided; 4
	 *               = 4 sided)
	 * @param n1     supershape parameter 1
	 * @param n2     supershape parameter 2
	 * @param n3     supershape parameter 3
	 * @return
	 */
	public static PShape createSuperShape(double centerX, double centerY, double radius, double m, double n1, double n2, double n3) {
		// http://paulbourke.net/geometry/supershape/
		PShape shape = new PShape(PShape.PATH);
		shape.setFill(true);
		shape.setFill(RGB.WHITE);
		shape.beginShape();

		final int points = 180;
		final double angleInc = Math.PI * 2 / points;
		double angle = 0;
		while (angle < Math.PI * 2) {
			double r;
			double t1, t2;

			t1 = Math.cos(m * angle / 4);
			t1 = Math.abs(t1);
			t1 = Math.pow(t1, n2);

			t2 = Math.sin(m * angle / 4);
			t2 = Math.abs(t2);
			t2 = Math.pow(t2, n3);

			r = Math.pow(t1 + t2, 1 / n1);
			if (Math.abs(r) != 0) {
				r *= radius; // multiply r (0...1) by (max) radius
//				r = radius/r;
				shape.vertex((float) (centerX + r * Math.cos(angle)), (float) (centerY + r * Math.sin(angle)));
			}

			angle += angleInc;
		}

		shape.endShape();
		return shape;

	}

	/**
	 * Creates an elliptical arc polygon. The polygon is formed from the specified
	 * arc of an ellipse and the two radii connecting the endpoints to the centre of
	 * the ellipse.
	 * 
	 * @param x           centre point X
	 * @param y           centre point Y
	 * @param width
	 * @param height
	 * @param orientation start angle/orientation in radians (where 0 is 12 o'clock)
	 * @param angle       size of the arc angle in radians
	 * @return
	 */
	public static PShape createArc(double x, double y, double width, double height, double orientation, double angle) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(PGS.SHAPE_SAMPLES);
		shapeFactory.setCentre(new Coordinate(x, y));
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createArcPolygon(-Math.PI / 2 + orientation, angle));
	}

	/**
	 * Creates a star shape.
	 * 
	 * @param x           The x coordinate of the center
	 * @param y           The y coordinate of the center
	 * @param numRays     The number of rays that the star should have
	 * @param innerRadius The inner radius of the star
	 * @param outerRadius The outer radius of the star
	 * @param roundness   A roundness value between 0.0 and 1.0, for the inner and
	 *                    outer corners of the star.
	 * @return The star shape
	 */
	public static PShape createStar(double x, double y, int numRays, double innerRadius, double outerRadius, double roundness) {
		roundness = Math.max(Math.min(1, roundness), 0);
		final PShape shape = Star.createStarShape(x, y, innerRadius, outerRadius, numRays, roundness);
		shape.setFill(true);
		shape.setFill(255);
		return shape;
	}

	/**
	 * Creates a heart shape.
	 * 
	 * @param x     The x coordinate of the center of the heart
	 * @param y     The y coordinate of the center of the heart
	 * @param width Maximum width of the widest part of the heart
	 * @return
	 * @since 1.1.0
	 */
	public static PShape createHeart(final double x, final double y, final double width) {
		// https://mathworld.wolfram.com/HeartCurve.html
		PShape heart = new PShape(PShape.PATH);
		heart.setFill(true);
		heart.setFill(RGB.WHITE);
		heart.beginShape();

		final int points = 180;
		final double angleInc = Math.PI * 2 / points;
		double angle = 0;
		while (angle < Math.PI * 2) {
			final double s = Math.sin(angle);
			double vx = s * s * s;
			double vy = 13 * Math.cos(angle) - 5 * Math.cos(2 * angle) - 2 * Math.cos(3 * angle) - Math.cos(4 * angle);
			vy /= 17; // normalise to 1
			heart.vertex((float) (x + vx * width / 2), (float) (y - vy * width / 2));
			angle += angleInc;
		}

		heart.endShape(PConstants.CLOSE);
		return heart;
	}

	/**
	 * Creates a joined ring (a "donut") shape.
	 * 
	 * @param x           the x coordinate of the center
	 * @param y           the y coordinate of the center
	 * @param outerRadius radius of ring exterior
	 * @param innerRadius radius of ring hole
	 * @return the ring shape
	 * @since 1.1.3
	 */
	public static PShape createRing(double x, double y, double outerRadius, double innerRadius) {
		return createRing(x, y, outerRadius, innerRadius, 0, PConstants.TWO_PI);
	}

	/**
	 * Creates an (un)joined ring shape.
	 * 
	 * @param x           the x coordinate of the center
	 * @param y           the y coordinate of the center
	 * @param outerRadius radius of ring exterior
	 * @param innerRadius radius of ring hole
	 * @param orientation start angle/orientation in radians (where 0 is 12 o'clock)
	 * @param angle       size of the ring arc angle in radians
	 * @return the ring shape
	 * @since 1.1.3
	 */
	public static PShape createRing(double x, double y, double outerRadius, double innerRadius, double orientation, double angle) {
		final double outerR = Math.max(outerRadius, innerRadius);
		final double innerR = Math.min(outerRadius, innerRadius);

		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(48);
		shapeFactory.setCentre(new Coordinate(x, y));
		shapeFactory.setWidth(outerR * 2);
		shapeFactory.setHeight(outerR * 2);

		final Geometry outer;
		if (angle > PConstants.TWO_PI - 0.001) {
			outer = shapeFactory.createCircle();
		} else {
			outer = shapeFactory.createArcPolygon(-Math.PI / 2 + orientation, angle);
		}

		shapeFactory.setWidth(innerR * 2);
		shapeFactory.setHeight(innerR * 2);
		final Geometry inner = shapeFactory.createCircle();

		return toPShape(outer.difference(inner));
	}

	/**
	 * Creates a closed Sierpiński curve (a recursive space-filling curve).
	 * 
	 * @param centerX    the x coordinate of the curve center point
	 * @param centerY    the y coordinate of the curve center point
	 * @param width      length (the maximum width and height) of the curve (the
	 *                   curve will approach this width as the curve order is
	 *                   increased and more space is filled)
	 * @param curveOrder the order of the curve (the number of recursive
	 *                   subdivisions). Must be 1 or greater.
	 * @return a Sierpiński curve of the specified order
	 * @since 1.2.0
	 */
	public static PShape createSierpinskiCurve(double centerX, double centerY, double width, int curveOrder) {
		// https://piratefsh.github.io/2020/08/08/sierpinski-curve.html
		// when centerX/Y = 0, the actual center of the curve is width/2 (in both
		// dimensions).
		curveOrder = Math.max(1, curveOrder);
		centerX -= width / 2;
		centerY -= width / 2;
		// create the two triangles from the first subdivision
		final double[][] tri1 = new double[][] { { 0 + centerX, width + centerY }, { centerX, centerY }, { width + centerX, centerY } };
		final double[][] tri2 = new double[][] { { width + centerX, centerY }, { width + centerX, width + centerY },
				{ centerX, width + centerY } };

		// get points for each half of square recursively
		final List<double[]> half1 = subdivide(tri1, curveOrder);
		final List<double[]> half2 = subdivide(tri2, curveOrder);

		half1.addAll(half2); // combine points

		final PShape curve = new PShape(PShape.PATH);
		curve.setFill(true);
		curve.setFill(RGB.WHITE);
		curve.beginShape();
		half1.forEach(p -> curve.vertex((float) p[0], (float) p[1]));
		curve.endShape(PConstants.CLOSE);

		return curve;
	}

	/**
	 * Sierpinski curve subdivide.
	 */
	private static List<double[]> subdivide(double[][] triangle, int iters) {
		final List<double[]> points = new ArrayList<>();
		final double[] centroid = centroid(triangle); // find center of current triangle

		if (iters == 0) {
			// if recursed all the way down, add point
			points.add(centroid);
		} else {
			// else, subdivide triangle into two right angle triangle
			// and add the points for each
			double[][] sub1 = new double[][] { triangle[0], midpoint(triangle[0], triangle[2]), triangle[1] };
			double[][] sub2 = new double[][] { triangle[1], midpoint(triangle[0], triangle[2]), triangle[2] };
			points.addAll(subdivide(sub1, iters - 1));
			points.addAll(subdivide(sub2, iters - 1));
		}

		return points;
	}

	private static double[] midpoint(double[] p1, double[] p2) {
		return new double[] { (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2 };
	}

	private static double[] centroid(double[][] triangle) {
		double x = triangle[0][0] + triangle[1][0] + triangle[2][0];
		double y = triangle[0][1] + triangle[1][1] + triangle[2][1];
		return new double[] { x / 3, y / 3 };
	}

}
