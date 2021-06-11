package micycle.pgs;

import static micycle.pgs.PGS_Conversion.toPShape;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pgs.color.RGB;
import micycle.pgs.utility.RandomPolygon;
import micycle.pgs.utility.Star;
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
		PShape shape = new PShape(PShape.GEOMETRY);
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
		PShape heart = new PShape(PShape.GEOMETRY);
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

		heart.endShape();
		return heart;
	}

}
