package micycle.pgs;

import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.linemerge.LineMerger;
import org.locationtech.jts.shape.GeometricShapeBuilder;
import org.locationtech.jts.shape.fractal.HilbertCurveBuilder;
import org.locationtech.jts.shape.fractal.KochSnowflakeBuilder;
import org.locationtech.jts.shape.fractal.SierpinskiCarpetBuilder;
import org.locationtech.jts.util.GeometricShapeFactory;

import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import micycle.pgs.PGS_Contour.OffsetStyle;
import micycle.pgs.color.Colors;
import micycle.pgs.commons.BezierShapeGenerator;
import micycle.pgs.commons.PEdge;
import micycle.pgs.commons.RandomPolygon;
import micycle.pgs.commons.RandomSpaceFillingCurve;
import micycle.pgs.commons.Star;
import micycle.spacefillingcurves.SierpinskiFiveSteps;
import micycle.spacefillingcurves.SierpinskiFourSteps;
import micycle.spacefillingcurves.SierpinskiTenSteps;
import micycle.spacefillingcurves.SierpinskiThreeSteps;
import micycle.spacefillingcurves.SpaceFillingCurve;
import micycle.srpg.SRPolygonGenerator;
import net.jafama.FastMath;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Construct uncommon/interesting 2D geometries (beyond those offered in
 * Processing).
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Construction {

	private PGS_Construction() {
	}

	private static final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();

	static {
		shapeFactory.setNumPoints(PGS.SHAPE_SAMPLES);
	}

	/**
	 * Generates a random simple convex polygon (n-gon).
	 * 
	 * @param n         number of polygon vertices/sides
	 * @param maxWidth  maximum width of generated random polygon
	 * @param maxHeight maximum height of generated random polygon
	 * @return a PShape representing the generated polygon
	 * @see {@link #createRandomPolygonExact(int, double, double)}
	 *      createRandomPolygonExact()} to specify exact dimensions.
	 */
	public static PShape createRandomPolygon(int n, double maxWidth, double maxHeight) {
		return createRandomPolygon(n, maxWidth, maxHeight, System.nanoTime());
	}

	/**
	 * Generates a random simple convex polygon (n-gon), having a given random seed.
	 * 
	 * @param n         number of polygon vertices/sides
	 * @param maxWidth  maximum width of generated random polygon
	 * @param maxHeight maximum height of generated random polygon
	 * @param seed      a seed value used to generate the random polygon (optional)
	 * @return a PShape representing the generated polygon
	 * @see {@link #createRandomPolygonExact(int, double, double, long)
	 *      createRandomPolygonExact()} to specify exact dimensions
	 */
	public static PShape createRandomPolygon(int n, double maxWidth, double maxHeight, long seed) {
		return PGS_Transformation.translateEnvelopeTo(PGS_Conversion.fromPVector(RandomPolygon.generateRandomConvexPolygon(n, maxWidth, maxHeight, seed)),
				maxWidth / 2, maxHeight / 2);
	}

	/**
	 * Generates a random simple convex polygon (n-gon), where the output's bounding
	 * box has the dimensions of those specified.
	 * 
	 * @param n      number of polygon vertices/sides
	 * @param width  width of generated random polygon
	 * @param height height of generated random polygon
	 * @return a PShape representing the generated polygon
	 */
	public static PShape createRandomPolygonExact(int n, double width, double height) {
		return createRandomPolygonExact(n, width, height, System.nanoTime());
	}

	/**
	 * Generates a random simple convex polygon (n-gon), where the output's bounding
	 * box has the dimensions of those specified.
	 * 
	 * @param n      number of polygon vertices/sides
	 * @param width  width of generated random polygon
	 * @param height height of generated random polygon
	 * @param seed   a seed value used to generate the random polygon
	 * @return a PShape representing the generated polygon
	 */
	public static PShape createRandomPolygonExact(int n, double width, double height, long seed) {
		return PGS_Transformation.resize(PGS_Conversion.fromPVector(RandomPolygon.generateRandomConvexPolygon(n, width, height, seed)), width, height);
	}

	/**
	 * Generates an N-sided regular polygon.
	 * 
	 * @param n       number of sides
	 * @param centerX centre point X
	 * @param centerY centre point Y
	 * @param width   polygon width
	 * @since 2.0
	 */
	public static PShape createRegularPolyon(int n, double centerX, double centerY, double width) {
		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(n);
		shapeFactory.setCentre(new Coordinate(centerX, centerY));
		shapeFactory.setWidth(width * 2);
		shapeFactory.setHeight(width * 2);
		double ia = (Math.PI * (n - 2)) / n;
		// flat edge facing down
		shapeFactory.setRotation(ia / 2);

		return toPShape(shapeFactory.createCircle());
	}

	/**
	 * Generates a <b>supercircle</b> (also known as a superellipse or Lamé curve)
	 * shape centered at the specified coordinates.
	 * <p>
	 * The supercircle is defined by the equation |x/a|<sup>n</sup> +
	 * |y/a|<sup>n</sup> = 1, where the <code>power</code> (n) controls the shape's
	 * roundness:
	 * <ul>
	 * <li><b>power &lt; 1</b>: produces star-like (concave) forms</li>
	 * <li><b>power = 1</b>: produces a square (with rounded corners due to
	 * sampling)</li>
	 * <li><b>power &gt; 1</b>: increasingly resembles a circle as power
	 * increases</li>
	 * </ul>
	 * 
	 * @param centerX  the x-coordinate of the supercircle center
	 * @param centerY  the y-coordinate of the supercircle center
	 * @param diameter the diameter of the supercircle (distance from side to side)
	 * @param power    the exponent controlling the "circularity" (roundness or
	 *                 squareness) of the shape; values &lt; 1 yield star-like
	 *                 shapes, 1 is a square, values &gt; 1 become more circular
	 * @return a {@link PShape} representing the generated supercircle
	 */
	public static PShape createSupercircle(double centerX, double centerY, double diameter, double power) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setCentre(new Coordinate(centerX, centerY));
		shapeFactory.setWidth(diameter);
		shapeFactory.setHeight(diameter);

		final double outerC = Math.PI * 2 * diameter / 2;
		final int samples = (int) Math.max(PGS.SHAPE_SAMPLES, Math.ceil(outerC / BufferParameters.DEFAULT_QUADRANT_SEGMENTS));
		shapeFactory.setNumPoints(samples);

		return toPShape(shapeFactory.createSupercircle(power));
	}

	/**
	 * Generates a <b>supershape</b> using the superformula, centered at the
	 * specified coordinates.
	 * <p>
	 * The superformula is a versatile 2D mathematical equation capable of
	 * describing an enormous variety of shapes—ranging from rounded polygons and
	 * starfish to petals and snowflakes—depending on parameter choices. See
	 * <a href="https://en.wikipedia.org/wiki/Superformula">Superformula
	 * (Wikipedia)</a> or <a href="http://paulbourke.net/geometry/supershape/">Paul
	 * Bourke's page</a> for details.
	 * <p>
	 * Brief parameter effects:
	 * <ul>
	 * <li>Equal n-values pinch the form; reducing them increases pinching.</li>
	 * <li>n₁ larger than n₂/n₃ produces bloated forms.</li>
	 * <li>Very large n₁, n₂, and n₃ create polygonal shapes.</li>
	 * <li>Unequal n-values yield asymmetric forms.</li>
	 * <li>n₁ smaller than n₂/n₃ yields smooth starfish shapes.</li>
	 * </ul>
	 *
	 * @param centerX the x-coordinate for the center of the supershape
	 * @param centerY the y-coordinate for the center of the supershape
	 * @param radius  maximum radius of the supershape (outermost point)
	 * @param m       symmetry number (e.g. 3 for 3-pointed, 4 for 4-pointed shapes)
	 * @param n1      superformula parameter controlling general shape (pinching,
	 *                inflation)
	 * @param n2      superformula parameter affecting shape symmetry
	 * @param n3      superformula parameter affecting shape symmetry
	 * @return a {@link PShape} instance representing the generated supershape
	 */
	public static PShape createSuperShape(double centerX, double centerY, double radius, double m, double n1, double n2, double n3) {
		// http://paulbourke.net/geometry/supershape/
		PShape shape = new PShape(PShape.PATH);
		shape.setFill(true);
		shape.setFill(Colors.WHITE);
		shape.beginShape();

		final int points = 180;
		final double angleInc = Math.PI * 2 / points;
		double angle = 0;
		while (angle < Math.PI * 2) {
			double r;
			double t1, t2;

			t1 = FastMath.cos(m * angle / 4);
			t1 = Math.abs(t1);
			t1 = Math.pow(t1, n2);

			t2 = FastMath.sin(m * angle / 4);
			t2 = Math.abs(t2);
			t2 = Math.pow(t2, n3);

			r = Math.pow(t1 + t2, 1 / n1);
			if (Math.abs(r) != 0) {
				r *= radius; // multiply r (0...1) by (max) radius
//				r = radius/r;
				shape.vertex((float) (centerX + r * FastMath.cos(angle)), (float) (centerY + r * FastMath.sin(angle)));
			}

			angle += angleInc;
		}

		shape.endShape();
		return shape;
	}

	/**
	 * Creates an elliptical arc polygon—a filled "slice" of an ellipse defined by
	 * the specified arc and the two radii connecting its endpoints to the ellipse
	 * center.
	 * <p>
	 * The arc begins at the given orientation angle (measured from the 12 o'clock
	 * position, in radians) and extends counterclockwise by the specified angular
	 * size. If the arc angle is equal to or greater than 2&pi;, a full ellipse is
	 * generated.
	 * <p>
	 * The resulting {@link PShape} is suitable for rendering a pie-slice or sector
	 * from an ellipse or circle.
	 *
	 * @param centerX     the x-coordinate of the ellipse center
	 * @param centerY     the y-coordinate of the ellipse center
	 * @param width       the full width (diameter on the x-axis) of the ellipse
	 * @param height      the full height (diameter on the y-axis) of the ellipse
	 * @param orientation the starting angle (in radians), where 0 corresponds to
	 *                    "12 o'clock" (upwards)
	 * @param angle       the arc angle (in radians), defining the sweep of the arc;
	 *                    if equal or greater than 2&pi;, a full ellipse is produced
	 * @return a {@link PShape} representing the elliptical arc polygon
	 */
	public static PShape createArc(double centerX, double centerY, double width, double height, double orientation, double angle) {
		if (angle == 0) {
			return new PShape();
		}
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setCentre(new Coordinate(centerX, centerY));
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);

		final double outerC = Math.min(Math.PI * 2, angle) * Math.max(width, height) / 2;
		final int samples = (int) Math.max(PGS.SHAPE_SAMPLES, Math.ceil(outerC / BufferParameters.DEFAULT_QUADRANT_SEGMENTS));
		shapeFactory.setNumPoints(samples);

		PShape out;
		if (angle >= PConstants.TWO_PI) {
			out = toPShape(shapeFactory.createCircle());
		} else {
			out = toPShape(shapeFactory.createArcPolygon(-Math.PI / 2 + orientation, angle));
		}

		out.setStroke(false);
		return out;
	}

	/**
	 * 
	 * Creates a <i>Taijitu</i> shape (a geometric representation of the Taoist
	 * symbol of yin and yang).
	 * 
	 * @param centerX the x-coordinate of the center of the shape
	 * @param centerY the y-coordinate of the center of the shape
	 * @param radius  the radius of the shape
	 * @return a PShape representing the Taijitu shape
	 * @since 1.4.0
	 */
	public static PShape createTaijitu(double centerX, double centerY, double radius) {
		Coordinate center = new Coordinate(centerX, centerY);
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setCentre(center);
		shapeFactory.setWidth(radius * 2);
		shapeFactory.setHeight(radius * 2);

		Geometry a = shapeFactory.createArcPolygon(-Math.PI / 2, Math.PI);
		Geometry b = createCirclePoly(center.x, centerY + radius / 2, radius / 2);
		Geometry c = createCirclePoly(center.x, centerY - radius / 2, radius / 2);

		Geometry yinG = a.union(b).difference(c).getGeometryN(0);
		AffineTransformation t = AffineTransformation.rotationInstance(Math.PI, centerX, centerY);
		Geometry yangG = t.transform(yinG);

		PShape yin = toPShape(yinG);
		PShape yang = toPShape(yangG);
		yin.setFill(0);
		yang.setFill(255);

		return PGS_Conversion.flatten(yin, yang);
	}

	/**
	 * Creates an Arbelos shape, a mathematical figure bounded by three semicircles.
	 * <p>
	 * The position of the central notch is arbitrary and can be located anywhere
	 * along the diameter.
	 * 
	 * @param centerX       the x-coordinate of the center of the shape
	 * @param centerY       the y-coordinate of the center of the shape
	 * @param radius        radius of the largest (enclosing) circle
	 * @param notchPosition the fractional position, between 0 and 1, along the
	 *                      diameter where the notch will be
	 * @since 1.4.0
	 * @return a PShape representing the Arbelos shape
	 */
	public static PShape createArbelos(double centerX, double centerY, double radius, double notchPosition) {
		centerY += radius / 2;
		final double rA = radius * notchPosition;
		final double rB = radius * (1 - notchPosition);
		final double xA = (centerX - radius) + rA;
		final double xB = xA + rA + rB;

		final Coordinate center = new Coordinate(centerX, centerY);
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setCentre(center);
		shapeFactory.setWidth(radius * 2);
		shapeFactory.setHeight(radius * 2);

		Geometry a = shapeFactory.createArcPolygon(Math.PI, Math.PI);
		Geometry b = createCirclePoly(xA, centerY, rA);
		Geometry c = createCirclePoly(xB, centerY, rB);
		Geometry curve = a.difference(b).difference(c).buffer(1e-3);
		@SuppressWarnings("unchecked")
		List<Polygon> polygons = PolygonExtracter.getPolygons(curve);
		polygons.sort((m, n) -> Integer.compare(n.getNumPoints(), m.getNumPoints()));

		PShape arbelos = toPShape(polygons.get(0));
		arbelos.setStroke(false);
		return arbelos;
	}

	/**
	 * Creates a star shape, having a specified number of rays.
	 * 
	 * @param centerX     The x coordinate of the center
	 * @param centerY     The y coordinate of the center
	 * @param numRays     The number of rays the star has
	 * @param innerRadius The inner radius of the star
	 * @param outerRadius The outer radius of the star
	 * @param roundness   A roundness value between 0.0 and 1.0, for the inner and
	 *                    outer corners of the star.
	 * @return The star shape
	 */
	public static PShape createStar(double centerX, double centerY, int numRays, double innerRadius, double outerRadius, double roundness) {
		roundness = Math.max(Math.min(1, roundness), 0);
		final PShape shape = Star.createStarShape(centerX, centerY, innerRadius, outerRadius, numRays, roundness);
		shape.setFill(true);
		shape.setFill(255);
		return shape;
	}

	/**
	 * Generates a "blobbie" shape—a deformable, organic, blobby closed curve
	 * defined by four parameters.
	 * <p>
	 * The blobbie shape is based on functions involving cosine waves of different
	 * frequencies, allowing for smooth, natural-looking deformations. Typical
	 * results are lobed or petal-like closed shapes useful in generative graphics
	 * or modeling organic phenomena.
	 * <p>
	 * To avoid self-intersections, the sum of the parameters <code>a</code> and
	 * <code>b</code> should be less than 1. If self-intersection occurs, the result
	 * will attempt to be geometrically fixed but may yield unexpected forms. See
	 * <a href="http://paulbourke.net/geometry/blobbie/">Paul Bourke: Blobbie</a>
	 * for background.
	 *
	 * <pre>{@code
	 * r(theta) = (maxWidth/2) * [1 + a*cos(2θ + c) + b*cos(3θ + d)]
	 * }</pre>
	 *
	 * @param centerX  the x-coordinate of the blobbie center
	 * @param centerY  the y-coordinate of the blobbie center
	 * @param maxWidth the maximum width (diameter) of the blobbie
	 * @param a        2-lobed deformation parameter; together with b, controls the
	 *                 main shape undulations (<b>a + b < 1</b> for simple forms)
	 * @param b        3-lobed deformation parameter; see note above
	 * @param c        phase offset for the 2-lobed term
	 * @param d        phase offset for the 3-lobed term
	 * @return a {@link PShape} representing the generated blobbie shape
	 * @since 1.3.0
	 */
	public static PShape createBlobbie(double centerX, double centerY, double maxWidth, double a, double b, double c, double d) {
		// http://paulbourke.net/geometry/blobbie/
		final double cirumference = 2 * Math.PI * maxWidth / 2;
		final int samples = (int) (cirumference / 2); // 1 point every 2 distance
		double dt = Math.PI * 2 / samples;

		final CoordinateList blobbieCoords = new CoordinateList();

		for (int i = 0; i <= samples; i++) {
			final double theta = i * dt;
			final double r = maxWidth / 2 * (1 + a * FastMath.cos(2 * theta + c) + b * FastMath.cos(3 * theta + d));

			// polar to cartesian
			double x = r * FastMath.cos(theta) + centerX;
			double y = r * FastMath.sin(theta) + centerY;
			blobbieCoords.add(new Coordinate(x, y), false);
		}
		blobbieCoords.closeRing();

		if (blobbieCoords.size() > 1) {
			Geometry g = PGS.GEOM_FACTORY.createPolygon(blobbieCoords.toCoordinateArray());
			if (a + b > 1) {
				g = GeometryFixer.fix(g); // fix self intersection
			}
			PShape blob = PGS_Conversion.toPShape(g);
			blob.setStroke(false);
			return blob;
		} else {
			return new PShape();
		}
	}

	/**
	 * Creates a classic "heart" shape using a parametric curve, centered and scaled
	 * as specified.
	 * <p>
	 * The curve is based on the well-known heart equations, producing a symmetric
	 * heart that is widest at the center and comes to a point below.
	 * <p>
	 * The maximum width parameter controls the distance across the heart at its
	 * widest part.
	 * <p>
	 * See <a href="https://mathworld.wolfram.com/HeartCurve.html">Heart Curve
	 * (MathWorld)</a> for the mathematical reference.
	 *
	 * @param centerX the x-coordinate for the center of the heart shape
	 * @param centerY the y-coordinate for the center of the heart shape
	 * @param width   the maximum width (horizontal extent) of the heart
	 * @return a {@link PShape} representing the generated heart curve
	 * @since 1.1.0
	 */
	public static PShape createHeart(final double centerX, final double centerY, final double width) {
		// https://mathworld.wolfram.com/HeartCurve.html
		PShape heart = new PShape(PShape.PATH);
		heart.setFill(true);
		heart.setFill(Colors.WHITE);
		heart.beginShape();

		final double length = 6.3855 * width; // Arc length of parametric curve from wolfram alpha
		final int points = (int) length / 2; // sample every 2 units along curve (roughly)
		final double angleInc = Math.PI * 2 / points;
		double angle = 0;
		while (angle < Math.PI * 2) {
			final double s = FastMath.sin(angle);
			double vx = s * s * s;
			double vy = 13 * FastMath.cos(angle) - 5 * FastMath.cos(2 * angle) - 2 * FastMath.cos(3 * angle) - FastMath.cos(4 * angle);
			vy /= 16; // normalise to 1
			heart.vertex((float) (centerX + vx * width / 2), (float) (centerY - vy * width / 2));
			angle += angleInc;
		}

		heart.endShape(PConstants.CLOSE);
		PShape out = PGS_Processing.fix(heart); // fix pinch
		out.setStroke(false);
		return out;
	}

	/**
	 * Creates a teardrop shape using a parametric polar curve, centered and scaled
	 * as specified.
	 * <p>
	 * This method generates a classic teardrop or droplet outline, where the
	 * parameter {@code m} controls the taper and sharpness of the pointed end.
	 * Lower values for {@code m} (such as 2) create softer drops, while higher
	 * values (up to 5) yield sharper, pointier ends.
	 * <p>
	 * The generated {@link PShape} is suitable for use in generative design,
	 * infographics, and iconography.
	 * <p>
	 * See <a href="https://mathworld.wolfram.com/TeardropCurve.html">Teardrop Curve
	 * (MathWorld)</a> for mathematical background.
	 *
	 * @param centerX the x-coordinate of the center of the teardrop shape
	 * @param centerY the y-coordinate of the center of the teardrop shape
	 * @param height  the full vertical height of the teardrop, from its base to tip
	 * @param m       the order/tapering factor of the curve; recommended range is
	 *                2–5 for visually pleasing shapes
	 * @return a {@link PShape} representing the teardrop outline
	 * @since 1.4.0
	 */
	public static PShape createTeardrop(final double centerX, final double centerY, double height, final double m) {
		// https://mathworld.wolfram.com/TeardropCurve.html
		height /= 2; // get height in terms of radius
		PShape curve = new PShape(PShape.PATH);
		curve.setFill(true);
		curve.setFill(Colors.WHITE);
		curve.beginShape();
		final double angleInc = Math.PI * 2 / 360;
		double angle = 0;

		while (angle < Math.PI * 2) {
			double x = FastMath.cos(angle);
			double y = FastMath.sin(angle) * FastMath.pow(FastMath.sin(0.5 * angle), m);
			curve.vertex((float) (centerX + x * height), (float) (centerY + y * height * 1));
			angle += angleInc;
		}
		curve.endShape(PConstants.CLOSE);

		return PGS_Transformation.rotate(curve, new PVector((float) centerX, (float) centerY), -Math.PI / 2);
	}

	/**
	 * Creates a gear shape from a parametric gear curve.
	 * 
	 * @param centerX The x coordinate of the center of the gear
	 * @param centerY The y coordinate of the center of the gear
	 * @param radius  maximum radius of gear teeth
	 * @param n       number of gear teeth
	 * @return the gear shape
	 * @since 1.4.0
	 */
	public static PShape createGear(final double centerX, final double centerY, final double radius, final int n) {
		// https://mathworld.wolfram.com/GearCurve.html
		PShape curve = new PShape(PShape.PATH);
		curve.setFill(true);
		curve.setFill(Colors.WHITE);
		curve.beginShape();

		final double cirumference = 2 * Math.PI * radius;
		final int samples = (int) (cirumference / 5); // 1 point every 5 distance
		final double angleInc = Math.PI * 2 / samples;
		double angle = 0;

		final double a = 1; // wolfram default
		final double b = 10; // wolfram default
		while (angle < Math.PI * 2) {
			double r = a + (1 / b) * FastMath.tanh(b * FastMath.sin(n * angle));
			r *= radius;
			curve.vertex((float) (centerX + r * FastMath.cos(angle)), (float) (centerY + r * FastMath.sin(angle)));
			angle += angleInc;
		}
		curve.endShape(PConstants.CLOSE);

		return curve;
	}

	/**
	 * Creates a joined ring (a "donut") shape.
	 * 
	 * @param centerX     the x coordinate of the center
	 * @param centerY     the y coordinate of the center
	 * @param outerRadius radius of ring exterior
	 * @param innerRadius radius of ring hole
	 * @return the ring shape
	 * @since 1.1.3
	 */
	public static PShape createRing(double centerX, double centerY, double outerRadius, double innerRadius) {
		return createRing(centerX, centerY, outerRadius, innerRadius, 0, PConstants.TWO_PI);
	}

	/**
	 * Creates an (un)joined ring shape.
	 * 
	 * @param centerX     the x coordinate of the center
	 * @param centerY     the y coordinate of the center
	 * @param outerRadius radius of ring exterior
	 * @param innerRadius radius of ring hole
	 * @param orientation start angle/orientation in radians (where 0 is 12 o'clock)
	 * @param angle       size of the ring arc angle in radians
	 * @return the ring shape
	 * @since 1.1.3
	 */
	public static PShape createRing(double centerX, double centerY, double outerRadius, double innerRadius, double orientation, double angle) {
		final double outerR = Math.max(outerRadius, innerRadius);
		final double innerR = Math.min(outerRadius, innerRadius);

		final double outerC = Math.min(Math.PI * 2, angle) * outerR;
		final int outerSamples = (int) Math.max(PGS.SHAPE_SAMPLES, Math.ceil(outerC / BufferParameters.DEFAULT_QUADRANT_SEGMENTS));
		final double innerC = Math.min(Math.PI * 2, angle) * innerR;
		final int innerSamples = (int) Math.max(PGS.SHAPE_SAMPLES, Math.ceil(innerC / BufferParameters.DEFAULT_QUADRANT_SEGMENTS));

		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(outerSamples);
		shapeFactory.setCentre(new Coordinate(centerX, centerY));
		shapeFactory.setWidth(outerR * 2);
		shapeFactory.setHeight(outerR * 2);

		final Geometry outer;
		if (angle > PConstants.TWO_PI - 1e-5) {
			outer = shapeFactory.createCircle();
		} else {
			outer = shapeFactory.createArcPolygon(-Math.PI / 2 + orientation, angle);
		}

		shapeFactory.setWidth(innerR * 2);
		shapeFactory.setHeight(innerR * 2);
		shapeFactory.setNumPoints(innerSamples);
		final Geometry inner = shapeFactory.createCircle();

		PShape out = toPShape(outer.difference(inner));
		out.setStroke(false);
		return out;
	}

	/**
	 * Creates a sponge-like porous structure.
	 * <p>
	 * The sponge structure is formed by randomly merging adjacent cells of a
	 * Voronoi tessellation and then smoothing them; the final structure is obtained
	 * by subtracting that result from a rectangle.
	 *
	 * @param width      the width of the sponge bounds
	 * @param height     the height of the sponge bounds
	 * @param generators the number of generator points for the underlying Voronoi
	 *                   tessellation. Should be >5.
	 * @param thickness  thickness of sponge structure walls
	 * @param smoothing  the cell smoothing factor which determines how rounded the
	 *                   cells are. a value of 6 is a good starting point.
	 * @param classes    the number of classes to use for the cell merging process,
	 *                   where lower results in more merging (or larger "blob-like"
	 *                   shapes).
	 * @param seed       the seed for the random number generator
	 * @return the sponge shape
	 * @since 1.4.0
	 */
	public static PShape createSponge(double width, double height, int generators, double thickness, double smoothing, int classes, long seed) {
		// A Simple and Effective Geometric Representation for Irregular Porous
		// Structure Modeling
		List<PVector> points = PGS_PointSet.random(thickness, thickness / 2, width - thickness / 2, height - thickness / 2, generators, seed);
		if (points.size() < 6) {
			return new PShape();
		}
		PShape voro = PGS_Voronoi.innerVoronoi(points, 2);

		List<PShape> blobs = PGS_Conversion.getChildren(PGS_Meshing.stochasticMerge(voro, classes, seed)).stream().map(c -> {
			c = PGS_Morphology.buffer(c, -thickness / 2, OffsetStyle.MITER);
			if (smoothing != 0) {
				c = PGS_Morphology.smoothGaussian(c, smoothing);
			}
			return c;
		}).collect(Collectors.toList());

		/*
		 * Although faster, can't use .simpleSubtract() here because holes (cell
		 * islands) are *sometimes* nested.
		 */
		PShape s = PGS_ShapeBoolean.subtract(PGS.createRect(0, 0, width, height), PGS_Conversion.flatten(blobs));
		s.setStroke(false);
		return s;
	}

	/**
	 * Creates an linear/archimedean spiral shape, where the distance between any 2
	 * successive windings is constant.
	 * 
	 * @param centerX     the x coordinate of the spiral center point / origin
	 * @param centerY     the y coordinate of the spiral center point / origin
	 * @param coils       number of coils/rotations in the spiral
	 * @param outerRadius the outer-most radius of the spiral / final coil
	 * @return a stroked PATH PShape
	 * @since 1.2.0
	 */
	public static PShape createLinearSpiral(double centerX, double centerY, double coils, double outerRadius) {
		// from https://stackoverflow.com/a/13901170/9808792
		final int chord = 2; // distance between points to plot
		// spiral is rotated by this number of radians
		final double rotation = 0; // -Math.PI / 2;
		// value of theta corresponding to end of last coil
		final double thetaMax = coils * 2 * Math.PI;
		// How far to step away from center for each side.
		final double awayStep = Math.max(outerRadius, 1) / thetaMax;
		final int direction = 1; // either +1 (CW) or -1 (CCW)

		final CoordinateList coords = new CoordinateList();
		coords.add(new Coordinate(centerX, centerY), false);

		double theta = chord / awayStep;
		/*
		 * For every side, step around and away from center. Start at the angle
		 * corresponding to a distance of chord away from centre.
		 */
		while (theta <= thetaMax) {
			final double away = awayStep * theta; // How far away from center
			final double around = direction * theta + rotation; // How far around the center
			// Convert 'around' and 'away' to X and Y.
			double x = centerX + FastMath.cos(around) * away;
			double y = centerY + FastMath.sin(around) * away;
			coords.add(new Coordinate(x, y), false);

			double delta = (-2 * away + Math.sqrt(4 * away * away + 8 * awayStep * chord)) / (2 * awayStep);
			theta += delta;
		}

		if (coords.size() > 1) {
			Geometry lineString = PGS.GEOM_FACTORY.createLineString(coords.toCoordinateArray());
			PShape spiral = PGS_Conversion.toPShape(lineString);
			spiral.setStrokeWeight(10);
			spiral.setStroke(Colors.WHITE);
			spiral.setStrokeCap(PConstants.ROUND);
			return spiral;
		} else {
			return new PShape();
		}
	}

	/**
	 * Creates Fermat's spiral, a parabolic spiral which is symmetrical about the
	 * origin.
	 * 
	 * @param centerX     the x coordinate of the spiral center point / origin
	 * @param centerY     the y coordinate of the spiral center point / origin
	 * @param coils       number of coils/rotations in the spiral
	 * @param outerRadius the outer-most radius of the spiral / final coil
	 * @return a stroked PATH PShape
	 * @since 1.2.0
	 */
	public static PShape createFermatSpiral(double centerX, double centerY, double coils, double outerRadius) {
		double thetaEnd = (Math.PI * coils);
		final int samples = (int) (50 * coils);
		double dt = thetaEnd / samples;

		final CoordinateList yin = new CoordinateList();
		final CoordinateList yang = new CoordinateList();

		final double z = outerRadius / Math.sqrt(thetaEnd);

		for (int i = 0; i <= samples; i++) {
			final double theta = i * dt; // archimedean spiral
			final double r = Math.sqrt(theta) * z; // Specific to made a Fermat Spiral.

			// polar to cartesian
			double x = r * FastMath.cos(theta) + centerX;
			double y = r * FastMath.sin(theta) + centerY;
			yin.add(new Coordinate(x, y), false);

			x = -r * FastMath.cos(theta) + centerX;
			y = -r * FastMath.sin(theta) + centerY;
			yang.add(new Coordinate(x, y), false);
		}

		Collections.reverse(yin);
		yin.addAll(yang);

		LineString lineString = PGS.GEOM_FACTORY.createLineString(yin.toCoordinateArray());
		PShape spiral = PGS_Conversion.toPShape(lineString);
		spiral.setStrokeWeight(10);
		spiral.setStroke(Colors.WHITE);
		spiral.setStrokeCap(PConstants.ROUND);
		spiral.setFill(false);
		return spiral;
	}

	/**
	 * Creates a rectangular-shaped spiral shape.
	 * 
	 * @param x       x position of the top-left of spiral
	 * @param y       y position of the top-left of spiral
	 * @param width   width of outer-most coil
	 * @param height  height of outer-most coil
	 * @param spacing the distance between successive coils
	 * @return a stroked PATH PShape with SQUARE stroke cap and MITER joins
	 * @since 1.3.0
	 */
	public static PShape createRectangularSpiral(float x, float y, float width, float height, float spacing) {
		float xx = -width / 2;
		float yy = -height / 2;
		int count = 0;
		// below is used to dictate spiral orientation & starting point
//		if (count == 1) {
//			xx *= -1;
//		}
//		if (count == 2) {
//			xx *= -1;
//			yy *= -1;
//		}
//		if (count == 3) {
//			yy *= -1;
//		}
		final float offX = x + width / 2;
		final float offY = y + height / 2;
		final PShape spiral = new PShape(PShape.PATH);
		spiral.setFill(false);
		spiral.setStroke(true);
		spiral.setStrokeWeight(5);
//		spiral.setStrokeWeight(spacing * 0.333f);
		spiral.setStroke(Colors.WHITE);
		spiral.setStrokeJoin(PConstants.MITER);
		spiral.setStrokeCap(PConstants.SQUARE);
		spiral.beginShape();
		while (width > 0 && height > 0) {
			int dir = count % 4;
			spiral.vertex(xx + offX, yy + offY);
			if (dir == 0) {
				xx += width;
				width -= spacing;
			} else if (dir == 1) {
				yy += height;
				height -= spacing;
			} else if (dir == 2) {
				xx -= width;
				width -= spacing;
			} else if (dir == 3) {
				yy -= height;
				height -= spacing;
			}
			count++;
		}
		spiral.endShape();
		return spiral;
	}

	/**
	 * Creates a random space-filling curve.
	 * <p>
	 * A space-filling curve is a continuous curve that (in this case) traverses
	 * every cell of a grid exactly once.
	 * 
	 * @param nColumns   number of columns in the underlying grid
	 * @param nRows      number of rows in the underlying grid
	 * @param cellWidth  visual/pixel width of each cell
	 * @param cellHeight visual/pixel width of each cell
	 * @return a stroked PATH PShape
	 * @see #createRandomSFCurve(int, int, double, double, long)
	 * @since 1.4.0
	 */
	public static PShape createRandomSFCurve(int nColumns, int nRows, double cellWidth, double cellHeight) {
		return createRandomSFCurve(nColumns, nRows, cellWidth, cellHeight, System.nanoTime());
	}

	/**
	 * Creates a random space-filling curve, having a specific random seed.
	 * <p>
	 * A space-filling curve is a continuous curve that (in this case) traverses
	 * every cell of a grid exactly once.
	 * 
	 * @param nColumns   number of columns in the underlying grid
	 * @param nRows      number of rows in the underlying grid
	 * @param cellWidth  visual/pixel width of each cell
	 * @param cellHeight visual/pixel width of each cell
	 * @param seed       random seed
	 * @return a mitered stroked PATH PShape
	 * @see #createRandomSFCurve(int, int, double, double)
	 * @since 1.4.0
	 */
	public static PShape createRandomSFCurve(int nColumns, int nRows, double cellWidth, double cellHeight, long seed) {
		RandomSpaceFillingCurve factory = new RandomSpaceFillingCurve(nColumns, nRows, seed);
		PShape curve = factory.getCurve((float) cellWidth, (float) cellHeight);
		curve.setFill(255);
		curve.setStroke(0);
		curve.setStrokeWeight(3);

		return curve;
	}

	/**
	 * Generates a highly customisable random polygon based on a square grid of NxN
	 * cells.
	 * <p>
	 * The number of vertices of the polygon generated is not configurable, but
	 * depends on the size of <code>cells</code> and the percentage
	 * <code>markPercent</code>: for larger values of <code>cells</code>, the more
	 * vertices the polygon tends to have for a given markPercent.
	 * <p>
	 * Visually pleasing "random" polygons can be achieved by selecting fairly small
	 * values for markFraction, e.g., <code>markFraction=0.1</code> or even
	 * <code>markPercent=0.01</code>.
	 * 
	 * @param dimensions   pixel dimensions of the polygon in its longest axis.
	 * @param cells        the number of cells in the X and Y directions of the
	 *                     grid.
	 * @param markFraction The fraction of vertices marked on the grid. The larger
	 *                     the percentage, the more vertices the polygon tends to
	 *                     have. Should generally be between 0...0.5.
	 * @param holes        If true, generates a holes in the polygon.
	 * @param orthogonal   Whether the polygon vertices lie exactly on integer grid
	 *                     points and only form horizontal and vertical lines.
	 * @param smoothing    The number of rounds of corner cutting to apply to the
	 *                     polygon. A small positive integer value is recommended. A
	 *                     value of 3 is probably sufficient.
	 * @param depth        The number of rounds of recursive refinement to apply to
	 *                     the polygon generated. A small positive integer value is
	 *                     recommended. This is akin to increasing the depth of
	 *                     fractal curve.
	 * @param seed         the seed for the random number generator
	 * @since 1.4.0
	 */
	public static PShape createSuperRandomPolygon(double dimensions, int cells, double markFraction, int smoothing, int depth, boolean orthogonal,
			boolean holes, long seed) {
		Random r = new XoRoShiRo128PlusRandom(seed);
		SRPolygonGenerator generator = new SRPolygonGenerator(cells, cells, markFraction, holes, orthogonal, !orthogonal, smoothing, depth, !orthogonal, r);
		List<List<double[]>> rings = generator.getPolygon();
		PShape exterior = PGS_Conversion.fromArray(rings.get(0).toArray(new double[0][0]), true);
		List<PShape> interiorRings = rings.subList(1, rings.size()).stream().map(l -> PGS_Conversion.fromArray(l.toArray(new double[0][0]), true))
				.collect(Collectors.toList());

		PShape polygon = PGS_ShapeBoolean.simpleSubtract(exterior, PGS_Conversion.flatten(interiorRings));
		polygon = PGS_Transformation.resizeByMajorAxis(polygon, dimensions);
		polygon.setStroke(false);
		return PGS_Transformation.translateToOrigin(polygon);
	}

	/**
	 * Generates a smooth or spiky random polygon comprising Bezier curves.
	 * 
	 * @param nPoints   The number of bezier curves the polygon consists of.
	 * @param scale     Polygon scale. Determines the maximum width or height the
	 *                  polygon could have.
	 * @param radius    The radius relative to the distance between adjacent points.
	 *                  The radius is used to position the control points of the
	 *                  bezier curve, should be between 0 and 1. Larger values
	 *                  result in sharper features on the curve.
	 * @param spikiness A measure of the curve's smoothness. If 0, the curve's angle
	 *                  through each point will be the average between the direction
	 *                  to adjacent points. As it increases, the angle will be
	 *                  determined mostly by one adjacent point, making the curve
	 *                  more "spiky".
	 * @param seed      the seed for the random number generator
	 * @return the random polygon shape
	 * @since 1.4.0
	 */
	public static PShape createRandomBezierPolygon(int nPoints, double scale, double radius, double spikiness, long seed) {
		BezierShapeGenerator bsg = new BezierShapeGenerator(nPoints, 3, radius, spikiness);
		Geometry shape = PGS.GEOM_FACTORY.createPolygon(bsg.generate(false, false, scale, seed));
		PShape poly = toPShape(shape);
		poly.setStroke(false);
		return PGS_Transformation.translateToOrigin(poly);
	}

	/**
	 * Creates a closed Sierpiński curve (a recursive space-filling curve), having a
	 * user-defined degree/order.
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
		final double[][] tri2 = new double[][] { { width + centerX, centerY }, { width + centerX, width + centerY }, { centerX, width + centerY } };

		// get points for each half of square recursively
		final List<double[]> half1 = subdivide(tri1, curveOrder);
		final List<double[]> half2 = subdivide(tri2, curveOrder);

		half1.addAll(half2); // combine points

		final PShape curve = new PShape(PShape.PATH);
		curve.setFill(true);
		curve.setFill(Colors.WHITE);
		curve.beginShape();
		half1.forEach(p -> curve.vertex((float) p[0], (float) p[1]));
		curve.endShape(PConstants.CLOSE);

		return curve;
	}

	/**
	 * Creates a Hilbert Curve shape, a type of plane-filling curve.
	 * 
	 * @param width  pixel width of the curve
	 * @param height pixel height of the curve
	 * @param order  order of the hilbert curve. should be at least 0
	 * @return a stroked PATH PShape, anchored at (0, 0)
	 * @since 1.3.0
	 */
	public static PShape createHilbertCurve(double width, double height, int order) {
		final GeometricShapeBuilder builder = new HilbertCurveBuilder(PGS.GEOM_FACTORY);
		builder.setExtent(new Envelope(0, width, 0, height));
		final int points = (int) Math.ceil(Math.exp((order * Math.log(4))) * 3); // inverse of recursionLevelForSize()
		builder.setNumPoints(points);
		PShape out = toPShape(builder.getGeometry());
		out.setStroke(255);
		out.setStrokeWeight(4);
		return out;
	}

	/**
	 * Creates a Sierpiński Carpet shape, a type of plane fractal.
	 * 
	 * @param width  pixel width of the curve
	 * @param height pixel height of the curve
	 * @param order  the number of recursive subdivisions (at least 0, probably no
	 *               more than 5)
	 * @return carpet shape, anchored at (0, 0)
	 * @since 1.3.0
	 */
	public static PShape createSierpinskiCarpet(double width, double height, int order) {
		final GeometricShapeBuilder builder = new SierpinskiCarpetBuilder(PGS.GEOM_FACTORY);
		builder.setExtent(new Envelope(0, width, 0, height));
		final int points = (int) Math.ceil(Math.exp((order * Math.log(4))) * 3); // inverse of recursionLevelForSize()
		builder.setNumPoints(points);
		PShape out = toPShape(builder.getGeometry());
		out.setStroke(false);
		return out;
	}

	/**
	 * Creates a Koch Snowflake shape, a fractal curve.
	 * 
	 * @param width  pixel width of the curve's envelope
	 * @param height pixel width of the curve's envelope
	 * @param order  the number of recursive subdivisions (at least 0)
	 * @return snowflake shape, whose envelope is anchored at (0, 0)
	 * @since 1.3.0
	 */
	public static PShape createKochSnowflake(double width, double height, int order) {
		order++; //
		final GeometricShapeBuilder builder = new KochSnowflakeBuilder(PGS.GEOM_FACTORY);
		builder.setExtent(new Envelope(0, width, 0, height));
		final int points = (int) Math.ceil(Math.exp((order * Math.log(4))) * 3); // inverse of recursionLevelForSize()
		builder.setNumPoints(points);
		PShape out = toPShape(builder.getGeometry());
		out.setStroke(false);
		return out;
	}

	public enum SierpinskiTriCurveType {
		TRI, TETRA, PENTA, DECA;
	}

	/**
	 * Creates one of a family of trifurcating Sierpinski curves.
	 * 
	 * @param type  the type of tri-curve: {TRI, TETRA, PENTA, DECA}
	 * @param width pixel width of the curve's envelope
	 * @param order the number of recursive subdivisions (at least 1)
	 * @return curve shape, anchored at (0, 0)
	 * @since 1.3.0
	 */
	public static PShape createSierpinskiTriCurve(SierpinskiTriCurveType type, double width, int order) {
		SpaceFillingCurve fractal = switch (type) {
			case TRI -> new SierpinskiThreeSteps(width, width);
			case TETRA -> new SierpinskiFourSteps(width, width);
			case PENTA -> new SierpinskiFiveSteps(width, width);
			case DECA -> new SierpinskiTenSteps(width, width);
		};

		fractal.setN(order);

		fractal.start();

		final List<PEdge> edges = new ArrayList<>(fractal.getLineSegments().size());
		fractal.getLineSegments().forEach(l -> edges.add(new PEdge(l[0], l[1], l[2], l[3])));
		/*
		 * Fractal is comprised of unconnected line segments -- they need to be merged
		 * together. Use linemerger since PGS.fromEdges() does not support unclosed edge
		 * collections.
		 */
		final LineMerger lm = new LineMerger();
		edges.forEach(e -> {
			final LineString l = PGS.createLineString(e.a, e.b);
			lm.add(l);
		});
		LineString path = (LineString) lm.getMergedLineStrings().iterator().next();
		CoordinateList list = new CoordinateList(path.getCoordinates());
		list.closeRing();

		PShape out = toPShape(PGS.GEOM_FACTORY.createLinearRing(list.toCoordinateArray()));
		out.setStroke(false);
		out = PGS_Transformation.resizeByWidth(out, width);
		out = PGS_Transformation.translateToOrigin(out);
		return out;
	}

	/**
	 * Shortcut for creating a rectangle with uniformly rounded corners, using
	 * {@link PConstants#CORNER} mode.
	 * <p>
	 * This convenience method creates a rectangle where all four corners have the
	 * same radius {@code r}. The rectangle uses Processing's
	 * {@link PConstants#CORNER} mode coordinates: {@code (a, b)} specify the
	 * top-left corner, and {@code c} and {@code d} specify width and height,
	 * respectively.
	 *
	 * @param a the x-coordinate of the top-left corner of the rectangle
	 * @param b the y-coordinate of the top-left corner of the rectangle
	 * @param c the width of the rectangle
	 * @param d the height of the rectangle
	 * @param r the uniform radius to be applied to all four corners (<code>0</code>
	 *          gives a regular rectangle)
	 * @return a {@link PShape} representing the rounded rectangle
	 * @since 2.1
	 */
	public static PShape createRect(double a, double b, double c, double d, double r) {
		return createRect(PConstants.CORNER, a, b, c, d, r);
	}

	/**
	 * Creates a rectangle with specified corner radii, in any Processing-style
	 * rectangle mode.
	 * <p>
	 * The meaning of the {@code a}, {@code b}, {@code c}, and {@code d} parameters
	 * depends on the {@code rectMode}:
	 * <ul>
	 * <li>{@link PConstants#CORNER}: {@code (a, b)} is the top-left corner;
	 * {@code c} is width, {@code d} is height</li>
	 * <li>{@link PConstants#CORNERS}: {@code (a, b)} is the top-left corner;
	 * {@code (c, d)} is the bottom-right corner</li>
	 * <li>{@link PConstants#CENTER}: {@code (a, b)} is the rectangle center;
	 * {@code c} is width, {@code d} is height</li>
	 * <li>{@link PConstants#RADIUS}: {@code (a, b)} is the center; {@code c} is
	 * half width, {@code d} is half height</li>
	 * </ul>
	 * Each corner radius parameter refers to a specific corner:
	 * <ul>
	 * <li>{@code tl} – top-left corner radius</li>
	 * <li>{@code tr} – top-right corner radius</li>
	 * <li>{@code br} – bottom-right corner radius</li>
	 * <li>{@code bl} – bottom-left corner radius</li>
	 * </ul>
	 *
	 * @param rectMode rectangle mode as in Processing; one of
	 *                 {@link PConstants#CORNER}, {@link PConstants#CORNERS},
	 *                 {@link PConstants#CENTER}, or {@link PConstants#RADIUS}
	 * @param a        first coordinate: x (or center x, depending on mode)
	 * @param b        second coordinate: y (or center y, depending on mode)
	 * @param c        width, x2, or half-width (see mode above)
	 * @param d        height, y2, or half-height (see mode above)
	 * @param r        the uniform radius to be applied to all four corners
	 * @return a {@link PShape} representing the rounded rectangle
	 * @since 2.1
	 */
	public static PShape createRect(int rectMode, double a, double b, double c, double d, double r) {
		return rect(rectMode, a, b, c, d, r, r, r, r);
	}

	static PShape rect(int rectMode, double a, double b, double c, double d, double tl, double tr, double br, double bl) {
		double hradius, vradius;
		switch (rectMode) {
			case PConstants.CORNERS :
				break;
			case PConstants.CORNER :
				c += a;
				d += b;
				break;
			case PConstants.RADIUS :
				hradius = c;
				vradius = d;
				c = a + hradius;
				d = b + vradius;
				a -= hradius;
				b -= vradius;
				break;
			case PConstants.CENTER :
				hradius = c / 2.0;
				vradius = d / 2.0;
				c = a + hradius;
				d = b + vradius;
				a -= hradius;
				b -= vradius;
				break;
		}
		if (a > c) {
			double t = a;
			a = c;
			c = t;
		}
		if (b > d) {
			double t = b;
			b = d;
			d = t;
		}
		double maxRounding = Math.min((c - a) / 2, (d - b) / 2);
		tl = Math.min(tl, maxRounding);
		tr = Math.min(tr, maxRounding);
		br = Math.min(br, maxRounding);
		bl = Math.min(bl, maxRounding);
		return rectImpl((float) a, (float) b, (float) c, (float) d, (float) tl, (float) tr, (float) br, (float) bl);
	}

	private static PShape rectImpl(float x1, float y1, float x2, float y2, float tl, float tr, float br, float bl) {
		PShape sh = new PShape(PShape.PATH);
		sh.setFill(true);
		sh.setFill(Colors.WHITE);
		sh.beginShape();
		// Top edge and top-right corner
		if (tr != 0) {
			sh.vertex(x2 - tr, y1);
			sh.quadraticVertex(x2, y1, x2, y1 + tr);
		} else {
			sh.vertex(x2, y1);
		}
		// Right edge and bottom-right
		if (br != 0) {
			sh.vertex(x2, y2 - br);
			sh.quadraticVertex(x2, y2, x2 - br, y2);
		} else {
			sh.vertex(x2, y2);
		}
		// Bottom edge and bottom-left
		if (bl != 0) {
			sh.vertex(x1 + bl, y2);
			sh.quadraticVertex(x1, y2, x1, y2 - bl);
		} else {
			sh.vertex(x1, y2);
		}
		// Left edge and top-left
		if (tl != 0) {
			sh.vertex(x1, y1 + tl);
			sh.quadraticVertex(x1, y1, x1 + tl, y1);
		} else {
			sh.vertex(x1, y1);
		}
		sh.endShape(PConstants.CLOSE);
		return sh;
	}

	/**
	 * Creates a polygon finely approximating a circle.
	 * 
	 * @since 2.0
	 */
	static Polygon createCircle(Coordinate c, double r) {
		return createCircle(c.x, c.y, r, 0.5);
	}

	/**
	 * Creates a circle of radius r centered on (x,y).
	 * 
	 * @since 2.0
	 */
	public static PShape createCircle(double x, double y, double r) {
		return toPShape(createCirclePoly(x, y, r));
	}

	/**
	 * Creates a polygon finely approximating a circle.
	 * 
	 * @since 2.0
	 */
	static Polygon createCirclePoly(double x, double y, double r) {
		return createCircle(x, y, r, 0.5); // 0.5 still very generous
	}

	/**
	 * Creates a polygon approximating a circle.
	 *
	 * <p>
	 * The `maxDeviation` parameter controls the maximum acceptable deviation of any
	 * segment of the polygon from the true arc of the circle. Smaller values of
	 * `maxDeviation` result in smoother circles but require more vertices,
	 * potentially increasing computational cost.
	 * </p>
	 *
	 * @param x            The x-coordinate of the circle's center.
	 * @param y            The y-coordinate of the circle's center.
	 * @param r            The radius of the circle.
	 * @param maxDeviation The maximum acceptable deviation of a segment from the
	 *                     true arc of the circle.
	 * @return A polygon approximating the specified circle.
	 */
	static Polygon createCircle(double x, double y, double r, final double maxDeviation) {
		// Calculate the number of points based on the radius and maximum deviation.
		int nPts = (int) Math.ceil(2 * Math.PI / Math.acos(1 - maxDeviation / r));
		nPts = Math.max(nPts, 21); // min of 21 points for tiny circles
		final int circumference = (int) (Math.PI * r * 2);
		if (nPts > circumference * 2) {
			// AT MOST 1 point every half pixel
			nPts = circumference * 2;
		}

		Coordinate[] pts = new Coordinate[nPts + 1];
		for (int i = 0; i < nPts; i++) {
			double ang = i * (2 * Math.PI / nPts);
			double px = r * FastMath.cos(ang) + x;
			double py = r * FastMath.sin(ang) + y;
			pts[i] = new Coordinate(px, py);
		}
		pts[nPts] = new Coordinate(pts[0]); // Close the circle

		return PGS.GEOM_FACTORY.createPolygon(pts);
	}

	/**
	 * Sample a circular arc from startAngle→endAngle (radians), producing a List of
	 * PVectors so that no straight‐line chord deviates by more than maxDev pixels
	 * from the true circle.
	 *
	 * @param cx         center x
	 * @param cy         center y
	 * @param r          radius (>0)
	 * @param startAngle start angle in radians
	 * @param endAngle   end angle in radians
	 * @param maxDev     maximum sagitta error in pixels (>0)
	 * @return List of PVectors, inclusive of both endpoints.
	 * @since 2.1
	 */
	static List<PVector> arcPoints(double cx, double cy, double r, double startAngle, double endAngle, double maxDev) {
		List<PVector> pts = new ArrayList<>();

		// 1) Compute raw sweep and direction
		double rawRange = endAngle - startAngle;
		int dir = rawRange >= 0 ? +1 : -1;
		double range = Math.abs(rawRange);

		// 2) Cap at a full circle if the user overshoots
		if (range > 2.0 * Math.PI) {
			range = 2.0 * Math.PI;
		}

		// 3) Solve for maximum step so that sagitta s = R*(1–cos(dθ/2)) ≤ maxDev
		// ⇒ dθ ≤ 2 * acos(1 – maxDev / R)
		double cosArg = 1.0 - maxDev / r;
		// clamp to avoid NaNs
		cosArg = Math.max(-1.0, Math.min(1.0, cosArg));
		double maxStep = 2.0 * Math.acos(cosArg);

		// 4) If that step is wider than the entire arc, just use one segment
		if (maxStep >= range) {
			maxStep = range;
		}

		// 5) Determine how many segments are needed (at least 1)
		int nSeg = Math.max(1, (int) Math.ceil(range / maxStep));

		// 6) Compute the actual signed step
		double step = dir * (range / nSeg);

		// 7) Sample from i=0…nSeg (inclusive) so endpoints are exact
		for (int i = 0; i <= nSeg; i++) {
			double a = startAngle + i * step;
			double px = cx + r * FastMath.cos(a);
			double py = cy + r * FastMath.sin(a);
			pts.add(new PVector((float) px, (float) py));
		}

		return pts;
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
