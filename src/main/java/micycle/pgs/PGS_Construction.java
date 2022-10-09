package micycle.pgs;

import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.linemerge.LineMerger;
import org.locationtech.jts.shape.GeometricShapeBuilder;
import org.locationtech.jts.shape.fractal.HilbertCurveBuilder;
import org.locationtech.jts.shape.fractal.KochSnowflakeBuilder;
import org.locationtech.jts.shape.fractal.SierpinskiCarpetBuilder;
import org.locationtech.jts.util.GeometricShapeFactory;

import micycle.pgs.color.RGB;
import micycle.pgs.commons.PEdge;
import micycle.pgs.commons.RandomPolygon;
import micycle.pgs.commons.Star;
import micycle.spacefillingcurves.SierpinskiFiveSteps;
import micycle.spacefillingcurves.SierpinskiFourSteps;
import micycle.spacefillingcurves.SierpinskiTenSteps;
import micycle.spacefillingcurves.SierpinskiThreeSteps;
import micycle.spacefillingcurves.SpaceFillingCurve;
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
	 * @return
	 * @see {@link #createRandomPolygonExact(int, double, double)} to specify exact
	 *      dimensions
	 */
	public static PShape createRandomPolygon(int n, double maxWidth, double maxHeight) {
		return PGS_Transformation.translateEnvelopeTo(
				PGS_Conversion.fromPVector(RandomPolygon.generateRandomConvexPolygon(n, maxWidth, maxHeight)), maxWidth / 2, maxHeight / 2);
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
	 * Creates a supercircle shape.
	 * 
	 * @param centerX centre point X
	 * @param centerY centre point Y
	 * @param width
	 * @param height
	 * @param power   circularity of super circle. Values less than 1 create
	 *                star-like shapes; power=1 is a square;
	 * @return
	 */
	public static PShape createSupercircle(double centerX, double centerY, double width, double height, double power) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setCentre(new Coordinate(centerX, centerY));
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);

		final double outerC = Math.PI * 2 * Math.max(width, height) / 2;
		final int samples = (int) Math.max(PGS.SHAPE_SAMPLES, Math.ceil(outerC / BufferParameters.DEFAULT_QUADRANT_SEGMENTS));
		shapeFactory.setNumPoints(samples);

		return toPShape(shapeFactory.createSupercircle(power));
	}

	/**
	 * Creates a supershape PShape. The parameters feed into the superformula, which
	 * is a simple 2D analytical expression allowing to draw a wide variety of
	 * geometric and natural shapes (starfish, petals, snowflakes) by choosing
	 * suitable values relevant to few parameters.
	 * 
	 * @param centerX centre point X
	 * @param centerY centre point Y
	 * @param radius  maximum radius
	 * @param m       specifies the rotational symmetry of the shape (3 = 3 sided; 4
	 *                = 4 sided)
	 * @param n1      supershape parameter 1
	 * @param n2      supershape parameter 2
	 * @param n3      supershape parameter 3
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
	 * Creates an elliptical arc polygon (a slice of a circle). The polygon is
	 * formed from the specified arc of an ellipse and the two radii connecting the
	 * endpoints to the centre of the ellipse.
	 * 
	 * @param centerX     centre point X
	 * @param centerY     centre point Y
	 * @param width
	 * @param height
	 * @param orientation start angle/orientation in radians (where 0 is 12 o'clock)
	 * @param angle       size of the arc angle in radians
	 * @return
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
	 * Creates a heart shape.
	 * 
	 * @param centerX The x coordinate of the center of the heart
	 * @param centerY The y coordinate of the center of the heart
	 * @param width   Maximum width of the widest part of the heart
	 * @return
	 * @since 1.1.0
	 */
	public static PShape createHeart(final double centerX, final double centerY, final double width) {
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
			heart.vertex((float) (centerX + vx * width / 2), (float) (centerY - vy * width / 2));
			angle += angleInc;
		}

		heart.endShape(PConstants.CLOSE);
		return heart;
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
	public static PShape createRing(double centerX, double centerY, double outerRadius, double innerRadius, double orientation,
			double angle) {
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
			double x = centerX + Math.cos(around) * away;
			double y = centerY + Math.sin(around) * away;
			coords.add(new Coordinate(x, y), false);

			double delta = (-2 * away + Math.sqrt(4 * away * away + 8 * awayStep * chord)) / (2 * awayStep);
			theta += delta;
		}

		if (coords.size() > 1) {
			Geometry lineString = PGS.GEOM_FACTORY.createLineString(coords.toCoordinateArray());
			PShape spiral = PGS_Conversion.toPShape(lineString);
			spiral.setStrokeWeight(10);
			spiral.setStroke(RGB.WHITE);
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
			double x = r * Math.cos(theta) + centerX;
			double y = r * Math.sin(theta) + centerY;
			yin.add(new Coordinate(x, y), false);

			x = -r * Math.cos(theta) + centerX;
			y = -r * Math.sin(theta) + centerY;
			yang.add(new Coordinate(x, y), false);
		}

		Collections.reverse(yin);
		yin.addAll(yang);

		LineString lineString = PGS.GEOM_FACTORY.createLineString(yin.toCoordinateArray());
		PShape spiral = PGS_Conversion.toPShape(lineString);
		spiral.setStrokeWeight(10);
		spiral.setStroke(RGB.WHITE);
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
	 * @return a stroked PATH PShape
	 * @since 1.2.1
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
		spiral.setStrokeWeight(spacing * 0.333f);
		spiral.setStroke(RGB.WHITE);
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
	 * Creates a Hilbert Curve shape, a type of plane-filling curve.
	 * 
	 * @param width  pixel width of the curve
	 * @param height pixel height of the curve
	 * @param order  order of the hilbert curve. should be at least 0
	 * @return a stroked PATH PShape, anchored at (0, 0)
	 * @since 1.2.1
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
	 * @since 1.2.1
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
	 * @since 1.2.1
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
	 * @since 1.2.1
	 */
	public static PShape createSierpinskiTriCurve(SierpinskiTriCurveType type, double width, int order) {
		SpaceFillingCurve fractal;

		switch (type) {
			default :
			case TRI :
				fractal = new SierpinskiThreeSteps(width, width);
				break;
			case TETRA :
				fractal = new SierpinskiFourSteps(width, width);
				break;
			case PENTA :
				fractal = new SierpinskiFiveSteps(width, width);
				break;
			case DECA :
				fractal = new SierpinskiTenSteps(width, width);
		}
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

	static Polygon createEllipse(Coordinate center, double width, double height) {
		shapeFactory.setCentre(center);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return shapeFactory.createEllipse();
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
