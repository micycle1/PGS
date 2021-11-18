package micycle.pgs.utility;

import java.util.List;

import micycle.pgs.PGS_Conversion;
import micycle.pgs.color.RGB;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Corner rounding for PShape polygons. Replaces corners with arcs.
 * 
 * @author Michael Carleton
 *
 */
public final class CornerRounding {

	// Inspired by https://observablehq.com/@daformat/rounding-polygon-corners

	private CornerRounding() {
	}

	/**
	 * 
	 * @param shape
	 * @param extent 0...1
	 * @return
	 */
	public static PShape round(PShape shape, double extent) {
		if (shape.getChildCount() > 1) {
			PShape groupRound = new PShape(PConstants.GROUP);
			for (PShape child : shape.getChildren()) {
				groupRound.addChild(roundPolygon(child, extent));
			}
			return groupRound;
		}
		else {
			return roundPolygon(shape, extent);			
		}
	}
	
	private static PShape roundPolygon(PShape shape, double extent) {
		PShape rounded = new PShape(PShape.GEOMETRY);
		PGS_Conversion.setAllFillColor(rounded, RGB.PINK);
		rounded.beginShape();

		final List<PVector> l = PGS_Conversion.toPVector(shape);
		final int size = l.size();
		for (int i = 0; i < l.size(); i++) {
			final PVector p1 = l.get(Math.floorMod(i - 1, size));
			final PVector p2 = l.get(i);
			final PVector p3 = l.get((i + 1) % size);
			roundCorner(p1, p2, p3, extent, rounded);
		}
		rounded.endShape(PConstants.CLOSE);
		return rounded;
	}

	/**
	 * Round a triplet of points.
	 * 
	 * @param a
	 * @param b      middle/enclosed point
	 * @param c
	 * @param extent
	 * @param shape
	 */
	private static void roundCorner(PVector a, PVector b, PVector c, double extent, PShape shape) {
		if (clockwise(a, b, c)) {
			PVector temp = a;
			a = c;
			c = temp;
		}

		// line vectors
		PVector ab = PVector.sub(a, b);
		PVector cb = PVector.sub(c, b);

		float theta = PApplet.acos(ab.dot(cb) / (ab.mag() * cb.mag())); // same as a.angleBetween(a, c)

//		final float maxRadius = PApplet.min(ab.div(2).mag(), cb.div(2).mag());
//		extent = extent * maxRadius;
		final float extentF = (float) extent;

		final PVector A = PVector.add(b, ab.mult(extentF / ab.mag())); // where circle touches AB
		final PVector C = PVector.add(b, cb.mult(extentF / cb.mag())); // where circle touches CB

		PVector vBC = PVector.sub(C, b);

		final float lengthBC = vBC.mag();
		final float lengthBD = PApplet.cos(theta / 2); // length
		final float lengthFD = PApplet.sin(theta / 2) * lengthBD; // length
		final float lengthBF = lengthBD * lengthBD; // length
		final float r = lengthFD / (lengthBF / lengthBC); // radius of circle
		PVector tangent = ab.normalize().rotate(PConstants.HALF_PI).mult(r); // tangent to ab
		tangent.add(A); // find tangent to ab at point A -- this is circle center

		shape.vertex(C.x, C.y);
		final float a1 = angleBetween(C, tangent);
		final float a2 = angleBetween(A, tangent);
		sampleArc(tangent, r, a1, a2, shape);
		shape.vertex(A.x, A.y);
	}

	/**
	 * Sub-sample an arc into invididual vertices.
	 * 
	 * @param center
	 * @param r          arc radius
	 * @param startAngle
	 * @param endAngle
	 * @param shape
	 */
	private static void sampleArc(PVector center, float r, float startAngle, float endAngle, PShape shape) {
		if (startAngle > endAngle) {
			startAngle -= PConstants.TWO_PI;
		}
		final float n = 4; // every n degrees // TODO magic constant
		final float angleInc = (endAngle - startAngle) / (360 / n);
		float angle = startAngle;
		while (angle < endAngle) {
			shape.vertex(r * PApplet.cos(angle) + center.x, center.y + r * PApplet.sin(angle));
			angle += angleInc;
		}
	}

	/**
	 * Are the three given points in clockwise order?
	 * 
	 * @param p1
	 * @param p2 middle point
	 * @param p3
	 * @return true if the points consitute a "right turn" or clockwise orientation
	 */
	private static boolean clockwise(PVector p1, PVector p2, PVector p3) {
		float val = (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y);
		return val < 0;
	}

	/**
	 * East = 0; North = -1/2PI; West = -PI; South = -3/2PI | 1/2PI
	 * 
	 * @param tail PVector Coordinate 1.
	 * @param head PVector Coordinate 2.
	 * @return float θ in radians.
	 */
	private static float angleBetween(PVector tail, PVector head) {
		float a = PApplet.atan2(tail.y - head.y, tail.x - head.x);
		if (a < 0) {
			a += PConstants.TWO_PI;
		}
		return a;
	}

}
