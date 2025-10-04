package micycle.pgs.commons;

import micycle.pgs.color.Colors;
import processing.core.PConstants;
import processing.core.PShape;

/**
 * Random rectangle subdivision.
 * 
 * @author Java port of this
 *         <a href="https://openprocessing.org/sketch/1115849">sketch</a> by
 *         Michael Carleton
 *
 */
public class RectangularSubdivision {

	private int minDepth = 3;
	private int maxDepth = 10; // ... 14
	private double divProb = 0.7; // 0...1
	private double splitUniformity = 0.15; // 0...1
	private double balance = 0.5; // 0...1
	private double towardsCenter = 0.5; // 0...1

	private final double[] divPt;
	private final double cornerDist;
	private final double width, height;

	// Base seed for deterministic per-rectangle randomness
	private final long baseSeed;

	private PShape division;

	// Salts to separate different purposes
	private static final long SALT_SPLIT_H = 0xA0761D6478BD642FL;
	private static final long SALT_SPLIT_V = 0xE7037ED1A0B428DBL;
	private static final long SALT_ORIENT_GATE = 0x8EBC6AF09C88C6E3L;
	private static final long SALT_ORIENT_BOOL = 0x589965CC75374CC3L;
	private static final long SALT_DIVIDE = 0x1D8E4E27C47D124FL;

	public RectangularSubdivision(double width, double height, long seed) {
		this.width = width;
		this.height = height;
		this.baseSeed = seed;
		divPt = new double[] { width / 2, height / 2 };
		cornerDist = Math.sqrt(divPt[0] * divPt[0] + divPt[1] * divPt[1]);
	}

	public RectangularSubdivision(double width, double height, int maxDepth, long seed) {
		this.width = width;
		this.height = height;
		this.maxDepth = maxDepth;
		this.baseSeed = seed;
		divPt = new double[] { width / 2, height / 2 };
		cornerDist = Math.sqrt(divPt[0] * divPt[0] + divPt[1] * divPt[1]);
	}

	/**
	 * Produces a new rectangular subdivision using the configured parameters.
	 *
	 * @return a GROUP PShape where each child shape is one rectangle
	 */
	public PShape divide() {
		division = new PShape(PConstants.GROUP);
		sliceDivide(0, 0, 0, height, width, height, width, 0, maxDepth, 1, true);
		return division;
	}

	private void sliceDivide(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, int depth, int base, boolean horizontal) {

		if (depth == base) {
			rect(x1, y1, x2, y2, x3, y3, x4, y4);
			return;
		}

		// Current rectangle bounds (order-independent)
		double minX = Math.min(Math.min(x1, x2), Math.min(x3, x4));
		double maxX = Math.max(Math.max(x1, x2), Math.max(x3, x4));
		double minY = Math.min(Math.min(y1, y2), Math.min(y3, y4));
		double maxY = Math.max(Math.max(y1, y2), Math.max(y3, y4));

		// Deterministic split positions for this rectangle
		double uH = rectUniform(minX, minY, maxX, maxY, SALT_SPLIT_H);
		double uV = rectUniform(minX, minY, maxX, maxY, SALT_SPLIT_V);
		double randValHor = y1 + (y2 - y1) * uH * (1 - splitUniformity) + (y2 - y1) / 2.0 * splitUniformity;
		double randValVer = x1 + (x4 - x1) * uV * (1 - splitUniformity) + (x4 - x1) / 2.0 * splitUniformity;

		// horizontal split points
		double nx12 = (x1 + x2) / 2.0;
		double ny12 = randValHor;
		double nx34 = (x3 + x4) / 2.0;
		double ny34 = randValHor;
		// vertical split points
		double nx14 = randValVer;
		double ny14 = (y1 + y4) / 2.0;
		double nx23 = randValVer;
		double ny23 = (y2 + y3) / 2.0;

		// Ternary assignments (current split direction)
		double a1 = x1;
		double a2 = y1;
		double a3 = (horizontal ? nx12 : x2);
		double a4 = (horizontal ? ny12 : y2);
		double a5 = (horizontal ? nx34 : nx23);
		double a6 = (horizontal ? ny34 : ny23);
		double a7 = (horizontal ? x4 : nx14);
		double a8 = (horizontal ? y4 : ny14);

		double b1 = (horizontal ? x2 : nx14);
		double b2 = (horizontal ? y2 : ny14);
		double b3 = (horizontal ? nx12 : nx23);
		double b4 = (horizontal ? ny12 : ny23);
		double b5 = (horizontal ? nx34 : x3);
		double b6 = (horizontal ? ny34 : y3);
		double b7 = (horizontal ? x3 : x4);
		double b8 = (horizontal ? y3 : y4);

		// Centers and distances from division point
		double center_aX = (a1 + a3) / 2.0;
		double center_aY = (a2 + a4) / 2.0;
		double center_bX = (b1 + b3) / 2.0;
		double center_bY = (b2 + b4) / 2.0;
		double distA = Math.hypot(center_aX - divPt[0], center_aY - divPt[1]);
		double distB = Math.hypot(center_bX - divPt[0], center_bY - divPt[1]);

		// Child bounds (order-independent), used for stable randomness per child
		double[] boundsA = rectBounds(a1, a2, a3, a4, a5, a6, a7, a8);
		double[] boundsB = rectBounds(b1, b2, b3, b4, b5, b6, b7, b8);

		// Decide each child's next split orientation deterministically
		boolean r1horwider = (Math.abs(a5 - a1) <= Math.abs(a6 - a2)); // is wider horizontally
		boolean r2horwider = (Math.abs(b5 - b1) <= Math.abs(b6 - b2));
		boolean r1hor = chooseOrientation(boundsA, r1horwider);
		boolean r2hor = chooseOrientation(boundsB, r2horwider);

		// Recursive calls or make rectangles
		double probA = towardsCenter * distA / cornerDist + (1 - towardsCenter) * (1 - divProb);
		double probB = towardsCenter * distB / cornerDist + (1 - towardsCenter) * (1 - divProb);

		double toDivide1 = (depth < maxDepth - minDepth ? rectUniform(boundsA[0], boundsA[1], boundsA[2], boundsA[3], SALT_DIVIDE) : 1.0);
		if (toDivide1 > probA) {
			sliceDivide(a1, a2, a3, a4, a5, a6, a7, a8, depth - 1, base, r1hor);
		} else {
			rect(a1, a2, a3, a4, a5, a6, a7, a8);
		}

		double toDivide2 = (depth < maxDepth - minDepth ? rectUniform(boundsB[0], boundsB[1], boundsB[2], boundsB[3], SALT_DIVIDE) : 1.0);
		if (toDivide2 > probB) {
			sliceDivide(b1, b2, b3, b4, b5, b6, b7, b8, depth - 1, base, r2hor);
		} else {
			rect(b1, b2, b3, b4, b5, b6, b7, b8);
		}
	}

	private boolean chooseOrientation(double[] bounds, boolean widerHorizontal) {
		double minX = bounds[0], minY = bounds[1], maxX = bounds[2], maxY = bounds[3];
		double gate = rectUniform(minX, minY, maxX, maxY, SALT_ORIENT_GATE);
		if (gate > 1 - balance) {
			return widerHorizontal;
		}
		double coin = rectUniform(minX, minY, maxX, maxY, SALT_ORIENT_BOOL);
		return coin > 0.5;
	}

	private static double[] rectBounds(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
		double minX = Math.min(Math.min(x1, x2), Math.min(x3, x4));
		double maxX = Math.max(Math.max(x1, x2), Math.max(x3, x4));
		double minY = Math.min(Math.min(y1, y2), Math.min(y3, y4));
		double maxY = Math.max(Math.max(y1, y2), Math.max(y3, y4));
		return new double[] { minX, minY, maxX, maxY };
	}

	private double rectUniform(double minX, double minY, double maxX, double maxY, long salt) {
		long h = hashRect(minX, minY, maxX, maxY, salt);
		return toUnit(mix64(h));
	}

	private long hashRect(double minX, double minY, double maxX, double maxY, long salt) {
		long h = mix64(baseSeed ^ salt);
		h = mix64(h ^ Double.doubleToLongBits(minX));
		h = mix64(h ^ Double.doubleToLongBits(minY));
		h = mix64(h ^ Double.doubleToLongBits(maxX));
		h = mix64(h ^ Double.doubleToLongBits(maxY));
		return h;
	}

	private static long mix64(long z) {
		z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
		z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
		return z ^ (z >>> 31);
	}

	private static double toUnit(long bits) {
		return (bits >>> 11) * 0x1.0p-53; // [0,1)
	}

	private void rect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
		final PShape rect = new PShape(PShape.PATH);
		rect.setFill(true);
		rect.setFill(Colors.PINK);
		rect.setStroke(true);
		rect.setStroke(255);
		rect.beginShape();
		rect.vertex((float) x1, (float) y1);
		rect.vertex((float) x2, (float) y2);
		rect.vertex((float) x3, (float) y3);
		rect.vertex((float) x4, (float) y4);
		rect.endShape(PConstants.CLOSE);
		division.addChild(rect);
	}
}
