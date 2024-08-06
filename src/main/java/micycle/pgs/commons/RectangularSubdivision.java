package micycle.pgs.commons;

import java.util.SplittableRandom;

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

	private final SplittableRandom r;

	private PShape division;

	public RectangularSubdivision(double width, double height, long seed) {
		this.width = width;
		this.height = height;
		divPt = new double[] { width / 2, height / 2 };
		cornerDist = Math.sqrt(Math.pow(divPt[0], 2) + Math.pow(divPt[1], 2));
		r = new SplittableRandom(seed);
	}

	public RectangularSubdivision(double width, double height, int maxDepth, long seed) {
		this.width = width;
		this.height = height;
		this.maxDepth = maxDepth;
		divPt = new double[] { width / 2, height / 2 };
		cornerDist = Math.sqrt(Math.pow(divPt[0], 2) + Math.pow(divPt[1], 2));
		r = new SplittableRandom(seed);
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

	private void sliceDivide(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, int depth, int base,
			boolean horizontal) {
		if (depth == base) {
			rect(x1, y1, x2, y2, x3, y3, x4, y4);
		} else {
			double randValHor = y1 + (y2 - y1) * r.nextDouble() * (1 - splitUniformity) + (y2 - y1) / 2 * splitUniformity;
			double randValVer = x1 + (x4 - x1) * r.nextDouble() * (1 - splitUniformity) + (x4 - x1) / 2 * splitUniformity;
			// horizontal split points
			double nx12 = (x1 + x2) / 2;
			double ny12 = randValHor;
			double nx34 = (x3 + x4) / 2;
			double ny34 = randValHor;
			// vertical split points
			double nx14 = randValVer;
			double ny14 = (y1 + y4) / 2;
			double nx23 = randValVer;
			double ny23 = (y2 + y3) / 2;
			// ternary assignments
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

			// center pt
			double center_aX = (a1 + a3) / 2;
			double center_aY = (a2 + a4) / 2;
			double center_bX = (b1 + b3) / 2;
			double center_bY = (b2 + b4) / 2;
			double distA = Math.sqrt(square(center_aX - divPt[0]) + square(center_aY - divPt[1]));
			double distB = Math.sqrt(square(center_bX - divPt[0]) + square(center_bY - divPt[1]));

			// Determine which way to split
			boolean r1horwider = (Math.abs(a5 - a1) <= Math.abs(a6 - a2)); // boolean for "is wider horizontally"
			boolean r2horwider = (Math.abs(b5 - b1) <= Math.abs(b6 - b2));
			boolean r1hor = (r.nextDouble() > 1 - balance ? r1horwider : r.nextDouble() > 0.5); // divide by wider or random
			boolean r2hor = (r.nextDouble() > 1 - balance ? r2horwider : r.nextDouble() > 0.5);

			// Recursive calls or make and display rectangles
			double probA = towardsCenter * distA / cornerDist + (1 - towardsCenter) * (1 - divProb);
			double probB = towardsCenter * distB / cornerDist + (1 - towardsCenter) * (1 - divProb);
			double toDivide1 = (depth < maxDepth - minDepth ? r.nextDouble() : 1.0); // 1 if less than min depth
			if (toDivide1 > probA) { // 0.2
				sliceDivide(a1, a2, a3, a4, a5, a6, a7, a8, depth - 1, base, r1hor);
			} else {
				rect(a1, a2, a3, a4, a5, a6, a7, a8);
			}

			double toDivide2 = (depth < maxDepth - minDepth ? r.nextDouble() : 1.0);
			if (toDivide2 > probB) {
				sliceDivide(b1, b2, b3, b4, b5, b6, b7, b8, depth - 1, base, r2hor);
			} else {
				rect(b1, b2, b3, b4, b5, b6, b7, b8);
			}
		}
	}

	private void rect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
		// rect(x1, y1, x3 - x1, y3 - y1);
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

	private static double square(double x) {
		return x * x;
	}

}
