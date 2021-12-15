package micycle.pgs.commons;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import micycle.pgs.color.RGB;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Balanced triangle subdivision.
 * https://openprocessing.org/sketch/970715
 * 
 * @author Michael Carleton
 *
 */
public class TriangleSubdivision {

	private static final double VARIANCE = 0.1; // gaussian variance

	private final int maxDiv;
	/** Probability that a triangle will subdivide */
	private final double divProb = 0.85;

	final RandomGenerator random;

	private PShape division;
	private final double width, height;

	public TriangleSubdivision(double width, double height, int maxDepth, long seed) {
		this.width = width;
		this.height = height;
		maxDiv = maxDepth;
		random = new Well19937c(seed);
	}

	public PShape divide() {
		division = new PShape(PConstants.GROUP);
		divideTriangle(0, 0, width, 0, width, height, maxDiv, 1); // top right half
		divideTriangle(0, 0, 0, height, width, height, maxDiv, 1); // bottom left half
		return division;
	}

	private void divideTriangle(double x1, double y1, double x2, double y2, double x3, double y3, int depth, int base) {
		final Triangle tri = new Triangle(x1, y1, x2, y2, x3, y3);
		if (depth == base) {
			division.addChild(tri.getShape());
		} else {
			final double toDivide = randomGaussian(0.5, 0.25);
			if (toDivide < divProb) {
				tri.computeOppositePoint();
				divideTriangle(tri.d.x, tri.d.y, tri.l.x, tri.l.y, tri.n1.x, tri.n1.y, depth - 1, base);
				divideTriangle(tri.d.x, tri.d.y, tri.l.x, tri.l.y, tri.n2.x, tri.n2.y, depth - 1, base);
			} else {
				division.addChild(tri.getShape());
			}
		}
	}

	private double randomGaussian(double mean, double sd) {
		return random.nextGaussian() * sd + mean;
	}

	private class Triangle {

		private final PVector p1, p2, p3;
		private PVector d, l;
		private PVector n1, n2;

		private Triangle(double x1, double y1, double x2, double y2, double x3, double y3) {
			p1 = new PVector((float) x1, (float) y1);
			p2 = new PVector((float) x2, (float) y2);
			p3 = new PVector((float) x3, (float) y3);
		}

		private void computeOppositePoint() {
			double d12 = p1.dist(p2);
			double d23 = p2.dist(p3);
			double d13 = p3.dist(p1);
			double maxLength = Math.max(Math.max(d12, d23), d13);

			float randVal = PApplet.constrain((float) randomGaussian(0.5, VARIANCE), 0, 1);
			if (maxLength == d12) {
				d = p3.copy();
				n1 = p1.copy();
				n2 = p2.copy();
				l = PVector.lerp(p1, p2, randVal);
			} else if (maxLength == d23) {
				d = p1.copy();
				n1 = p2.copy();
				n2 = p3.copy();
				l = PVector.lerp(p2, p3, randVal);
			} else {
				d = p2.copy();
				n1 = p1.copy();
				n2 = p3.copy();
				l = PVector.lerp(p1, p3, randVal);
			}
		}

		private PShape getShape() {
			final PShape triangle = new PShape(PShape.PATH);
			triangle.setFill(true);
			triangle.setFill(RGB.composeColor(255, 0, 255, 80));
			triangle.setStroke(true);
			triangle.setStroke(255);
			triangle.beginShape();
			triangle.vertex(p1.x, p1.y);
			triangle.vertex(p2.x, p2.y);
			triangle.vertex(p3.x, p3.y);
			triangle.endShape(PConstants.CLOSE);
			return triangle;
		}
	}
}
