package micycle.pgs.commons;

import micycle.pgs.color.ColorUtils;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Balanced triangle subdivision.
 * 
 * @author Michael Carleton
 *
 */
public class TriangleSubdivision {

	// implements openprocessing.org/sketch/970715

	private static final double VARIANCE = 0.1; // gaussian variance

	private final int maxDiv;
	/** Probability that a triangle will subdivide */
	private final double divProb = 0.85;

	// Base seed we’ll derive deterministic values from
	private final long baseSeed;

	private PShape division;
	private final double width, height;

	// salts to separate different random purposes
	private static final long SALT_DIVIDE = 0xA0761D6478BD642FL;
	private static final long SALT_LERP = 0xE7037ED1A0B428DBL;

	public TriangleSubdivision(double width, double height, int maxDepth, long seed) {
		this.width = width;
		this.height = height;
		maxDiv = maxDepth;
		this.baseSeed = seed;
	}

	public PShape divide() {
		return divide(mix64(baseSeed) % 2 == 0);
	}

	public PShape divide(boolean flip) {
		// flip==false -> diagonal TL->BR (original)
		// flip==true -> diagonal TR->BL (the other diagonal)
		division = new PShape(PConstants.GROUP);

		if (!flip) {
			// diagonal from top-left (0,0) to bottom-right (width,height)
			divideTriangle(0, 0, width, 0, width, height, maxDiv, 1); // top-right triangle
			divideTriangle(0, 0, 0, height, width, height, maxDiv, 1); // bottom-left triangle
		} else {
			// diagonal from top-right (width,0) to bottom-left (0,height)
			divideTriangle(width, 0, width, height, 0, height, maxDiv, 1); // right-bottom triangle
			divideTriangle(0, 0, width, 0, 0, height, maxDiv, 1); // left-top triangle
		}

		return division;
	}

	private void divideTriangle(double x1, double y1, double x2, double y2, double x3, double y3, int depth, int base) {
		final Triangle tri = new Triangle(x1, y1, x2, y2, x3, y3);
		if (depth == base) {
			division.addChild(tri.getShape());
		} else {
			// Deterministic per-triangle Gaussian value for "should we divide?"
			final double toDivide = triGaussian(0.5, 0.25, tri.p1, tri.p2, tri.p3, SALT_DIVIDE);
			if (toDivide < divProb) {
				tri.computeOppositePoint();
				divideTriangle(tri.d.x, tri.d.y, tri.l.x, tri.l.y, tri.n1.x, tri.n1.y, depth - 1, base);
				divideTriangle(tri.d.x, tri.d.y, tri.l.x, tri.l.y, tri.n2.x, tri.n2.y, depth - 1, base);
			} else {
				division.addChild(tri.getShape());
			}
		}
	}

	// Deterministic Gaussian based on triangle vertices and a salt
	private double triGaussian(double mean, double sd, PVector a, PVector b, PVector c, long salt) {
		long x = hashTriangle(baseSeed, a, b, c, salt);
		// Two uniforms via splitmix-like stepping for Box-Muller
		long r1 = mix64(x + 0x9E3779B97F4A7C15L);
		long r2 = mix64(x + 2L * 0x9E3779B97F4A7C15L);
		double u1 = toUnit(r1);
		double u2 = toUnit(r2);
		// Guard against log(0)
		u1 = (u1 <= 1e-15) ? 1e-15 : (u1 >= 1.0) ? 1.0 - 1e-15 : u1;
		double z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
		return z * sd + mean;
	}

	// Convert 64-bit random bits to double in [0,1)
	private static double toUnit(long bits) {
		// 53 significant bits for double fraction
		return (bits >>> 11) * 0x1.0p-53;
	}

	// Hash the triangle in a stable way (order-sensitive, which is fine because
	// construction is deterministic)
	private static long hashTriangle(long seed, PVector p1, PVector p2, PVector p3, long salt) {
		long h = seed ^ salt;
		h = mix64(h ^ Float.floatToIntBits(p1.x));
		h = mix64(h ^ Float.floatToIntBits(p1.y));
		h = mix64(h ^ Float.floatToIntBits(p2.x));
		h = mix64(h ^ Float.floatToIntBits(p2.y));
		h = mix64(h ^ Float.floatToIntBits(p3.x));
		h = mix64(h ^ Float.floatToIntBits(p3.y));
		return h;
	}

	// Strong 64-bit mixer (SplitMix64 finalizer)
	private static long mix64(long z) {
		z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
		z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
		return z ^ (z >>> 31);
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

			// Deterministic per-triangle Gaussian for lerp t
			float randVal = PApplet.constrain((float) triGaussian(0.5, VARIANCE, p1, p2, p3, SALT_LERP), 0, 1);

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
			triangle.setFill(ColorUtils.composeColor(237, 50, 162));
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
