package micycle.pgs;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

import micycle.pgs.color.RGB;
import micycle.pgs.commons.DoyleSpiral;
import micycle.pgs.commons.RectangularSubdivision;
import micycle.pgs.commons.TriangleSubdivision;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Tiling, tessellation and subdivision of the plane using periodic or
 * non-periodic geometric shapes.
 * <p>
 * A tiling is created when a collection of plane figures (tiles) fills a plane
 * such that no gaps occur between the tiles and no two tiles overlap each
 * other.
 * 
 * @author Michael Carleton
 * @since 1.2.0
 */
public final class PGS_Tiling {

	// @formatter:off
	/**
	 * Rhombille Penrose tiling?
	 * Truchet tiling
	 * L-System too?
	 * Pentagonal tilings (upto 15)
	 * Gailiunas's Spiral Tilings
	 * https://openprocessing.org/browse/?time=anytime&type=all&q=tiling#
	 * https://openprocessing.org/browse/?time=anytime&type=all&q=tiling#
	 * overview https://www.nathaniel.ai/de-bruijn-projections/
	 */
	// @formatter:on

	// ideas @
	// https://demonstrations.wolfram.com/topic.html?topic=Art&start=21&limit=20&sortmethod=recent

	private PGS_Tiling() {
	}

	/**
	 * Produces a shape made by randomly subdividing a rectangular plane into
	 * rectangles.
	 * 
	 * @param width    width of the rectangular plane to subdivide
	 * @param height   height of the rectangular plane to subdivide
	 * @param maxDepth maximum number of subdivisions
	 * @return a GROUP PShape, where each child is a rectangle, that together
	 *         comprise the plane
	 * @see #rectSubdivsion(double, double, int, long) seeded rectSubdivsion()
	 */
	public static PShape rectSubdivsion(double width, double height, int maxDepth) {
		final RectangularSubdivision rectangularSubdivision = new RectangularSubdivision(width, height, maxDepth,
				System.currentTimeMillis());
		return rectangularSubdivision.divide();
	}

	/**
	 * 
	 * @param width
	 * @param height
	 * @param maxDepth
	 * @param seed
	 * @return
	 * @see #rectSubdivsion(double, double, int) non-seeded rectSubdivsion()
	 */
	public static PShape rectSubdivsion(double width, double height, int maxDepth, long seed) {
		final RectangularSubdivision rectangularSubdivision = new RectangularSubdivision(width, height, maxDepth, seed);
		return rectangularSubdivision.divide();
	}

	public static PShape triangleSubdivsion(double width, double height, int maxDepth) {
		final TriangleSubdivision subdivision = new TriangleSubdivision(width, height, maxDepth, System.currentTimeMillis());
		return subdivision.divide();
	}

	public static PShape triangleSubdivsion(double width, double height, int maxDepth, long seed) {
		final TriangleSubdivision subdivision = new TriangleSubdivision(width, height, maxDepth, seed);
		return subdivision.divide();
	}

	/**
	 * 
	 * @param width
	 * @param height
	 * @param depth
	 * @return
	 * @see #quadSubdivision(double, double, int, long) seeded quadSubdivision()
	 */
	public static PShape quadSubdivision(double width, double height, int depth) {
		return quadSubdivision(width, height, depth, System.currentTimeMillis());
	}

	/**
	 * Recursively and randomly subdivides the given/bounded plane into convex QUADS
	 * POLYGONS.
	 * 
	 * randomSubdivision() ?
	 * 
	 * @param width
	 * @param height
	 * @param depth
	 * @param seed
	 * @return
	 * @see #quadSubdivision(double, double, int) non-seeded quadSubdivision()
	 */
	public static PShape quadSubdivision(double width, double height, int depth, long seed) {
		// https://openprocessing.org/sketch/1045334
		final float w = (float) width;
		final float h = (float) height;
		final float off = 20;
		final SplittableRandom r = new SplittableRandom(seed);

		final PVector p1 = new PVector(off, off);
		final PVector p2 = new PVector(w - off, off);
		final PVector p3 = new PVector(w - off, h - off);
		final PVector p4 = new PVector(off, h - off);

		final PShape divisions = new PShape(PConstants.GROUP);
		divideRect(p1, p2, p3, p4, depth, divisions, r);
		PGS_Conversion.setAllFillColor(divisions, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(divisions, RGB.PINK, 2);
		return divisions;
	}

	/**
	 * Generates a Doyle spiral. A Doyle spiral fills the plane with closely packed
	 * circles, where the radius of each circle in a packing is proportional to the
	 * distance of its centre from a central point. Each circle is tangent to six
	 * others that surround it by a ring of tangent circles
	 * 
	 * @param centerX   x coordinate of the center of the spiral
	 * @param centerY   y coordinate of the center of the spiral
	 * @param p         at least 2
	 * @param q         at least p + 1
	 * @param maxRadius the maximum radius of the packing arrangement (the maximum
	 *                  distance a circle centroid can be from the center of the
	 *                  arrangement)
	 * @return A list of PVectors, each representing one circle in the spiral: (.x,
	 *         .y) represent the center point and .z represents radius.
	 */
	public static List<PVector> doyleSpiral(double centerX, double centerY, int p, int q, double maxRadius) {
		// A closed-form solution for a single p, q (now deprecated).
		/*
		 * double start = 0; // starting circle n double sr, ang, cr;
		 * 
		 * for (int i = 0; i < nCircles; i++) { sr = Math.exp((start + i) * 0.06101); //
		 * spiral radius ang = (start + i) * 0.656; // spiral angle cr = 0.3215 *
		 * Math.exp((start + i) * 0.06101); // circle radius circles.add(new
		 * PVector((float) (sr * Math.cos(ang) + centerX), (float) (sr * Math.sin(ang) +
		 * centerY), (float) cr)); }
		 */
		final DoyleSpiral doyleSpiral = new DoyleSpiral(p, q, maxRadius);
		doyleSpiral.getCircles().forEach(c -> c.add((float) centerX, (float) centerY));
		return new ArrayList<>(doyleSpiral.getCircles());
	}

	public static PShape islamicTiling(double width, double height, double w, double h) {
		// adapted from https://openprocessing.org/sketch/320133
		final double[] vector = { -w, 0, w, -h, w, 0, -w, h };
		final ArrayList<PVector> segments = new ArrayList<>();
		for (int x = 0; x < width; x += w * 2) {
			for (int y = 0; y < height; y += h * 2) {
				for (int i = 0; i <= vector.length; i++) {
					segments.add(
							new PVector((float) (vector[i % vector.length] + x + w), (float) (vector[(i + 6) % vector.length] + y + h)));
					segments.add(new PVector((float) (vector[(i + 1) % vector.length] + x + w),
							(float) (vector[(i + 1 + 6) % vector.length] + y + h)));
				}
			}
		}
		return PGS_Processing.polygonizeLines(segments);
	}

	private static void divideRect(PVector p1, PVector p2, PVector p3, PVector p4, int n, PShape parent, SplittableRandom r) {
		n--;
		if (n == 0) {
			PShape division = new PShape(PShape.PATH);
			division.setStrokeJoin(PConstants.MITER);
			division.beginShape();
			division.vertex(p1.x, p1.y);
			division.vertex(p2.x, p2.y);
			division.vertex(p3.x, p3.y);
			division.vertex(p4.x, p4.y);
			division.endShape(PConstants.CLOSE);
			parent.addChild(division);
		} else if (n > 0) {
			float w = PVector.dist(p1, p2) + PVector.dist(p3, p4);
			float h = PVector.dist(p1, p4) + PVector.dist(p2, p3);
			int t = 3;
			float r1 = (1f / t) * r.nextInt(1, t);
			float r2 = (1f / t) * r.nextInt(1, t);
			if (w < h) {
				PVector v1 = PVector.lerp(p1, p4, r1);
				PVector v2 = PVector.lerp(p2, p3, r2);
				divideRect(p1, p2, v2, v1, n, parent, r);
				divideRect(v1, v2, p3, p4, n, parent, r);
			} else {
				PVector v1 = PVector.lerp(p1, p2, r1);
				PVector v2 = PVector.lerp(p3, p4, r2);
				divideRect(p1, v1, v2, p4, n, parent, r);
				divideRect(v1, p2, p3, v2, n, parent, r);
			}
		}
	}

}
