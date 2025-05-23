package micycle.pgs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.SplittableRandom;

import micycle.pgs.color.Colors;
import micycle.pgs.commons.DoyleSpiral;
import micycle.pgs.commons.HatchTiling;
import micycle.pgs.commons.PenroseTiling;
import micycle.pgs.commons.RectangularSubdivision;
import micycle.pgs.commons.SquareTriangleTiling;
import micycle.pgs.commons.TriangleSubdivision;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Tiling, tessellation and subdivision of the plane using periodic or
 * non-periodic geometric shapes.
 * <p>
 * A tiling is created when a collection of tiles fills a plane such that no
 * gaps occur between the tiles and no two tiles overlap each other.
 * 
 * @author Michael Carleton
 * @since 1.2.0
 */
public final class PGS_Tiling {

	private static final double ROOT3 = Math.sqrt(3);

	private PGS_Tiling() {
	}

	/**
	 * Recursively and randomly subdivides the given/bounded plane into rectangles.
	 * 
	 * @param width    width of the quad subdivision plane
	 * @param height   height of the quad subdivision plane
	 * @param maxDepth maximum number of subdivisions (recursion depth)
	 * @return a GROUP PShape, where each child shape is a face of the subdivision
	 * @see #rectSubdivision(double, double, int, long) seeded rectSubdivsion()
	 */
	public static PShape rectSubdivision(double width, double height, int maxDepth) {
		return rectSubdivision(width, height, maxDepth, System.nanoTime());
	}

	/**
	 * Recursively and randomly subdivides the given/bounded plane into rectangles.
	 * 
	 * @param width    width of the quad subdivision plane
	 * @param height   height of the quad subdivision plane
	 * @param maxDepth maximum number of subdivisions (recursion depth)
	 * @param seed     the random seed
	 * @return a GROUP PShape, where each child shape is a face of the subdivision
	 * @see #rectSubdivision(double, double, int) non-seeded rectSubdivsion()
	 */
	public static PShape rectSubdivision(double width, double height, int maxDepth, long seed) {
		maxDepth++; // so that given depth==0 returns non-divided square
		final RectangularSubdivision rectangularSubdivision = new RectangularSubdivision(width, height, maxDepth, seed);
		return rectangularSubdivision.divide();
	}

	/**
	 * Recursively and randomly subdivides the given/bounded plane into triangles.
	 * 
	 * @param width    width of the subdivision plane
	 * @param height   height of the subdivision plane
	 * @param maxDepth maximum number of subdivisions (recursion depth)
	 * @return a GROUP PShape, where each child shape is a face of the subdivision
	 * @see #triangleSubdivision(double, double, int, long) seeded
	 *      triangleSubdivsion()
	 */
	public static PShape triangleSubdivision(double width, double height, int maxDepth) {
		return triangleSubdivision(width, height, maxDepth, System.nanoTime());
	}

	/**
	 * Recursively and randomly subdivides the given/bounded plane into triangles.
	 * 
	 * @param width    width of the subdivision plane
	 * @param height   height of the subdivision plane
	 * @param maxDepth maximum number of subdivisions (recursion depth)
	 * @param seed     the random seed
	 * @return a GROUP PShape, where each child shape is a face of the subdivision
	 * @see PGS_Tiling#triangleSubdivision(double, double, int) non-seeded
	 *      triangleSubdivision()
	 */
	public static PShape triangleSubdivision(double width, double height, int maxDepth, long seed) {
		maxDepth++; // so that given depth==0 returns non-divided triangle
		final TriangleSubdivision subdivision = new TriangleSubdivision(width, height, maxDepth, seed);
		return subdivision.divide();
	}

	/**
	 * Recursively and randomly subdivides the given/bounded plane into convex quad
	 * polygons.
	 * 
	 * @param width  width of the plane that is subdivided
	 * @param height height of the plane that is subdivided
	 * @param depth  number of subdivisions (recursion depth)
	 * @return a GROUP PShape, where each child shape is a face of the subdivision
	 * @see #quadSubdivision(double, double, int, long) seeded quadSubdivision()
	 */
	public static PShape quadSubdivision(double width, double height, int depth) {
		return quadSubdivision(width, height, depth, System.nanoTime());
	}

	/**
	 * Recursively and randomly subdivides the given/bounded plane into convex quad
	 * polygons.
	 * 
	 * @param width  width of the quad subdivision plane
	 * @param height height of the quad subdivision plane
	 * @param depth  number of subdivisions (recursion depth)
	 * @param seed   the random seed
	 * @return a GROUP PShape, where each child shape is a face of the subdivision
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
		PGS_Conversion.setAllFillColor(divisions, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(divisions, Colors.PINK, 2);
		return divisions;
	}

	/**
	 * Randomly subdivides the plane into equal-width strips having varying lengths.
	 * 
	 * @param width      width of the subdivision plane
	 * @param height     height of the subdivision plane
	 * @param gridCountX horizontal grid count
	 * @param gridCountY vertical grid count
	 * @param seed       the random seed
	 * @return a GROUP PShape, where each child shape is a face of the subdivision
	 * @since 1.3.0
	 */
	public static PShape hatchSubdivision(double width, double height, int gridCountX, int gridCountY, long seed) {
		final HatchTiling ht = new HatchTiling((int) width, (int) height, gridCountX, gridCountY);
		PShape tiling = ht.getTiling(seed);
		PGS_Conversion.setAllStrokeColor(tiling, Colors.PINK, 4);
		return tiling;
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

	/**
	 * Generates a hexagonal tiling of the plane.
	 * 
	 * @param width      width of the tiling plane
	 * @param height     height of the tiling plane
	 * @param sideLength side length of each hexagon
	 * @param flat       determines the orientation of the hexagons -- whether the
	 *                   top is flat, or pointy
	 * @return a GROUP PShape, where each child shape is a hexagon of the tiling
	 */
	public static PShape hexTiling(double width, double height, double sideLength, boolean flat) {
		final double span = sideLength * ROOT3;
		final PShape tiling = new PShape(PConstants.GROUP);
		double x = 0;
		double y = 0;

		int row = 0;
		while (y < height) {
			int column = 0;

			if (row % 2 == 0 || flat) {
				x = 0;
			} else {
				x = span / 2;
			}

			while (x < width) {
				if (flat) {
					y = span * row;
					if (column % 2 != 0) {
						y += span / 2;
					}
				}
				tiling.addChild(hexagon(x, y, sideLength, flat));
				x += flat ? sideLength * 1.5 : span;
				column++;
			}
			y += flat ? span / 2 : sideLength * 1.5;
			row++;
		}

		return tiling;
	}

	/**
	 * Generates an "islamic-style" (Girih) tiling of the plane.
	 * 
	 * @param width  width of the tiling plane
	 * @param height height of the tiling plane
	 * @param w
	 * @param h
	 * @return a GROUP PShape, where each child shape is a tile of the tiling
	 */
	public static PShape islamicTiling(double width, double height, double w, double h) {
		// adapted from https://openprocessing.org/sketch/320133
		final double[] vector = { -w, 0, w, -h, w, 0, -w, h };
		final ArrayList<PVector> segments = new ArrayList<>();
		for (int x = 0; x < width; x += w * 2) {
			for (int y = 0; y < height; y += h * 2) {
				for (int i = 0; i <= vector.length; i++) {
					segments.add(new PVector((float) (vector[i % vector.length] + x + w), (float) (vector[(i + 6) % vector.length] + y + h)));
					segments.add(new PVector((float) (vector[(i + 1) % vector.length] + x + w), (float) (vector[(i + 1 + 6) % vector.length] + y + h)));
				}
			}
		}
		return PGS_Processing.polygonizeLines(segments);
	}

	/**
	 * Generates a Penrose Tiling (consisting of rhombi).
	 * 
	 * @param centerX x coordinate of the center/origin of the tiling
	 * @param centerY y coordinate of the center/origin of the tiling
	 * @param radius  maximum radius of the tiling (measured from the center)
	 * @param steps   number of tiling subdivisions
	 * @return a GROUP PShape, where each child shape is a face of the tiling
	 */
	public static PShape penroseTiling(double centerX, double centerY, final double radius, final int steps) {
		final PenroseTiling pr = new PenroseTiling(centerX, centerY, radius, steps);
		return PGS.polygonizeEdges(pr.getEdges());
	}

	/**
	 * Generates a non-periodic tiling, comprising squares and equilateral
	 * triangles.
	 * 
	 * @param width    width of the tiling plane
	 * @param height   height of the tiling plane
	 * @param tileSize diameter of each tile
	 * @return a GROUP PShape, where each child shape is a tile of the tiling
	 * @since 1.3.0
	 */
	public static PShape squareTriangleTiling(double width, double height, double tileSize) {
		return squareTriangleTiling(width, height, tileSize, System.nanoTime());
	}

	/**
	 * Generates a non-periodic tiling, comprising squares and equilateral
	 * triangles, having a given seed.
	 * 
	 * @param width    width of the tiling plane
	 * @param height   height of the tiling plane
	 * @param tileSize diameter of each tile
	 * @param seed     the random seed
	 * @return a GROUP PShape, where each child shape is a tile of the tiling
	 * @since 1.3.0
	 */
	public static PShape squareTriangleTiling(double width, double height, double tileSize, long seed) {
		final SquareTriangleTiling stt = new SquareTriangleTiling(width, height, tileSize);
		return stt.getTiling(seed);
	}

	/**
	 * Generates a geometric arrangement composed of annular-sector bricks arranged
	 * in concentric circular rings. Rings progressively expand from the inside out
	 * based on the growth rates provided. Brick sizes (arc length) adapt radially
	 * according to growth factors.
	 *
	 * @param nRings      Number of annular rings to generate.
	 * @param cx          The x-coordinate of the center of the generated pattern.
	 * @param cy          The y-coordinate of the center of the generated pattern.
	 * @param innerRadius The radius of the innermost ring in pixels.
	 * @param ringGrowth  The growth factor of each successive ring radius; values
	 *                    greater than 1.0 cause rings to expand radially outward.
	 * @param segGrowth   The growth factor controlling segment (brick) arc length
	 *                    adjustment along each ring; values greater than 1.0
	 *                    increase brick length progressively outward.
	 * @return A single flattened {@code PShape} consisting of annular-sector bricks
	 *         forming concentric rings around the specified center point
	 *         ({@code cx}, {@code cy}).
	 * @since 2.1
	 */
	public static PShape annularBricks(int nRings, double cx, double cy, double innerRadius, double ringGrowth, double segGrowth) {
		segGrowth = Math.max(segGrowth, 1e-6);
		ringGrowth = Math.max(ringGrowth, 1e-6);
		List<PShape> bricks = new ArrayList<>();
		final double GOLDEN_ANGLE = 2.0 * Math.PI * 0.38196601125;
		double segSize = innerRadius;
		final double maxDev = 0.25; // max chord sagitta in pixels
	
		for (int i = 0; i < nRings; i++) {
			// ring geometry in double
			double r = innerRadius * Math.pow(ringGrowth, i);
			double dr = (r * (ringGrowth - 1.0)) * 0.90;
			double perim = 2.0 * Math.PI * r;
			int count = (int) Math.max(1.0, Math.floor(perim / segSize));
			double dang = 2.0 * Math.PI / count;
			double gapAng = (dr / 6.0) / perim * 2.0 * Math.PI;
			double offset = (i * GOLDEN_ANGLE) % (2.0 * Math.PI);
	
			// tile this ring
			for (int k = 0; k < count; k++) {
				double a0 = offset + k * dang + 0.5 * gapAng;
				double a1 = offset + (k + 1) * dang - 0.5 * gapAng;
	
				List<PVector> outer = PGS_Construction.arcPoints(cx, cy, r + dr, a0, a1, maxDev);
	
				// generate inner in the same direction and reverse it
				List<PVector> inner = PGS_Construction.arcPoints(cx, cy, r, a0, a1, maxDev);
				Collections.reverse(inner);
	
				outer.addAll(inner);
				outer.add(outer.get(0)); // close the loop
	
				PShape brick = PGS_Conversion.fromPVector(outer);
				brick.setStroke(false);
				bricks.add(brick);
			}
	
			segSize *= segGrowth;
			segSize = Math.max(segSize, 1);
		}
	
		return PGS_Conversion.flatten(bricks);
	}

	/**
	 * Generates a hexagon shape.
	 * 
	 * @param x          x-position of hexagon envelope's top left corner
	 * @param y          y-position of hexagon envelope's top left corner
	 * @param sideLength hexagon side length
	 * @param flat       the orientation of the hexagon (whether the top is flat, or
	 *                   pointy)
	 * @return a PATH PShape
	 */
	private static PShape hexagon(double x, double y, double sideLength, boolean flat) {
		final double span = sideLength * ROOT3;
		final PShape hexagon = new PShape(PShape.PATH);
		hexagon.setStroke(true);
		hexagon.setStroke(Colors.PINK);
		hexagon.setStrokeWeight(2);
		hexagon.setFill(true);
		hexagon.setFill(Colors.WHITE);

		hexagon.beginShape();
		if (flat) {
			hexagon.vertex((float) (x + 0.5 * sideLength), (float) (y + span));
			hexagon.vertex((float) (x + 1.5 * sideLength), (float) (y + span));
			hexagon.vertex((float) (x + 2.0 * sideLength), (float) (y + span / 2.0));
			hexagon.vertex((float) (x + 1.5 * sideLength), (float) y);
			hexagon.vertex((float) (x + 0.5 * sideLength), (float) y);
			hexagon.vertex((float) x, (float) (y + span / 2.0));
		} else {
			hexagon.vertex((float) (x + 0.5 * span), (float) (y + 2.0 * sideLength));
			hexagon.vertex((float) (x + span), (float) (y + 1.5 * sideLength));
			hexagon.vertex((float) (x + span), (float) (y + 0.5 * sideLength));
			hexagon.vertex((float) (x + 0.5 * span), (float) y);
			hexagon.vertex((float) x, (float) (y + 0.5 * sideLength));
			hexagon.vertex((float) x, (float) (y + 1.5 * sideLength));
		}
		hexagon.endShape(PConstants.CLOSE);

		return hexagon;
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
