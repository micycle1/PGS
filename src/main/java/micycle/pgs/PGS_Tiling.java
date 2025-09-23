package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.SplittableRandom;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.overlayng.RingClipper;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;

import micycle.pgs.color.Colors;
import micycle.pgs.commons.DoyleSpiral;
import micycle.pgs.commons.HatchTiling;
import micycle.pgs.commons.PEdge;
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
 * </p>
 * <p>
 * <b>Naming convention in this class:</b><br>
 * Methods ending with "subdivision" recursively break the plane down into
 * smaller and smaller shapes, usually of the same geometric type at each step
 * (for example, rectangles are split into smaller rectangles).<br>
 * Methods ending with "division" split the plane all at once, in a single step,
 * often by making several cuts or slices, with no recursion.
 * </p>
 *
 * @author Michael Carleton
 * @since 1.2.0
 */
public final class PGS_Tiling {

	private static final double ROOT3 = Math.sqrt(3); // for hex

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
	public static PShape rectSubdivision(final double width, final double height, final int maxDepth) {
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
	public static PShape rectSubdivision(final double width, final double height, int maxDepth, final long seed) {
		maxDepth++; // so that given depth==0 returns non-divided square
		final RectangularSubdivision rectangularSubdivision = new RectangularSubdivision(width, height, maxDepth, seed);
		var division = rectangularSubdivision.divide();
		division = PGS_Conversion.setAllFillColor(division, Colors.WHITE);
		division = PGS_Conversion.setAllStrokeColor(division, Colors.PINK, 2);
		return division;
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
	public static PShape triangleSubdivision(final double width, final double height, final int maxDepth) {
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
	public static PShape triangleSubdivision(final double width, final double height, int maxDepth, final long seed) {
		maxDepth++; // so that given depth==0 returns non-divided triangle
		final TriangleSubdivision subdivision = new TriangleSubdivision(width, height, maxDepth, seed);
		var division = subdivision.divide();
		division = PGS_Conversion.setAllFillColor(division, Colors.WHITE);
		division = PGS_Conversion.setAllStrokeColor(division, Colors.PINK, 2);
		return division;
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
	public static PShape quadSubdivision(final double width, final double height, final int depth) {
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
	public static PShape quadSubdivision(final double width, final double height, final int depth, final long seed) {
		// https://openprocessing.org/sketch/1045334
		final float w = (float) width;
		final float h = (float) height;
		final float off = 0;
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
	public static PShape hatchSubdivision(final double width, final double height, final int gridCountX, final int gridCountY, final long seed) {
		final HatchTiling ht = new HatchTiling((int) width, (int) height, gridCountX, gridCountY);
		final PShape tiling = ht.getTiling(seed);
		PGS_Conversion.setAllStrokeColor(tiling, Colors.PINK, 4);
		return tiling;
	}

	/**
	 * Divides the plane into randomly “sliced” polygonal regions.
	 * <p>
	 * {@code slices} random cuts are generated across the plane (dimensions w×h, at
	 * (0,0)). Each cut connects a random point on one side of the plane to a random
	 * point on another side. If {@code forceOpposite} is true, each cut always
	 * connects opposite sides; otherwise the two sides are chosen at random (but
	 * never the same side).
	 * </p>
	 * <p>
	 * <strong>In practice:</strong>
	 * <ul>
	 * <li>forceOpposite == <code>true</code> → mostly long, quadrilateral strips
	 * that span the full width or height of the rectangle.</li>
	 * <li>forceOpposite == <code>false</code> → a richer variety of cell shapes
	 * (triangles, trapezoids, L-shapes, etc.) that may not stretch all the way
	 * across.</li>
	 * </ul>
	 * </p>
	 *
	 * @param width         the width of the plane
	 * @param height        the height of the plane
	 * @param slices        the number of random interior cuts to perform
	 * @param forceOpposite if true, each cut connects opposite sides of the
	 *                      rectangle; if false, cuts connect any two distinct sides
	 * @param seed          the random seed for reproducible slice placement
	 * @return a GROUP PShape containing the subdivided polygonal regions
	 * @since 2.1
	 */
	public static PShape sliceDivision(final double width, final double height, final int slices, final boolean forceOpposite, final long seed) {
		final List<PEdge> cuts = new ArrayList<>(slices + 4);
		final double x = 0, y = 0;
		final PVector A = new PVector((float) x, (float) y);
		final PVector B = new PVector((float) (x + width), (float) y);
		final PVector C = new PVector((float) (x + width), (float) (y + height));
		final PVector D = new PVector((float) x, (float) (y + height));

		cuts.add(new PEdge(A, B));
		cuts.add(new PEdge(B, C));
		cuts.add(new PEdge(C, D));
		cuts.add(new PEdge(D, A));

		final SplittableRandom r = new SplittableRandom(seed);
		for (int i = 0; i < slices; i++) {
			final int s1 = r.nextInt(4);
			int s2;
			if (forceOpposite) {
				// always pick the side directly opposite s1
				s2 = (s1 + 2) % 4;
			} else {
				// pick any side except s1
				do {
					s2 = r.nextInt(4);
				} while (s2 == s1);
			}

			final var p1 = PGS.toPVector(pointOnSide(s1, r, x, y, width, height));
			final var p2 = PGS.toPVector(pointOnSide(s2, r, x, y, width, height));
			final PEdge cut = new PEdge(p1, p2);
			cuts.add(cut);
		}

		return PGS.polygonizeEdges(cuts);
	}

	/**
	 * Creates a cellular partition of the plane using arcs formed by circles seeded
	 * along its boundary. Each circle’s radius is chosen large enough to guarantee
	 * it intersects at least two distinct sides of the plane, forming an arc cut.
	 * <p>
	 * <strong>In practice:</strong>
	 * <ul>
	 * <li>You get a “cellular” subdivision where each circle carves out roughly
	 * circular holes or bulges against the rectangle boundary.</li>
	 * <li>Because every radius ≥ the minimum distance to a second side, each circle
	 * always spans at least two sides (forming an arc).</li>
	 * </ul>
	 * </p>
	 * 
	 * @param width        the width of the plane
	 * @param height       the height of the plane
	 * @param arcs         the number of circles to seed along the boundary
	 * @param circlePoints the number of linear segments used to approximate each
	 *                     circle
	 * @param seed         the random seed for reproducible circle placement & radii
	 * @return a GROUP PShape containing the polygonal faces formed by the plane
	 *         plus the seeded circles
	 * @since 2.1
	 */
	public static PShape arcDivision(final double width, final double height, final int arcs, final long seed) {
		final SplittableRandom rnd = new SplittableRandom(seed);
		final double maxDim = Math.max(width, height);
		final double x = 0, y = 0;
		final var e = new Envelope(x, x + width, y, y + height);
		final RingClipper clipper = new RingClipper(e);

		final List<Geometry> polys = new ArrayList<>(arcs);
		polys.add(((Polygon) GEOM_FACTORY.toGeometry(e)).getExteriorRing());

		// 2) Seed 'seeds' circles along random sides
		for (int i = 0; i < arcs; i++) {
			// pick a side 0=top,1=right,2=bottom,3=left
			final int side = rnd.nextInt(4);
			// pick a random point along that side
			final Coordinate center = pointOnSide(side, rnd, x, y, width, height);

			// compute min distance from center to any *other* side
			double minReq = Double.POSITIVE_INFINITY;
			// top edge y=y, bottom y=y+h, left x=x, right x=x+w
			final double cx = center.x, cy = center.y;
			if (side != 0) {
				minReq = Math.min(minReq, Math.abs(cy - y));
			}
			if (side != 2) {
				minReq = Math.min(minReq, Math.abs(cy - (y + height)));
			}
			if (side != 3) {
				minReq = Math.min(minReq, Math.abs(cx - x));
			}
			if (side != 1) {
				minReq = Math.min(minReq, Math.abs(cx - (x + width)));
			}

			final double radius = minReq + rnd.nextDouble() * (maxDim - minReq);

			final var circle = PGS_Construction.createCirclePoly(cx, cy, radius);
			final var clip = clipper.clip(circle.getCoordinates());
			polys.add(GEOM_FACTORY.createLineString(clip));

		}

		final Polygonizer polygonizer = new Polygonizer();
		polygonizer.setCheckRingsValid(false);
		polygonizer.add(UnaryUnionOp.union(polys));
		return PGS_Conversion.toPShape(polygonizer.getPolygons());
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
	public static List<PVector> doyleSpiral(final double centerX, final double centerY, final int p, final int q, final double maxRadius) {
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
	public static PShape hexTiling(final double width, final double height, final double sideLength, final boolean flat) {
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
	public static PShape islamicTiling(final double width, final double height, final double w, final double h) {
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
	public static PShape penroseTiling(final double centerX, final double centerY, final double radius, final int steps) {
		final PenroseTiling pr = new PenroseTiling(centerX, centerY, radius, steps);
		return PGS.polygonizeNodedEdges(pr.getEdges());
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
	public static PShape squareTriangleTiling(final double width, final double height, final double tileSize) {
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
	public static PShape squareTriangleTiling(final double width, final double height, final double tileSize, final long seed) {
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
	public static PShape annularBricks(final int nRings, final double cx, final double cy, final double innerRadius, double ringGrowth, double segGrowth) {
		segGrowth = Math.max(segGrowth, 1e-6);
		ringGrowth = Math.max(ringGrowth, 1e-6);
		final List<PShape> bricks = new ArrayList<>();
		final double GOLDEN_ANGLE = 2.0 * Math.PI * 0.38196601125;
		double segSize = innerRadius;
		final double maxDev = 0.25; // max chord sagitta in pixels

		for (int i = 0; i < nRings; i++) {
			// ring geometry in double
			final double r = innerRadius * Math.pow(ringGrowth, i);
			final double dr = (r * (ringGrowth - 1.0)) * 0.90;
			final double perim = 2.0 * Math.PI * r;
			final int count = (int) Math.max(1.0, Math.floor(perim / segSize));
			final double dang = 2.0 * Math.PI / count;
			final double gapAng = (dr / 6.0) / perim * 2.0 * Math.PI;
			final double offset = (i * GOLDEN_ANGLE) % (2.0 * Math.PI);

			// tile this ring
			for (int k = 0; k < count; k++) {
				final double a0 = offset + k * dang + 0.5 * gapAng;
				final double a1 = offset + (k + 1) * dang - 0.5 * gapAng;

				final List<PVector> outer = PGS_Construction.arcPoints(cx, cy, r + dr, a0, a1, maxDev);

				// generate inner in the same direction and reverse it
				final List<PVector> inner = PGS_Construction.arcPoints(cx, cy, r, a0, a1, maxDev);
				Collections.reverse(inner);

				outer.addAll(inner);
				outer.add(outer.get(0)); // close the loop

				final PShape brick = PGS_Conversion.fromPVector(outer);
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
	private static PShape hexagon(final double x, final double y, final double sideLength, final boolean flat) {
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

	private static void divideRect(final PVector p1, final PVector p2, final PVector p3, final PVector p4, int n, final PShape parent,
			final SplittableRandom r) {
		n--;
		if (n == 0) {
			final PShape division = new PShape(PShape.PATH);
			division.setStrokeJoin(PConstants.MITER);
			division.beginShape();
			division.vertex(p1.x, p1.y);
			division.vertex(p2.x, p2.y);
			division.vertex(p3.x, p3.y);
			division.vertex(p4.x, p4.y);
			division.endShape(PConstants.CLOSE);
			parent.addChild(division);
		} else if (n > 0) {
			final float w = PVector.dist(p1, p2) + PVector.dist(p3, p4);
			final float h = PVector.dist(p1, p4) + PVector.dist(p2, p3);
			final int t = 3;
			final float r1 = (1f / t) * r.nextInt(1, t);
			final float r2 = (1f / t) * r.nextInt(1, t);
			if (w < h) {
				final PVector v1 = PVector.lerp(p1, p4, r1);
				final PVector v2 = PVector.lerp(p2, p3, r2);
				divideRect(p1, p2, v2, v1, n, parent, r);
				divideRect(v1, v2, p3, p4, n, parent, r);
			} else {
				final PVector v1 = PVector.lerp(p1, p2, r1);
				final PVector v2 = PVector.lerp(p3, p4, r2);
				divideRect(p1, v1, v2, p4, n, parent, r);
				divideRect(v1, p2, p3, v2, n, parent, r);
			}
		}
	}

	/**
	 * Returns a random point on one of the four sides of the rectangle. side =
	 * 0→top, 1→right, 2→bottom, 3→left.
	 */
	private static Coordinate pointOnSide(final int side, final SplittableRandom r, final double x, final double y, final double w, final double h) {
		double px, py;
		switch (side) {
			case 0 : // top
				px = x + r.nextDouble() * w;
				py = y;
				break;
			case 1 : // right
				px = x + w;
				py = y + r.nextDouble() * h;
				break;
			case 2 : // bottom
				px = x + r.nextDouble() * w;
				py = y + h;
				break;
			default : // left
				px = x;
				py = y + r.nextDouble() * h;
				break;
		}
		return new Coordinate(px, py);
	}

}
