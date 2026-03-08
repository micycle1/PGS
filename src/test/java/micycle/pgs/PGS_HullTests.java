package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

class PGS_HullTests {

	private static final double EPSILON = 1E-5;

	static PShape square;

	@BeforeAll
	static void initShapes() {
		square = new PShape(PShape.GEOMETRY); // 10x10 square at origin
		square.beginShape();
		square.vertex(0, 0);
		square.vertex(10, 0);
		square.vertex(10, 10);
		square.vertex(0, 10);
		square.endShape(PConstants.CLOSE);
	}

	@Test
	void testBoundingBoxSquare() {
		// A square's bounding box is itself
		PShape bb = PGS_Hull.boundingBox(square);
		assertEquals(100, PGS_ShapePredicates.area(bb), EPSILON);
		assertEquals(10, PGS_ShapePredicates.width(bb), EPSILON);
		assertEquals(10, PGS_ShapePredicates.height(bb), EPSILON);
	}

	@Test
	void testBoundingBoxTriangle() {
		// Right triangle: (0,0), (10,0), (0,10). Bounding box should be 10x10.
		PShape tri = new PShape(PShape.GEOMETRY);
		tri.beginShape();
		tri.vertex(0, 0);
		tri.vertex(10, 0);
		tri.vertex(0, 10);
		tri.endShape(PConstants.CLOSE);

		PShape bb = PGS_Hull.boundingBox(tri);
		assertEquals(100, PGS_ShapePredicates.area(bb), EPSILON);
	}

	@Test
	void testBoundingBoxWithOut() {
		double[] out = new double[4];
		PGS_Hull.boundingBox(square, out);
		assertEquals(0, out[0], EPSILON); // minX
		assertEquals(0, out[1], EPSILON); // minY
		assertEquals(10, out[2], EPSILON); // maxX
		assertEquals(10, out[3], EPSILON); // maxY
	}

	@Test
	void testConvexHullAlreadyConvex() {
		// Convex hull of a square's vertices is the square itself
		List<PVector> pts = Arrays.asList(new PVector(0, 0), new PVector(10, 0), new PVector(10, 10), new PVector(0, 10));
		PShape hull = PGS_Hull.convexHull(pts);
		assertEquals(100, PGS_ShapePredicates.area(hull), EPSILON);
	}

	@Test
	void testConvexHullWithInteriorPoint() {
		// Add an interior point — hull should still be the square
		List<PVector> pts = Arrays.asList(new PVector(0, 0), new PVector(10, 0), new PVector(10, 10), new PVector(0, 10), new PVector(5, 5)); // interior
		PShape hull = PGS_Hull.convexHull(pts);
		assertEquals(100, PGS_ShapePredicates.area(hull), EPSILON);
		assertEquals(4, PGS_ShapePredicates.vertexCount(hull));
	}

	@Test
	void testConvexHullCollinearPoints() {
		// All points on a line — hull should be a degenerate shape (line) with area 0
		List<PVector> pts = Arrays.asList(new PVector(0, 0), new PVector(5, 0), new PVector(10, 0));
		PShape hull = PGS_Hull.convexHull(pts);
		assertEquals(0, PGS_ShapePredicates.area(hull), EPSILON);
	}

	@Test
	void testConvexHullShapeConvex() {
		// Convex hull of a square should have area 100
		PShape hull = PGS_Hull.convexHull(square);
		assertEquals(100, PGS_ShapePredicates.area(hull), EPSILON);
	}

	@Test
	void testConvexHullShapeConcave() {
		// L-shaped polygon (concave) — its hull should be larger
		PShape lshape = new PShape(PShape.GEOMETRY);
		lshape.beginShape();
		lshape.vertex(0, 0);
		lshape.vertex(10, 0);
		lshape.vertex(10, 5);
		lshape.vertex(5, 5);
		lshape.vertex(5, 10);
		lshape.vertex(0, 10);
		lshape.endShape(PConstants.CLOSE);

		double lArea = PGS_ShapePredicates.area(lshape); // 75
		PShape hull = PGS_Hull.convexHull(lshape);
		double hullArea = PGS_ShapePredicates.area(hull); // 87.5
		assertTrue(hullArea > lArea, "Hull area should exceed concave shape area");
		assertEquals(87.5, hullArea, EPSILON);
	}

	@Test
	void testSnapHullConvexityOne() {
		// convexity=1 should produce the convex hull
		PShape hull = PGS_Hull.snapHull(square, 1.0);
		assertEquals(100, PGS_ShapePredicates.area(hull), EPSILON);
	}

	@Test
	void testSnapHullConvexityZero() {
		// convexity=0 should reproduce the original shape (area matches exactly)
		PShape hull = PGS_Hull.snapHull(square, 0.0);
		assertEquals(100, PGS_ShapePredicates.area(hull), EPSILON);
	}

}
