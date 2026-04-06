package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;

class PGS_ShapeBooleanTests {

	@Test
	void testPolygonPolygonUnion() {
		PShape a = createSquare(0, 0, 10);
		PShape b = createSquare(5, 0, 10);
		PShape union = PGS_ShapeBoolean.union(a, b);

		// Expected area: (10*10) + (10*10) - (5*10) = 150
		assertEquals(150.0, PGS_ShapePredicates.area(union), 1e-6);
		assertEquals(1, union.getChildCount() == 0 ? 1 : union.getChildCount());
	}

	@Test
	void testSelfUnion() {
		PShape a = createSquare(0, 0, 10);
		PShape b = createSquare(5, 0, 10);
		var s = PGS_Conversion.flatten(a, b);
		PShape union = PGS_ShapeBoolean.union(s);

		// Expected area: (10*10) + (10*10) - (5*10) = 150
		assertEquals(150.0, PGS_ShapePredicates.area(union), 1e-6);
		assertEquals(1, union.getChildCount() == 0 ? 1 : union.getChildCount());
	}

	@Test
	void testPolygonPolygonIntersection() {
		PShape a = createSquare(0, 0, 10);
		PShape b = createSquare(5, 5, 10);
		PShape intersection = PGS_ShapeBoolean.intersect(a, b);

		// Expected area: 5*5 = 25
		assertEquals(25.0, PGS_ShapePredicates.area(intersection), 1e-6);
	}

	@Test
	void testPolygonPolygonSubtraction() {
		PShape a = createSquare(0, 0, 10);
		PShape b = createSquare(5, 0, 10);
		PShape difference = PGS_ShapeBoolean.subtract(a, b);

		// Expected area: 100 - 50 = 50
		assertEquals(50.0, PGS_ShapePredicates.area(difference), 1e-6);
	}

	@Test
	void testPolygonPolygonSymDifference() {
		PShape a = createSquare(0, 0, 10);
		PShape b = createSquare(5, 0, 10);
		PShape symDiff = PGS_ShapeBoolean.symDifference(a, b);

		// (A union B) - (A intersect B) = 150 - 50 = 100
		assertEquals(100.0, PGS_ShapePredicates.area(symDiff), 1e-6);
	}

	@Test
	void testLineLineIntersection() {
		PShape line1 = createLine(0, 5, 10, 5);
		PShape line2 = createLine(5, 0, 5, 10);
		PShape intersection = PGS_ShapeBoolean.intersect(line1, line2);

		// Intersection of two lines is a point
		assertEquals(1, intersection.getVertexCount());
		assertEquals(5, intersection.getVertexX(0));
		assertEquals(5, intersection.getVertexY(0));
	}

	@Test
	void testLineLineSubtraction() {
		PShape line1 = createLine(0, 5, 10, 5);
		PShape line2 = createLine(5, 5, 15, 5);
		PShape diff = PGS_ShapeBoolean.subtract(line1, line2);

		// (0,5 -> 10,5) minus (5,5 -> 15,5) should be (0,5 -> 5,5)
		assertEquals(2, diff.getVertexCount());
		assertEquals(0, diff.getVertexX(0));
		assertEquals(5, diff.getVertexX(1));
	}

	@Test
	void testPolygonLineIntersection() {
		PShape square = createSquare(0, 0, 10);
		PShape line = createLine(-5, 5, 15, 5);

		PShape intersection = PGS_ShapeBoolean.intersect(square, line);

		assertEquals(2, intersection.getVertexCount());
		assertEquals(0, intersection.getVertexX(0));
		assertEquals(10, intersection.getVertexX(1));
	}

	@Test
	void testPolygonLineDifference() {
		PShape square = createSquare(0, 0, 10);
		PShape line = createLine(-5, 5, 15, 5);

		PShape difference = PGS_ShapeBoolean.subtract(square, line);
		// Difference between poly and line is usually the poly itself (with points
		// added at intersection)
		// unless it's subtraction logic specifically handles it.
		assertTrue(PGS_ShapePredicates.area(difference) > 99.9);
	}

	@Test
	void testMultiUnion() {
		PShape s1 = createSquare(0, 0, 10);
		PShape s2 = createSquare(5, 0, 10);
		PShape s3 = createSquare(0, 5, 10);

		PShape union = PGS_ShapeBoolean.union(s1, s2, s3);
		// Area: 3 * 100 - (overlap 1&2 = 50) - (overlap 1&3 = 50) - (overlap 2&3 = 25)
		// + (overlap 1&2&3 = 25)
		// 300 - 50 - 50 - 25 + 25 = 200
		assertEquals(200.0, PGS_ShapePredicates.area(union), 1e-6);
	}

	@Test
	void testMultiIntersection() {
		PShape s1 = createSquare(0, 0, 10);
		PShape s2 = createSquare(5, 0, 10);
		PShape s3 = createSquare(5, 5, 10);

		PShape intersect = PGS_ShapeBoolean.intersect(s1, s2, s3);
		// Intersection: [5,10]x[5,10] -> Area 25
		assertEquals(25.0, PGS_ShapePredicates.area(intersect), 1e-6);
	}

	@Test
	void testOverlapRegions() {
		PShape s1 = createSquare(0, 0, 10);
		PShape s2 = createSquare(5, 0, 10);
		PShape s3 = createSquare(0, 5, 10);

		// Overlap regions are s1∩s2, s1∩s3, and s2∩s3∩s1
		PShape overlaps = PGS_ShapeBoolean.overlapRegions(List.of(s1, s2, s3), true);
		// s1∩s2 is [5,10]x[0,10], area 50
		// s1∩s3 is [0,10]x[5,10], area 50
		// They overlap in [5,10]x[5,10], area 25
		// Total area: 50 + 50 - 25 = 75
		assertEquals(75.0, PGS_ShapePredicates.area(overlaps), 1e-6);
	}

	@Test
	void testUnionMesh() {
		// Create two squares that share an edge
		PShape s1 = createSquare(0, 0, 10);
		PShape s2 = createSquare(10, 0, 10);

		PShape mesh = new PShape(PConstants.GROUP);
		mesh.addChild(s1);
		mesh.addChild(s2);

		PShape combined = PGS_ShapeBoolean.unionMesh(mesh);
		// Should be a 20x10 rectangle
		assertEquals(200.0, PGS_ShapePredicates.area(combined), 1e-6);
		// Should be a single polygon (not a group)
		assertTrue(combined.getChildCount() == 0);
	}

	private static PShape createSquare(float x, float y, float size) {
		PShape s = new PShape(PShape.GEOMETRY);
		s.beginShape();
		s.vertex(x, y);
		s.vertex(x + size, y);
		s.vertex(x + size, y + size);
		s.vertex(x, y + size);
		s.endShape(PConstants.CLOSE);
		return s;
	}

	private static PShape createLine(float x1, float y1, float x2, float y2) {
		PShape s = new PShape(PShape.PATH);
		s.beginShape();
		s.vertex(x1, y1);
		s.vertex(x2, y2);
		s.endShape();
		return s;
	}

}
