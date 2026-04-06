package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

class PGS_ShapePredicatesTests {

	private static final double EPSILON = 1E-5;

	static PShape square, triangle, rect;

	@BeforeAll
	static void initShapes() {
		square = new PShape(PShape.GEOMETRY); // 10x10 square
		square.beginShape();
		square.vertex(0, 0);
		square.vertex(10, 0);
		square.vertex(10, 10);
		square.vertex(0, 10);
		square.endShape(PConstants.CLOSE); // close affects rendering only -- does not append another vertex

		rect = new PShape(PShape.GEOMETRY); // 10x20 rect
		rect.beginShape();
		rect.vertex(0, 0);
		rect.vertex(10, 0);
		rect.vertex(10, 20);
		rect.vertex(0, 20);
		rect.endShape(PConstants.CLOSE); // close affects rendering only -- does not append another vertex

		float[] centroid = new float[] { 0, 0 };
		float side_length = 10;
		triangle = new PShape(PShape.GEOMETRY); // equilateral triangle
		triangle.beginShape();
		triangle.vertex(centroid[0], centroid[1] + ((float) Math.sqrt(3) / 3) * side_length); // Top vertex
		triangle.vertex(centroid[0] - (side_length / 2), centroid[1] - ((float) Math.sqrt(3) / 6) * side_length); // Bottom
																													// left
																													// vertex
		triangle.vertex(centroid[0] + (side_length / 2), centroid[1] - ((float) Math.sqrt(3) / 6) * side_length); // Bottom
																													// right
																													// vertex
		triangle.endShape(PConstants.CLOSE);
	}

	/** Helper: creates a closed square PShape at given origin with given size. */
	private static PShape makeSquare(float x, float y, float size) {
		PShape s = new PShape(PShape.GEOMETRY);
		s.beginShape();
		s.vertex(x, y);
		s.vertex(x + size, y);
		s.vertex(x + size, y + size);
		s.vertex(x, y + size);
		s.endShape(PConstants.CLOSE);
		return s;
	}

	@Test
	void testArea() {
		assertEquals(100, PGS_ShapePredicates.area(square));
		assertEquals(Math.sqrt(3) / 4 * 10 * 10, PGS_ShapePredicates.area(triangle), EPSILON);
	}

	@Test
	void testCentroid() {
		assertEquals(new PVector(5, 5), PGS_ShapePredicates.centroid(square));
		assertEquals(new PVector(0, 0), PGS_ShapePredicates.centroid(triangle));
	}

	@Test
	void testDiameter() {
		assertEquals(10 * Math.sqrt(2), PGS_ShapePredicates.diameter(square), EPSILON);
	}

	@Test
	void testMaximumInteriorAngle() {
		assertEquals(Math.PI / 2, PGS_ShapePredicates.maximumInteriorAngle(square));
		assertEquals(Math.PI / 3, PGS_ShapePredicates.maximumInteriorAngle(triangle), EPSILON);
	}

	@Test
	void testInteriorAnglesSquare() {
		var angles = PGS_ShapePredicates.interiorAngles(square);
		assertEquals(4, angles.size(), "Square should have 4 angles");

		double expectedAngleRadians = Math.PI / 2.0; // 90 degrees in radians
		double expectedAngleSumRadians = Math.PI * 2; // 360 degrees for a square
		double actualAngleSumRadians = 0;

		for (double angle : angles.values()) {
			assertEquals(expectedAngleRadians, angle, 1e-6, "Interior angle should be approximately 90 degrees");
			actualAngleSumRadians += angle;
		}
		assertEquals(expectedAngleSumRadians, actualAngleSumRadians, 1e-6, "Sum of square interior angles should be approximately 360 degrees");

	}

	@Test
	void testInteriorAnglesTriangle() {
		var angles = PGS_ShapePredicates.interiorAngles(triangle);

		assertEquals(3, angles.size(), "Triangle should have 3 angles");

		double expectedAngleSumRadians = Math.PI; // 180 degrees for a triangle
		double actualAngleSumRadians = 0;
		for (double angle : angles.values()) {
			actualAngleSumRadians += angle;
		}
		assertEquals(expectedAngleSumRadians, actualAngleSumRadians, 1e-6, "Sum of triangle interior angles should be approximately 180 degrees");
	}

	@Test
	void testHoles() {
		assertEquals(0, PGS_ShapePredicates.holes(square));
		PShape withHole = PGS_ShapeBoolean.subtract(square, PGS_Transformation.scale(square, 0.5));
		assertEquals(1, PGS_ShapePredicates.holes(withHole));
		PShape groupHoles = PGS_Conversion.flatten(withHole, withHole);
		assertEquals(2, PGS_ShapePredicates.holes(groupHoles));

		PShape coverage = PGS_Processing.split(withHole); // test coverage; no face has a hole, but together they form a
															// mesh with a hole
		assertEquals(1, PGS_ShapePredicates.holes(coverage));

		coverage.removeChild(0); // remove a mesh face; mesh no longer forms a hole
		assertEquals(0, PGS_ShapePredicates.holes(coverage));
	}

	@Test
	void testIsClockwise() {
		assertTrue(PGS_ShapePredicates.isClockwise(square));
		List<PVector> ccw = PGS_Conversion.toPVector(square);
		Collections.reverse(ccw);
		ccw.add(ccw.get(0)); // close
		assertFalse(PGS_ShapePredicates.isClockwise(PGS_Conversion.fromPVector(ccw)));

	}

	@Test
	void testElongation() {
		assertEquals(0, PGS_ShapePredicates.elongation(square));
		assertEquals(0.5, PGS_ShapePredicates.elongation(rect));
	}

	@Test
	void testMinimumInteriorAngle() {
		assertEquals(Math.PI / 2, PGS_ShapePredicates.minimumInteriorAngle(rect), EPSILON);
		assertEquals(Math.PI / 3, PGS_ShapePredicates.minimumInteriorAngle(triangle), EPSILON);
	}

	@Test
	void testMedian() {
		PVector median = PGS_ShapePredicates.median(square);
		assertEquals(5, median.x, EPSILON);
		assertEquals(5, median.y, EPSILON);
	}

	@Test
	void testContains() {
		PShape inner = new PShape(PShape.GEOMETRY);
		inner.beginShape();
		inner.vertex(2, 2);
		inner.vertex(8, 2);
		inner.vertex(8, 8);
		inner.vertex(2, 8);
		inner.endShape(PConstants.CLOSE);
		assertTrue(PGS_ShapePredicates.contains(square, inner));
		assertFalse(PGS_ShapePredicates.contains(inner, square));
	}

	@Test
	void testContainsSelf() {
		// A shape contains itself
		assertTrue(PGS_ShapePredicates.contains(square, square));
	}

	@Test
	void testContainsPoint() {
		assertTrue(PGS_ShapePredicates.containsPoint(square, new PVector(5, 5)));
		assertFalse(PGS_ShapePredicates.containsPoint(square, new PVector(15, 15)));
	}

	@Test
	void testContainsPointOnBoundary() {
		// Points on the boundary should be considered contained
		assertTrue(PGS_ShapePredicates.containsPoint(square, new PVector(0, 0)));
		assertTrue(PGS_ShapePredicates.containsPoint(square, new PVector(10, 5)));
		assertTrue(PGS_ShapePredicates.containsPoint(square, new PVector(5, 0)));
	}

	@Test
	void testContainsAllPoints() {
		assertTrue(PGS_ShapePredicates.containsAllPoints(square, Arrays.asList(new PVector(5, 5), new PVector(2, 2))));
		assertFalse(PGS_ShapePredicates.containsAllPoints(square, Arrays.asList(new PVector(5, 5), new PVector(15, 15))));
	}

	@Test
	void testContainsPoints() {
		List<Boolean> contained = PGS_ShapePredicates.containsPoints(square, Arrays.asList(new PVector(5, 5), new PVector(15, 15)));
		assertTrue(contained.get(0));
		assertFalse(contained.get(1));
	}

	@Test
	void testFindContainedPoints() {
		List<PVector> contained = PGS_ShapePredicates.findContainedPoints(square, Arrays.asList(new PVector(5, 5), new PVector(15, 15)));
		assertEquals(1, contained.size());
		assertEquals(new PVector(5, 5), contained.get(0));
	}

	@Test
	void testBoundsCenter() {
		PVector center = PGS_ShapePredicates.boundsCenter(square);
		assertEquals(5, center.x, EPSILON);
		assertEquals(5, center.y, EPSILON);
	}

	@Test
	void testBoundsCenterRect() {
		PVector center = PGS_ShapePredicates.boundsCenter(rect);
		assertEquals(5, center.x, EPSILON);
		assertEquals(10, center.y, EPSILON);
	}

	@Test
	void testIntersect() {
		PShape notIntersecting = new PShape(PShape.GEOMETRY);
		notIntersecting.beginShape();
		notIntersecting.vertex(20, 20);
		notIntersecting.vertex(30, 20);
		notIntersecting.vertex(30, 30);
		notIntersecting.vertex(20, 30);
		notIntersecting.endShape(PConstants.CLOSE);

		assertTrue(PGS_ShapePredicates.intersect(square, square));
		assertFalse(PGS_ShapePredicates.intersect(square, notIntersecting));
	}

	@Test
	void testIntersectPartialOverlap() {
		// Overlapping square shifted right by 5
		PShape shifted = makeSquare(5, 0, 10);
		assertTrue(PGS_ShapePredicates.intersect(square, shifted));
	}

	@Test
	void testTouch() {
		PShape adjacent = new PShape(PShape.GEOMETRY);
		adjacent.beginShape();
		adjacent.vertex(10, 0);
		adjacent.vertex(20, 0);
		adjacent.vertex(20, 10);
		adjacent.vertex(10, 10);
		adjacent.endShape(PConstants.CLOSE);

		assertTrue(PGS_ShapePredicates.touch(square, adjacent));
		assertFalse(PGS_ShapePredicates.touch(square, square)); // Covering, not just touching
	}

	@Test
	void testDistance() {
		PShape far = new PShape(PShape.GEOMETRY);
		far.beginShape();
		far.vertex(20, 0);
		far.vertex(30, 0);
		far.vertex(30, 10);
		far.vertex(20, 10);
		far.endShape(PConstants.CLOSE);

		assertEquals(10, PGS_ShapePredicates.distance(square, far), EPSILON);
		assertEquals(0, PGS_ShapePredicates.distance(square, square), EPSILON);
	}

	@Test
	void testDistanceAdjacent() {
		// Two touching squares should have distance 0
		PShape adjacent = makeSquare(10, 0, 10);
		assertEquals(0, PGS_ShapePredicates.distance(square, adjacent), EPSILON);
	}

	@Test
	void testDensity() {
		assertEquals(1.0, PGS_ShapePredicates.density(square), EPSILON);
		assertEquals(0.5, PGS_ShapePredicates.density(triangle), EPSILON); // Triangle bounding box is double its area
	}

	@Test
	void testDensityRect() {
		// An axis-aligned rectangle has density 1 (fills its envelope perfectly)
		assertEquals(1.0, PGS_ShapePredicates.density(rect), EPSILON);
	}

	@Test
	void testWidthHeightLength() {
		assertEquals(10, PGS_ShapePredicates.width(square), EPSILON);
		assertEquals(10, PGS_ShapePredicates.height(square), EPSILON);
		assertEquals(40, PGS_ShapePredicates.length(square), EPSILON); // Perimeter

		assertEquals(10, PGS_ShapePredicates.width(rect), EPSILON);
		assertEquals(20, PGS_ShapePredicates.height(rect), EPSILON);
		assertEquals(60, PGS_ShapePredicates.length(rect), EPSILON);
	}

	@Test
	void testCircularity() {
		double squareArea = 100;
		double squarePerim = 40;
		double expectedCircularity = (4 * Math.PI * squareArea) / (squarePerim * squarePerim);
		assertEquals(expectedCircularity, PGS_ShapePredicates.circularity(square), EPSILON);
	}

	@Test
	void testSphericity() {
		// Minimum bounding circle radius for 10x10 square is 5*sqrt(2) = 7.071
		// Maximum inscribed circle radius is 5
		// Sphericity = 5 / (5*sqrt(2)) = 1/sqrt(2) = ~0.707
		assertEquals(1.0 / Math.sqrt(2), PGS_ShapePredicates.sphericity(square), EPSILON);
	}

	@Test
	void testConvexity() {
		assertEquals(1.0, PGS_ShapePredicates.convexity(square), EPSILON);

		// Create a pacman shape (concave)
		PShape pacman = new PShape(PShape.GEOMETRY);
		pacman.beginShape();
		pacman.vertex(0, 0);
		pacman.vertex(10, 0);
		pacman.vertex(10, 10);
		pacman.vertex(0, 10);
		pacman.vertex(5, 5); // Concave part
		pacman.endShape(PConstants.CLOSE);

		// Area of 10x10 square = 100. The triangular bite is base 10, height 5 -> area
		// 25.
		// Pacman area = 75
		// Convex hull area = 100
		// Convexity = 75 / 100 = 0.75
		assertEquals(0.75, PGS_ShapePredicates.convexity(pacman), EPSILON);
	}

	@Test
	void testVertexCount() {
		assertEquals(4, PGS_ShapePredicates.vertexCount(square));
		assertEquals(4, PGS_ShapePredicates.vertexCount(rect));
		assertEquals(3, PGS_ShapePredicates.vertexCount(triangle));
	}

	@Test
	void testIsSimple() {
		assertTrue(PGS_ShapePredicates.isSimple(square));

		PShape selfIntersecting = new PShape(PShape.GEOMETRY);
		selfIntersecting.beginShape();
		selfIntersecting.vertex(0, 0);
		selfIntersecting.vertex(10, 10);
		selfIntersecting.vertex(0, 10);
		selfIntersecting.vertex(10, 0);
		selfIntersecting.endShape(PConstants.CLOSE);

		assertFalse(PGS_ShapePredicates.isSimple(selfIntersecting));
	}

	@Test
	void testIsConvex() {
		assertTrue(PGS_ShapePredicates.isConvex(square));
		assertTrue(PGS_ShapePredicates.isConvex(triangle));

		PShape pacman = new PShape(PShape.GEOMETRY);
		pacman.beginShape();
		pacman.vertex(0, 0);
		pacman.vertex(10, 0);
		pacman.vertex(10, 10);
		pacman.vertex(0, 10);
		pacman.vertex(5, 5);
		pacman.endShape(PConstants.CLOSE);

		assertFalse(PGS_ShapePredicates.isConvex(pacman));
	}

	@Test
	void testIsValid() {
		assertTrue(PGS_ShapePredicates.isValid(square));
	}

	@Test
	void testIsValidBowtie() {
		// A bowtie/figure-eight shape is self-intersecting and therefore invalid
		PShape bowtie = new PShape(PShape.GEOMETRY);
		bowtie.beginShape();
		bowtie.vertex(0, 0);
		bowtie.vertex(10, 10);
		bowtie.vertex(10, 0);
		bowtie.vertex(0, 10);
		bowtie.endShape(PConstants.CLOSE);

		assertFalse(PGS_ShapePredicates.isValid(bowtie));
	}

	@Test
	void testEqualsMethods() {
		assertTrue(PGS_ShapePredicates.equalsExact(square, square));
		assertTrue(PGS_ShapePredicates.equalsNorm(square, square));
		assertTrue(PGS_ShapePredicates.equalsTopo(square, square));

		assertFalse(PGS_ShapePredicates.equalsExact(square, rect));
	}

	@Test
	void testEqualsNormRotatedOrder() {
		// Same square but starting from a different vertex — equalsNorm should still
		// match because it normalises vertex order
		PShape rotated = new PShape(PShape.GEOMETRY);
		rotated.beginShape();
		rotated.vertex(10, 0);
		rotated.vertex(10, 10);
		rotated.vertex(0, 10);
		rotated.vertex(0, 0);
		rotated.endShape(PConstants.CLOSE);

		assertFalse(PGS_ShapePredicates.equalsExact(square, rotated));
		assertTrue(PGS_ShapePredicates.equalsNorm(square, rotated));
		assertTrue(PGS_ShapePredicates.equalsTopo(square, rotated));
	}

	@Test
	void testSimilarityIdentical() {
		// Identical shapes should have similarity 1.0
		assertEquals(1.0, PGS_ShapePredicates.similarity(square, square), EPSILON);
	}

	@Test
	void testSimilarityDifferent() {
		// Very different shapes should have low similarity
		PShape tiny = makeSquare(100, 100, 1);
		double sim = PGS_ShapePredicates.similarity(square, tiny);
		assertTrue(sim < 0.5, "Distant/different shapes should have low similarity, got " + sim);
	}

	@Test
	void testOverlapIdentical() {
		// Identical shapes overlap completely → 1.0
		assertEquals(1.0, PGS_ShapePredicates.overlap(square, square), EPSILON);
	}

	@Test
	void testOverlapDisjoint() {
		// Disjoint shapes have no overlap → 0.0
		PShape far = makeSquare(50, 50, 10);
		assertEquals(0.0, PGS_ShapePredicates.overlap(square, far), EPSILON);
	}

	@Test
	void testOverlapHalf() {
		// Two 10x10 squares, one shifted right by 5 → overlap region is 5x10 = 50
		// a1 = 100, a2 = 100, total = 200, w1 = w2 = 0.5
		// overlap = 0.5*(50/100) + 0.5*(50/100) = 0.5
		PShape shifted = makeSquare(5, 0, 10);
		assertEquals(0.5, PGS_ShapePredicates.overlap(square, shifted), EPSILON);
	}

	@Test
	void testEfdSimilarityIdentical() {
		// EFD distance between identical shapes should be 0
		assertEquals(0.0, PGS_ShapePredicates.efdSimilarity(square, square), EPSILON);
	}

	@Test
	void testEfdSimilarityDifferent() {
		// A very elongated rectangle should differ from a square
		PShape longRect = new PShape(PShape.GEOMETRY);
		longRect.beginShape();
		longRect.vertex(0, 0);
		longRect.vertex(100, 0);
		longRect.vertex(100, 1);
		longRect.vertex(0, 1);
		longRect.endShape(PConstants.CLOSE);

		double dist = PGS_ShapePredicates.efdSimilarity(square, longRect);
		assertTrue(dist > 0, "EFD distance between square and elongated rect should be > 0, got " + dist);
	}

	@Test
	void testFindContainingShape() {
		// Build a GROUP with two side-by-side children: left [0,0]-[10,10] and right
		// [10,0]-[20,10]
		PShape left = makeSquare(0, 0, 10);
		PShape right = makeSquare(10, 0, 10);
		PShape group = PGS_Conversion.flatten(left, right);

		// Point inside left child
		PShape found = PGS_ShapePredicates.findContainingShape(group, new PVector(5, 5));
		assertNotNull(found);

		// Point inside right child
		PShape found2 = PGS_ShapePredicates.findContainingShape(group, new PVector(15, 5));
		assertNotNull(found2);

		// Point outside both children
		PShape found3 = PGS_ShapePredicates.findContainingShape(group, new PVector(50, 50));
		assertNull(found3);
	}

	@Test
	void testIsConformingMesh() {
		// Two edge-adjacent squares form a conforming mesh
		PShape left = makeSquare(0, 0, 10);
		PShape right = makeSquare(10, 0, 10);
		PShape mesh = PGS_Conversion.flatten(left, right);
		assertTrue(PGS_ShapePredicates.isConformingMesh(mesh));
	}

	@Test
	void testIsConformingMeshOverlapping() {
		// Two overlapping squares do NOT form a conforming mesh
		PShape a = makeSquare(0, 0, 10);
		PShape b = makeSquare(5, 0, 10);
		PShape mesh = PGS_Conversion.flatten(a, b);
		assertFalse(PGS_ShapePredicates.isConformingMesh(mesh));
	}

}
