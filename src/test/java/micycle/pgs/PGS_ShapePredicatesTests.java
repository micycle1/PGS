package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

class PGS_ShapePredicatesTests {

	private static final double EPSILON = 1E-4;

	static PShape square, triangle;

	@BeforeAll
	static void initShapes() {
		square = new PShape(PShape.GEOMETRY); // 10x10 square
		square.beginShape();
		square.vertex(0, 0);
		square.vertex(10, 0);
		square.vertex(10, 10);
		square.vertex(0, 10);
		square.endShape(PConstants.CLOSE); // close affects rendering only -- does not append another vertex

		float[] centroid = new float[] { 0, 0 };
		float side_length = 10;
		triangle = new PShape(PShape.GEOMETRY); // equilateral triangle
		triangle.beginShape();
		triangle.vertex(centroid[0], centroid[1] + ((float) Math.sqrt(3) / 3) * side_length); // Top vertex
		triangle.vertex(centroid[0] - (side_length / 2), centroid[1] - ((float) Math.sqrt(3) / 6) * side_length); // Bottom left vertex
		triangle.vertex(centroid[0] + (side_length / 2), centroid[1] - ((float) Math.sqrt(3) / 6) * side_length); // Bottom right vertex
		triangle.endShape(PConstants.CLOSE);
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

		PShape coverage = PGS_Processing.split(withHole); // test coverage; no face has a hole, but together they form a mesh with a hole
		assertEquals(1, PGS_ShapePredicates.holes(coverage));

		coverage.removeChild(0); // remove a mesh face; mesh no longer forms a hole
		assertEquals(0, PGS_ShapePredicates.holes(coverage));
	}

	@Test
	void testIsClockwise() {
		assertTrue(PGS_ShapePredicates.isClockwise(square));
		List<PVector> ccw = PGS_Conversion.toPVector(square);// .reversed();
		Collections.reverse(ccw);
		ccw.add(ccw.get(0)); // close
		assertFalse(PGS_ShapePredicates.isClockwise(PGS_Conversion.fromPVector(ccw)));

	}

}
