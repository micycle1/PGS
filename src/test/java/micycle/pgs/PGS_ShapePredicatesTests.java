package micycle.pgs;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

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
		square.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex

		float[] centroid = new float[] { 0, 0 };
		float side_length = 10;
		triangle = new PShape(PShape.GEOMETRY); // equilateral triangle
		triangle.beginShape();
		triangle.vertex(centroid[0], centroid[1] + ((float) Math.sqrt(3) / 3) * side_length); // Top vertex
		triangle.vertex(centroid[0] - (side_length / 2), centroid[1] - ((float) Math.sqrt(3) / 6) * side_length); // Bottom left vertex
		triangle.vertex(centroid[0] + (side_length / 2), centroid[1] - ((float) Math.sqrt(3) / 6) * side_length); // Bottom right vertex
		triangle.endShape(PShape.CLOSE);
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
		assertEquals(10*Math.sqrt(2), PGS_ShapePredicates.diameter(square), EPSILON);
	}

	@Test
	void testMaximumInteriorAngle() {
		assertEquals(Math.PI / 2, PGS_ShapePredicates.maximumInteriorAngle(square));
		assertEquals(Math.PI / 3, PGS_ShapePredicates.maximumInteriorAngle(triangle), EPSILON);
	}

}
