package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static micycle.pgs.PGS_ShapePredicates.area;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import processing.core.PShape;

class PGS_TransformationTests {

	private static final double EPSILON = 1e-6;

	static PShape square;

	@BeforeAll
	static void initShapes() {
		square = new PShape(PShape.GEOMETRY); // 10x10 square
		square.beginShape();
		square.vertex(0, 0);
		square.vertex(10, 0);
		square.vertex(10, 10);
		square.vertex(0, 10);
		square.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex
		assertEquals(100, area(square));
	}

	@Test
	void testScaleAreaBy() {
		assertEquals(50, area(PGS_Transformation.scaleArea(square, 0.5)), EPSILON);
		assertEquals(125, area(PGS_Transformation.scaleArea(square, 1.25)), EPSILON);
	}

	@Test
	void testScaleAreaTo() {
		assertEquals(50, area(PGS_Transformation.scaleAreaTo(square, 50)), EPSILON);
		assertEquals(125, area(PGS_Transformation.scaleAreaTo(square, 125)), EPSILON);
	}

	@Test
	void testScale() {
		assertEquals(100 * 1.5 * 1.5, area(PGS_Transformation.scale(square, 1.5)), EPSILON);
	}

}
