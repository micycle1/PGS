package micycle.pgs;

import static micycle.pgs.PGS_ShapePredicates.area;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import processing.core.PConstants;
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
		square.endShape(PConstants.CLOSE); // close affects rendering only -- does not append another vertex
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
		// Test polygon scaling
		assertEquals(100 * 1.5 * 1.5, area(PGS_Transformation.scale(square, 1.5)), EPSILON);
		
		// Test empty geometry - should not throw error
		PShape emptyGeom = PGS_Conversion.fromWKT("POLYGON EMPTY");
		PShape scaledEmpty = PGS_Transformation.scale(emptyGeom, 1.5);
		assertNotNull(scaledEmpty);
		assertEquals(0, scaledEmpty.getVertexCount());
		
		// Test point (0-dimensional/pointal) - should not throw error
		PShape point = PGS_Conversion.fromWKT("POINT (10 10)");
		PShape scaledPoint = PGS_Transformation.scale(point, 2.0);
		assertNotNull(scaledPoint);
		
		// Test linestring (1-dimensional/lineal) - should not throw error
		PShape line = PGS_Conversion.fromWKT("LINESTRING (0 0, 10 10)");
		PShape scaledLine = PGS_Transformation.scale(line, 2.0);
		assertNotNull(scaledLine);
		
		// Test empty point - should not throw error
		PShape emptyPoint = PGS_Conversion.fromWKT("POINT EMPTY");
		PShape scaledEmptyPoint = PGS_Transformation.scale(emptyPoint, 1.5);
		assertNotNull(scaledEmptyPoint);
		
		// Test empty linestring - should not throw error
		PShape emptyLine = PGS_Conversion.fromWKT("LINESTRING EMPTY");
		PShape scaledEmptyLine = PGS_Transformation.scale(emptyLine, 1.5);
		assertNotNull(scaledEmptyLine);
	}

}
