package micycle.pgs;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;

class PGS_ShapeBooleanTests {

	@Test
	void testPolygonLineIntersection() { // a.k.a clipping
		PShape square = new PShape(PShape.GEOMETRY); // 10x10 square
		square.beginShape();
		square.vertex(0, 0);
		square.vertex(10, 0);
		square.vertex(10, 10);
		square.vertex(0, 10);
		square.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex

		PShape line = new PShape(PShape.PATH);
		line.beginShape(PConstants.LINES);
		line.vertex(-20, 5);
		line.vertex(20, 5);
		line.endShape();

		PShape intersection = PGS_ShapeBoolean.intersect(square, line);

		assertEquals(2, intersection.getVertexCount());
		assertEquals(0, intersection.getVertexX(0));
		assertEquals(10, intersection.getVertexX(1));
	}

	@Test
	void testPolygonLineDifference() {
		PShape square = new PShape(PShape.GEOMETRY); // 10x10 square
		square.beginShape();
		square.vertex(0, 0);
		square.vertex(10, 0);
		square.vertex(10, 10);
		square.vertex(0, 10);
		square.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex

		PShape line = new PShape(PShape.PATH);
		line.beginShape(PConstants.LINES);
		line.vertex(-20, 5);
		line.vertex(20, 5);
		line.endShape();

		PShape difference = PGS_ShapeBoolean.subtract(square, line);
		assertEquals(4 + 2, difference.getVertexCount());
		assertTrue(PGS_ShapePredicates.equalsTopo(square, difference));
	}
	
	@Test
	void testPolygonLineUnion() {
		PShape square = new PShape(PShape.GEOMETRY); // 10x10 square
		square.beginShape();
		square.vertex(0, 0);
		square.vertex(10, 0);
		square.vertex(10, 10);
		square.vertex(0, 10);
		square.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex
		
		PShape line = new PShape(PShape.PATH);
		line.beginShape(PConstants.LINES);
		line.vertex(-20, 5);
		line.vertex(20, 5);
		line.endShape();
		
		PShape union = PGS_ShapeBoolean.union(square, line);
		assertEquals(3, union.getChildCount());
		assertTrue(PGS_ShapePredicates.equalsTopo(square, union.getChild(0)));
	}

}
