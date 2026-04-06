package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import processing.core.PShape;

class PGS_ConstructionTests {

	private static final double EPSILON = 1E-5;

	@Test
	void testCreateRegularPolygonVertexCount() {
		// A hexagon should have 6 vertices
		PShape hex = PGS_Construction.createRegularPolygon(6, 50, 50, 20);
		assertEquals(6, PGS_ShapePredicates.vertexCount(hex));
	}

	@Test
	void testCreateRegularPolygonIsConvex() {
		PShape hex = PGS_Construction.createRegularPolygon(6, 50, 50, 20);
		assertTrue(PGS_ShapePredicates.isConvex(hex));
		assertTrue(PGS_ShapePredicates.isValid(hex));
	}

	@Test
	void testCreateRegularPolygonSquare() {
		// A 4-sided "regular polygon" should have 4 vertices
		PShape sq = PGS_Construction.createRegularPolygon(4, 0, 0, 10);
		assertEquals(4, PGS_ShapePredicates.vertexCount(sq));
		assertTrue(PGS_ShapePredicates.isConvex(sq));
	}

	@Test
	void testCreateRegularPolygonTriangle() {
		PShape tri = PGS_Construction.createRegularPolygon(3, 0, 0, 10);
		assertEquals(3, PGS_ShapePredicates.vertexCount(tri));
		assertTrue(PGS_ShapePredicates.isConvex(tri));
	}

	@Test
	void testCreateRandomPolygonVertexCount() {
		PShape poly = PGS_Construction.createRandomPolygon(7, 100, 100, 42);
		assertEquals(7, PGS_ShapePredicates.vertexCount(poly));
		assertTrue(PGS_ShapePredicates.isConvex(poly)); // random polygon is convex
		assertTrue(PGS_ShapePredicates.isValid(poly));
	}

	@Test
	void testCreateRandomPolygonDeterministic() {
		// Two calls with the same seed should produce the same polygon
		PShape a = PGS_Construction.createRandomPolygon(5, 100, 100, 42);
		PShape b = PGS_Construction.createRandomPolygon(5, 100, 100, 42);
		assertTrue(PGS_ShapePredicates.equalsExact(a, b));
	}

	@Test
	void testCreateRandomPolygonExactDimensions() {
		PShape poly = PGS_Construction.createRandomPolygonExact(6, 50, 80, 99);
		assertEquals(50, PGS_ShapePredicates.width(poly), EPSILON);
		assertEquals(80, PGS_ShapePredicates.height(poly), EPSILON);
	}

	@Test
	void testCreateStarNotConvex() {
		// Stars are concave
		PShape star = PGS_Construction.createStar(50, 50, 5, 10, 30, 0);
		assertTrue(!PGS_ShapePredicates.isConvex(star));
	}

	@Test
	void testCreateStarValid() {
		PShape star = PGS_Construction.createStar(50, 50, 5, 10, 30, 0);
		assertTrue(PGS_ShapePredicates.isValid(star));
	}

	@Test
	void testCreateRingHasHole() {
		PShape ring = PGS_Construction.createRing(50, 50, 30, 10);
		assertNotNull(ring);
		assertEquals(1, PGS_ShapePredicates.holes(ring));
	}

	@Test
	void testCreateRingArea() {
		// Area of ring = pi*(R² - r²) ≈ pi*(30² - 10²) = pi*800 ≈ 2513.27
		PShape ring = PGS_Construction.createRing(50, 50, 30, 10);
		double expectedArea = Math.PI * (30 * 30 - 10 * 10);
		assertEquals(expectedArea, PGS_ShapePredicates.area(ring), expectedArea * 0.02); // within 2%
	}

	@Test
	void testCreateHeartValid() {
		PShape heart = PGS_Construction.createHeart(50, 50, 40);
		assertNotNull(heart);
		assertTrue(PGS_ShapePredicates.area(heart) > 0);
		assertTrue(PGS_ShapePredicates.isValid(heart));
	}

	@Test
	void testCreateGearValid() {
		PShape gear = PGS_Construction.createGear(50, 50, 30, 8);
		assertNotNull(gear);
		assertTrue(PGS_ShapePredicates.area(gear) > 0);
		assertTrue(PGS_ShapePredicates.isValid(gear));
		assertTrue(!PGS_ShapePredicates.isConvex(gear)); // gears are concave
	}

	@Test
	void testCreateSupercircleSquarish() {
		// power=1 produces a diamond/rhombus shape (Lamé with p=1)
		PShape sc = PGS_Construction.createSupercircle(0, 0, 20, 1);
		assertNotNull(sc);
		assertTrue(PGS_ShapePredicates.area(sc) > 0);
	}

}
