package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Collection;
import java.util.function.Function;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import processing.core.PShape;
import processing.core.PVector;

class PGS_PolygonisationTests {

	private Collection<PVector> points;

	@BeforeEach
	void setup() {
		points = PGS_PointSet.random(0, 0, 1, 1, 1000, 1337L);
	}

	private void assertValidPolygonisation(Function<Collection<PVector>, PShape> polygoniserMethod) {
		final PShape out = polygoniserMethod.apply(points);

		assertNotNull(out, "Polygonisation returned null.");
		assertTrue(PGS_ShapePredicates.isSimple(out), "Polygonisation output is not a simple polygon.");
		assertEquals(points.size(), out.getVertexCount(), "Unexpected vertex count: polygon does not appear to use exactly all input points.");
	}

	@Test
	void testMinArea() {
		assertValidPolygonisation(PGS_Polygonisation::minArea);
	}

	@Test
	void testMaxArea() {
		assertValidPolygonisation(PGS_Polygonisation::maxArea);
	}

	@Test
	void testMinPerimeter() {
		assertValidPolygonisation(PGS_Polygonisation::minPerimeter);
	}

	@Test
	void testHorizontal() {
		assertValidPolygonisation(PGS_Polygonisation::horizontal);
	}

	@Test
	void testVertical() {
		assertValidPolygonisation(PGS_Polygonisation::vertical);
	}

	@Test
	void testHilbert() {
		assertValidPolygonisation(PGS_Polygonisation::hilbert);
	}

	@Test
	void testCircular() {
		assertValidPolygonisation(PGS_Polygonisation::circular);
	}

	@Test
	void testAngular() {
		assertValidPolygonisation(PGS_Polygonisation::angular);
	}

	@Test
	void testOnion() {
		assertValidPolygonisation(PGS_Polygonisation::onion);
	}
}