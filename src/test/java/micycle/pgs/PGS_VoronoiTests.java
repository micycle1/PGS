package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.Collection;
import java.util.List;
import java.util.function.BiFunction;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import processing.core.PShape;
import processing.core.PVector;

class PGS_VoronoiTests {

	static final int N = 100;
	static List<PVector> sites;
	static double[] bounds;

	@BeforeAll
	static void initSites() {
		sites = PGS_PointSet.random(200, 200, 800, 800, N, 0);
		// assign slightly different weights, but not enough to collapse cells
		// (so we can check voronoi faces against N)
		sites = PGS_PointSet.applyRandomWeights(sites, 1, 3, 0);
		bounds = new double[] { 0, 0, 1000, 1000 };
	}

	@Test
	void testVoronoi() {
		assertValidVoronoi(PGS_Voronoi::innerVoronoi);
	}

	@Test
	void testManhattenVoronoi() {
		assertValidVoronoi(PGS_Voronoi::manhattanVoronoi);
	}

	@Test
	void testPowerDiagram() {
		assertValidVoronoi(PGS_Voronoi::powerDiagram);
	}

	@Test
	void testAdditivelyWeightedVoronoi() {
		assertValidVoronoi(PGS_Voronoi::additivelyWeightedVoronoi);
	}

	@Test
	void testMultiplicativelyWeightedVoronoi() {
		assertValidVoronoi(PGS_Voronoi::multiplicativelyWeightedVoronoi, true);
	}

	private void assertValidVoronoi(BiFunction<Collection<PVector>, double[], PShape> op) {
		PShape vd = op.apply(sites, bounds);
		assertTrue(PGS_ShapePredicates.isConformingMesh(vd));
		assertEquals(N, vd.getChildCount());
	}

	private void assertValidVoronoi(TriFunction<Collection<PVector>, double[], Boolean, PShape> op, boolean flag) {
		PShape vd = op.apply(sites, bounds, flag);
		assertTrue(PGS_ShapePredicates.isConformingMesh(vd), "Voronoi output is non-conforming");
		assertEquals(N, vd.getChildCount());
	}

	@FunctionalInterface
	public interface TriFunction<T, U, V, R> {
		R apply(T t, U u, V v);
	}

}
