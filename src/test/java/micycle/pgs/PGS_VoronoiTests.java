package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import processing.core.PVector;

class PGS_VoronoiTests {

	static final int N = 50;
	static List<PVector> sites;
	static double[] bounds;

	@BeforeAll
	static void initSites() {
		sites = PGS_PointSet.random(200, 200, 800, 800, N, 0);
		bounds = new double[] { 0, 0, 1000, 1000 };
	}

	@Test
	void testVoronoi() {
		assertEquals(N, PGS_Voronoi.innerVoronoi(sites, bounds).getChildCount());
	}
	
	@Test
	void testManhattenVoronoi() {
		assertEquals(N, PGS_Voronoi.manhattanVoronoi(sites, bounds).getChildCount());
	}
	
	@Test
	void testPowerDiagram() {
		assertEquals(N, PGS_Voronoi.powerDiagram(sites, bounds).getChildCount());
	}
	
	@Test
	void testMultiplicativelyWeightedVoronoi() {
		assertEquals(N, PGS_Voronoi.multiplicativelyWeightedVoronoi(sites, bounds).getChildCount());
	}

}
