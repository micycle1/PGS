package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import micycle.pgs.commons.PEdge;
import processing.core.PVector;

class PGSTests {

	@Test
	void testFromEdgesSimple() {
		PEdge a = new PEdge(0, 0, 1, 1);
		PEdge b = new PEdge(1, 1, 1, 0);
		PEdge c = new PEdge(1, 0, 0, 0);

		List<PEdge> edges = Arrays.asList(a, c, b); // a, c, b

		List<PVector> orderedVertices = PGS.fromEdges(edges);
		assertEquals(3, orderedVertices.size());
	}

	@Test
	void testFromEdges() {
		List<PEdge> edges = new ArrayList<>();
		for (int i = 0; i < 15; i++) {
			edges.add(new PEdge(i, i, i + 1, i + 1));
		}
		edges.add(new PEdge(15, 15, 0, 0)); // close

		Collections.shuffle(edges);

		List<PVector> orderedVertices = PGS.fromEdges(edges);
		PGS.fromEdges(edges).forEach(q -> System.out.println(q));
		assertEquals(16, orderedVertices.size());
	}

	@Test
	void testOrientation() {
		/*
		 * NOTE the isClockwise method tests for orientation in a y-axis-up coordinate
		 * system. The lists below are defined in terms of visual orientation in
		 * Processing (which uses y-axis-down orientation). Hence the results are
		 * inverted to check the method is geometrically correct.
		 */

		// @formatter:off
        //		(0,0) --------> (1,0)  (x-axis increasing to the right)
        //		   |              |
        //		   |              |
        //		   V              V
        //		(0,1) <--------- (1,1)  (y-axis increasing downwards)
		// @formatter:on
		List<PVector> clockwisePoints = List.of(new PVector(0, 0), new PVector(1, 0), new PVector(1, 1), new PVector(0, 1));
		assertTrue(!PGS.isClockwise(clockwisePoints)); // NOTE inverted

		List<PVector> counterClockwisePoints = List.of(new PVector(0, 0), new PVector(0, 1), new PVector(1, 1), new PVector(1, 0));
		assertTrue(PGS.isClockwise(counterClockwisePoints)); // NOTE inverted

		List<PVector> clockwisePointsClosed = new ArrayList<>(
				List.of(new PVector(0, 0), new PVector(1, 0), new PVector(1, 1), new PVector(0, 1), new PVector(0, 0)));
		assertTrue(!PGS.isClockwise(clockwisePointsClosed)); // NOTE inverted

		List<PVector> counterClockwisePointsClosed = new ArrayList<>(
				List.of(new PVector(0, 0), new PVector(0, 1), new PVector(1, 1), new PVector(1, 0), new PVector(0, 0)));
		assertTrue(PGS.isClockwise(counterClockwisePointsClosed)); // NOTE inverted
	}

}
