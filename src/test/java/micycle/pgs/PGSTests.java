package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;

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
		edges.add(new PEdge(15, 15, 0,0)); // close

		Collections.shuffle(edges);

		List<PVector> orderedVertices = PGS.fromEdges(edges);
		PGS.fromEdges(edges).forEach(q -> System.out.println(q));
		assertEquals(16, orderedVertices.size());
	}

}
