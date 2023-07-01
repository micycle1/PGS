package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import micycle.pgs.commons.PEdge;
import processing.core.PShape;
import processing.core.PVector;

class PGS_SegmentSetTests {

	@Test
	void testDissolve() {
		PVector p1 = new PVector(0, 0);
		PVector p2 = new PVector(10, 10);
		PVector p3 = new PVector(20, 20);
		PVector p4 = new PVector(30, 30);

		PEdge a = new PEdge(p1, p2);
		PEdge b = new PEdge(p2, p3);
		List<PEdge> edges = new ArrayList<>();
		edges.add(a);
		edges.add(b);
		edges.add(b); // note test duplicate edge

		PShape dissolved = PGS_SegmentSet.dissolve(edges);
		assertEquals(0, dissolved.getChildCount());
		assertEquals(3, dissolved.getVertexCount());
		
		PEdge c = new PEdge(p3, p4);
		edges.add(c);
		dissolved = PGS_SegmentSet.dissolve(edges);
		assertEquals(0, dissolved.getChildCount());
		assertEquals(4, dissolved.getVertexCount()); // a merged line of p1,p2,p3,p4
		
		PEdge d = new PEdge(p2, p4);
		edges.add(d); // add a branch to existing linestring
		dissolved = PGS_SegmentSet.dissolve(edges);
		assertEquals(2, dissolved.getChildCount()); // it's now 2 distinct merged linestrings
	}

	@Test
	void testToPShape() {
		PVector p1 = new PVector(0, 0);
		PVector p2 = new PVector(10, 10);
		PVector p3 = new PVector(20, 20);

		PEdge a = new PEdge(p1, p2);
		PEdge b = new PEdge(p2, p3);
		List<PEdge> edges = new ArrayList<>();
		edges.add(a);
		edges.add(b);

		PShape lines = PGS_SegmentSet.toPShape(edges);
		assertEquals(0, lines.getChildCount());
		assertEquals(4, lines.getVertexCount());

	}

}
