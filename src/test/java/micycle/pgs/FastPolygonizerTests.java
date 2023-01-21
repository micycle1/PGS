package micycle.pgs;

import static org.junit.jupiter.api.Assertions.*;
import static micycle.pgs.commons.FastPolygonizer.polygonize;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.junit.Ignore;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.operation.polygonize.Polygonizer;

import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import micycle.pgs.commons.PEdge;
import processing.core.PShape;

class FastPolygonizerTests {

	@Test
	void testSimpleTriangle() {
		PEdge a, b, c, d, e;
		a = new PEdge(0, 0, 10, 0);
		b = new PEdge(10, 0, 5, 5);
		c = new PEdge(5, 5, 0, 0);
		d = new PEdge(10, 0, -5, -5);
		e = new PEdge(-5, -5, 0, 0);
		PShape out = polygonize(Arrays.asList(a, b, c, d, e));
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testSimpleDangling() {
		PEdge a, b, c, d;
		a = new PEdge(0, 0, 10, 0);
		b = new PEdge(10, 0, 5, 5);
		c = new PEdge(5, 5, 0, 0);
		d = new PEdge(10, 0, -5, -5);
		List<PEdge> edges = new ArrayList<>();
		edges.addAll(Arrays.asList(a, b, c, d));
		for (int i = 11; i < 25; i++) { // add many dangling edges
			edges.add(new PEdge(i * 2, i * 2, i * 2 + 1, i * 2 + 1));
		}

		PShape out = polygonize(edges);
		assertEquals(1, out.getChildCount());
	}

	@Test
	void testComplexHalfDangling() {
		PEdge r, l, u, d;
		r = new PEdge(0, 0, 10, 0);
		l = new PEdge(0, 0, -10, 0);
		u = new PEdge(0, 0, 0, 10);
		d = new PEdge(0, 0, -10, 0);
		List<PEdge> edges = new ArrayList<>();
		edges.addAll(Arrays.asList(l, r, u, d));

		PShape out = polygonize(edges);
		assertEquals(0, out.getChildCount());

		PEdge join1 = new PEdge(0, 10, 10, 10);
		PEdge join2 = new PEdge(10, 0, 10, 10);
		edges.add(join1);
		edges.add(join2);

		out = polygonize(edges);
		assertEquals(1, out.getChildCount());

		PEdge join3 = new PEdge(0, 10, -10, 10);
		PEdge join4 = new PEdge(-10, 0, -10, 10);
		edges.add(join3);
		edges.add(join4);

		out = polygonize(edges);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testSimpleHalfDangling() {
		PEdge a, b, c, d, e;
		a = new PEdge(0, 0, 10, 0);
		b = new PEdge(10, 0, 5, 5);
		c = new PEdge(5, 5, 0, 0);
		d = new PEdge(10, 0, -5, -5);
		e = new PEdge(-5, -5, 1, 1); // half dangling

		PShape out = polygonize(Arrays.asList(a, b, c, d, e));
		assertEquals(1, out.getChildCount());
	}

	@Test
	@Ignore
	void testRobustnessRandomly() {
		Random r = new XoRoShiRo128PlusRandom(0);
		for (int k = 0; k < 100; k++) {
			List<SegmentString> segmentStrings = new ArrayList<>(111 + k);
			for (int i = 0; i < 111 + k; i++) {
				segmentStrings.add(new NodedSegmentString(new Coordinate[] { new Coordinate(r.nextDouble() * 10000, r.nextDouble() * 10000),
						new Coordinate(r.nextDouble() * 10000, r.nextDouble() * 13337) }, null));
			}
			Collection<SegmentString> nodedSS = PGS.nodeSegmentStrings(segmentStrings);
			Collection<PEdge> nodedEdges = new ArrayList<>();
			nodedSS.forEach(ss -> nodedEdges.add(new PEdge(PGS.toPVector(ss.getCoordinate(0)), PGS.toPVector(ss.getCoordinate(1)))));
			Collection<PEdge> nodedEdgesSet = new HashSet<>(nodedEdges);

			final Polygonizer polygonizerJTS = new Polygonizer();
			polygonizerJTS.setCheckRingsValid(false);
			nodedEdgesSet.forEach(e -> {
				final LineString l = PGS.GEOM_FACTORY
						.createLineString(new Coordinate[] { PGS.coordFromPVector(e.a), PGS.coordFromPVector(e.b) });
				polygonizerJTS.add(l);
			});

			long t1 = System.currentTimeMillis();
			@SuppressWarnings("unchecked")
			Collection<Polygon> JTS = polygonizerJTS.getPolygons();
			long t2 = System.currentTimeMillis();
			long timeJTS = t2 - t1;

			t1 = System.currentTimeMillis();
			PShape FP = polygonize(nodedEdgesSet);
			t2 = System.currentTimeMillis();
			long timeFP = t2 - t1;
//			System.out.println(timeJTS + " " + timeFP + " " + JTS.size());
			assertEquals(JTS.size(), FP.getChildCount());
		}
	}

}
