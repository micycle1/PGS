package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.noding.BasicSegmentString;
import org.locationtech.jts.noding.SegmentString;

import processing.core.PShape;
import processing.core.PVector;

class PGS_ProcessingTests {

	@Test
	void extractPerimeter() {
		var r = PGS_Construction.createRect(0, 0, 1, 1, 0);
		var b1 = PGS_Processing.extractPerimeter(r, 0, 0.5);
		assertEquals(2, PGS_ShapePredicates.length(b1), 1e-6);

		var b2 = PGS_Processing.extractPerimeter(r, 0, 2);
		assertEquals(4, PGS_ShapePredicates.length(b2), 1e-6);

		// todo
//		var b3 = PGS_Processing.extractPerimeter(r, 1, 0);
//		assertEquals(4, PGS_ShapePredicates.length(b3), 1e-6);
//		assertFalse(boundary.isClosed());
	}

	@Test
	void intersectionPoints() {
		// Bow / X-shaped polyline in a single PATH - cross at (5,5)
		PShape path = new PShape(PShape.PATH);
		path.beginShape();
		path.vertex(0, 0);
		path.vertex(10, 10);
		path.vertex(0, 10);
		path.vertex(10, 0);
		path.endShape();

		List<PVector> hits = PGS_Processing.intersectionPoints(path);

		assertEquals(1, hits.size());
		assertContainsPoint(hits, 5, 5);
	}

	@Test
	void intersections() {
		// 1) Proper crossing (interior-interior) -> always included
		{
			SegmentString a = seg(0, 0, 10, 10);
			SegmentString b = seg(0, 10, 10, 0);

			List<PVector> hits = PGS_Processing.intersections(List.of(a, b), false);

			assertEquals(1, hits.size());
			assertContainsPoint(hits, 5, 5);
		}

		// 2) Endpoint touch (T-junction): endpoint of one hits interior of other
		// excluded when countEndpointTouches=false; included when true
		{
			SegmentString a = seg(0, 0, 10, 0); // horizontal
			SegmentString b = seg(5, 0, 5, 10); // vertical starting at (5,0) (endpoint touch)

			List<PVector> hitsNoTouches = PGS_Processing.intersections(List.of(a, b), false);
			assertEquals(0, hitsNoTouches.size());

			List<PVector> hitsWithTouches = PGS_Processing.intersections(List.of(a, b), true);
			assertEquals(1, hitsWithTouches.size());
			assertContainsPoint(hitsWithTouches, 5, 0);
		}

		// 3) Collinear overlap -> should return overlap endpoints
		{
			SegmentString a = seg(0, 0, 10, 0);
			SegmentString b = seg(5, 0, 15, 0);

			List<PVector> hits = PGS_Processing.intersections(List.of(a, b), false);

			assertEquals(2, hits.size());
			assertContainsPoint(hits, 5, 0);
			assertContainsPoint(hits, 10, 0);
		}

		// 4) De-duplication: multiple segments crossing at same point -> only one
		// output point
		{
			SegmentString a = seg(0, 0, 10, 10);
			SegmentString b = seg(0, 10, 10, 0);
			SegmentString c = seg(5, -10, 5, 20); // also passes through (5,5)

			List<PVector> hits = PGS_Processing.intersections(List.of(a, b, c), false);

			assertEquals(1, hits.size());
			assertContainsPoint(hits, 5, 5);
		}

		// 5) Endpoint touch (standalone)
		{
			SegmentString a = seg(0, 0, 1, 1);
			SegmentString b = seg(1, 1, 2, 2);

			List<PVector> hits = PGS_Processing.intersections(List.of(a, b), false);

			assertEquals(0, hits.size());

			hits = PGS_Processing.intersections(List.of(a, b), true);

			assertEquals(1, hits.size());
			assertContainsPoint(hits, 1, 1);
		}
	}

	private static SegmentString seg(double x0, double y0, double x1, double y1) {
		return new BasicSegmentString(new Coordinate[] { new Coordinate(x0, y0), new Coordinate(x1, y1) }, null);
	}

	private static void assertContainsPoint(List<PVector> pts, float x, float y) {
		final float eps = 1e-6f;
		for (PVector p : pts) {
			if (Math.abs(p.x - x) < eps && Math.abs(p.y - y) < eps) {
				return;
			}
		}
		assertTrue(false, "Expected point (" + x + ", " + y + ") not found in " + pts);
	}

}
