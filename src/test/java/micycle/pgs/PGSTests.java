package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.UnaryOperator;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;

import micycle.pgs.commons.PEdge;
import processing.core.PConstants;
import processing.core.PShape;
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

	@Test
	void testApplyToLinealGeometries() {
		GeometryFactory gf = new GeometryFactory();

		// 1) Single LineString -> keep
		LineString ls = gf.createLineString(new Coordinate[] { new Coordinate(0, 0), new Coordinate(1, 1) });
		PShape lineShape = PGS_Conversion.toPShape(ls);
		UnaryOperator<LineString> keepAll = (LineString in) -> in; // identity
		PShape outLineShape = PGS.applyToLinealGeometries(lineShape, keepAll);
		assertNotNull(outLineShape, "LineString should be kept when function returns non-null");
		Geometry outGeom = PGS_Conversion.fromPShape(outLineShape);
		assertTrue(outGeom instanceof LineString, "Result should be a LineString");
		assertArrayEquals(ls.getCoordinates(), ((LineString) outGeom).getCoordinates(), "Coordinates should be unchanged");

		// 2) Single LineString -> drop (function returns null)
		UnaryOperator<LineString> dropAll = (lsIn) -> null;
		PShape dropped = PGS.applyToLinealGeometries(lineShape, dropAll);
		assertEquals(0, dropped.getChildCount());
		assertEquals(0, dropped.getVertexCount());

		// 3) Polygon with exterior + one hole -> drop hole only
		// exterior: square (0,0)-(4,0)-(4,4)-(0,4)-(0,0)
		LinearRing exterior = gf.createLinearRing(
				new Coordinate[] { new Coordinate(0, 0), new Coordinate(4, 0), new Coordinate(4, 4), new Coordinate(0, 4), new Coordinate(0, 0) });
		// hole: square (1,1)-(3,1)-(3,3)-(1,3)-(1,1)
		LinearRing hole = gf.createLinearRing(
				new Coordinate[] { new Coordinate(1, 1), new Coordinate(3, 1), new Coordinate(3, 3), new Coordinate(1, 3), new Coordinate(1, 1) });
		Polygon polyWithHole = gf.createPolygon(exterior, new LinearRing[] { hole });
		PShape polyShape = PGS_Conversion.toPShape(polyWithHole);

		// function that drops any ring whose first coordinate x == 1 (i.e., the hole)
		UnaryOperator<LineString> dropHoleIfStartsAt1 = (LineString in) -> {
			Coordinate c0 = in.getCoordinateN(0);
			if (Double.compare(c0.x, 1.0) == 0) {
				return null; // drop this ring (hole)
			}
			return in;
		};

		PShape polyProcessed = PGS.applyToLinealGeometries(polyShape, dropHoleIfStartsAt1);
		assertNotNull(polyProcessed, "Polygon with hole should remain when only hole is dropped");
		Geometry polyProcessedGeom = PGS_Conversion.fromPShape(polyProcessed);
		assertTrue(polyProcessedGeom instanceof Polygon, "Result should be a Polygon");
		Polygon pRes = (Polygon) polyProcessedGeom;
		assertEquals(0, pRes.getNumInteriorRing(), "Hole should have been removed");

		// 4) Polygon -> drop exterior => entire polygon is dropped
		UnaryOperator<LineString> dropExteriorIfStartsAt0 = (LineString in) -> {
			Coordinate c0 = in.getCoordinateN(0);
			if (Double.compare(c0.x, 0.0) == 0 && Double.compare(c0.y, 0.0) == 0) {
				return null; // drop exterior -> polygon should be dropped entirely
			}
			return in;
		};
		PShape polyDropped = PGS.applyToLinealGeometries(polyShape, dropExteriorIfStartsAt0);
		assertTrue(polyDropped.getChildCount() == 0 && polyDropped.getVertexCount() == 0,
				"If exterior ring is dropped, the whole polygon should be dropped (null returned)");

		// 5) MultiPolygon where one child is dropped and one kept
		// Polygon A (kept): square at origin without hole
		Polygon polyA = gf.createPolygon(
				gf.createLinearRing(
						new Coordinate[] { new Coordinate(0, 0), new Coordinate(2, 0), new Coordinate(2, 2), new Coordinate(0, 2), new Coordinate(0, 0) }),
				null);

		// Polygon B (to be dropped): square starting at x==10
		Polygon polyB = gf.createPolygon(gf.createLinearRing(
				new Coordinate[] { new Coordinate(10, 10), new Coordinate(12, 10), new Coordinate(12, 12), new Coordinate(10, 12), new Coordinate(10, 10) }),
				null);

		MultiPolygon multi = gf.createMultiPolygon(new Polygon[] { polyA, polyB });
		PShape multiShape = PGS_Conversion.toPShape(multi);

		// function that drops any ring starting at x >= 10 (so polygon B dropped)
		UnaryOperator<LineString> dropXge10 = (LineString in) -> {
			Coordinate c0 = in.getCoordinateN(0);
			if (c0.x >= 10.0) {
				return null;
			}
			return in;
		};

		PShape multiProcessed = PGS.applyToLinealGeometries(multiShape, dropXge10);
		assertNotNull(multiProcessed, "MultiPolygon with one surviving child should not be null");
		assertEquals(PConstants.GROUP, multiProcessed.getKind(), "Resulting PShape should be a GROUP");
		assertEquals(1, multiProcessed.getChildCount(), "GROUP should have exactly one child after dropping one polygon");

		Geometry multiProcGeom = PGS_Conversion.fromPShape(multiProcessed);
		// After transformation, should be a MultiPolygon or a Polygon depending on
		// builder; accept either but verify one child/polygon remains
		if (multiProcGeom instanceof MultiPolygon) {
			MultiPolygon mp = (MultiPolygon) multiProcGeom;
			assertEquals(1, mp.getNumGeometries(), "One polygon should remain in the MultiPolygon");
			assertTrue(mp.getGeometryN(0) instanceof Polygon, "Remaining geometry should be a Polygon");
		} else if (multiProcGeom instanceof Polygon) {
			// Possible that transformer collapses to single Polygon; check it's polyA
			// coordinates
			Polygon p = (Polygon) multiProcGeom;
			assertEquals(0, p.getNumInteriorRing(), "Remaining polygon should have no holes");

			assertTrue(polyA.getExteriorRing().equalsTopo(p.getExteriorRing()), "Remaining polygon should match polyA");
		} else {
			fail("Unexpected geometry type after processing MultiPolygon: " + multiProcGeom.getGeometryType());
		}
	}

}
