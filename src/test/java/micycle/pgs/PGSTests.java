package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.util.ArrayList;
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

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

class PGSTests {

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
//		assertEquals(PConstants.GROUP, multiProcessed.getKind(), "Resulting PShape should be a GROUP");
//		assertEquals(1, multiProcessed.getChildCount(), "GROUP should have exactly one child after dropping one polygon");

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

	@Test
	void testApplyToLinealGeometriesProcessingOrder() {
		GeometryFactory gf = new GeometryFactory();

		// Polygon with exterior + 2 holes so order is unambiguous
		LinearRing exterior = gf.createLinearRing(
				new Coordinate[] { new Coordinate(0, 0), new Coordinate(4, 0), new Coordinate(4, 4), new Coordinate(0, 4), new Coordinate(0, 0) });
		LinearRing hole1 = gf.createLinearRing(
				new Coordinate[] { new Coordinate(1, 1), new Coordinate(2, 1), new Coordinate(2, 2), new Coordinate(1, 2), new Coordinate(1, 1) });
		LinearRing hole2 = gf.createLinearRing(
				new Coordinate[] { new Coordinate(3, 3), new Coordinate(3.5, 3), new Coordinate(3.5, 3.5), new Coordinate(3, 3.5), new Coordinate(3, 3) });
		Polygon poly = gf.createPolygon(exterior, new LinearRing[] { hole1, hole2 });
		PShape polyShape = PGS_Conversion.toPShape(poly);

		// MultiPolygon with 3 polygons; middle one will be dropped, so we can verify
		// survivor order too
		Polygon polyA = gf.createPolygon(
				gf.createLinearRing(
						new Coordinate[] { new Coordinate(0, 0), new Coordinate(2, 0), new Coordinate(2, 2), new Coordinate(0, 2), new Coordinate(0, 0) }),
				null);
		Polygon polyB = gf.createPolygon(gf.createLinearRing(
				new Coordinate[] { new Coordinate(10, 10), new Coordinate(12, 10), new Coordinate(12, 12), new Coordinate(10, 12), new Coordinate(10, 10) }),
				null);
		Polygon polyC = gf.createPolygon(gf.createLinearRing(
				new Coordinate[] { new Coordinate(20, 20), new Coordinate(22, 20), new Coordinate(22, 22), new Coordinate(20, 22), new Coordinate(20, 20) }),
				null);

		MultiPolygon mp = gf.createMultiPolygon(new Polygon[] { polyA, polyB, polyC });
		PShape mpShape = PGS_Conversion.toPShape(mp);

		// (A) Verify CALLING order for Polygon rings
		List<String> polyCallOrder = new ArrayList<>();
		UnaryOperator<LineString> recordPolyCalls = (LineString in) -> {
			Coordinate c0 = in.getCoordinateN(0);
			polyCallOrder.add(c0.x + "," + c0.y);
			return in;
		};

		PGS.applyToLinealGeometries(polyShape, recordPolyCalls);

		assertEquals(List.of("0.0,0.0", "1.0,1.0", "3.0,3.0"), polyCallOrder,
				"Polygon ring processing order should be: exterior, then holes in interior-ring index order");

		// Verify CALLING order for MultiPolygon children
		List<String> mpCallOrder = new ArrayList<>();
		UnaryOperator<LineString> recordMpCalls = (LineString in) -> {
			Coordinate c0 = in.getCoordinateN(0);
			mpCallOrder.add(c0.x + "," + c0.y);
			return in;
		};

		PGS.applyToLinealGeometries(mpShape, recordMpCalls);

		assertEquals(List.of("0.0,0.0", "10.0,10.0", "20.0,20.0"), mpCallOrder,
				"MultiPolygon processing order should follow geometry index order (A, then B, then C)");

		// Verify OUTPUT order of surviving geometries is preserved after dropping B
		UnaryOperator<LineString> dropB = (LineString in) -> {
			double x0 = in.getCoordinateN(0).x;
			return (x0 == 10.0) ? null : in; // drop polygon B's exterior ring => polygon B removed
		};

		PShape outShape = PGS.applyToLinealGeometries(mpShape, dropB);
		Geometry outGeom = PGS_Conversion.fromPShape(outShape);

		if (outGeom instanceof MultiPolygon outMp) {
			assertEquals(2, outMp.getNumGeometries(), "After dropping B, exactly 2 polygons should remain");

			Polygon first = (Polygon) outMp.getGeometryN(0);
			Polygon second = (Polygon) outMp.getGeometryN(1);

			assertEquals(0.0, first.getExteriorRing().getCoordinateN(0).x, 0.0, "First survivor should be A");
			assertEquals(20.0, second.getExteriorRing().getCoordinateN(0).x, 0.0, "Second survivor should be C");
		} else if (outGeom instanceof Polygon) {
			fail("Expected MultiPolygon with survivors A and C, but got single Polygon (ordering cannot be verified)");
		} else {
			fail("Unexpected geometry type after dropping B: " + outGeom.getGeometryType());
		}
	}

}
