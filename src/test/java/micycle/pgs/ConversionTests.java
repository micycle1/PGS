package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;

import processing.core.PShape;
import processing.core.PVector;

class ConversionTests {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	@Test
	void testFromPShapeSimple() {
		final PShape shape = new PShape(PShape.GEOMETRY);
		shape.beginShape();
		shape.vertex(0, 0);
		shape.vertex(10, 0);
		shape.vertex(0, 10);
		shape.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex

		final Geometry g = fromPShape(shape);

		assertEquals(shape.getVertexCount() + 1, g.getCoordinates().length); // geometry is closed so has one more point
		assertEquals(1, g.getNumGeometries());
		for (int i = 0; i < g.getCoordinates().length; i++) {
			assertTrue(pointsAreEqual(g.getCoordinates()[i], shape.getVertex(i)));
		}
	}

	@Test
	void testFromPShapeHoles() {
		final PShape shape = new PShape(PShape.GEOMETRY);
		shape.beginShape();
		shape.vertex(0, 0);
		shape.vertex(10, 0);
		shape.vertex(0, 10);
		shape.vertex(10, 10);
		shape.beginContour();
		shape.vertex(2, 2);
		shape.vertex(8, 2);
		shape.vertex(2, 8);
		shape.vertex(8, 8);
		shape.endContour();
		shape.endShape(PShape.CLOSE);

		final Geometry g = fromPShape(shape);
		assertEquals(Geometry.TYPENAME_POLYGON, g.getGeometryType());

		final Polygon p = (Polygon) g;

		assertEquals(1, p.getNumInteriorRing());
		assertEquals(shape.getVertexCount() + 2, g.getNumPoints());
	}

	@Test
	void testToPShapeSimple() {
		Coordinate c1 = new Coordinate(0, 0);
		Coordinate c2 = new Coordinate(10, 0);
		Coordinate c3 = new Coordinate(0, 10);
		Coordinate c4 = new Coordinate(0, 0); // closed

		Coordinate[] coords = new Coordinate[] { c1, c2, c3, c4 };
		final Geometry g = GEOM_FACTORY.createPolygon(coords);

		final PShape shape = toPShape(g);

		assertEquals(g.getCoordinates().length - 1, shape.getVertexCount());
		assertEquals(0, shape.getChildCount());
		for (int i = 0; i < shape.getVertexCount(); i++) {
			assertTrue(pointsAreEqual(g.getCoordinates()[i], shape.getVertex(i)));
		}
	}

	@Test
	void testToPShapeHoles() {
		Coordinate c1 = new Coordinate(0, 0);
		Coordinate c2 = new Coordinate(10, 0);
		Coordinate c3 = new Coordinate(0, 10);
		Coordinate c4 = new Coordinate(0, 0); // closed
		Coordinate[] coords = new Coordinate[] { c1, c2, c3, c4 };

		// hole
		Coordinate h1 = new Coordinate(2, 2);
		Coordinate h2 = new Coordinate(8, 2);
		Coordinate h3 = new Coordinate(2, 8);
		Coordinate h4 = new Coordinate(2, 2); // closed
		Coordinate[] hole = new Coordinate[] { h1, h2, h3, h4 };
		LinearRing[] holes = new LinearRing[] { GEOM_FACTORY.createLinearRing(hole) };

		final Geometry g = GEOM_FACTORY.createPolygon(GEOM_FACTORY.createLinearRing(coords), holes);

		final PShape shape = toPShape(g);
		assertEquals(PShape.BREAK, shape.getVertexCode(3)); // contour break after third vertex (marks new hole)
		assertEquals(hole.length - 1 + coords.length - 1, shape.getVertexCount()); // contour break after third vertex
	}

	@Test
	void testPathToLinestring() {
		final PShape shape = new PShape(PShape.PATH);
		shape.beginShape();
		shape.vertex(0, 0);
		shape.vertex(10, 0);
		shape.vertex(0, 10);
		shape.vertex(10, 10);
		shape.endShape(); // unclosed

		final Geometry g = fromPShape(shape);
		assertEquals(Geometry.TYPENAME_LINESTRING, g.getGeometryType());

		assertEquals(shape.getVertexCount(), g.getCoordinates().length); // geometry is closed so has one more point
		for (int i = 0; i < g.getCoordinates().length; i++) {
			assertTrue(pointsAreEqual(g.getCoordinates()[i], shape.getVertex(i)));
		}
	}

	@Test
	void testLineStringToPath() {
		Coordinate c1 = new Coordinate(0, 0);
		Coordinate c2 = new Coordinate(10, 0);
		Coordinate c3 = new Coordinate(0, 10);
		Coordinate c4 = new Coordinate(10, 10); // unclosed

		Coordinate[] coords = new Coordinate[] { c1, c2, c3, c4 };
		final Geometry g = GEOM_FACTORY.createLineString(coords);

		final PShape shape = toPShape(g);
		assertEquals(PShape.PATH, shape.getFamily());

		assertEquals(g.getCoordinates().length, shape.getVertexCount());
		assertEquals(0, shape.getChildCount());
		for (int i = 0; i < shape.getVertexCount(); i++) {
			assertTrue(pointsAreEqual(g.getCoordinates()[i], shape.getVertex(i)));
		}
	}

	@Test
	void testPathsToMultiLinestring() {
		final PShape path1 = new PShape(PShape.PATH);
		path1.beginShape();
		path1.vertex(0, 0);
		path1.vertex(10, 0);
		path1.vertex(0, 10);
		path1.vertex(10, 10);
		path1.endShape(); // unclosed

		final PShape path2 = new PShape(PShape.PATH);
		path2.beginShape();
		path2.vertex(70, 70);
		path2.vertex(710, 70);
		path2.vertex(70, 710);
		path2.vertex(710, 710);
		path2.endShape(); // unclosed

		final PShape shape = new PShape(PShape.GROUP);
		shape.addChild(path1);
		shape.addChild(path2);

		final Geometry g = fromPShape(shape);
		assertEquals(Geometry.TYPENAME_MULTILINESTRING, g.getGeometryType());
		assertEquals(shape.getChildCount(), g.getNumGeometries());

		assertEquals(path1.getVertexCount() + path2.getVertexCount(), g.getCoordinates().length);
		for (int k = 0; k < g.getNumGeometries(); k++) {
			for (int i = 0; i < g.getGeometryN(k).getCoordinates().length; i++) {
				assertTrue(pointsAreEqual(g.getGeometryN(k).getCoordinates()[i], shape.getChild(k).getVertex(i)));
			}
		}
	}

	@Test
	void testMultiLinestringToPaths() {
		Coordinate c1 = new Coordinate(0, 0);
		Coordinate c2 = new Coordinate(10, 0);
		Coordinate c3 = new Coordinate(0, 10);
		Coordinate c4 = new Coordinate(10, 10); // unclosed

		Coordinate[] coords = new Coordinate[] { c1, c2, c3, c4 };
		final LineString path1 = GEOM_FACTORY.createLineString(coords);
		coords = new Coordinate[] { c4, c2, c1, c3 };
		final LineString path2 = GEOM_FACTORY.createLineString(coords);

		final Geometry g = GEOM_FACTORY.createMultiLineString(new LineString[] { path1, path2 });

		final PShape shape = toPShape(g);
		assertEquals(PShape.GROUP, shape.getFamily());
		assertEquals(g.getNumGeometries(), shape.getChildCount());

		assertEquals(g.getCoordinates().length, shape.getChild(0).getVertexCount() + shape.getChild(1).getVertexCount());
		for (int k = 0; k < g.getNumGeometries(); k++) {
			for (int i = 0; i < g.getGeometryN(k).getNumPoints(); i++) {
				assertTrue(pointsAreEqual(g.getGeometryN(k).getCoordinates()[i], shape.getChild(k).getVertex(i)));
			}
		}
	}

	@Test
	void testPointsToMultipoint() {
		final PShape points = new PShape(PShape.GEOMETRY);
		points.beginShape(PShape.POINTS);
		points.vertex(0, 0);
		points.vertex(10, 0);
		points.vertex(0, 10);
		points.vertex(10, 10);
		points.endShape();

		Geometry g = fromPShape(points);
		assertEquals(Geometry.TYPENAME_MULTIPOINT, g.getGeometryType());

		assertEquals(points.getVertexCount(), g.getCoordinates().length); // geometry is closed so has one more point
		for (int i = 0; i < g.getCoordinates().length; i++) {
			assertTrue(pointsAreEqual(g.getCoordinates()[i], points.getVertex(i)));
		}
	}

	@Test
	void testMultipointToPoints() {
		Coordinate c1 = new Coordinate(0, 0);
		Coordinate c2 = new Coordinate(10, 0);
		Coordinate c3 = new Coordinate(0, 10);
		Coordinate c4 = new Coordinate(10, 10);

		Coordinate[] coords = new Coordinate[] { c1, c2, c3, c4 };
		final Geometry g = GEOM_FACTORY.createMultiPointFromCoords(coords);

		final PShape shape = toPShape(g);
		assertEquals(PShape.POINTS, shape.getKind());

		assertEquals(g.getCoordinates().length, shape.getVertexCount());
		for (int i = 0; i < g.getCoordinates().length; i++) {
			assertTrue(pointsAreEqual(g.getCoordinates()[i], shape.getVertex(i)));
		}
	}

	@Test
	void testVertexRounding() {
		final PShape shape = new PShape(PShape.GEOMETRY);
		shape.beginShape();
		shape.vertex(12.4985f, -97.234f);
		shape.vertex(10, -10);
		shape.vertex(999.99f, 0.0001f);
		shape.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex

		PGS_Conversion.roundVertexCoords(shape);

		assertEquals(12, shape.getVertex(0).x);
		assertEquals(-97, shape.getVertex(0).y);
		assertEquals(10, shape.getVertex(1).x);
		assertEquals(-10, shape.getVertex(1).y);
		assertEquals(1000, shape.getVertex(2).x);
		assertEquals(0, shape.getVertex(2).y);
	}

	@Test
	void testDuplicateVertices() {
		final PShape shape = new PShape(PShape.GEOMETRY);
		shape.beginShape();
		shape.vertex(12.4985f, -97.234f);
		shape.vertex(12.4985f, -97.234f);
		shape.vertex(12.4985f, -97.234f);
		shape.vertex(10, -10);
		shape.vertex(10, -10);
		shape.vertex(12.4985f, -97.234f);
		shape.vertex(999.99f, 0.0001f);
		shape.vertex(999.99f, 0.0001f);
		shape.vertex(.99f, 0.0001f);
		shape.endShape(PShape.CLOSE); // close affects rendering only -- does not append another vertex

		final PShape processed = toPShape(fromPShape(shape));
		assertEquals(shape.getVertex(0), processed.getVertex(0));
		assertEquals(shape.getVertex(3), processed.getVertex(1));
		assertEquals(shape.getVertex(5), processed.getVertex(2));
		assertEquals(shape.getVertex(6), processed.getVertex(3));
		assertEquals(shape.getVertex(8), processed.getVertex(4));
		assertEquals(5, processed.getVertexCount());
	}

	private static boolean pointsAreEqual(Coordinate c, PVector p) {
		return (c.x == p.x && c.y == p.y);
	}

}
