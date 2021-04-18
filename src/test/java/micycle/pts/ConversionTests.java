package micycle.pts;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;

import processing.core.PShape;
import processing.core.PVector;

public class ConversionTests {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(
			new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	@Test
	public void testFromPShapeSimple() {
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
	public void testFromPShapeHoles() {
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
	public void testToPShapeSimple() {
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
	public void testToPShapeHoles() {
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

	private static boolean pointsAreEqual(Coordinate c, PVector p) {
		return (c.x == p.x && c.y == p.y);
	}

}
