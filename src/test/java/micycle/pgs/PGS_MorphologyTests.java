package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;

import micycle.pgs.commons.DiscreteCurveEvolution.DCETerminationCallback;
import processing.core.PShape;

public class PGS_MorphologyTests {

	private static final GeometryFactory GF = new GeometryFactory();

	private PShape inShape;
	private Geometry inGeom;
	private double[] originalAreas;
	private int[] originalHoleCounts;

	@BeforeEach
	public void setUpPolygons() {
		// Build two polygons with holes and outward/inward spikes (same as in your
		// smoothGaussian test)
		Polygon p1 = buildPolygonWithHoles(spikyRectRing(0, 0, 10, 8, 0.8, true), new LinearRing[] { spikyRectRing(2, 2, 4.5, 4.5, 0.6, false) });

		Polygon p2 = buildPolygonWithHoles(spikyRectRing(20, 0, 34, 12, 1.0, true),
				new LinearRing[] { spikyRectRing(22.5, 2.5, 27.5, 5.5, 0.7, false), spikyRectRing(28.0, 7.0, 32.0, 10.0, 0.5, false) });

		MultiPolygon mp = GF.createMultiPolygon(new Polygon[] { p1, p2 });
		inShape = PGS_Conversion.toPShape(mp);

		// Basic sanity
		assertNotNull(inShape);
		assertTrue(inShape.getChildCount() >= 2);

		inGeom = PGS_Conversion.fromPShape(inShape);
		assertEquals(2, inGeom.getNumGeometries());

		originalAreas = new double[2];
		originalHoleCounts = new int[2];
		for (int i = 0; i < 2; i++) {
			Geometry gi = inGeom.getGeometryN(i);
			assertTrue(gi instanceof Polygon);
			originalHoleCounts[i] = ((Polygon) gi).getNumInteriorRing();
			originalAreas[i] = PGS_ShapePredicates.area(inShape.getChild(i));
		}
	}

	@Test
	public void testHobbySimplify() {
		PShape outShape = PGS_Morphology.simplifyHobby(inShape, 1);
		Geometry outGeom = getOutputGeom(outShape);
		assertPolygonsAndHoleCounts(outGeom);
		assertAreasDecreased(outShape, "after elliptic Fourier smoothing");
	}

	@Test
	public void testRounding() {
		PShape outShape = PGS_Morphology.round(inShape, 20);
		Geometry outGeom = getOutputGeom(outShape);
		assertPolygonsAndHoleCounts(outGeom);
		assertAreasDecreased(outShape, "after rounding");
	}

	@Test
	public void testChaikinCut() {
		PShape outShape = PGS_Morphology.chaikinCut(inShape, 0.5, 3);
		Geometry outGeom = getOutputGeom(outShape);
		assertPolygonsAndHoleCounts(outGeom);
		assertAreasDecreased(outShape, "after Chaikin cutting");
	}

	@Test
	public void testSmoothEllipticFourier() {
		// choose descriptors moderately large to preserve smoothing effect
		int descriptors = 10;

		PShape outShape = PGS_Morphology.smoothEllipticFourier(inShape, descriptors);
		Geometry outGeom = getOutputGeom(outShape);
		assertPolygonsAndHoleCounts(outGeom);
		assertAreasDecreased(outShape, "after elliptic Fourier smoothing");
	}

	@Test
	public void testSimplifyDCE() {
		// Create a simple termination callback: stop when remaining vertices <=
		// threshold
		final int targetVertices = 8;
		DCETerminationCallback callback = new DCETerminationCallback() {
			@Override
			public boolean shouldTerminate(Coordinate currentVertex, double relevance, int verticesRemaining) {
				return verticesRemaining <= targetVertices;
			}
		};

		// Apply DCE simplification
		PShape outShape = PGS_Morphology.simplifyDCE(inShape, callback);

		// Output should be a non-null group with two polygonal children
		Geometry outGeom = getOutputGeom(outShape);
		// additional DCE-specific assertions (ensure output geometries are polygons)
		for (int i = 0; i < 2; i++) {
			assertTrue(outGeom.getGeometryN(i) instanceof Polygon);
		}

		// Structure: preserve polygon count and hole counts
		assertPolygonsAndHoleCounts(outGeom);

		// Areas: each polygon's area should be reduced after simplification
		assertAreasDecreased(outShape, "after DCE simplification");
	}

	@Test
	public void testSmoothGaussian() {
		// Apply smoothing
		double sigma = 1.5;
		PShape outShape = PGS_Morphology.smoothGaussian(inShape, sigma);

		Geometry outGeom = getOutputGeom(outShape);
		assertPolygonsAndHoleCounts(outGeom);
		assertAreasDecreased(outShape, "after smoothing");
	}

	/* Helper factories and geometry builders */

	private static Polygon buildPolygonWithHoles(LinearRing exterior, LinearRing[] holes) {
		return GF.createPolygon(exterior, holes);
	}

// Build a rectangular ring with two spike points per edge.
// For spikesOutward = true: spikes point outside the rectangle bounds.
// For spikesOutward = false (holes): spikes point toward the rectangle center.
	private static LinearRing spikyRectRing(double minX, double minY, double maxX, double maxY, double amplitude, boolean spikesOutward) {

		List<Coordinate> coords = new ArrayList<>();

		// Bottom edge: (minX,minY) -> (maxX,minY)
		coords.add(new Coordinate(minX, minY));
		coords.add(new Coordinate(lerp(minX, maxX, 1.0 / 3.0), minY + (spikesOutward ? -amplitude : +amplitude)));
		coords.add(new Coordinate(lerp(minX, maxX, 2.0 / 3.0), minY + (spikesOutward ? -amplitude : +amplitude)));
		coords.add(new Coordinate(maxX, minY));

		// Right edge: (maxX,minY) -> (maxX,maxY)
		coords.add(new Coordinate(maxX + (spikesOutward ? +amplitude : -amplitude), lerp(minY, maxY, 1.0 / 3.0)));
		coords.add(new Coordinate(maxX + (spikesOutward ? +amplitude : -amplitude), lerp(minY, maxY, 2.0 / 3.0)));
		coords.add(new Coordinate(maxX, maxY));

		// Top edge: (maxX,maxY) -> (minX,maxY)
		coords.add(new Coordinate(lerp(maxX, minX, 2.0 / 3.0), maxY + (spikesOutward ? +amplitude : -amplitude)));
		coords.add(new Coordinate(lerp(maxX, minX, 1.0 / 3.0), maxY + (spikesOutward ? +amplitude : -amplitude)));
		coords.add(new Coordinate(minX, maxY));

		// Left edge: (minX,maxY) -> (minX,minY)
		coords.add(new Coordinate(minX + (spikesOutward ? -amplitude : +amplitude), lerp(maxY, minY, 2.0 / 3.0)));
		coords.add(new Coordinate(minX + (spikesOutward ? -amplitude : +amplitude), lerp(maxY, minY, 1.0 / 3.0)));

		// Close ring
		coords.add(new Coordinate(minX, minY));

		return GF.createLinearRing(coords.toArray(new Coordinate[0]));
	}

	private static double lerp(double a, double b, double t) {
		return a + (b - a) * t;
	}

	/* New helper assertion methods to remove duplication */

	private Geometry getOutputGeom(PShape outShape) {
		assertNotNull(outShape, "Output shape must not be null");
		Geometry outGeom = PGS_Conversion.fromPShape(outShape);
		assertEquals(2, outGeom.getNumGeometries(), "Output geometry must contain two geometries");
		assertEquals(2, outShape.getChildCount(), "Output PShape must have two children");
		return outGeom;
	}

	private void assertPolygonsAndHoleCounts(Geometry outGeom) {
		for (int i = 0; i < 2; i++) {
			assertTrue(outGeom.getGeometryN(i) instanceof Polygon, "Geometry " + i + " must be a Polygon");
			Polygon op = (Polygon) outGeom.getGeometryN(i);
			assertEquals(originalHoleCounts[i], op.getNumInteriorRing(), "Hole count should be preserved for polygon " + i);
		}
	}

	private void assertAreasDecreased(PShape outShape, String messagePrefix) {
		double outArea0 = PGS_ShapePredicates.area(outShape.getChild(0));
		double outArea1 = PGS_ShapePredicates.area(outShape.getChild(1));
		assertTrue(outArea0 < originalAreas[0], "Polygon 0 area should decrease " + messagePrefix);
		assertTrue(outArea1 < originalAreas[1], "Polygon 1 area should decrease " + messagePrefix);
	}
}