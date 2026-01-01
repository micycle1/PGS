package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.tinfour.common.IIncrementalTin;

import processing.core.PShape;
import processing.core.PVector;

class PGS_TriangulationTests {

	@Test
	void testFromPoints() {
		List<PVector> points = PGS_PointSet.random(0, 0, 1000, 1000, 1000, 1337);
		IIncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(points);
		assertEquals(points.size(), tin.getVertices().size());
		int h = PGS_Hull.convexHull(points).getVertexCount(); // points on convex hull
		assertEquals(2 * points.size() - h - 2, tin.countTriangles().getCount());

		PShape triangulation = PGS_Triangulation.delaunayTriangulation(points);
		assertEquals(tin.countTriangles().getCount(), triangulation.getChildCount());
	}

	@Test
	void testEarCutTriangulation() {
		// Build a square via convexHull of 4 corners
		List<PVector> corners = new ArrayList<>();
		corners.add(new PVector(0, 0));
		corners.add(new PVector(100, 0));
		corners.add(new PVector(100, 100));
		corners.add(new PVector(0, 100));

		PShape square = PGS_Hull.convexHull(corners);
		assertNotNull(square);
		assertTrue(square.getVertexCount() >= 4);

		PShape triangles = PGS_Triangulation.earCutTriangulation(square);
		assertNotNull(triangles);

		// A simple convex n-gon triangulates into (n-2) triangles.
		assertEquals(2, triangles.getChildCount());
	}

	@Test
	void testRefine() {
		List<PVector> points = PGS_PointSet.random(0, 0, 1500, 1500, 500, 321);
		IIncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(points);

		int before = tin.countTriangles().getCount();

		PGS_Triangulation.refine(tin, 20);

		int after = tin.countTriangles().getCount();
		assertTrue(after >= before);
		assertTrue(tin.getVertices().size() >= points.size());
	}

}