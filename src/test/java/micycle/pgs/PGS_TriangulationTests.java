package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;

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

}