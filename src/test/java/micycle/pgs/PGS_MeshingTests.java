package micycle.pgs;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import processing.core.PShape;

public class PGS_MeshingTests {

	@Test
	void testAreaMerge() {
		PShape mesh = PGS_Triangulation.delaunayTriangulation(PGS_PointSet.random(0, 0, 1000, 1000, 1111, 0));
		List<PShape> faces = PGS_Conversion.getChildren(mesh);
		faces.sort((a, b) -> Double.compare(PGS_ShapePredicates.area(a), PGS_ShapePredicates.area(b)));
		double areaThreshold = PGS_ShapePredicates.area(faces.get(faces.size() / 2));

		PShape mergedMesh = PGS_Meshing.areaMerge(mesh, areaThreshold);
		assertTrue(PGS_Conversion.getChildren(mergedMesh).stream().allMatch(f -> PGS_ShapePredicates.area(f) >= areaThreshold));
		assertTrue(faces.size() >= mergedMesh.getChildCount());
		assertEquals(PGS_ShapePredicates.area(mesh), PGS_ShapePredicates.area(mergedMesh), 1e-6);
	}

}
