package micycle.pgs;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.Stream;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import micycle.pgs.PGS_Coloring.ColoringAlgorithm;
import processing.core.PConstants;
import processing.core.PShape;

/**
 * Tests whether {@link micycle.pgs.PGS_Coloring PGS_Coloring} colors conforming
 * and non-conforming meshes as expected.
 * 
 * @author Michael Carleton
 *
 */
class PGS_ColoringTests {

	private PShape GROUP_SHAPE;

	@BeforeEach
	/**
	 * Recreate the test shape (two adjacent triangles) before each test case in
	 * case some methods mutate the shape.
	 */
	void prepareGroupShape() {
		final PShape a = new PShape(PShape.GEOMETRY);
		a.beginShape();
		a.vertex(0, 0);
		a.vertex(10, 0);
		a.vertex(0, 10);
		a.endShape(PConstants.CLOSE);

		final PShape b = new PShape(PShape.GEOMETRY);
		b.beginShape();
		b.vertex(0, 10);
		b.vertex(10, 10);
		b.vertex(10, 0);
		b.endShape(PConstants.CLOSE);

		GROUP_SHAPE = new PShape(PShape.GROUP);
		GROUP_SHAPE.setKind(PShape.GROUP);
		GROUP_SHAPE.addChild(a);
		GROUP_SHAPE.addChild(b);
	}

	@Test
	void testMeshColoring() {
		Map<PShape, Integer> coloring = PGS_Coloring.colorMesh(GROUP_SHAPE, ColoringAlgorithm.RLF);
		assertEquals(2, coloring.size());
		List<Integer> colorClasses = new ArrayList<>(coloring.values());
		assertNotSame(colorClasses.get(0), colorClasses.get(1));
	}

	@Test
	void testNonMeshColoring() {
		/*
		 * Add a non-conforming face (this face is flush with one triangle but the flush
		 * / overlapping edge is not identical).
		 */
		final PShape c = new PShape(PShape.GEOMETRY);
		c.beginShape();
		c.vertex(3, 10);
		c.vertex(6, 10);
		c.vertex(6, 15);
		c.vertex(3, 15);
		c.endShape(PConstants.CLOSE);

		GROUP_SHAPE.addChild(c);

		Map<PShape, Integer> coloring = PGS_Coloring.colorNonMesh(GROUP_SHAPE, ColoringAlgorithm.RLF);
		assertEquals(3, coloring.size());

		Supplier<Stream<PShape>> faces = () -> coloring.keySet().stream();
		// find the faces by vertex count
		Integer aC = coloring.get(faces.get().filter(f -> f.getVertexCount() == 3).findFirst().get());
		Integer bC = coloring.get(faces.get().filter(f -> f.getVertexCount() > 4).findFirst().get());
		Integer cC = coloring.get(faces.get().filter(f -> f.getVertexCount() == 4).findFirst().get());

		assertNotNull(aC);
		assertNotNull(bC);
		assertNotNull(cC);
		assertNotSame(aC, bC);
		assertNotSame(bC, cC);
	}

}
