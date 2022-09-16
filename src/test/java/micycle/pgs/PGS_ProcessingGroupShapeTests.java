package micycle.pgs;

import static org.junit.Assume.assumeTrue;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;

/**
 * Tests to determine which methods from {@link micycle.pgs.PGS_Processing
 * PGS_Processing} natively support GROUP PShape (MultiPolygon) inputs (or at
 * least handle them without failing). These tests do not validate the
 * particular output of the method in question.
 */
class PGS_ProcessingGroupShapeTests {

	private PShape GROUP_SHAPE;

	@BeforeEach
	/**
	 * Recreate the test shape before each test case in case some methods mutate the
	 * shape.
	 */
	void prepareGroupShape() {
		final PShape a = new PShape(PShape.GEOMETRY);
		a.beginShape();
		a.vertex(0, 0);
		a.vertex(10, 0);
		a.vertex(0, 10);
		a.vertex(10, 10);
		a.endShape(PConstants.CLOSE);

		final PShape b = new PShape(PShape.GEOMETRY);
		b.beginShape();
		b.vertex(70, 70);
		b.vertex(710, 70);
		b.vertex(70, 710);
		b.vertex(710, 710);
		b.endShape(PConstants.CLOSE);

		GROUP_SHAPE = new PShape(PShape.GROUP);
		GROUP_SHAPE.setKind(PShape.GROUP);
		GROUP_SHAPE.addChild(a);
		GROUP_SHAPE.addChild(b);
	}

	@Test
	void test_PGS_Processing_densify() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Processing.densify(GROUP_SHAPE, 1);
		assertEquals(2, out.getChildCount());
	}
	
	@Test
	void test_PGS_Processing_removeSmallHoles() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Processing.removeSmallHoles(GROUP_SHAPE, 10);
		assertEquals(2, out.getChildCount());
	}
	
	@Test
	void test_PGS_Processing_convexPartition() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Processing.convexPartition(GROUP_SHAPE);
		assertEquals(2, out.getChildCount());
	}

}
