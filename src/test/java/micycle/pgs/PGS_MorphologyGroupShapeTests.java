package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;

/**
 * Tests to determine which methods from {@link micycle.pgs.PGS_Morphology
 * PGS_Morphology} natively support GROUP PShape (MultiPolygon) inputs (or at
 * least handle them without failing). These tests do not validate the
 * particular output of the method in question.
 */
class PGS_MorphologyGroupShapeTests {

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
		a.vertex(10, 10);
		a.vertex(0, 10);
		a.endShape(PConstants.CLOSE);

		final PShape b = new PShape(PShape.GEOMETRY);
		b.beginShape();
		b.vertex(70, 70);
		b.vertex(710, 70);
		b.vertex(710, 710);
		b.vertex(70, 710);
		b.endShape(PConstants.CLOSE);

		GROUP_SHAPE = new PShape(PConstants.GROUP);
		GROUP_SHAPE.setKind(PConstants.GROUP);
		GROUP_SHAPE.addChild(a);
		GROUP_SHAPE.addChild(b);
	}

	@Test
	void testBuffer() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.buffer(GROUP_SHAPE, -1);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testChaikinCut() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.chaikinCut(GROUP_SHAPE, 0.5, 2);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testErosionDilation() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.erosionDilation(GROUP_SHAPE, 0);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testFieldWarp() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.fieldWarp(GROUP_SHAPE, 10, 1, false);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testMinkDifference() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		final PShape mink = new PShape(PShape.PATH);
		mink.beginShape();
		mink.vertex(0, 0);
		mink.vertex(5, 0);
		mink.vertex(5, 5);
		mink.vertex(0, 5);
		mink.endShape(PConstants.CLOSE);

		PShape out = PGS_Morphology.minkDifference(GROUP_SHAPE, mink);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testMinkSum() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		final PShape mink = new PShape(PShape.PATH);
		mink.beginShape();
		mink.vertex(0, 0);
		mink.vertex(5, 0);
		mink.vertex(5, 5);
		mink.vertex(0, 5);
		mink.endShape(PConstants.CLOSE);

		PShape out = PGS_Morphology.minkSum(GROUP_SHAPE, mink);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testRadialWarp() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.radialWarp(GROUP_SHAPE, 10, 1, false);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testRound() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.round(GROUP_SHAPE, 0.5);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testSimplify() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.simplify(GROUP_SHAPE, 1);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testSimplifyTopology() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.simplifyTopology(GROUP_SHAPE, 1);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testSimplifyVW() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.simplifyVW(GROUP_SHAPE, 1);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testSmooth() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.smooth(GROUP_SHAPE, 0.5);
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testSmoothGaussian() {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);
		PShape out = PGS_Morphology.smoothGaussian(GROUP_SHAPE, 10);
		assertEquals(2, out.getChildCount());
	}

}
