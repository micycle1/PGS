package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.util.List;
import java.util.function.Function;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Tests to determine which methods from {@link micycle.pgs.PGS_Morphology
 * PGS_Morphology} natively support GROUP PShape (MultiPolygon) inputs (or at
 * least handle them without failing). These tests do not validate the
 * particular output of the method in question.
 */
class PGS_MorphologyGroupShapeTests {

	private PShape GROUP_SHAPE;

	// Reused auxiliary inputs for methods that need additional arguments
	private PShape minkShape;
	private PShape toShape;
	private PVector pinchPoint;
	private List<PVector> arapFrom;
	private List<PVector> arapTo;

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

		// Minkowski operand shape
		minkShape = new PShape(PShape.PATH);
		minkShape.beginShape();
		minkShape.vertex(0, 0);
		minkShape.vertex(5, 0);
		minkShape.vertex(5, 5);
		minkShape.vertex(0, 5);
		minkShape.endShape(PConstants.CLOSE);

		// "to" shape for interpolate() (same structure: GROUP with 2 children)
		final PShape a2 = new PShape(PShape.GEOMETRY);
		a2.beginShape();
		a2.vertex(2, 2);
		a2.vertex(12, 2);
		a2.vertex(12, 12);
		a2.vertex(2, 12);
		a2.endShape(PConstants.CLOSE);

		final PShape b2 = new PShape(PShape.GEOMETRY);
		b2.beginShape();
		b2.vertex(80, 80);
		b2.vertex(720, 80);
		b2.vertex(720, 720);
		b2.vertex(80, 720);
		b2.endShape(PConstants.CLOSE);

		toShape = new PShape(PConstants.GROUP);
		toShape.setKind(PConstants.GROUP);
		toShape.addChild(a2);
		toShape.addChild(b2);

		pinchPoint = new PVector(5, 5);

		// ARAP control points (simple, deterministic)
		arapFrom = List.of(new PVector(0, 0), new PVector(10, 0), new PVector(10, 10), new PVector(0, 10));
		arapTo = List.of(new PVector(0, 0), new PVector(12, -2), new PVector(9, 13), new PVector(-1, 11));
	}

	private void assertGroupInGroupOut(Function<PShape, PShape> op) {
		assumeTrue(GROUP_SHAPE.getChildCount() == 2);

		PShape out = op.apply(GROUP_SHAPE);

		// For "supports GROUP" tests, we expect the output to preserve
		// multipolygon-ness.
		assertEquals(2, out.getChildCount());
	}

	@Test
	void testBuffer() {
		assertGroupInGroupOut(s -> PGS_Morphology.buffer(s, -1));
	}

	@Test
	void testVariableBuffer() {
		assertGroupInGroupOut(s -> PGS_Morphology.variableBuffer(s, -2, 2));
	}

	@Test
	void testVariableBufferCallback() {
		assertGroupInGroupOut(s -> PGS_Morphology.variableBuffer(s, (coord, t) -> t * 10 + 1));
	}

	@Test
	void testNormalisedErosion() {
		assertGroupInGroupOut(s -> PGS_Morphology.normalisedErosion(s, 0.1));
	}

	@Test
	void testErosionDilation() {
		assertGroupInGroupOut(s -> PGS_Morphology.erosionDilation(s, 1));
	}

	@Test
	void testDilationErosion() {
		assertGroupInGroupOut(s -> PGS_Morphology.dilationErosion(s, 1));
	}

	@Test
	void testSimplify() {
		assertGroupInGroupOut(s -> PGS_Morphology.simplify(s, 1));
	}

	@Test
	void testSimplifyVW() {
		assertGroupInGroupOut(s -> PGS_Morphology.simplifyVW(s, 1));
	}

	@Test
	void testSimplifyTopology() {
		assertGroupInGroupOut(s -> PGS_Morphology.simplifyTopology(s, 1));
	}

	@Test
	void testSimplifyDCE() {
		assertGroupInGroupOut(s -> PGS_Morphology.simplifyDCE(s, 0.25));
	}

	@Test
	void testSimplifyHobby() {
		assertGroupInGroupOut(s -> PGS_Morphology.simplifyHobby(s, 1));
	}

	@Test
	void testMinkSum() {
		assertGroupInGroupOut(s -> PGS_Morphology.minkSum(s, minkShape));
	}

	@Test
	void testMinkDifference() {
		assertGroupInGroupOut(s -> PGS_Morphology.minkDifference(s, minkShape));
	}

	@Test
	void testSmooth() {
		assertGroupInGroupOut(s -> PGS_Morphology.smooth(s, 0.5));
	}

	@Test
	void testSmoothGaussian() {
		assertGroupInGroupOut(s -> PGS_Morphology.smoothGaussian(s, 10));
	}

	@Test
	void testSmoothGaussianNormalised() {
		assertGroupInGroupOut(s -> PGS_Morphology.smoothGaussianNormalised(s, 0.25));
	}

	@Test
	void testSmoothEllipticFourier() {
		assertGroupInGroupOut(s -> PGS_Morphology.smoothEllipticFourier(s, 12));
	}

	@Test
	void testSmoothLaneRiesenfeld() {
		assertGroupInGroupOut(s -> PGS_Morphology.smoothLaneRiesenfeld(s, 3, 2, 0.5));
	}

	@Test
	void testRound() {
		assertGroupInGroupOut(s -> PGS_Morphology.round(s, 0.5));
	}

	@Test
	void testChaikinCut() {
		assertGroupInGroupOut(s -> PGS_Morphology.chaikinCut(s, 0.5, 2));
	}

	@Test
	void testRadialWarp() {
		assertGroupInGroupOut(s -> PGS_Morphology.radialWarp(s, 10, 1, false));
	}

	@Test
	void testSineWarp() {
		assertGroupInGroupOut(s -> PGS_Morphology.sineWarp(s, 5, 2, 0));
	}

	@Test
	void testFieldWarp() {
		assertGroupInGroupOut(s -> PGS_Morphology.fieldWarp(s, 10, 1, 0.0, false, 1337L));
	}

	@Test
	void testPinchWarp() {
		assertGroupInGroupOut(s -> PGS_Morphology.pinchWarp(s, pinchPoint, 0.75));
	}

	@Test
	@Disabled // returns the input unchanged if not polygonal
	void testInterpolate() {
		assertGroupInGroupOut(s -> PGS_Morphology.interpolate(s, toShape, 0.5));
	}

	@Test
	void testArapDeform() {
		// NOTE doesn't support GROUP
		assertThrows(Exception.class, () -> PGS_Morphology.arapDeform(GROUP_SHAPE, arapFrom, arapTo));
	}

	@Test
	void testRegularise() {
		assertGroupInGroupOut(s -> PGS_Morphology.regularise(s, 0.5));
	}
	
	@Test
	void testSmoothBezierFit() {
		assertGroupInGroupOut(s -> PGS_Morphology.smoothBezierFit(s, 1));
	}
	
	@Test
	void testReducePrecision() {
		assertGroupInGroupOut(s -> PGS_Morphology.reducePrecision(s, 1));
	}
}