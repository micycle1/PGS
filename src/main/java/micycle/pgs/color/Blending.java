package micycle.pgs.color;

import static micycle.pgs.color.RGB.composeclr;
import static micycle.pgs.color.RGB.decomposeclr;

/**
 * Color blending for Processing colors (32bit ARGB integers).
 * <p>
 * Blend modes take two multicomponent colors, namely a source color and a
 * destination color, and blend them to create a new color.
 * 
 * <p>
 * Methods operate component-wise in the RGB color space.
 * 
 * @author Michael Carleton
 *
 */
public class Blending {
	
	private Blending() {
	}

	public static int subtract(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(subtract(decomposedA, decomposedB));
	}

	private static float[] subtract(float[] src, float[] dst) {
		return new float[] { Math.max(0, src[0] - dst[0]), Math.max(0, src[1] - dst[1]), Math.max(0, src[2] - dst[2]),
				Math.max(0, src[3] - dst[3]) };
	}

	public static int add(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(add(decomposedA, decomposedB));
	}

	private static float[] add(float[] src, float[] dst) {
		return new float[] { Math.min(1, src[0] + dst[0]), Math.min(1, src[1] + dst[1]), Math.min(1, src[2] + dst[2]),
				Math.min(1, src[3] + dst[3]) };
	}

	public static int difference(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(difference(decomposedA, decomposedB));
	}

	private static float[] difference(float[] src, float[] dst) {
		return new float[] { Math.abs(src[0] - dst[0]), Math.abs(src[1] - dst[1]), Math.abs(src[2] - dst[2]),
				Math.abs(src[3] - dst[3]) };
	}

	public static int multiply(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(multiply(decomposedA, decomposedB));
	}

	private static float[] multiply(float[] src, float[] dst) {
		return new float[] { src[0] * dst[0], src[1] * dst[1], src[2] * dst[2], src[3] * dst[3] };
	}

	public static int screen(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(screen(decomposedA, decomposedB));
	}

	private static float[] screen(float[] src, float[] dst) {
		return new float[] { 1 - (1 - dst[0]) * (1 - src[0]), 1 - (1 - dst[1]) * (1 - src[1]),
				1 - (1 - dst[2]) * (1 - src[2]), 1 - (1 - dst[3]) * (1 - src[3]) };
	}

	public static int average(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(average(decomposedA, decomposedB));
	}

	private static float[] average(float[] src, float[] dst) {
		return new float[] { src[0] + (dst[0] - src[0]) * 0.5f, src[1] + (dst[1] - src[1]) * 0.5f,
				src[2] + (dst[2] - src[2]) * 0.5f, src[3] + (dst[3] - src[3]) * 0.5f };
	}

	public static int exclusion(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(exclusion(decomposedA, decomposedB));
	}

	private static float[] exclusion(float[] src, float[] dst) {
		return new float[] { src[0] - 2 * src[0] * dst[0] + dst[1], src[1] - 2 * src[1] * dst[1] + dst[1],
				src[2] - 2 * src[2] * dst[2] + dst[2], src[3] - 2 * src[3] * dst[3] + dst[3] };
	}

}
