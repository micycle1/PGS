package micycle.pts.color;

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

	private static final float INV_255 = 1f / 255f; // used to normalise RGB values to 0...1

	// Normal: src.
	// Lighten: max(src, dst).
	// Darken: min(src, dst).

	// Add: min(1.0, src + dst).
	// Subtract: max(0.0, src - dst).
	// Multiply: (src * dst).
	// Screen: 1 - (1 - dst) * (1 - src).
	// Average: src + (dst - src) * 0.5.
	// Difference: abs(src - dst).
	// Exclusion: src - 2 * src * dst + dst.

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

	/**
	 * Decompose and pre-multiply alpha
	 * 
	 * @param clr
	 * @return RGBA
	 */
	private static float[] decomposeclr(int clr) {
		final float alpha = (clr >> 24 & 0xff) == 255 ? 1 : (clr >> 24 & 0xff) * INV_255;
		return new float[] { (clr >> 16 & 0xff) * INV_255 * alpha, (clr >> 8 & 0xff) * INV_255 * alpha,
				(clr & 0xff) * INV_255 * alpha, alpha };
	}

	/**
	 * Compose a 32 bit sARGB int from float[] 0...1
	 * 
	 * @param in RGBA
	 * @return
	 */
	private static int composeclr(float[] RGBA) {
		return (int) (RGBA[3] * 255) << 24 | (int) (RGBA[0] * 255) << 16 | (int) (RGBA[1] * 255) << 8
				| (int) (RGBA[2] * 255);
	}

}
