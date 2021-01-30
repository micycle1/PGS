package micycle.pts.color;

/**
 * Various methods to composite two RGB colors.
 * 
 * Porter and Duff (1984) define twelve formulas for combining (compositing) two
 * RGBA colors(30). In the formulas below, it is assumed that the two colors and
 * the output are in the 0-1 format and have been premultiplied (that is, their
 * red, green, and blue components have been multiplied beforehand by their
 * alpha component). Given src, the source RGBA color, and dst, the destination
 * RGBA color, the Porter–Duff formulas are as follows.
 * 
 * @author MCarleton
 *
 */
public class Composite {

	private static final float INV_255 = 1f / 255f; // used to normalise RGB values

//	Source Over: [src[0]-dst[0]*(src[3] - 1), src[1]-dst[1]*(src[3] - 1), src[2]-dst[2]*(src[3] - 1), src[3]-dst[3]*(src[3] - 1)].
//	Source In: [dst[3]*src[0], dst[3]*src[1], dst[3]*src[2], dst[3]*src[3]].
//	Source Held Out: [src[0]*(1 - dst[3]), src[1]*(1 - dst[3]), src[2]*(1 - dst[3]), src[3]*(1 - dst[3])].
//	Source Atop: [dst[0]*src[3] - src[0]*(dst[3] - 1), dst[1]*src[3] - src[1]*(dst[3] - 1), dst[2]*src[3] - src[2]*(dst[3] - 1), src[3]].
//	Destination Over: [dst[0] - src[0]*(dst[3] - 1), dst[1] - src[1]*(dst[3] - 1), dst[2] - src[2]*(dst[3] - 1), dst[3] - src[3]*(dst[3] - 1)].
//	Destination In: [dst[0]*src[3], dst[1]*src[3], dst[2]*src[3], dst[3]*src[3]]. Uses the destination color/alpha with the source alpha as the "mask".
//	Destination Held Out: [dst[0]*(1 - src[3]), dst[1]*(1 - src[3]), dst[2]*(1 - src[3]), dst[3]*(1 - src[3])].
//	Destination Atop: [dst[3]*src[0] - dst[0]*(src[3] - 1), dst[3]*src[1] - dst[1]*(src[3] - 1), dst[3]*src[2] - dst[2]*(src[3] - 1), dst[3]].
//	Xor: [-dst[3]*src[0] - dst[0]*src[3] + dst[0] + src[0], -dst[3]*src[1] - dst[1]*src[3] + dst[1] + src[1], -dst[3]*src[2] - dst[2]*src[3] + dst[2] + src[2], -2*dst[3]*src[3] + dst[3] + src[3]].

	public static int xor(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(xor(decomposedA, decomposedB));
	}

	public static int sourceOver(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(sourceOver(decomposedA, decomposedB));
	}

	public static int sourceHeldOut(int colorA, int colorB) {
		float[] decomposedA = decomposeclr(colorA);
		float[] decomposedB = decomposeclr(colorB);
		return composeclr(sourceHeldOut(decomposedA, decomposedB));
	}

	private static float[] xor(float[] src, float[] dst) {
		return new float[] { -dst[3] * src[0] - dst[0] * src[3] + dst[0] + src[0],
				-dst[3] * src[1] - dst[1] * src[3] + dst[1] + src[1],
				-dst[3] * src[2] - dst[2] * src[3] + dst[2] + src[2], -2 * dst[3] * src[3] + dst[3] + src[3] };
	}

	private static float[] sourceOver(float[] src, float[] dst) {
		return new float[] { src[0] - dst[0] * (src[3] - 1), src[1] - dst[1] * (src[3] - 1),
				src[2] - dst[2] * (src[3] - 1), src[3] - dst[3] * (src[3] - 1) };
	}

	private static float[] sourceHeldOut(float[] src, float[] dst) {
		return new float[] { src[0] * (1 - dst[3]), src[1] * (1 - dst[3]), src[2] * (1 - dst[3]),
				src[3] * (1 - dst[3]) };
	}

	/**
	 * Decompose and pre-multiply alpha
	 * 
	 * @param clr
	 * @return premultiplied clr
	 */
	public static float[] decomposeclr(int clr) {
		final float alpha = (clr >> 24 & 0xff) == 255 ? 1 : (clr >> 24 & 0xff) * INV_255;
		return new float[] { (clr >> 16 & 0xff) * INV_255 * alpha, (clr >> 8 & 0xff) * INV_255 * alpha,
				(clr & 0xff) * INV_255 * alpha, alpha };
	}

	/**
	 * Compose a 32 bit sARGB int from float[] 0...1
	 * 
	 * @param in
	 * @return
	 */
	public static int composeclr(float[] RGBA) {
		return (int) (RGBA[3] * 255) << 24 | (int) (RGBA[0] * 255) << 16 | (int) (RGBA[1] * 255) << 8
				| (int) (RGBA[2] * 255);
	}

}
