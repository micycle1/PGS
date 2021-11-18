package micycle.pgs.color;

/**
 * Provides static methods for Processing color generation (that doesn't require
 * a PApplet instance).
 * 
 * @author Michael Carleton
 *
 */
public class RGB {

	public static final int BLACK = composeColor(0, 0, 0);
	public static final int WHITE = composeColor(255, 255, 255);
	public static final int PINK = composeColor(237, 50, 162);

	private static final float INV_255 = 1f / 255f; // used to normalise RGB values to 0...1

	private RGB() {
	}

	/**
	 * @param red   ∈[0, 255]
	 * @param green ∈[0, 255]
	 * @param blue  ∈[0, 255]
	 * @param alpha ∈[0, 255] (where 0 is transparent; 255 is opaque)
	 * @return
	 */
	public static int composeColor(final int red, final int green, final int blue, final int alpha) {
		return alpha << 24 | red << 16 | green << 8 | blue;
	}

	/**
	 * @param red   ∈[0, 255]
	 * @param green ∈[0, 255]
	 * @param blue  ∈[0, 255]
	 * @return
	 */
	public static int composeColor(final int red, final int green, final int blue) {
		return -16777216 | red << 16 | green << 8 | blue;
	}

	/**
	 * @param grey  ∈[0, 255]
	 * @param alpha ∈[0, 255] (where 0 is transparent; 255 is opaque)
	 * @return
	 */
	public static int composeColor(final int grey, final int alpha) {
		return alpha << 24 | grey << 16 | grey << 8 | grey;
	}

	/**
	 * 
	 * @param color
	 * @param alpha ∈[0, 255] (where 0 is transparent; 255 is opaque)
	 * @return
	 */
	public static int setAlpha(int color, int alpha) {
		return (color & 16777215) | alpha << 24;
	}

	/**
	 * Compose a 32 bit sARGB int from float[] 0...1
	 * 
	 * @param in
	 * @return
	 */
	static int composeclr(float[] RGBA) {
		return (int) (RGBA[3] * 255) << 24 | (int) (RGBA[0] * 255) << 16 | (int) (RGBA[1] * 255) << 8 | (int) (RGBA[2] * 255);
	}

	/**
	 * Decompose and pre-multiply alpha
	 * 
	 * @param clr
	 * @return premultiplied clr
	 */
	static float[] decomposeclr(int clr) {
		final float alpha = (clr >> 24 & 0xff) == 255 ? 1 : (clr >> 24 & 0xff) * INV_255;
		return new float[] { (clr >> 16 & 0xff) * INV_255 * alpha, (clr >> 8 & 0xff) * INV_255 * alpha, (clr & 0xff) * INV_255 * alpha,
				alpha };
	}

}
