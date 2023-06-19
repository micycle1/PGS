package micycle.pgs.color;

/**
 * Provides static methods for Processing color generation (that doesn't require
 * a PApplet instance).
 * 
 * @author Michael Carleton
 *
 */
public class RGB {

	/** RGB (0, 0, 0) */
	public static final int BLACK = composeColor(0, 0, 0);
	/** RGB (255, 255, 255) */
	public static final int WHITE = composeColor(255, 255, 255);
	/** RGB (237, 50, 162) */
	public static final int PINK = composeColor(237, 50, 162);
	/** RGBA (237, 50, 162, 128) */
	public static final int HALF_PINK = micycle.pgs.color.RGB.setAlpha(PINK, 128);
	/** RGB (255, 0, 0) */
	public static final int RED = composeColor(255, 0, 0);
	/** RGB (0, 255, 0) */
	public static final int GREEN = composeColor(0, 255, 0);
	/** RGB (0, 0, 255) */
	public static final int BLUE = composeColor(0, 0, 255);
	/** RGB (255, 255, 0) */
	public static final int YELLOW = composeColor(255, 255, 0);
	/** RGB (0, 255, 255) */
	public static final int CYAN = composeColor(0, 255, 255);
	/** RGB (255, 0, 255) */
	public static final int FUCHSIA = composeColor(255, 0, 255);

	private static final float INV_255 = 1f / 255f; // used to normalise RGB values to 0...1

	private RGB() {
	}

	/**
	 * Composes an integer value that represents a color in RGB format.
	 * 
	 * @param red   ∈[0, 255]
	 * @param green ∈[0, 255]
	 * @param blue  ∈[0, 255]
	 * @param alpha ∈[0, 255] (where 0 is transparent; 255 is opaque)
	 * @return the integer representation of the color in RGB format
	 */
	public static int composeColor(final int red, final int green, final int blue, final int alpha) {
		return alpha << 24 | red << 16 | green << 8 | blue;
	}

	/**
	 * Composes an integer value that represents a color in RGB format.
	 * 
	 * @param red   ∈[0, 255]
	 * @param green ∈[0, 255]
	 * @param blue  ∈[0, 255]
	 * @return the integer representation of the color in RGB format
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
	 * Decompose color integer (ARGB) into its 3 separate RGB components (0...255)
	 * 
	 * @param clr
	 * @return [R,G,B] 0...255
	 */
	public static int[] decomposeclrRGB(int clr) {
		int[] out = new int[3];
		out[0] = (clr >> 16 & 0xff);
		out[1] = (clr >> 8 & 0xff);
		out[2] = (clr & 0xff);
		return out;
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

	/**
	 * Converts hex color strings to Processing integer colors (RRGGBB).
	 * 
	 * @param colors list of hex e.g "021d34" or "#FF6F91"
	 * @return
	 */
	public static int[] hexToColor(String[] colors) {
		int[] out = new int[colors.length];
		for (int i = 0; i < colors.length; i++) {
			if (colors[i].charAt(0) == '#') { // handle leading '#'
				colors[i] = colors[i].substring(1);
			}
			out[i] = -16777216 + (int) (Long.parseLong(colors[i], 16));
		}
		return out;
	}

}
