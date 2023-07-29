package micycle.pgs.color;

import com.scrtwpns.Mixbox;

import net.jafama.FastMath;

/**
 * Provides static methods for Processing color generation (that doesn't require
 * a PApplet instance).
 * 
 * @author Michael Carleton
 *
 */
public class ColorUtils {

	private static final float INV_255 = 1f / 255f; // used to normalise ColorUtils values to 0...1

	private ColorUtils() {
	}

	/**
	 * Composes an integer value that represents a color in ColorUtils format.
	 * 
	 * @param red   ∈[0, 255]
	 * @param green ∈[0, 255]
	 * @param blue  ∈[0, 255]
	 * @param alpha ∈[0, 255] (where 0 is transparent; 255 is opaque)
	 * @return the integer representation of the color in ColorUtils format
	 */
	public static int composeColor(final int red, final int green, final int blue, final int alpha) {
		return alpha << 24 | red << 16 | green << 8 | blue;
	}

	/**
	 * Composes an integer value that represents a color in ColorUtils format.
	 * 
	 * @param red   ∈[0, 255]
	 * @param green ∈[0, 255]
	 * @param blue  ∈[0, 255]
	 * @return the integer representation of the color in ColorUtils format
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
	 * Decompose color integer (ARGB) into its 3 separate ColorUtils components
	 * (0...255)
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

	static int composeclr(double[] RGBA) {
		return (int) (RGBA[3] * 255) << 24 | (int) (RGBA[0] * 255) << 16 | (int) (RGBA[1] * 255) << 8 | (int) (RGBA[2] * 255);
	}
	
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

	/**
	 * Mixes/blends two colors using natural color mixing.
	 * <p>
	 * It produces saturated gradients with hue shifts and natural secondary colors
	 * during blending. For instance, yellow and blue make green.
	 * 
	 * @return the new mixed color
	 */
	public static int pigmentMix(int colorA, int colorB, float t) {
		return Mixbox.lerp(colorA, colorB, t);
	}

	/**
	 * Produces a smooth hue-cycling rainbow.
	 * 
	 * @param t 0...1]
	 * @return RGB color integer
	 */
	public static int sinebow(double t) {
		if (t > 1) {
			t %= 1;
		}
		t = 0.5f - t;
		double[] cols = new double[] { sin2(t), sin2(t + (1 / 3d)), sin2(t + (2 / 3d)), 1 };
		return composeclr(cols);
	}

	private static double sin2(double t) {
		double z = FastMath.sin(Math.PI * t);
		return z * z;
	}

}
