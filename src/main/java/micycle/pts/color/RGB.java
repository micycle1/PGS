package micycle.pts.color;

public class RGB {

	public static final int BLACK = composeColor(0, 0, 0);
	public static final int WHITE = composeColor(255, 255, 255);
	public static final int PINK = composeColor(237, 50, 162);

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

}
