package micycle.pts.color;

public class RGB {

	/**
	 * @param red   ∈[0, 255]
	 * @param green ∈[0, 255]
	 * @param blue  ∈[0, 255]
	 * @param alpha ∈[0, 255]
	 * @return
	 */
	public static int composeclr(int red, int green, int blue, int alpha) {
		return alpha << 24 | red << 16 | green << 8 | blue;
	}

}
