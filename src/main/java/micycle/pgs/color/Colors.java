package micycle.pgs.color;

public final class Colors {

	/** ColorUtils (0, 0, 0) */
	public static final int BLACK = ColorUtils.composeColor(0, 0, 0);
	/** ColorUtils (255, 255, 255) */
	public static final int WHITE = ColorUtils.composeColor(255, 255, 255);
	/** ColorUtils (237, 50, 162) */
	public static final int PINK = ColorUtils.composeColor(237, 50, 162);
	/** RGBA (237, 50, 162, 128) */
	public static final int HALF_PINK = ColorUtils.setAlpha(PINK, 128);
	/** ColorUtils (255, 0, 0) */
	public static final int RED = ColorUtils.composeColor(255, 0, 0);
	/** ColorUtils (0, 255, 0) */
	public static final int GREEN = ColorUtils.composeColor(0, 255, 0);
	/** ColorUtils (0, 0, 255) */
	public static final int BLUE = ColorUtils.composeColor(0, 0, 255);
	/** ColorUtils (255, 255, 0) */
	public static final int YELLOW = ColorUtils.composeColor(255, 255, 0);
	/** ColorUtils (0, 255, 255) */
	public static final int CYAN = ColorUtils.composeColor(0, 255, 255);
	/** ColorUtils (255, 0, 255) */
	public static final int FUCHSIA = ColorUtils.composeColor(255, 0, 255);

	private Colors() {
	}

}
