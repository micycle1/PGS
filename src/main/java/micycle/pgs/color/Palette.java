package micycle.pgs.color;

/**
 * A selection of color palettes.
 * <p>
 * All palettes were taken from (and are visible at)
 * <a href="https://ronikaufman.github.io/color_pals/">Roni Kaufman Color
 * Pals</a>.
 * 
 * @author Michael Carleton
 *
 */
public enum Palette {

	// @formatter:off
		_1115(new String[] { "#ffe03d", "#fe4830", "#d33033", "#6d358a", "#1c509e", "#00953c" }),
		ARCADE(new String[] { "#021d34", "#228fca", "#dcedf0" }),
		ARCS(new String[] { "#ff3232", "#ff9932", "#ffff32", "#32ff32", "#32ffff", "#3232ff", "#ff32ff" }),
		AVANTGARDE(new String[] { "#f398c3", "#cf3895", "#a0d28d", "#06b4b0", "#fef45f", "#fed000", "#0c163f" }),
		BASE(new String[] { "#f2eb8a", "#fed000", "#fc8405", "#ed361a", "#e2f0f3", "#b3dce0", "#4464a1", "#203051", "#ffc5c7", "#f398c3", "#cf3895", "#6d358a", "#06b4b0", "#4b8a5f" }),
		BRAIN(new String[] { "#ee726b", "#ffc5c7", "#fef9c6" }),
		CHEDDAR(new String[] { "#ff7b00", "#ff8800", "#ff9500", "#ffa200", "#ffaa00", "#ffb700", "#ffc300" }),
		CLOUDY(new String[] { "#044e9e", "#6190d3", "#fcf7ed", "#fcd494", "#f4b804" }),
		COMMIT(new String[] { "#563d7c", "#0096d8", "#f4e361", "#f24679" }),
		CONNECTION(new String[] { "#f24358", "#f2a643", "#f2e343", "#43f278", "#43a0f2", "#c343f2" }),
		CUBE(new String[] { "#fde200", "#ef2626", "#5600ae", "#713df5" }),
		DATA(new String[] { "#f9f9f9", "#ff0000", "#0000ff" }),
		EDIT(new String[] { "#1c1c1c", "#f47a9d", "#f4ea7a", "#f2f2f2" }),
		ESCAPE(new String[] { "#f3e17e", "#dd483c", "#4b8a5f", "#0d150b", "#faf8e2" }),
		EVENING(new String[] { "#ff4f19", "#15084d", "#5ce6e6" }),
		FANS(new String[] { "#000000", "#83a6bc", "#faece1", "#bab1a8" }),
		FLAT(new String[] { "#203051", "#4464a1", "#62b6de", "#b3dce0", "#e2f0f3" }),
		FREAK(new String[] { "#07224f", "#ed361a", "#fc8405", "#f7c72a" }),
		HER(new String[] { "#cd1440", "#fafafa" }),
		HOMAGE(new String[] { "#fef9c6", "#ffcc4d", "#f5b800", "#56a1c4", "#4464a1", "#ee726b", "#df5f50", "#5a3034" }),
		JELLY(new String[] { "#102340", "#fe01ec", "#8a07da" }),
		LEFT(new String[] { "#070213", "#1c0541", "#300560", "#3e137f", "#412287", "#3b3d8c", "#e6f41c" }),
		LESSON(new String[] { "#ed225d", "#3caf65", "#0d40bf", "#f5b800", "#2a2a2a" }),
		LIGHT(new String[] { "#00b2b4", "#fdd35b", "#3b4696", "#f4737d", "#333333" }),
		LIZARD(new String[] { "#fff247", "#47e9ff", "#8447ff", "#ff47ce" }),
		MARBLE(new String[] { "#218ad4", "#76df55", "#0b1435", "#ffffff" }),
		META(new String[] { "#226699", "#dd2e44", "#ffcc4d" }),
		MONDRIAN(new String[] { "#0a0a0a", "#f7f3f2", "#0077e1", "#f5d216", "#fc3503" }),
		MORNING(new String[] { "#ffd919", "#262104", "#fffbe6" }),
		NORTH(new String[] { "#dc060e", "#ffd400", "#0064b0", "#001a5b", "#ffffff" }),
		OPTICAL(new String[] { "#f2eb8a", "#fed000", "#fc8405", "#ed361a", "#e2f0f3", "#b3dce0", "#4464a1", "#203051", "#ffc5c7", "#f398c3", "#cf3895", "#6d358a", "#06b4b0", "#4b8a5f" }),
		PAINT(new String[] { "#b30000", "#e6cf00", "#1283b3", "#fafafa", "#0a0a0a" }),
		PIE(new String[] { "#fdcc01", "#f60001", "#001ea6", "#007d45", "#040a07", "#f1ebf7" }),
		PLAY(new String[] { "#e20404", "#fcd202", "#ffffff", "#000000" }),
		REPETITION(new String[] { "#2c2060", "#4bd3e5", "#fffbe6", "#ffd919", "#ff4f19" }),
		SCIENCE(new String[] { "#ffe819", "#000000" }),
		SHEETS(new String[] { "#025760", "#7fd0d3", "#b3ead7", "#eff7F5", "#fae9c1", "#f8ca75", "#e88d44" }),
		SPLASH(new String[] { "#32312e", "#795330", "#c7ae82", "#f5f2e3" }),
		STARS(new String[] { "#0a1966", "#ffef0d", "#fafafa" }),
		TEK(new String[] { "#fcf7ed", "#c6c3ba", "#a5a29b", "#fcde97", "#fccf64" }),
		TEN(new String[] { "#fffbe6", "#050505", "#abcd5e", "#29ac9f", "#14976b", "#b3dce0", "#62b6de", "#2b67af", "#f589a3", "#ef562f", "#fc8405", "#f9d531" }),
		TROLL(new String[] { "#294984", "#6ca0a7", "#ffc789", "#df5f50", "#5a3034", "#fff1dd" }),
		TROPICAL(new String[] { "#3f7373", "#4d8c8c", "#5ba6a6", "#69bfbf", "#77d9d9", "#f0de84", "#c5d419" }),
		WAVE(new String[] { "#008cff", "#0099ff", "#00a5ff", "#00b2ff", "#00bfff", "#00cbff", "#00d8ff", "#00e5ff", "#00f2ff", "#00ffff" }),
		WHEEL(new String[] { "#ffe140", "#ffa922", "#1bc0c6", "#2484ae", "#134e6e" });
		// @formatter:on

	private final String[] value;
	private final int[] intValue;

	private Palette(String[] value) {
		this.value = value;
		this.intValue = RGB.hexToColor(value);
	}

	public String[] stringValue() {
		return value;
	}

	public String[] stringValue(int rotation) {
		return rotate(stringValue(), rotation);
	}

	public int[] intValue() {
		return intValue;
	}

	public int[] intValue(int rotation) {
		return rotate(intValue, rotation);
	}

	private static int[] rotate(int[] data, int rotation) {
		int[] result = new int[data.length];
		for (int i = 0; i < data.length; i++) {
			result[(i + (data.length - rotation)) % data.length] = data[i];
		}
		return result;
	}

	private static String[] rotate(String[] data, int rotation) {
		String[] result = new String[data.length];
		for (int i = 0; i < data.length; i++) {
			result[(i + (data.length - rotation)) % data.length] = data[i];
		}
		return result;
	}

}
