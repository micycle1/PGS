package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

import processing.core.PConstants;
import processing.core.PShape;

/**
 * Port of https://openprocessing.org/sketch/1523350/
 * 
 * @author Michael Carleton
 *
 */
public class HatchTiling {

	private final int width, height;
	private final int gridCountX, gridCountY;
	private PShape tiling;
	private SplittableRandom random;

	public HatchTiling(int width, int height) {
		this(width, height, 40, 40);
	}

	public HatchTiling(int width, int height, int gridCountX, int gridCountY) {
		this.width = width;
		this.height = height;
		this.gridCountX = gridCountX;
		this.gridCountY = gridCountY;
	}

	public PShape getTiling(long seed) {
		random = new SplittableRandom(seed);
		tiling = new PShape(PConstants.GROUP);
		tiling();
		return tiling;
	}

	public PShape getTiling() {
		return getTiling(System.nanoTime());
	}

	private void tiling() {
		int offset = 0;
		int gridCountW = gridCountX;
		int gridCountH = gridCountY;
		int gridW = (width - (offset * 2)) / gridCountW;
		int gridH = (height - (offset * 2)) / gridCountH;
		int emp = gridCountW * gridCountH;
		boolean[][] grids = new boolean[gridCountW][gridCountH];
		List<int[]> rects = new ArrayList<>();
		for (int j = 0; j < gridCountW; j++) {
			boolean[] arr = new boolean[gridCountH];
			for (int i = 0; i < gridCountH; i++) {
				arr[i] = false;
			}
			grids[j] = arr;
		}

		while (emp > 0) {
			int w = (random.nextInt(1, gridCountW));
			int h = 1;
			if (random.nextDouble() < 0.5) {
				w = 1;
				h = (random.nextInt(1, gridCountH));
			}
			int x = (random.nextInt(gridCountW - w + 1));
			int y = (random.nextInt(gridCountH - h + 1));
			boolean lap = true;
			for (int j = 0; j < h; j++) {
				for (int i = 0; i < w; i++) {
					if (grids[x + i][y + j]) {
						lap = false;
						break;
					}
				}
			}

			if (lap) {
				for (int j = 0; j < h; j++) {
					for (int i = 0; i < w; i++) {
						grids[x + i][y + j] = true;
					}
				}
				int xx = offset + x * gridW;
				int yy = offset + y * gridH;
				int ww = w * gridW;
				int hh = h * gridH;
				rects.add(new int[] { xx, yy, ww, hh });
				emp -= w * h;
			}
		}
		for (int[] rect : rects) {
			final int off = 0; // border between adjacent rects
			rect(rect[0] + off, rect[1] + off, rect[2] - off * 2, rect[3] - off * 2);
		}
	}

	private void rect(int x, int y, int w, int h) {
		final PShape quad = new PShape(PShape.PATH);
		quad.setFill(true);
		quad.setFill(255);
		quad.setStroke(true);
		quad.setStroke(0);
		quad.setStrokeWeight(1);
		quad.beginShape();
		quad.vertex(x, y);
		quad.vertex(x + w, y);
		quad.vertex(x + w, y + h);
		quad.vertex(x, y + h);
		quad.endShape(PConstants.CLOSE);
		tiling.addChild(quad);
	}

}
