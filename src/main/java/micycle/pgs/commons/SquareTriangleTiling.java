package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.SplittableRandom;

import micycle.pgs.color.ColorUtils;
import micycle.pgs.color.Colors;
import processing.core.PConstants;
import processing.core.PShape;

/**
 * A non-periodic (quasiperiodic?) tiling, comprising squares and equilateral
 * triangles.
 * 
 * @author uila
 * @author Java port by Michael Carleton
 */
public class SquareTriangleTiling {

	// SOURCE: https://editor.p5js.org/uila/sketches/DrmqpckUK

	final double w, h;
	final int maxTiles;
	final double maxRadius;
	int colorMode = 0; // 0 = alternate, 1=fixed, 2=petal, 3=random
	float strokeWeight = 2; // strokeweight
	int strokeColor = 0;
	List<Vert> verts;
	List<Double> scos = new ArrayList<>();
	List<Double> ssin = new ArrayList<>();
	int tileCount = 0;
	List<Integer> next = new ArrayList<>();
	List<Integer> colors = new ArrayList<>();
	int activeColor = 0;
	boolean circles = false; // whether to draw circles in some midpoints
	private SplittableRandom random;
	private float xOff;
	private float yOff;
	private float tileSize;

	private PShape tiling;

	/**
	 * Creates a tiling. Tiles begin at the center point and wrap clockwise.
	 */
	SquareTriangleTiling(double maxRadius, double tileSize) {
		this(maxRadius * 2, maxRadius * 2, maxRadius, tileSize, Integer.MAX_VALUE);
	}

	/**
	 * Creates a tiling. Tiles begin at the center point (given by x and y offset)
	 * and wrap clockwise.
	 */
	SquareTriangleTiling(double maxRadius, double tileSize, double xOffset, double yOffset) {
		this(maxRadius, tileSize);
		xOff = (float) xOffset;
		yOff = (float) yOffset;
	}

	public SquareTriangleTiling(double width, double height, double tileSize) {
		this(width, height, Math.max(width / 2, height / 2), tileSize, Integer.MAX_VALUE);
		xOff = (float) width / 2;
		yOff = (float) height / 2;
	}

	/**
	 * Creates a tiling.
	 * 
	 * @param maxWidth  maximum width of the tiling
	 * @param maxHeight maximum height of the tiling
	 * @param maxRadius maximum radius of the tiling
	 * @param tileSize
	 * @param maxTiles  maximum number of tiles
	 */
	SquareTriangleTiling(double maxWidth, double maxHeight, double maxRadius, double tileSize, int maxTiles) {
		this.w = maxWidth;
		this.h = maxHeight;
		this.maxRadius = maxRadius;
		this.maxTiles = maxTiles;
		this.tileSize = (float) tileSize;

		colors.add(ColorUtils.composeColor(200, 30, 99));
		colors.add(Colors.PINK);
		colors.add(Colors.WHITE);

		// store frequent calculations
		for (int i = 0; i < 12; i++) {
			scos.add(tileSize * Math.cos((i * Math.PI) / 6));
			ssin.add(tileSize * Math.sin((i * Math.PI) / 6));
		}
		// initialise the vertices
		verts = new ArrayList<>();
		verts.add(new Vert(0, 0, 0, new ArrayList<>(), 1, 1, 12));
		verts.add(new Vert(tileSize, 0, 6, new ArrayList<>(), 0, 0, 12));

		xOff = (float) Math.max((w / 2), maxRadius); // x-offset (centered on 0,0 by default)
		yOff = (float) Math.max((h / 2), maxRadius); // y-offset (centered on 0,0 by default)
	}

	public PShape getTiling() {
		return getTiling(System.nanoTime());
	}

	public PShape getTiling(long seed) {
		random = new SplittableRandom(seed);
		tiling = new PShape(PConstants.GROUP);

		while (maxTiles > tileCount) {

			int sweep = 0;
			Vert vert = null;
			while (sweep == 0 && !next.isEmpty()) {
				int n = next.get(0);
				vert = verts.get(n);
				sweep = vert.getSweep();
				if (sweep == 0) {
					next.remove(0);
					if (colorMode == 2) {
						activeColor = switch (colorCount(vert.colors)) {
							case 0 -> random.nextInt(0, 3);
							case 1 -> getSecondColor(activeColor);
							default -> getLastColor(vert.colors);
						};
					}
				}
			}

			if (vert == null || dist(0, 0, vert.x, vert.y) > w) {
				return tiling;
			}

			switch (sweep) {
				case 0 :
					return tiling;
				case 2 :
				case 4 :
					constructTriangle(vert);
					break;
				case 3 :
					constructSquare(vert);
					break;
				default :
					if (random.nextDouble() > 0.5) {
						constructSquare(vert);
					} else {
						constructTriangle(vert);
					}
			}
			tileCount++;
		}
		return tiling;
	}

	private void constructTriangle(Vert root) {
		final int ang = (root.startAngle + 2) % 12;
		if (root.getSweep() > 2) {
			verts.add(new Vert(root.x + scos.get(ang), root.y + ssin.get(ang), (ang + 8) % 12, new ArrayList<>(), root.id, root.last, 12));
			drawTriangle(root.id, root.last, verts.size() - 1);
		} else {
			drawTriangle(root.id, root.last, root.first);
		}
	}

	private void constructSquare(Vert root) {
		Vert base = verts.get(root.last);
		Vert third, fourth;
		int ang = (root.startAngle + 3) % 12;

		if (root.getSweep() == 3) {
			fourth = verts.get(root.first);
			if (fourth.getSweep() == 3) {
				third = verts.get(fourth.first);
			} else {
				third = new Vert(base.x + scos.get(ang), base.y + ssin.get(ang), (ang + 3) % 12, new ArrayList<>(), root.first, root.last,
						12);
				verts.add(third);
			}
		} else { // create third and fourth vert
			third = new Vert(base.x + scos.get(ang), base.y + ssin.get(ang), (ang + 3) % 12, new ArrayList<>(), verts.size() + 1, root.last,
					12);
			verts.add(third);
			fourth = new Vert(root.x + scos.get(ang), root.y + ssin.get(ang), (root.startAngle + 9) % 12, new ArrayList<>(), root.id,
					third.id, 12);
			verts.add(fourth);
		}
		drawSquare(root.id, base.id, third.id, fourth.id);
	}

	class Vert {

		final int id;
		final double x, y;
		int startAngle;
		List<Integer> colors;
		int first, last;
		final List<Integer> sweeps;

		Vert(double x, double y, int startAngle, List<Integer> colors, int first, int last, int sweep) {
			this.id = verts.size();
			this.x = x;
			this.y = y;
			this.startAngle = startAngle;
			this.colors = colors;
			this.first = first;
			this.last = last;
			sweeps = new ArrayList<>();
			this.sweeps.add(sweep);
			if (Math.abs(x) < (w - tileSize) / 2 && Math.abs(y) < (h - tileSize) / 2) {
				next.add(this.id);
			}
//			if (dist(0, 0, x, y) < maxRadius) { // NOTE uses radius
//			}
		}

		int getSweep() {
			if (!this.sweeps.isEmpty()) {
				return this.sweeps.get(0);
			} else {
				return 0;
			}
		}

		void reduceSweep(int val) {
			reduceSweep(val, -1);
		}

		void reduceSweep(int val, int pos) {
			int index = 0;
			if (pos == 1) {
				index = this.sweeps.size() - 1;
			}
			this.sweeps.set(index, this.sweeps.get(index) - val);
			if (this.sweeps.get(index) < 0) {
				return;
			}
			if (this.sweeps.get(index) == 0) {
				if (pos == 1) {
					this.sweeps.remove(this.sweeps.size() - 1); // pop
				} else if (this.sweeps.size() > 1) {
					this.sweeps.remove(0);
				}
			}
			int sum = 0;
			for (Integer element : this.sweeps) {
				sum += element;
			}
//			if (sum < 5 && dist(0, 0, this.x, this.y) < maxRadius) { // NOTE uses radius
//				next.add(0, this.id);
//			}
			if (sum < 5 && this.x < (w - tileSize) && this.y < (h - tileSize)) {
				next.add(0, this.id);
			}
		}

		void splitSweep(int val) {
			int i = this.sweeps.size() - 1;
			this.sweeps.set(i, this.sweeps.get(i) - val);
			this.sweeps.add(val);
		}
	}

	private int getColor(List<Integer> colors, int sweep, int tile) {
		switch (colorMode) {
			case 1 :
				return tile;
			case 2 :
				return activeColor;
			case 3 :
				return random.nextInt(0, 3);
		}

		if ((sweep == 2) || (sweep == 3)) {
			return getLastColor(colors);
		}
		return switch (colors.size()) {
			case 0 -> random.nextInt(0, 3);
			case 1 -> getSecondColor(colors.get(0));
			case 2 -> colors.get(0);
			default -> getNextColor(colors);
		};
	}

	private int getSecondColor(int c1) {
		return switch (c1) {
			case 0 -> random.nextInt(1, 3);
			case 1 -> random.nextInt(0, 3);
			case 2 -> random.nextInt(0, 2);
			default -> 0; // shouldn't be hit
		};
	}

	private int getThirdColor(int c1, int c2) {
		return switch (c1 + c2) {
			case 1 -> 2;
			case 2 -> 1;
			case 3 -> 0;
			default -> 0; // shouldn't be hit
		};
	}

	private int getLastColor(List<Integer> colors) {
		return getLastColor(colors, 0);
	}

	private int getLastColor(List<Integer> colors, int tile) {
		if (colorMode == 1) {
			return tile;
		}
		if (colorMode == 2) {
			return activeColor;
		}
		int c1 = colors.get(0);
		int c2 = colors.get(colors.size() - 1);
		if (c1 == c2) {
			return colors.get(1);
		} else {
			return getThirdColor(c1, c2);
		}
	}

	private int getNextColor(List<Integer> colors) {
		if (areAllColorsUsed(colors)) {
			return getLastColor(colors);
		} else {
			return colors.get(colors.size() - 2);
		}
	}

	private static boolean areAllColorsUsed(List<Integer> colors) {
		boolean c0 = false;
		boolean c1 = false;
		boolean c2 = false;
		for (Integer color : colors) {
			switch (color) {
				case 0 :
					c0 = true;
					break;
				case 1 :
					c1 = true;
					break;
				case 2 :
					c2 = true;
					break;
			}
		}
		return c0 && c1 && c2;
	}

	private static int colorCount(List<Integer> colors) {
		return new HashSet<>(colors).size();
	}

	private void adjustSweeps(int currentSweep, int last_id) {
		while (currentSweep < 5) {
			final Vert last = verts.get(last_id);
			switch (currentSweep) {
				case 2 :
				case 4 :
					last.splitSweep(2);
					break;
				case 3 :
					last.splitSweep(3);
			}
			currentSweep = last.getSweep();
			last_id = last.last;
			if (last.id == last.last) {
				break;
			}
		}
	}

	private void drawTriangle(int root_id, int base_id, int last_id) {
		Vert root = verts.get(root_id);
		Vert base = verts.get(base_id);
		Vert last = verts.get(last_id);

		int newColor = getColor(root.colors, root.getSweep(), 1);

		final PShape triangle = new PShape(PShape.PATH);
		triangle.setFill(true);
		triangle.setFill(colors.get(newColor));
		triangle.setStroke(true);
		triangle.setStroke(strokeColor);
		triangle.setStrokeWeight(strokeWeight);
		triangle.beginShape();
		triangle.vertex((float) root.x + xOff, (float) root.y + yOff);
		triangle.vertex((float) base.x + xOff, (float) base.y + yOff);
		triangle.vertex((float) last.x + xOff, (float) last.y + yOff);
		triangle.endShape(PConstants.CLOSE);
		tiling.addChild(triangle);

		root.reduceSweep(2);
		root.colors.add(newColor);
		root.last = last_id;
		root.startAngle = (root.startAngle + 2) % 12;

		base.reduceSweep(2, 1);
		base.colors.add(0, newColor);
		base.first = last_id;

		last.reduceSweep(2, 1);

		if (root.getSweep() == 0) {
			// just closed a gap
//			drawCircle(root);
			last.last = base.id;
			last.colors.add(newColor);
			last.startAngle = (last.startAngle + 2) % 12;
		} else {
			last.first = root_id;
			last.colors.add(0, newColor);
		}
	}

	private void drawSquare(int root_id, int base_id, int third_id, int fourth_id) {
		Vert root = verts.get(root_id);
		Vert base = verts.get(base_id);
		Vert third = verts.get(third_id);
		Vert fourth = verts.get(fourth_id);

		int newColor = getColor(root.colors, root.getSweep(), 2);

		root.colors.add(newColor);
		root.startAngle = (root.startAngle + 3) % 12;
		root.reduceSweep(3);
		root.last = fourth.id;

		base.colors.add(0, newColor);
		base.first = third.id;
		base.reduceSweep(3, 1);

		third.colors.add(newColor);
		third.reduceSweep(3);
		third.startAngle = (third.startAngle + 3) % 12;

		fourth.colors.add(newColor);
		fourth.startAngle = (fourth.startAngle + 3) % 12;
		fourth.reduceSweep(3, 0);
		fourth.last = third.id;

		final PShape quad = new PShape(PShape.PATH);
		quad.setFill(true);
		quad.setFill(colors.get(newColor));
		quad.setStroke(true);
		quad.setStroke(strokeColor);
		quad.setStrokeWeight(strokeWeight);
		quad.beginShape();
		quad.vertex((float) root.x + xOff, (float) root.y + yOff);
		quad.vertex((float) base.x + xOff, (float) base.y + yOff);
		quad.vertex((float) third.x + xOff, (float) third.y + yOff);
		quad.vertex((float) fourth.x + xOff, (float) fourth.y + yOff);
		quad.endShape(PConstants.CLOSE);
		tiling.addChild(quad);

		adjustSweeps(fourth.getSweep(), fourth.last);

		if (root.getSweep() == 0) {
//			drawCircle(root);
		}
	}

	private static final double dist(double x1, double y1, double x2, double y2) {
		return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
	}

}
