package micycle.pgs.commons;

import static java.lang.Math.PI;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import micycle.pgs.color.ColorUtils;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Generates a Penrose tiling of the plane.
 * <p>
 * A Java port of https://openprocessing.org/sketch/183715
 * 
 * @author Michael Carleton
 *
 */
public class PenroseTiling {

	private static final float PHI = (float) ((1 + Math.sqrt(5)) / 2);

//	private double width, height;
	private List<Triangle> triangles = new ArrayList<>();

	public PenroseTiling(final double centerX, final double centerY, final double radius, final int divisions) {
		makeWheel(new PVector((float) centerX, (float) centerY), radius);
		for (int i = 0; i < divisions; i++) {
			triangles = subdivide(triangles, 1);
		}
	}

	/**
	 * Returns the edge work of the tiling.
	 */
	public Set<PEdge> getEdges() {
		final Set<PEdge> edges = new HashSet<>(triangles.size() * 3);
		triangles.forEach(t -> {
			PEdge a = new PEdge(t.v[1], t.v[0]);
			PEdge c = new PEdge(t.v[0], t.v[2]);
			edges.add(a.round());
			edges.add(c.round());
		});
		return edges;
	}

	/**
	 * Returns raw triangles from the penrose tiling. Note these are not collapsed
	 * (the tiling should consist of quadrangles).
	 */
	public PShape getTriangles() {
		final PShape tiling = new PShape(PConstants.GROUP);

		triangles.forEach(t -> {
			final PShape triangle = new PShape(PShape.PATH);

			triangle.setFill(true);
			if (t.alt) {
				triangle.setFill(ColorUtils.composeColor(255, 0, 255));
			} else {
				triangle.setFill(ColorUtils.composeColor(255, 255, 0));
			}
			triangle.setStroke(true);
			triangle.setStrokeWeight(2);

			triangle.beginShape();
			triangle.vertex(t.v[0].x, t.v[0].y);
			triangle.vertex(t.v[1].x, t.v[1].y);
			triangle.vertex(t.v[2].x, t.v[2].y);
			triangle.endShape(PConstants.CLOSE);
			tiling.addChild(triangle);
		});

		return tiling;
	}

	/**
	 * Create the initial triangle wheel/fan.
	 */
	private void makeWheel(PVector origin, double r) {
		for (int i = 0; i < 10; i++) {
			PVector b = new PVector((float) (origin.x + r * Math.cos((2 * i - 1) * PI / 10)),
					(float) (origin.y + r * Math.sin((2 * i - 1) * PI / 10)));
			PVector c = new PVector((float) (origin.x + r * Math.cos((2 * i + 1) * PI / 10)),
					(float) (origin.y + r * Math.sin((2 * i + 1) * PI / 10)));
			if (i % 2 == 0) { // mirror every second triangle
				PVector temp = b;
				b = c;
				c = temp;
			}
			triangles.add(new Triangle(true, origin, b, c));
		}
	}

	private List<Triangle> subdivide(List<Triangle> tList, float s) {
		List<Triangle> result = new ArrayList<>();
		for (int i = 0; i < tList.size(); i++) {
			Triangle t = tList.get(i);
			if (t.alt) {
				PVector p = new PVector(t.v[0].x + s * (t.v[1].x - t.v[0].x) / PHI, t.v[0].y + s * (t.v[1].y - t.v[0].y) / PHI);
				result.add(new Triangle(true, t.v[2], p, t.v[1]));
				result.add(new Triangle(true, t.v[2], p, t.v[1]));
				result.add(new Triangle(false, p, t.v[2], t.v[0]));
			} else {
				PVector q = new PVector(t.v[0].x - (1 - 1 / PHI) * s * (t.v[0].x - t.v[1].x),
						t.v[0].y - (1 - 1 / PHI) * s * (t.v[0].y - t.v[1].y));
				PVector r = new PVector(t.v[1].x + (t.v[2].x - t.v[1].x) / PHI, t.v[1].y + (t.v[2].y - t.v[1].y) / PHI);
				result.add(new Triangle(false, r, t.v[2], t.v[0]));
				result.add(new Triangle(false, q, r, t.v[1]));
				result.add(new Triangle(true, r, q, t.v[0]));
			}
		}
		return result;
	}

	private class Triangle {

		/** Whether triangle is alternative layer/color */
		boolean alt;
		PVector[] v = new PVector[3];

		private Triangle(boolean red, PVector a, PVector b, PVector c) {
			this.alt = red;
			this.v[0] = a;
			this.v[1] = b;
			this.v[2] = c;
		}
	}

}
