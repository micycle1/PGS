package micycle.pgs.commons;

import static processing.core.PApplet.lerp;

import java.util.ArrayList;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Chaikin’s corner-cutting algorithm, which is used to perform polygon
 * smoothing. The algorithm involves repeatedly generating a new curve by
 * cutting corners off the original (based on some fixed ratio).
 * 
 * @author Manohar Vanga
 */
public class ChaikinCut {

	// https://sighack.com/post/chaikin-curves

	/*
	 * TODO implement 'Multi-step Subdivision Algorithm for Chaikin Curves'; a
	 * faster approach without recursion
	 */

	private ChaikinCut() {
	}
	
	public static PShape chaikin(PShape shape, float ratio, int iterations) {
		if (shape.getChildCount() > 1) {
			PShape groupCut = new PShape(PConstants.GROUP);
			for (PShape child : shape.getChildren()) {
				groupCut.addChild(cut(child, ratio, iterations));
			}
			return groupCut;
		}
		else {
			return cut(shape, ratio, iterations);
		}
	}

	/**
	 * 
	 * @param shape
	 * @param ratio      value between zero and one. For example, a ratio of 0.25
	 *                   would mean cutting each edge one quarter and three quarters
	 *                   in.
	 * @param iterations number of cutting iterations to perform. 1 will cut
	 *                   corners; higher values produces smoother cuts. 10 should be
	 *                   enough
	 * @param close      whether shape is closed
	 * @return
	 */
	private static PShape cut(PShape shape, float ratio, int iterations) {
		// If the number of iterations is zero, return shape as is
		if (iterations < 1) {
			return shape;
		}

		final boolean close = shape.isClosed();

		PShape next = new PShape(PShape.GEOMETRY);
		next.setFill(iterations);
		next.beginShape();

		/*
		 * Step 1: Figure out how many corners the shape has depending on whether it's
		 * open or closed.
		 */
		final int num_corners = shape.getVertexCount() - (close ? 0 : 1);

		/*
		 * Step 2: Since we don't have access to edges directly with a PShape object, do
		 * a pairwise iteration over vertices instead. Same thing.
		 */
		for (int i = 0; i < num_corners; i++) {

			// Get the i'th and (i+1)'th vertex to work on that edge.
			PVector a = shape.getVertex(i);
			PVector b = shape.getVertex((i + 1) % shape.getVertexCount());

			if (cut(a, b)) { // don't cut small edges

				// Step 3: Break it using our chaikin_break() function
				ArrayList<PVector> n = chaikinCut(a, b, ratio);

				/*
				 * Now we have to deal with one corner case. In the case of open shapes, the
				 * first and last endpoints shouldn't be moved. However, in the case of closed
				 * shapes, we cut all edges on both ends.
				 */
				if (!close && i == 0) {
					// For the first point of open shapes, ignore vertex A
					next.vertex(a.x, a.y);
					next.vertex(n.get(1).x, n.get(1).y);
				} else if (!close && i == num_corners - 1) {
					// For the last point of open shapes, ignore vertex B
					next.vertex(n.get(0).x, n.get(0).y);
					next.vertex(b.x, b.y);
				} else {
					// For all other cases (i.e. interior edges of open
					// shapes or edges of closed shapes), add both vertices
					// returned by our chaikin_break() method
					next.vertex(n.get(0).x, n.get(0).y);
					next.vertex(n.get(1).x, n.get(1).y);
				}
			} else {
				next.vertex(a.x, a.y);
			}
		}

		if (close) {
			next.endShape(PConstants.CLOSE);
		} else {
			next.endShape();
		}

		return cut(next, ratio, iterations - 1);
	}

	/**
	 * Cuts an edge, represented by two vertices.
	 * 
	 * @param a
	 * @param b
	 * @param ratio determines where along the edge to make the cut
	 * @return
	 */
	private static ArrayList<PVector> chaikinCut(PVector a, PVector b, float ratio) {
		float x, y;
		ArrayList<PVector> n = new ArrayList<>();

		/*
		 * If ratio is greater than 0.5 flip it so we avoid cutting across the midpoint
		 * of the line.
		 */
		if (ratio > 0.5) {
			ratio = 1 - ratio;
		}

		/* Find point at a given ratio going from A to B */
		x = lerp(a.x, b.x, ratio);
		y = lerp(a.y, b.y, ratio);
		n.add(new PVector(x, y));

		/* Find point at a given ratio going from B to A */
		x = lerp(b.x, a.x, ratio);
		y = lerp(b.y, a.y, ratio);
		n.add(new PVector(x, y));

		return n;
	}

	/**
	 * Determines whether to cut an edge. Returns false for edges with a euclidean
	 * distance less than 1.0.
	 */
	private static boolean cut(PVector a, PVector b) {
		// TODO expand to exclude almost coincident edge pairs
		final float dx = b.x - a.x;
		final float dy = b.y - a.y;
		final float d2 = dx * dx + dy * dy;
		return d2 > 1;
	}

}
