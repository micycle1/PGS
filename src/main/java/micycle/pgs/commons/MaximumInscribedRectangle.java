package micycle.pgs.commons;

import java.awt.Point;

import java.util.ArrayList;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;

/**
 * Computes the largest inscribed axis-aligned rectangle within a convex hull.
 * 
 * Credits: Daniel Sud for the Inscribed Rectangle algorithm
 * http://cgm.cs.mcgill.ca/~athens/cs507/Projects/2003/DanielSud/
 *
 * @author Peter (pborissow)
 * @author Adapted by Michael Carleton
 *
 ******************************************************************************/
public class MaximumInscribedRectangle {

	// https://github.com/pborissow/Poly2Rect

	private final ArrayList<Point> points; // use Point class for native integer x,y values
	private int[] rectangle;
	private int xmin, xmax, ymin, ymax; // position of hull

	public MaximumInscribedRectangle(Geometry polygon, double f) {
		points = new ArrayList<>();
		Coordinate[] coords = polygon.convexHull().getCoordinates();
		for (int i = 0; i < coords.length - 1; i++) {
			int x = (int) (coords[i].x / f);
			int y = (int) (coords[i].y / f);

			if (i > 0) {
				Point last = points.get(points.size() - 1); // points must be unclosed
				if (x == last.x && y == last.y) {
					continue;
				}
			}
			points.add(new Point(x, y));
		}
	}

	/**
	 * Returns the largest rectangle that will fit inside a convex hull
	 * 
	 * @return int[x, y, w, h]
	 */
	public int[] getInscribedRectangle() {
		if (rectangle == null) {
			rectangle = computeLargestRectangle();
		}
		return rectangle;
	}

	/**
	 * Used to find the largest rectangle that will fit inside a convex hull. This
	 * method uses a brute force algorithm to perform an exhaustive search for a
	 * solution.
	 */
	private int[] computeLargestRectangle() {
		if (points.size() < 3) {
			return new int[4];
		}
		int[] r = new int[] { 0, 0, 0, 0 };

		// Get list of edges that form the convex hull
		ArrayList<Edge> edgeList = computeEdgeList();

		// Find first top and bottom edges of the convex hull. The top edge forms
		// a 0-90 deg angle from minx. The bottom edge forms a 90-180 degree angle
		// from minx
		Edge topEdge = findEdge(xmin, true, edgeList);
		Edge bottomEdge = findEdge(xmin, false, edgeList);

		// Precompute a list of x-intercepts for every possible y coordinate value
		ArrayList<Point> xIntercepts = new ArrayList<>();
		for (int y = 0; y < ymax; y++) {
			int x = xIntersect(y, edgeList);
			xIntercepts.add(new Point(x, y));
		}

		// Scan for rectangle starting from the left-most position of the convex hull
		for (int x = xmin; x < xmax; x++) {

			// Find y-intercept for the top and bottom edges
			Integer top = yIntersect(x, topEdge);
			Integer bottom = yIntersect(x, bottomEdge);
			if (bottom == null) { // bottomEdge is vertical
				bottom = bottomEdge.ymax;
			}

			// Step through the y-intercepts
			for (int y = bottom; y >= top; y--) {// y = y-intercept from bottom to top
				for (int y1 = top; y1 < bottom; y1++) {// y1 = y-intercept from top to bottom
					if (y1 > y) {

						// Find right side (x-intercept of an edge to the right of the current position)
						int x1 = Math.max(xIntercepts.get(y).x, 0);
						int x2 = Math.max(xIntercepts.get(y1).x, 0);
						int right = Math.min(x1, x2);

						// Update rectangle
						if (right > 0) {
							int height = y1 - y;
							int width = right - x;
							int area = width * height;
							if (area > r[2] * r[3]) {
								r[0] = x;
								r[1] = y;
								r[2] = width;
								r[3] = height;
							}
						}
					}
				}
			}

			if (x == topEdge.xmax) {
				topEdge = findEdge(x, true, edgeList);
			}

			if (x == bottomEdge.xmax) {
				bottomEdge = findEdge(x, false, edgeList);
			}
		}
		
		return r;
	}

	private ArrayList<Edge> computeEdgeList() {
		ArrayList<Edge> edgeList = new ArrayList<>();
		Point a = points.get(points.size() - 1);
		for (int i = 0; i < points.size(); i++) {
			Point b = points.get(i);

			if (i == 0) {
				this.xmin = a.x;
				this.xmax = a.x;
				this.ymin = a.y;
				this.ymax = a.y;
			} else {
				if (a.x < this.xmin) {
					this.xmin = a.x;
				}
				if (a.x > this.xmax) {
					this.xmax = a.x;
				}
				if (a.y < this.ymin) {
					this.ymin = a.y;
				}
				if (a.y > this.ymax) {
					this.ymax = a.y;
				}
			}

			edgeList.add(new Edge(a, b));
			a = b;
		}
		return edgeList;
	}

	/**
	 * Returns the y-intercept of an edge for a given x coordinate
	 */
	private Integer yIntersect(int x, Edge e) {

		if (e.m == null) { // horizonal line
			return null;
		} else {
			double yfirst = e.m * (x - 0.5) + e.b;
			double ylast = e.m * (x + 0.5) + e.b;

			if (!e.isTop) {
				return (int) Math.floor(Math.min(yfirst, ylast));
			} else {
				return (int) Math.ceil(Math.max(yfirst, ylast));
			}
		}
	}

	/**
	 * Returns the x-intercept of an edge for a given y coordinate
	 */
	private int xIntersect(int y, ArrayList<Edge> edgeList) {
		int x = 0;
		for (Edge e : edgeList) {
			if (e.isRight && e.ymin <= y && e.ymax >= y) {
				if (e.m == null) {
					x = e.xmin;
				} else {
					double x0 = (y + 0.5 - e.b) / e.m;
					double x1 = (y - 0.5 - e.b) / e.m;
					x = (int) Math.floor(Math.min(x0, x1));
				}
			}
		}
		return x;
	}

	private Edge findEdge(int x, boolean isTop, ArrayList<Edge> edgeList) {
		ArrayList<Edge> edges = new ArrayList<>();
		for (Edge e : edgeList) {
			if (e.xmin == x) {
				if (e.isTop && isTop || !e.isTop && !isTop) {
					edges.add(e);
				}
			}
		}
		if (edges.size() == 1) {
			return edges.get(0);
		}
		for (Edge e : edges) {
			if (e.xmax == e.xmin) {
				return e;
			}
		}
		return edges.get(edges.size() - 1);
	}

	private class Edge {
		
		int xmin, xmax; /* horiz, +x is right */
		int ymin, ymax; /* vertical, +y is down */
		Double m, b; /* y = mx + b */
		boolean isTop, isRight; /* position of edge w.r.t. hull */
		Point p, q;

		public Edge(Point p, Point q) {
			this.xmin = Math.min(p.x, q.x);
			this.xmax = Math.max(p.x, q.x);
			this.ymin = Math.min(p.y, q.y);
			this.ymax = Math.max(p.y, q.y);

			if (p.x != q.x) {
				m = (double) (q.y - p.y) / (double) (q.x - p.x); // slope
				b = p.y - m * p.x; // y-intercept
			} else {
				// vertical line so no slope and no y intercept
			}

			this.isTop = p.x > q.x; // edge from right to left (ccw)
			this.isRight = p.y > q.y; // edge from bottom to top (ccw)
			this.p = p;
			this.q = q;
		}

		@Override
		public String toString() {
			return p + "->" + q;
		}
	}

}