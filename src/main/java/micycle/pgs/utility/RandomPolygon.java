package micycle.pgs.utility;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import processing.core.PVector;

/**
 * Random Convex Polygons
 */
public class RandomPolygon {

	private static final Random RAND = ThreadLocalRandom.current();

	/**
	 * 
	 * @param n number of vertices
	 * @param
	 * @return
	 */
	public static List<PVector> generateRandomConvexPolygon(int n, double xMax, double yMax) {
		// Generate two lists of random X and Y coordinates

		List<Double> xPool = new ArrayList<>(n);
		List<Double> yPool = new ArrayList<>(n);

		for (int i = 0; i < n; i++) {
			xPool.add(RAND.nextDouble() * xMax);
			yPool.add(RAND.nextDouble() * yMax);
		}

		// Sort them
		Collections.sort(xPool);
		Collections.sort(yPool);

		// Isolate the extreme points
		double minX = xPool.get(0);
		double maxX = xPool.get(n - 1);
		double minY = yPool.get(0);
		double maxY = yPool.get(n - 1);

		// Divide the interior points into two chains & Extract the vector components
		List<Double> xVec = new ArrayList<>(n);
		List<Double> yVec = new ArrayList<>(n);

		double lastTop = minX, lastBot = minX;

		for (int i = 1; i < n - 1; i++) {
			double x = xPool.get(i);

			if (RAND.nextBoolean()) {
				xVec.add(x - lastTop);
				lastTop = x;
			} else {
				xVec.add(lastBot - x);
				lastBot = x;
			}
		}

		xVec.add(maxX - lastTop);
		xVec.add(lastBot - maxX);

		double lastLeft = minY, lastRight = minY;

		for (int i = 1; i < n - 1; i++) {
			double y = yPool.get(i);

			if (RAND.nextBoolean()) {
				yVec.add(y - lastLeft);
				lastLeft = y;
			} else {
				yVec.add(lastRight - y);
				lastRight = y;
			}
		}

		yVec.add(maxY - lastLeft);
		yVec.add(lastRight - maxY);

		// Randomly pair up the X- and Y-components
		Collections.shuffle(yVec);

		// Combine the paired up components into vectors
		List<PVector> vec = new ArrayList<>(n);

		for (int i = 0; i < n; i++) {
			vec.add(new PVector(xVec.get(i).floatValue(), yVec.get(i).floatValue()));
		}

		// Sort the vectors by angle
		Collections.sort(vec, Comparator.comparingDouble(v -> Math.atan2(v.y, v.x)));

		// Lay them end-to-end
		double x = 0, y = 0;
		double minPolygonX = 0;
		double minPolygonY = 0;
		List<PVector> points = new ArrayList<>(n);

		for (int i = 0; i < n; i++) {
			points.add(new PVector((float) x, (float) y));

			x += vec.get(i).x;
			y += vec.get(i).y;

			minPolygonX = Math.min(minPolygonX, x);
			minPolygonY = Math.min(minPolygonY, y);
		}

		// Move the polygon to the original min and max coordinates
		double xShift = minX - minPolygonX;
		double yShift = minY - minPolygonY;

		for (int i = 0; i < n; i++) {
			PVector p = points.get(i);
			points.set(i, new PVector((float) (p.x + xShift), (float) (p.y + yShift)));
		}

		return points;
	}
}