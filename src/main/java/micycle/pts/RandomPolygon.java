package micycle.pts;

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
class RandomPolygon {

	private static final Random RAND = ThreadLocalRandom.current();

	/**
	 * 
	 * @param n number of vertices
	 * @param
	 * @return
	 */
	public static List<PVector> generateRandomConvexPolygon(int n, float xMax, float yMax) {
		// Generate two lists of random X and Y coordinates

		List<Float> xPool = new ArrayList<>(n);
		List<Float> yPool = new ArrayList<>(n);

		for (int i = 0; i < n; i++) {
			xPool.add(RAND.nextFloat() * yMax);
			yPool.add(RAND.nextFloat() * yMax);
		}

		// Sort them
		Collections.sort(xPool);
		Collections.sort(yPool);

		// Isolate the extreme points
		float minX = xPool.get(0);
		float maxX = xPool.get(n - 1);
		float minY = yPool.get(0);
		float maxY = yPool.get(n - 1);

		// Divide the interior points into two chains & Extract the vector components
		List<Float> xVec = new ArrayList<>(n);
		List<Float> yVec = new ArrayList<>(n);

		float lastTop = minX, lastBot = minX;

		for (int i = 1; i < n - 1; i++) {
			float x = xPool.get(i);

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

		float lastLeft = minY, lastRight = minY;

		for (int i = 1; i < n - 1; i++) {
			float y = yPool.get(i);

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
			vec.add(new PVector(xVec.get(i), yVec.get(i)));
		}

		// Sort the vectors by angle
		Collections.sort(vec, Comparator.comparingDouble(v -> Math.atan2(v.y, v.x)));

		// Lay them end-to-end
		float x = 0, y = 0;
		float minPolygonX = 0;
		float minPolygonY = 0;
		List<PVector> points = new ArrayList<>(n);

		for (int i = 0; i < n; i++) {
			points.add(new PVector(x, y));

			x += vec.get(i).x;
			y += vec.get(i).y;

			minPolygonX = Math.min(minPolygonX, x);
			minPolygonY = Math.min(minPolygonY, y);
		}

		// Move the polygon to the original min and max coordinates
		float xShift = minX - minPolygonX;
		float yShift = minY - minPolygonY;

		for (int i = 0; i < n; i++) {
			PVector p = points.get(i);
			points.set(i, new PVector(p.x + xShift, p.y + yShift));
		}

		return points;
	}
}