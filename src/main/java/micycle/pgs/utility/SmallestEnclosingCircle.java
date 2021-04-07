package micycle.pgs.utility;

/* 
 * Smallest enclosing circle - Library (Java)
 * 
 * Copyright (c) 2020 Project Nayuki
 * https://www.nayuki.io/page/smallest-enclosing-circle
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program (see COPYING.txt and COPYING.LESSER.txt).
 * If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import processing.core.PVector;

/**
 * Computes the smallest circle that encloses an arbitrary set of points. Runs
 * in linear time. Adapts Nayuki's original implementation for easy use with
 * Processing.
 * 
 * @author Michael Carleton
 *
 */
public final class SmallestEnclosingCircle {

	/*
	 * Returns the smallest circle that encloses all the given PVectors. Runs in
	 * expected O(n) time, randomized. Note: If 0 PVectors are given, null is
	 * returned. If 1 PVector is given, a circle of radius 0 is returned.
	 */
	public static Circle smallestEnclosingCircle(List<PVector> PVectors) {
		// Initially: No boundary PVectors known
		// Clone list to preserve the caller's data, randomize order
		List<PVector> shuffled = new ArrayList<>(PVectors);
		Collections.shuffle(shuffled, new Random());

		// Progressively add PVectors to circle or recompute circle
		Circle c = null;
		for (int i = 0; i < shuffled.size(); i++) {
			PVector p = shuffled.get(i);
			if (c == null || !c.contains(p))
				c = makeCircleOnePVector(shuffled.subList(0, i + 1), p);
		}
		return c;
	}

	public static Circle sec(float x1, float y1, float x2, float y2, float x3, float y3) {
		List<PVector> points = new ArrayList<PVector>(3);
		points.add(new PVector(x1, y1));
		points.add(new PVector(x2, y2));
		points.add(new PVector(x3, y3));
		return smallestEnclosingCircle(points);
	}

	public static Circle sec(PVector a, PVector b, PVector c) {
		List<PVector> points = new ArrayList<PVector>(3);
		points.add(a);
		points.add(b);
		points.add(c);
		return smallestEnclosingCircle(points);
	}

	public static Circle makeCircle(PVector... PVectors) {
		return smallestEnclosingCircle(Arrays.asList(PVectors));
	}

	// One boundary PVector known
	private static Circle makeCircleOnePVector(List<PVector> PVectors, PVector p) {
		Circle c = new Circle(p, 0);
		for (int i = 0; i < PVectors.size(); i++) {
			PVector q = PVectors.get(i);
			if (!c.contains(q)) {
				if (c.r == 0)
					c = makeDiameter(p, q);
				else
					c = makeCircleTwoPVectors(PVectors.subList(0, i + 1), p, q);
			}
		}
		return c;
	}

	// Two boundary PVectors known
	private static Circle makeCircleTwoPVectors(List<PVector> PVectors, PVector p, PVector q) {
		Circle circ = makeDiameter(p, q);
		Circle left = null;
		Circle right = null;

		// For each PVector not in the two-PVector circle
		PVector pq = q.sub(p);
		for (PVector r : PVectors) {
			if (circ.contains(r))
				continue;

			// Form a circumcircle and classify it on left or right side
			float cross = cross(pq, r.sub(p));
			Circle c = makeCircumcircle(p, q, r);
			if (c == null)
				continue;
			else if (cross > 0 && (left == null || cross(pq, c.c.sub(p)) > cross(pq, left.c.sub(p))))
				left = c;
			else if (cross < 0 && (right == null || cross(pq, c.c.sub(p)) < cross(pq, right.c.sub(p))))
				right = c;
		}

		// Select which circle to return
		if (left == null && right == null)
			return circ;
		else if (left == null)
			return right;
		else if (right == null)
			return left;
		else
			return left.r <= right.r ? left : right;
	}

	private static Circle makeDiameter(PVector a, PVector b) {
		PVector c = new PVector((a.x + b.x) / 2, (a.y + b.y) / 2);
		return new Circle(c, Math.max(c.dist(a), c.dist(b)));
	}

	private static Circle makeCircumcircle(PVector a, PVector b, PVector c) {
		// Mathematical algorithm from Wikipedia: Circumscribed circle
		float ox = (Math.min(Math.min(a.x, b.x), c.x) + Math.max(Math.max(a.x, b.x), c.x)) / 2;
		float oy = (Math.min(Math.min(a.y, b.y), c.y) + Math.max(Math.max(a.y, b.y), c.y)) / 2;
		float ax = a.x - ox, ay = a.y - oy;
		float bx = b.x - ox, by = b.y - oy;
		float cx = c.x - ox, cy = c.y - oy;
		float d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2;
		if (d == 0)
			return null;
		float x = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by))
				/ d;
		float y = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax))
				/ d;
		PVector p = new PVector(ox + x, oy + y);
		float r = Math.max(Math.max(p.dist(a), p.dist(b)), p.dist(c));
		return new Circle(p, r);
	}

	private static float cross(PVector a, PVector b) {
		// Signed area / determinant thing
		return a.x * b.y - a.y * b.x;
	}

}