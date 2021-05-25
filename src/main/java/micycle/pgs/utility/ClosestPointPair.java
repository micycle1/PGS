package micycle.pgs.utility;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import processing.core.PVector;

/* Copyright (c) 2012 Kevin L. Stern
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * An implementation of a divide-and-conquer algorithm for computing the closest
 * pair of elements of a set of points.
 * <p>
 * The algorithm consists of constructing an ordered list of points, then
 * recursively dividing the list into a left and right sublist towards finding
 * the closest point pairs for each sublist. The two sub-results are merged by
 * selecting the optimal among them and all closer point pairs that cross the
 * boundary of separation. Happily, only a linear amount of work is required to
 * find all closer point pairs that cross the boundary, giving a total runtime
 * of O(n*log(n)) for the algorithm.
 * 
 * @author Kevin L. Stern
 * @author Adapted for Processing by Michael Carleton
 */
public class ClosestPointPair {

	// https://github.com/KevinStern/software-and-algorithms

	private List<PVector> pointsOrderedByXCoordinate;
	private List<PVector> pointsOrderedByYCoordinate;
	
	/**
	 * Construct an instance of the algorithm for the specified point Collection.
	 * 
	 * @param points the Collection of points through which to search for the
	 *               closest pair.
	 */
	public ClosestPointPair(Collection<PVector> points) {
		if (points == null) {
			throw new NullPointerException("Point set is null.");
		}
		if (points.size() < 2) {
			throw new IllegalArgumentException("Point set is too small.");
		}
		pointsOrderedByXCoordinate = new ArrayList<>(points);
		Collections.sort(pointsOrderedByXCoordinate, (o1, o2) -> {
			float delta = o1.x - o2.x;
			if (delta == 0.0) {
				delta = o1.y - o2.y;
			}
			return delta < 0 ? -1 : delta > 0 ? 1 : 0;
		});
		pointsOrderedByYCoordinate = new ArrayList<>(points);
		Collections.sort(pointsOrderedByYCoordinate, (o1, o2) -> {
			float delta = o1.y - o2.y;
			if (delta == 0.0) {
				delta = o1.x - o2.x;
			}
			return delta < 0 ? -1 : delta > 0 ? 1 : 0;
		});
	}

	/**
	 * Execute the algorithm.
	 * 
	 * @return a List<PVector> containing exactly two elements which are the closest
	 *         pair of points among those in the collection used to construct this
	 *         instance.
	 */
	public List<PVector> execute() {
		PairStructure result = closestPair(0, pointsOrderedByXCoordinate.size(), pointsOrderedByYCoordinate);
		List<PVector> out = new ArrayList<>();
		out.add(result.p1);
		out.add(result.p2);
		return out;
	}

	/**
	 * Internal helper method which implements the closest point pair algorithm.
	 * 
	 * @param low                             the starting index, inclusive, of the
	 *                                        sublist in which to search for the
	 *                                        closest point pair.
	 * @param high                            the ending index, exclusive, of the
	 *                                        sublist in which to search for the
	 *                                        closest point pair.
	 * @param localPointsOrderedByYCoordinate the points from the target sublist,
	 *                                        ordered by y coordinate.
	 * @return a PairStructure containing the closest point pair among elements of
	 *         the target sublist.
	 */
	private PairStructure closestPair(int low, int high, List<PVector> localPointsOrderedByYCoordinate) {
		int size = high - low;
		if (size == 3) {
			return closestPair(pointsOrderedByXCoordinate.get(low), pointsOrderedByXCoordinate.get(low + 1),
					pointsOrderedByXCoordinate.get(low + 2));
		} else if (size == 2) {
			PVector p1 = pointsOrderedByXCoordinate.get(low);
			PVector p2 = pointsOrderedByXCoordinate.get(low + 1);
			return new PairStructure(p1, p2, distanceSq(p1, p2));
		}

		int mid = (low + high) >> 1;  // (low + high) / 2
		Set<PVector> leftSubtreeMemberSet = new HashSet<>(mid - low);
		for (int j = low; j < mid; j++) {
			leftSubtreeMemberSet.add(pointsOrderedByXCoordinate.get(j));
		}

		/*
		 * Construct the lists of points ordered by y coordinate for the left and right
		 * subtrees in linear time by drawing upon the master list of points ordered by
		 * y coordinate.
		 */
		List<PVector> leftPointsOrderedByYCoordinate = new ArrayList<>(mid - low);
		List<PVector> rightPointsOrderedByYCoordinate = new ArrayList<>(high - mid);
		for (PVector next : localPointsOrderedByYCoordinate) {
			if (leftSubtreeMemberSet.contains(next)) {
				leftPointsOrderedByYCoordinate.add(next);
			} else {
				rightPointsOrderedByYCoordinate.add(next);
			}
		}

		PairStructure leftSubtreeResult = closestPair(low, mid, leftPointsOrderedByYCoordinate);
		PairStructure rightSubtreeResult = closestPair(mid, high, rightPointsOrderedByYCoordinate);
		PairStructure result = leftSubtreeResult.distanceSq < rightSubtreeResult.distanceSq ? leftSubtreeResult : rightSubtreeResult;

		List<PVector> boundaryPointsOrderedByYCoordinate = new ArrayList<>();
		float midXCoordinate = pointsOrderedByXCoordinate.get(mid).x;
		for (PVector next : localPointsOrderedByYCoordinate) {
			float v = next.x - midXCoordinate;
			if (v * v < result.distanceSq) {
				boundaryPointsOrderedByYCoordinate.add(next);
			}
		}
		for (int i = 0; i < boundaryPointsOrderedByYCoordinate.size(); ++i) {
			PVector currentPoint = boundaryPointsOrderedByYCoordinate.get(i);
			int index;
			for (int j = 1; (index = i + j) < boundaryPointsOrderedByYCoordinate.size(); ++j) {
				PVector testPoint = boundaryPointsOrderedByYCoordinate.get(index);
				/*
				 * The number of points that can be situated within the boundary so that their y
				 * coordinate is within the minimum of the result distances for the left and
				 * right subtrees from currentPoint.y is bounded by a constant, since that
				 * distance value spatially limits the number of points that can be packed near
				 * one another on each side of the boundary.
				 */
				float v = testPoint.y - currentPoint.y;
				if (v * v >= result.distanceSq) {
					break;
				}
				float testDistance = distanceSq(currentPoint, testPoint);
				if (testDistance < result.distanceSq) {
					result = new PairStructure(currentPoint, testPoint, testDistance);
				}
			}
		}

		return result;
	}

	/**
	 * Find the closest pair of points among p1, p2 and p3.
	 */
	private static PairStructure closestPair(PVector p1, PVector p2, PVector p3) {
		float d1 = distanceSq(p1, p2);
		float d2 = distanceSq(p2, p3);
		float d3 = distanceSq(p3, p1);
		if (d1 < d2) {
			if (d1 < d3) {
				return new PairStructure(p1, p2, d1);
			} else {
				return new PairStructure(p1, p3, d3);
			}
		} else {
			if (d2 < d3) {
				return new PairStructure(p2, p3, d2);
			} else {
				return new PairStructure(p1, p3, d3);
			}
		}
	}

	/**
	 * Computes the squared distance between two PVectors.
	 */
	private static float distanceSq(PVector a, PVector b) {
		float dx = a.x - b.x;
		float dy = a.y - b.y;
		return (dx * dx + dy * dy);
	}

	/**
	 * Convenience data structure to hold a pair of points along with their distance
	 * from one another.
	 */
	private static class PairStructure {

		private PVector p1, p2;
		private float distanceSq;

		/**
		 * Constructor.
		 * 
		 * @param p1         the first point.
		 * @param p2         the second point.
		 * @param distanceSq the distance between p1 and p2, squared.
		 */
		public PairStructure(PVector p1, PVector p2, float distanceSq) {
			this.p1 = p1;
			this.p2 = p2;
			this.distanceSq = distanceSq;
		}
	}
}