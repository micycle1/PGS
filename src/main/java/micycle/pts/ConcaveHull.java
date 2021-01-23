package micycle.pts;

import java.util.*;

import org.locationtech.jts.geom.Coordinate;

/**
 * ConcaveHull.java - 14/10/16
 * 
 * https://github.com/Merowech/java-concave-hull
 *
 * @author Udo Schlegel - Udo.3.Schlegel(at)uni-konstanz.de
 * @version 1.0
 *
 *          This is an implementation of the algorithm described by Adriano
 *          Moreira and Maribel Yasmina Santos: CONCAVE HULL: A K-NEAREST
 *          NEIGHBOURS APPROACH FOR THE COMPUTATION OF THE REGION OCCUPIED BY A
 *          SET OF POINTS. GRAPP 2007 - International Conference on Computer
 *          Graphics Theory and Applications; pp 61-68.
 *
 *          https://repositorium.sdum.uminho.pt/bitstream/1822/6429/1/ConcaveHull_ACM_MYS.pdf
 *
 *          With help from
 *          https://github.com/detlevn/QGIS-ConcaveHull-Plugin/blob/master/concavehull.py
 * 
 */
public class ConcaveHull {

	// try
	// https://github.com/atolcd/pentaho-gis-plugins/blob/master/concave-hull/src/main/java/org/opensphere/geometry/algorithm/ConcaveHull.java

	// TODO JAFAMA OPTIMSE?
	// TODO make static?

	public ConcaveHull() {
	}

	/**
	 * 
	 * @param pointArrayList
	 * @param k              k-nearest points of the current point (point C) are
	 *                       selected as candidates to be the next point of the
	 *                       polygon
	 * @return
	 */
	public ArrayList<Coordinate> calculateConcaveHull(ArrayList<Coordinate> pointArrayList, Integer k) {

		// the resulting concave hull
		ArrayList<Coordinate> concaveHull = new ArrayList<>();

		// optional remove duplicates
		HashSet<Coordinate> set = new HashSet<>(pointArrayList);
		ArrayList<Coordinate> pointArraySet = new ArrayList<>(set);

		// k has to be greater than 3 to execute the algorithm
		int kk = Math.max(k, 3);

		// return Coordinates if already Concave Hull
		if (pointArraySet.size() < 3) {
			return pointArraySet;
		}

		// make sure that k neighbors can be found
		kk = Math.min(kk, pointArraySet.size() - 1);

		// find first point and remove from point list
		Coordinate firstCoordinate = findMinYCoordinate(pointArraySet);
		concaveHull.add(firstCoordinate);
		Coordinate currentCoordinate = firstCoordinate;
		pointArraySet.remove(firstCoordinate);

		double previousAngle = 0.0;
		int step = 2;

		while ((currentCoordinate != firstCoordinate || step == 2) && pointArraySet.size() > 0) {

			// after 3 steps add first point to dataset, otherwise hull cannot be closed
			if (step == 5) {
				pointArraySet.add(firstCoordinate);
			}

			// get k nearest neighbors of current point
			ArrayList<Coordinate> kNearestCoordinates = kNearestNeighbors(pointArraySet, currentCoordinate, kk);

			// sort points by angle clockwise
			ArrayList<Coordinate> clockwiseCoordinates = sortByAngle(kNearestCoordinates, currentCoordinate,
					previousAngle);

			// check if clockwise angle nearest neighbors are candidates for concave hull
			boolean its = true;
			int i = -1;
			while (its && i < clockwiseCoordinates.size() - 1) {
				i++;

				int lastCoordinate = 0;
				if (clockwiseCoordinates.get(i) == firstCoordinate) {
					lastCoordinate = 1;
				}

				// check if possible new concave hull point intersects with others
				int j = 2;
				its = false;
				while (!its && j < concaveHull.size() - lastCoordinate) {
					its = intersect(concaveHull.get(step - 2), clockwiseCoordinates.get(i),
							concaveHull.get(step - 2 - j), concaveHull.get(step - 1 - j));
					j++;
				}
			}

			// if there is no candidate increase k - try again
			if (its) {
				return calculateConcaveHull(pointArrayList, k + 1);
			}

			// add candidate to concave hull and remove from dataset
			currentCoordinate = clockwiseCoordinates.get(i);
			concaveHull.add(currentCoordinate);
			pointArraySet.remove(currentCoordinate);

			// calculate last angle of the concave hull line
			previousAngle = calculateAngle(concaveHull.get(step - 1), concaveHull.get(step - 2));

			step++;

		}

		// Check if all points are contained in the concave hull
		boolean insideCheck = true;
		int i = pointArraySet.size() - 1;

		while (insideCheck && i > 0) {
			insideCheck = pointInPolygon(pointArraySet.get(i), concaveHull);
			i--;
		}

		// if not all points inside - try again
		if (!insideCheck) {
			return calculateConcaveHull(pointArrayList, k + 1);
		} else {
			return concaveHull;
		}

	}

	private static double euclideanDistance(Coordinate a, Coordinate b) {
		return ((a.getX() - b.getX()) * (a.getX() - b.getX())) + ((a.getY() - b.getY()) * (a.getY() - b.getY()));
//		return Math.abs(Math.pow(a.getX() - b.getX(), 2) + Math.pow(a.getY() - b.getY(), 2));
	}

	private ArrayList<Coordinate> kNearestNeighbors(ArrayList<Coordinate> l, Coordinate q, Integer k) {
		ArrayList<Pair<Double, Coordinate>> nearestList = new ArrayList<>();
		for (Coordinate o : l) {
			nearestList.add(new Pair<>(euclideanDistance(q, o), o));
		}

		Collections.sort(nearestList, new Comparator<Pair<Double, Coordinate>>() {
			@Override
			public int compare(Pair<Double, Coordinate> o1, Pair<Double, Coordinate> o2) {
				return o1.getKey().compareTo(o2.getKey());
			}
		});

		ArrayList<Coordinate> result = new ArrayList<>();

		for (int i = 0; i < Math.min(k, nearestList.size()); i++) {
			result.add(nearestList.get(i).getValue());
		}

		return result;
	}

	private static Coordinate findMinYCoordinate(ArrayList<Coordinate> l) {
		Collections.sort(l, new Comparator<Coordinate>() {
			@Override
			public int compare(Coordinate o1, Coordinate o2) {
//				return o1.getY().compareTo(o2.getY());
				if (o1.getY() == o2.getY()) {
					return 0;
				} else {
					if (o1.getY() > o2.getY()) {
						return 1;
					} else {
						return -1;
					}
				}
			}
		});
		return l.get(0);
	}

	private static double calculateAngle(Coordinate o1, Coordinate o2) {
		return fastAtan2(o2.getY() - o1.getY(), o2.getX() - o1.getX());
	}

	private static double angleDifference(double a1, double a2) {
		// calculate angle difference in clockwise directions as radians
		if ((a1 > 0 && a2 >= 0) && a1 > a2) {
			return Math.abs(a1 - a2);
		} else if ((a1 >= 0 && a2 > 0) && a1 < a2) {
			return 2 * Math.PI + a1 - a2;
		} else if ((a1 < 0 && a2 <= 0) && a1 < a2) {
			return 2 * Math.PI + a1 + Math.abs(a2);
		} else if ((a1 <= 0 && a2 < 0) && a1 > a2) {
			return Math.abs(a1 - a2);
		} else if (a1 <= 0 && 0 < a2) {
			return 2 * Math.PI + a1 - a2;
		} else if (a1 >= 0 && 0 >= a2) {
			return a1 + Math.abs(a2);
		} else {
			return 0.0;
		}
	}

	private static ArrayList<Coordinate> sortByAngle(ArrayList<Coordinate> l, Coordinate q, double a) {
		// Sort by angle descending
		Collections.sort(l, new Comparator<Coordinate>() {
			@Override
			public int compare(final Coordinate o1, final Coordinate o2) {
				double a1 = angleDifference(a, calculateAngle(q, o1));
				double a2 = angleDifference(a, calculateAngle(q, o2));
				if (a2 == a1) {
					return 0;
				} else {
					if (a2 > a1) {
						return 1;
					} else {
						return -1;
					}
				}
//				return a2.compareTo(a1);
			}
		});
		return l;
	}

	private boolean intersect(Coordinate l1p1, Coordinate l1p2, Coordinate l2p1, Coordinate l2p2) {
		// calculate part equations for line-line intersection
		double a1 = l1p2.getY() - l1p1.getY();
		double b1 = l1p1.getX() - l1p2.getX();
		double c1 = a1 * l1p1.getX() + b1 * l1p1.getY();
		double a2 = l2p2.getY() - l2p1.getY();
		double b2 = l2p1.getX() - l2p2.getX();
		double c2 = a2 * l2p1.getX() + b2 * l2p1.getY();
		// calculate the divisor
		double tmp = (a1 * b2 - a2 * b1);

		// calculate intersection point x coordinate
		double pX = (c1 * b2 - c2 * b1) / tmp;

		// check if intersection x coordinate lies in line line segment
		if ((pX > l1p1.getX() && pX > l1p2.getX()) || (pX > l2p1.getX() && pX > l2p2.getX())
				|| (pX < l1p1.getX() && pX < l1p2.getX()) || (pX < l2p1.getX() && pX < l2p2.getX())) {
			return false;
		}

		// calculate intersection point y coordinate
		double pY = (a1 * c2 - a2 * c1) / tmp;

		// check if intersection y coordinate lies in line line segment
		if ((pY > l1p1.getY() && pY > l1p2.getY()) || (pY > l2p1.getY() && pY > l2p2.getY())
				|| (pY < l1p1.getY() && pY < l1p2.getY()) || (pY < l2p1.getY() && pY < l2p2.getY())) {
			return false;
		}

		return true;
	}

	private boolean pointInPolygon(Coordinate p, ArrayList<Coordinate> pp) {
		boolean result = false;
		for (int i = 0, j = pp.size() - 1; i < pp.size(); j = i++) {
			if ((pp.get(i).getY() > p.getY()) != (pp.get(j).getY() > p.getY())
					&& (p.getX() < (pp.get(j).getX() - pp.get(i).getX()) * (p.getY() - pp.get(i).getY())
							/ (pp.get(j).getY() - pp.get(i).getY()) + pp.get(i).getX())) {
				result = !result;
			}
		}
		return result;
	}

	private static final double QRTR_PI = (0.25f * Math.PI);
	private static final double THREE_QRTR_PI = (0.75f * Math.PI);

	private static double fastAtan2(final double y, final double x) {

		double r, angle;
		final double abs_y = Math.abs(y) + 1e-10f; // kludge to prevent 0/0 condition

		if (x < 0.0f) {
			r = (x + abs_y) / (abs_y - x); // (3)
			angle = THREE_QRTR_PI; // (4)
		} else {
			r = (x - abs_y) / (x + abs_y); // (1)
			angle = QRTR_PI; // (2)
		}
		angle += (0.1963f * r * r - 0.9817f) * r; // (2 | 4)
		if (y < 0.0f)
			return (-angle); // negate if in quad III or IV
		else
			return (angle);
	}

	class Pair<T, U> {
		private final T key;
		private final U value;

		public Pair(T key, U value) {
			this.key = key;
			this.value = value;
		}

		public T getKey() {
			return this.key;
		}

		public U getValue() {
			return this.value;
		}
	}

}
