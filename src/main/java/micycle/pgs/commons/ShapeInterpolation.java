package micycle.pgs.commons;

import java.util.Collections;
import java.util.SplittableRandom;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;

/**
 * Best-guess interpolation between any two linear rings.
 * 
 * @author Michael Carleton
 *
 */
public class ShapeInterpolation {

	/*-
	 * TODO:
	 * Cut shape with holes (hairline gap) to ensure each geometry has genus=1.
	 * Cuts are discussed in 'Polygon Vertex Set Matching Algorithm for Shapefile Tweening'
	 * Break-up single polygon to interpolate with a multipolygon.
	 * See https://github.com/veltman/openvis/blob/master/README.md
	 * See 'Guaranteed intersection-free polygon morphing'
	 * See CGAl Shape Deformation: https://doc.cgal.org/latest/Barycentric_coordinates_2/index.html#title10
	 * https://homepages.inf.ed.ac.uk/tkomura/cav/presentation14_2018.pdf
	 * RAP C++ : https://github.com/catherinetaylor2/Shape_Interpolation/blob/master/rigid_interp.cpp
	 * and https://github.com/deliagander/ARAPShapeInterpolation
	 */

	private final CoordinateList from, to;

	public ShapeInterpolation(Geometry from, Geometry to) {
		this(((Polygon) from).getExteriorRing(), ((Polygon) to).getExteriorRing());
	}

	public ShapeInterpolation(LinearRing from, LinearRing to) {
		if (!Orientation.isCCW(from.getCoordinates())) {
			from = from.reverse();
		}
		if (!Orientation.isCCW(to.getCoordinates())) {
			to = to.reverse();
		}

		// find the "smaller" ring (as measured by number of vertices)
		// NOTE use Kabsch algorithm?
		CoordinateList smaller, bigger;
		boolean smallerIsTo = false;
		if (from.getNumPoints() > to.getNumPoints()) {
			bigger = new CoordinateList(from.getCoordinates(), false);
			smaller = new CoordinateList(to.getCoordinates(), false);
			smallerIsTo = true;
		} else {
			bigger = new CoordinateList(to.getCoordinates(), false);
			smaller = new CoordinateList(from.getCoordinates(), false);
		}

		bigger.closeRing(); // ensure closed
		smaller.closeRing(); // ensure closed
		smaller.remove(smaller.size() - 1); // unclose (to be closed later, after array rotation)
		bigger.remove(bigger.size() - 1); // unclose (to be closed later, after array rotation)

		// densify smaller list
		final int diff = bigger.size() - smaller.size();
		SplittableRandom r = new SplittableRandom(1337);
		for (int i = 0; i < diff; i++) {
			int index = r.nextInt(smaller.size() - 1);
			Coordinate a = smaller.get(index);
			Coordinate c = smaller.get(index + 1);
			Coordinate b = new Coordinate((a.x + c.x) / 2, (a.y + c.y) / 2); // midpoint
			smaller.add(index + 1, b); // insert b between a and c (shift c onwards right)
		}

		/*
		 * The densified shape is rotated until the squared distance between point pairs
		 * of the two shapes is minimised.
		 */
		int bestOffset = findBestRotation(smaller, bigger);

		if (bestOffset != 0) {
			Collections.rotate(smaller, -bestOffset);
		}

		smaller.closeRing();
		bigger.closeRing();
		if (smallerIsTo) {
			this.to = smaller;
			this.from = bigger;
		} else {
			this.from = smaller;
			this.to = bigger;
		}
	}

	public Coordinate[] tween(double t) {
		if (t == 0) {
			return from.toCoordinateArray();
		} else if (t == 1) {
			return to.toCoordinateArray();
		}
		t %= 1;
		CoordinateList morph = new CoordinateList();
		for (int i = 0; i < from.size(); i++) {
			morph.add(lerp(from.get(i), to.get(i), t), true);
		}
		return morph.toCoordinateArray();
	}

	/**
	 * @return a rotation offset for list <code>a</code> that minimises the squared
	 *         distance between all point pairs
	 */
	private static int findBestRotation(CoordinateList a, CoordinateList b) {
		/*
		 * Ternary search optimisation. Will converge to the global best rotation if the
		 * rotation-distance "function" is unimodal. Not sure if it is, but seems so in
		 * practice.
		 */
		final int n = a.size();
		int low = 0;
		int high = n - 1;

		while (low < high) {
			int mid1 = low + (high - low) / 3;
			int mid2 = high - (high - low) / 3;
			double dist1 = calculateSumOfSquares(a, b, mid1, n);
			double dist2 = calculateSumOfSquares(a, b, mid2, n);

			if (dist1 < dist2) {
				high = mid2 - 1;
			} else {
				low = mid1 + 1;
			}
		}

		return low;
	}

	private static double calculateSumOfSquares(CoordinateList a, CoordinateList b, int offset, int n) {
		double sumOfSquares = 0;
		for (int i = 0; i < n; i++) {
			sumOfSquares += a.get((offset + i) % n).distanceSq(b.get(i));
		}
		return sumOfSquares;
	}

	private static Coordinate lerp(Coordinate from, Coordinate to, double t) {
		return new Coordinate(from.x + (to.x - from.x) * t, from.y + (to.y - from.y) * t);
	}

}
