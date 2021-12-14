package micycle.pgs.commons;

import java.util.Collections;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.linearref.LengthIndexedLine;

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
	 * RAP C++ : https://github.com/catherinetaylor2/Shape_Interpolation/blob/master/rigid_interp.cpp
	 */

	private final CoordinateList from, to;

	public ShapeInterpolation(Geometry from, Geometry to) {
		this(((Polygon) from).getExteriorRing(), ((Polygon) to).getExteriorRing());
	}

	public ShapeInterpolation(LinearRing from, LinearRing to) {
		if (!Orientation.isCCW(from.getCoordinates())) {
			System.out.println("sdas");
			from = from.reverse();
		}
		if (!Orientation.isCCW(to.getCoordinates())) {
			to = to.reverse();
		}

		// find the "smaller" ring (as measured by number of vertices)
		LinearRing smaller, bigger;
		boolean smallerIsTo = false;
		if (from.getNumPoints() > to.getNumPoints()) {
			bigger = from;
			smaller = to;
			smallerIsTo = true;
		} else {
			bigger = to;
			smaller = from;
		}

		CoordinateList smallerLine, biggerLine;
		if (from.getNumPoints() == to.getNumPoints()) {
			// don't attempt to densify if same number of points
			smallerLine = new CoordinateList(smaller.getCoordinates());
			smallerLine.remove(smallerLine.size() - 1); // remove closing vertex
		} else {
			/*
			 * Sample the ring that has fewer vertices N times (in this case each new edge
			 * is the same length). An alternative approach would be to split existing edges
			 * (starting with the longest) until it has N vertices.
			 */
			final LengthIndexedLine l = new LengthIndexedLine(smaller);
			smallerLine = new CoordinateList(); // densified version of smaller

			int n = bigger.getNumPoints() - 1; // don't count closed vertex
			final double length = l.getEndIndex(); // perimeter length
			for (int i = 0; i < n; i++) {
				double index = i / (n - 1d);
				smallerLine.add(l.extractPoint(index * length));
			}
		}

		biggerLine = new CoordinateList(bigger.getCoordinates(), true);
		biggerLine.remove(biggerLine.size() - 1); // remove closing vertex
		int n = biggerLine.size(); // number of coords (unclosed)

		int bestOffset = 0;
		double min = Double.MAX_VALUE;

		/*
		 * The densified shape is rotated until the squared distance between point pairs
		 * of the two shapes is minimised.
		 */
		for (int offset = 0; offset < n; offset += 2) { // NOTE +=2
			double sumOfSquares = 0;

			for (int i = 0; i < n; i++) {
				sumOfSquares += distSq(smallerLine.get((offset + i) % n), biggerLine.get(i));
			}

			if (sumOfSquares < min) {
				min = sumOfSquares;
				bestOffset = offset;
				if (sumOfSquares == 0) {
					break; // break early when shapes are identical
				}
			}
		}

		if (bestOffset != 0) {
			Collections.rotate(smallerLine, -bestOffset);
		}

		smallerLine.closeRing();
		biggerLine.closeRing();
		if (smallerIsTo) {
			this.to = smallerLine;
			this.from = biggerLine;
		} else {
			this.from = smallerLine;
			this.to = biggerLine;
		}
	}

	public Coordinate[] tween(double t) {
		CoordinateList morph = new CoordinateList();
		for (int i = 0; i < from.size(); i++) {
			morph.add(lerp(from.get(i), to.get(i), t), true);
		}
		return morph.toCoordinateArray();
	}

	private static Coordinate lerp(Coordinate from, Coordinate to, double t) {
		return new Coordinate(from.x + (to.x - from.x) * t, from.y + (to.y - from.y) * t);
	}

	private static double distSq(Coordinate a, Coordinate b) {
		double dx = a.x - b.x;
		double dy = a.y - b.y;
		return (dx * dx + dy * dy);
	}

}
