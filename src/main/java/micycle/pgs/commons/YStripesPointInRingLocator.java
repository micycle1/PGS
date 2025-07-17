package micycle.pgs.commons;

import java.util.Arrays;

import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;

/**
 * An implementation of the "YStripes" spatial index for fast point-in-polygon
 * tests.
 * <p>
 * The YStripes structure partitions a polygon into horizontal "stripes". Each
 * stripe stores a list of polygon segment indices that intersect it. This
 * allows for O(1) lookup of potentially intersecting segments for any given
 * Y-coordinate, drastically reducing the number of segments that need to be
 * checked for a point-in-polygon test.
 * <p>
 * For more details, see: <a href=
 * "https://github.com/tidwall/tg/blob/main/docs/POLYGON_INDEXING.md#ystripes">
 * YStripes: Polygon Indexing in 'tg' Library</a>
 * 
 * @author Michael Carleton
 */
final class YStripesPointInRingLocator implements PointOnGeometryLocator {

	// immutable, SoA-style geometry data
	private final double[] x1, y1, x2, y2; // segment end points
	private final double[] xmin, xmax, ymin, ymax; // per-segment bounding boxes

	// stripe index
	private final int[] stripeOfs; // length = nStripes + 1 (prefix sums)
	private final int[] stripeSeg; // flat list with all segment indices

	private final double minX, minY, maxX, maxY; // cached envelope
	private final int nStripes;
	private final double invHeight; // nStripes / (maxY - minY)

	YStripesPointInRingLocator(final LinearRing ring) {
		final Coordinate[] c = ring.getCoordinates();
		final int nSeg = c.length - 1;

		x1 = new double[nSeg];
		y1 = new double[nSeg];
		x2 = new double[nSeg];
		y2 = new double[nSeg];

		xmin = new double[nSeg];
		xmax = new double[nSeg];
		ymin = new double[nSeg];
		ymax = new double[nSeg];

		for (int i = 0; i < nSeg; i++) {
			final Coordinate a = c[i], b = c[i + 1];

			final double ax = a.x, ay = a.y;
			final double bx = b.x, by = b.y;

			x1[i] = ax;
			y1[i] = ay;
			x2[i] = bx;
			y2[i] = by;

			if (ax < bx) {
				xmin[i] = ax;
				xmax[i] = bx;
			} else {
				xmin[i] = bx;
				xmax[i] = ax;
			}
			if (ay < by) {
				ymin[i] = ay;
				ymax[i] = by;
			} else {
				ymin[i] = by;
				ymax[i] = ay;
			}
		}

		final Envelope env = ring.getEnvelopeInternal();
		minX = env.getMinX();
		maxX = env.getMaxX();
		minY = env.getMinY();
		maxY = env.getMaxY();

		// stripe count heuristic (Polsby–Popper)
		final double area = ring.getArea();
		final double length = ring.getLength();
		final double pp = area == 0 ? 0 : (4.0 * Math.PI * area) / (length * length);

		final int stripesWanted = Math.max(32, (int) (nSeg * pp));
		// round up to next power-of-two → cheaper maths later on
		nStripes = Integer.highestOneBit(stripesWanted - 1) << 1;

		final double height = maxY - minY;
		invHeight = height == 0 ? 0.0 : nStripes / height;

		// FIRST PASS: count how many segs per stripe
		final int[] cnt = new int[nStripes];
		for (int i = 0; i < nSeg; i++) {
			final int sMin = stripe(ymin[i]);
			final int sMax = stripe(ymax[i]);
			for (int s = sMin; s <= sMax; s++) {
				cnt[s]++;
			}
		}

		// build CSR arrays
		stripeOfs = new int[nStripes + 1];
		int total = 0;
		for (int s = 0; s < nStripes; s++) {
			stripeOfs[s] = total;
			total += cnt[s];
		}
		stripeOfs[nStripes] = total;
		stripeSeg = new int[total];

		// reuse cnt[] as cursors
		Arrays.fill(cnt, 0);
		for (int i = 0; i < nSeg; i++) {
			final int sMin = stripe(ymin[i]);
			final int sMax = stripe(ymax[i]);
			for (int s = sMin; s <= sMax; s++) {
				stripeSeg[stripeOfs[s] + cnt[s]++] = i;
			}
		}
	}

	@Override
	public int locate(final Coordinate p) {
		// envelope quick-reject
		if (p.x < minX || p.x > maxX || p.y < minY || p.y > maxY) {
			return Location.EXTERIOR;
		}

		boolean in = false;
		boolean onEdge = false;

		final int s = stripe(p.y);

		for (int k = stripeOfs[s], end = stripeOfs[s + 1]; k < end && !onEdge; k++) {

			final int i = stripeSeg[k];

			if (p.y < ymin[i] || p.y > ymax[i]) {
				continue;
			}

			if (p.x < xmin[i]) { // ray surely crosses
				if (p.y != ymin[i] && p.y != ymax[i]) {
					in ^= true;
					continue;
				}
			} else if (p.x > xmax[i]) { // ray surely misses
				if (ymin[i] != ymax[i] && xmin[i] != xmax[i]) {
					continue;
				}
			}

			switch (rayTest(i, p)) {
				case 0 :
					in ^= true;
					break; // IN (ray crosses segment)
				case 2 :
					onEdge = true;
					break; // ON (point sits on segment)
			}
		}
		return onEdge ? Location.BOUNDARY : (in ? Location.INTERIOR : Location.EXTERIOR);
	}

	/** Maps a Y value to its stripe index (clamped). */
	private int stripe(final double y) {
		final int s = (int) ((y - minY) * invHeight);
		return (s < 0) ? 0 : (s >= nStripes ? nStripes - 1 : s);
	}

	/**
	 * Very small, branch-light ray/segment test.
	 *
	 * return 0 = crosses (toggle), 1 = miss, 2 = point on segment.
	 */
	private int rayTest(final int i, final Coordinate p) {
		double ax = x1[i], ay = y1[i];
		double bx = x2[i], by = y2[i];

		if (ay > by) { // order by Y
			double t;
			t = ax;
			ax = bx;
			bx = t;
			t = ay;
			ay = by;
			by = t;
		}

		// vertex hit
		if ((p.x == ax && p.y == ay) || (p.x == bx && p.y == by)) {
			return 2;
		}

		// outside vertical range
		if (p.y < ay || p.y > by) {
			return 1;
		}

		// horizontal segment
		if (ay == by) {
			return (p.x >= xmin[i] && p.x <= xmax[i]) ? 2 : 1;
		}

		double y = p.y;
		if (y == ay || y == by) {
			y = Math.nextUp(y); // nudge off vertices
		}

		if (p.x >= Math.max(ax, bx)) {
			return 1;
		}
		if (p.x < Math.min(ax, bx)) {
			return 0;
		}

		// cross-product sign
		final double cross = (y - ay) * (bx - ax) - (by - ay) * (p.x - ax);
		if (cross == 0) {
			return 2;
		}
		return cross > 0 ? 0 : 1;
	}
}