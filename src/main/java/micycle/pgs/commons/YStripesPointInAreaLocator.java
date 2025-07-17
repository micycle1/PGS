package micycle.pgs.commons;

import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;

/**
 * A fast point-in-polygon locator for {@link Polygon} geometries using the
 * YStripes spatial index.
 * <p>
 * This locator efficiently determines the topological location (interior,
 * exterior, or boundary) of a point relative to a polygon, including its holes,
 * by leveraging the <a href=
 * "https://github.com/tidwall/tg/blob/main/docs/POLYGON_INDEXING.md#ystripes">YStripes
 * algorithm</a>.
 * <p>
 * The YStripes index is constructed separately for the exterior ring and each
 * interior ring (hole) of the polygon. For query points:
 * <ul>
 * <li>If the point lies outside the exterior ring or on its boundary, the
 * result is immediately returned.</li>
 * <li>Otherwise, each hole is checked: if the point is inside a hole, the point
 * is considered exterior to the polygon; if on a hole's boundary, the result is
 * {@link Location#BOUNDARY}.</li>
 * <li>If outside all holes, the point is interior to the polygon.</li>
 * </ul>
 * <p>
 * For more details, see the <a href=
 * "https://github.com/tidwall/tg/blob/main/docs/POLYGON_INDEXING.md#ystripes">
 * YStripes: Polygon Indexing in 'tg' Library</a>.
 *
 * @author Michael Carleton
 * @see YStripesPointInRingLocator
 * @see PointOnGeometryLocator
 */
public final class YStripesPointInAreaLocator implements PointOnGeometryLocator {

	private final YStripesPointInRingLocator exterior;
	private final YStripesPointInRingLocator[] holes;

	/**
	 * Constructs a YStripesPointInAreaLocator for the given polygon.
	 *
	 * @param poly the polygon geometry to index; must not be {@code null} or empty
	 * @throws IllegalArgumentException if the polygon is {@code null} or empty
	 */
	public YStripesPointInAreaLocator(final Polygon poly) {
		if (poly == null || poly.isEmpty()) {
			throw new IllegalArgumentException("polygon null / empty");
		}

		exterior = new YStripesPointInRingLocator(poly.getExteriorRing()); // full-blown heuristic
		final int nh = poly.getNumInteriorRing();
		holes = new YStripesPointInRingLocator[nh];
		for (int i = 0; i < nh; i++) {
			final LinearRing ring = poly.getInteriorRingN(i);
			holes[i] = new YStripesPointInRingLocator(ring);
		}
	}

	/**
	 * Determines the location of a point relative to the indexed polygon.
	 *
	 * @param p the coordinate to locate
	 * @return {@link Location#INTERIOR} if the point is inside the polygon (but not
	 *         in a hole), {@link Location#EXTERIOR} if outside the polygon or
	 *         inside a hole, {@link Location#BOUNDARY} if on the boundary of the
	 *         exterior ring or any hole
	 */
	@Override
	public int locate(final Coordinate p) {
		final int ext = exterior.locate(p);
		if (ext == Location.EXTERIOR || ext == Location.BOUNDARY) {
			return ext;
		}

		for (final YStripesPointInRingLocator h : holes) {
			final int hl = h.locate(p);
			if (hl == Location.EXTERIOR) {
				continue; // not in hole
			}

			if (hl == Location.BOUNDARY) {
				return Location.BOUNDARY;
			}
			/* inside a hole → outside the polygon */
			return Location.EXTERIOR;
		}
		return Location.INTERIOR;
	}
}