package micycle.pgs.commons;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.overlayng.OverlayNGRobust;
import org.locationtech.jts.operation.union.CascadedPolygonUnion;

import com.github.micycle1.geoblitz.HPRtreeX;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Determines what remains visible for each geometry in a draw-ordered list by
 * removing any polygonal areas covered by geometries drawn above it.
 * <p>
 * For each element {@code i}, returns
 * {@code geom[i] \ union(polygonalPart(geom[j]) for j>i)}. Non-polygonal
 * geometries can be targets, but only polygonal parts act as occluders.
 * 
 * @author Michael Carleton
 */
public final class OcclusionSubtract {

	private final List<? extends Geometry> drawOrder;
	private final int n;

	private final int smallK;
	private final boolean robust;

	private final HPRtreeX<OccItem> index;
	private final OccItem[] occItems; // keyed by draw-order idx (null if non-polygonal/empty)
	private final PreparedGeometry[] preparedOcc; // prepared occluders (polygonal only)
	private final PreparedGeometryFactory prepFactory = new PreparedGeometryFactory();

	/**
	 * Builds an occlusion subtractor for the given bottom-to-top draw order. Uses
	 * {@code smallK=8} and non-robust overlay.
	 *
	 * @param drawOrder geometries in draw order (index 0 is bottom, n-1 is top)
	 */
	public OcclusionSubtract(List<? extends Geometry> drawOrder) {
		this(drawOrder, 8, false);
	}

	/**
	 * Builds an occlusion subtractor with tuning parameters.
	 *
	 * @param drawOrder geometries in draw order (index 0 is bottom, n-1 is top)
	 * @param smallK    threshold for iterative subtraction vs. unioning occluders
	 * @param robust    whether to use {@link OverlayNGRobust} for differences
	 */
	public OcclusionSubtract(List<? extends Geometry> drawOrder, int smallK, boolean robust) {
		this.drawOrder = drawOrder;
		this.n = drawOrder.size();
		this.smallK = smallK;
		this.robust = robust;

		this.index = new HPRtreeX<>();
		this.occItems = new OccItem[n];
		this.preparedOcc = new PreparedGeometry[n];

		// Build tree over polygonal occluders only (areal masking).
		for (int i = 0; i < n; i++) {
			Geometry g = drawOrder.get(i);
			if (g == null || g.isEmpty()) {
				continue;
			}

			Geometry area = polygonalPart(g);
			if (area == null || area.isEmpty()) {
				continue;
			}

			OccItem it = new OccItem(i, area);
			occItems[i] = it;
			index.insert(it.env, it);
		}
		index.build();

		// Eagerly prepare occluders so parallel runs are race-free.
		// (Only polygonal occluders were inserted into occItems.)
		for (int i = 0; i < n; i++) {
			OccItem it = occItems[i];
			if (it != null) {
				preparedOcc[i] = prepFactory.create(it.area);
			}
		}
	}

	/**
	 * Builds an occlusion subtractor from a geometry container.
	 * <p>
	 * If {@code drawOrderGeom} is a GeometryCollection/Multi*, its components are
	 * used in internal order: {@code geometryN(0)} is bottom,
	 * {@code geometryN(n-1)} is top. If it is a single geometry, it becomes a
	 * 1-element draw list.
	 *
	 * @param drawOrderGeom single geometry or a collection of geometries
	 */
	public OcclusionSubtract(Geometry drawOrderGeom) {
		this(drawOrderGeom, 8, false);
	}

	/**
	 * Geometry constructor with tuning.
	 *
	 * @param drawOrderGeom single geometry or a collection of geometries
	 * @param smallK        threshold for iterative subtraction vs. unioning
	 *                      occluders
	 * @param robust        whether to use {@link OverlayNGRobust} for differences
	 */
	public OcclusionSubtract(Geometry drawOrderGeom, int smallK, boolean robust) {
		this(toDrawList(drawOrderGeom), smallK, robust);
	}

	/**
	 * Flattens a geometry into a draw list using {@code getGeometryN(i)} order.
	 *
	 * @param g a single geometry or a geometry collection (may be null)
	 * @return list of component geometries (possibly empty)
	 */
	private static List<Geometry> toDrawList(Geometry g) {
		if (g == null) {
			return List.of();
		}
		final int m = g.getNumGeometries(); // works for single geometry too (returns 1)
		List<Geometry> list = new ArrayList<>(m);
		for (int i = 0; i < m; i++) {
			list.add(g.getGeometryN(i));
		}
		return list;
	}

	/**
	 * Computes visible geometries by subtracting polygonal overlap of higher items
	 * from lower items.
	 * <p>
	 * This method is parallel: each output element is computed independently.
	 * Returns a list aligned with the input draw order (same size; elements may be
	 * null).
	 *
	 * @return visible geometries after areal occlusion subtraction
	 */
	public List<Geometry> subtractArealOcclusion() {
		final Geometry[] out = new Geometry[n];

		IntStream.range(0, n).parallel().forEach(i -> {
			Geometry target = drawOrder.get(i);
			if (target == null) {
				out[i] = null;
				return;
			}
			if (target.isEmpty()) {
				out[i] = target;
				return;
			}

			// Per-task buffers (no sharing)
			final List<OccItem> candidatesLocal = new ArrayList<>(32);
			final List<Geometry> unionLocal = new ArrayList<>(32);

			final Envelope env = target.getEnvelopeInternal();

			index.query(env).forEach(it -> {
				if (it.idx > i) {
					candidatesLocal.add(it);
				}
			});

			if (candidatesLocal.isEmpty()) {
				out[i] = target;
				return;
			}

			Geometry visible;
			if (candidatesLocal.size() == 1) {
				visible = subtractSinglePrepared(target, candidatesLocal.get(0));
			} else if (candidatesLocal.size() <= smallK) {
				visible = subtractIterativePrepared(target, candidatesLocal);
			} else {
				for (int k = 0; k < candidatesLocal.size(); k++) {
					unionLocal.add(candidatesLocal.get(k).area);
				}
				final Geometry union = (unionLocal.size() == 1) ? unionLocal.get(0) : CascadedPolygonUnion.union(unionLocal);
				visible = difference(target, union);
			}

			out[i] = visible;
		});

		return Arrays.stream(out).filter(java.util.Objects::nonNull).toList();
	}

	/**
	 * Subtracts a single prepared occluder from a target, with a cheap full-cover
	 * short-circuit.
	 */
	private Geometry subtractSinglePrepared(Geometry target, OccItem occ) {
		// Cheap full-cover short-circuit
		if (occ.env.contains(target.getEnvelopeInternal())) {
			if (preparedOcc[occ.idx].covers(target)) {
				return target.getFactory().createGeometryCollection(); // empty
			}
		}
		return difference(target, occ.area);
	}

	/**
	 * Iteratively subtracts prepared occluders, skipping candidates that cannot
	 * intersect.
	 */
	private Geometry subtractIterativePrepared(Geometry target, List<OccItem> candidates) {
		Geometry rem = target;
		for (int i = 0; i < candidates.size(); i++) {
			OccItem occ = candidates.get(i);
			if (!occ.env.intersects(rem.getEnvelopeInternal())) {
				continue;
			}
			if (!preparedOcc[occ.idx].intersects(rem)) {
				continue;
			}
			rem = difference(rem, occ.area);
			if (rem.isEmpty()) {
				return rem;
			}
		}
		return rem;
	}

	/**
	 * Computes {@code a \ b} using OverlayNG, optionally via the robust variant.
	 */
	private Geometry difference(Geometry a, Geometry b) {
		return robust ? OverlayNGRobust.overlay(a, b, OverlayNG.DIFFERENCE) : OverlayNG.overlay(a, b, OverlayNG.DIFFERENCE);
	}

	/**
	 * Extracts the polygonal (areal) part of a geometry; returns an empty polygon
	 * if none.
	 *
	 * @param g input geometry
	 * @return polygonal geometry comprising all polygon components
	 */
	private static Geometry polygonalPart(Geometry g) {
		if (g instanceof Polygonal) {
			return g;
		}

		GeometryFactory gf = g.getFactory();
		List<Polygon> polys = new ArrayList<>();
		collectPolygons(g, polys);

		if (polys.isEmpty()) {
			return gf.createPolygon(); // empty polygonal
		}
		if (polys.size() == 1) {
			return polys.get(0);
		}
		return gf.createMultiPolygon(polys.toArray(new Polygon[0]));
	}

	/**
	 * Recursively collects {@link Polygon} components from a (possibly nested)
	 * geometry.
	 */
	private static void collectPolygons(Geometry g, List<Polygon> out) {
		if (g == null || g.isEmpty()) {
			return;
		}
		if (g instanceof Polygon p) {
			out.add(p);
			return;
		}
		for (int i = 0; i < g.getNumGeometries(); i++) {
			collectPolygons(g.getGeometryN(i), out);
		}
	}

	/**
	 * Indexed occluder record (polygonal part only), keyed by draw-order index.
	 */
	private static final class OccItem {
		final int idx;
		final Geometry area; // polygonal part only
		final Envelope env;

		OccItem(int idx, Geometry area) {
			this.idx = idx;
			this.area = area;
			this.env = area.getEnvelopeInternal();
		}
	}
}