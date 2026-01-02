package micycle.pgs;

import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.List;

import micycle.pgs.commons.AreaOptimalPolygonizer;
import micycle.pgs.commons.AreaOptimalPolygonizer.AreaObjective;
import micycle.pgs.commons.Uncrossing2Opt;
import net.jafama.FastMath;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Generates simple polygonisations of point sets.
 * <p>
 * A polygonisation is a simple polygon whose vertex set is exactly the given
 * point set, i.e. a non-self-intersecting Hamiltonian cycle through all points.
 * Different algorithms may produce different polygonisations of the same point
 * set.
 * <p>
 * Polygonisations are distinct from geometric hulls: hulls may select a
 * <b>subset</b> of extreme points to form an enclosing boundary, whereas
 * polygonisations must use all points as vertices.
 *
 * @author Michael Carleton
 * @since 2.2
 */
public class PGS_Polygonisation {

	/**
	 * Produces a simple polygonisation that attempts to minimise the polygon area
	 * while using every point in the supplied set as a vertex.
	 * 
	 * @param points the input point set (must not be {@code null}) containing >2
	 *               points.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon that polygonises the input points and (attempts to) minimise
	 *         area.
	 * @see {@link #maxArea(Collection)}
	 * @since 2.2
	 */
	public static PShape minArea(Collection<PVector> points) {
		var coords = points.stream().map(p -> PGS.coordFromPVector(p)).toList();
		var g = AreaOptimalPolygonizer.polygonize(coords, AreaObjective.MINIMIZE);
		return toPShape(g);
	}

	/**
	 * Produces a simple polygonisation that attempts to maximise the polygon area
	 * while using every point in the supplied set as a vertex.
	 *
	 * @param points the input point set (must not be {@code null}) containing >2
	 *               points.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon that polygonises the input points and (attempts to) maximise
	 *         area.
	 * @see #minArea(Collection)
	 * @since 2.2
	 */
	public static PShape maxArea(Collection<PVector> points) {
		var coords = points.stream().map(p -> PGS.coordFromPVector(p)).toList();
		var g = AreaOptimalPolygonizer.polygonize(coords, AreaObjective.MAXIMIZE);
		return toPShape(g);
	}

	/**
	 * Computes a polygonisation that approximates a shortest closed tour visiting
	 * every point exactly once (a Hamiltonian cycle with small perimeter).
	 * <p>
	 * This method is effectively a TSP-style polygonisation: it returns a simple
	 * polygon whose total edge length is minimised (or approximated by the
	 * underlying shortest-tour routine).
	 *
	 * @param points the input point set (must not be {@code null}). If the set
	 *               contains fewer than three distinct points an appropriate
	 *               degenerate {@link processing.core.PShape PShape} containing the
	 *               input points will be returned.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon that attempts to minimise perimeter.
	 * @since 2.2
	 */
	public static PShape minPerimeter(Collection<PVector> points) {
		return PGS_PointSet.findShortestTour(points);
	}

	/**
	 * Builds a polygonisation by scanning points horizontally (primary sort by Y,
	 * secondary by X) and then removing edge crossings via a 2-opt (segment
	 * reversal) uncrossing routine.
	 * <p>
	 * The produced polygon has "horizontal scanline" characteristics.
	 *
	 * @param points the input point set (may be {@code null}). If {@code null} an
	 *               empty {@link processing.core.PShape PShape} is returned. If the
	 *               set contains fewer than three distinct points a degenerate
	 *               {@link processing.core.PShape PShape} containing the input
	 *               points is returned.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon constructed by horizontal scan and uncrossing.
	 * @since 2.2
	 */
	public static PShape horizontal(Collection<PVector> points) {
		return scanAndResolve(points, true);
	}

	/**
	 * Builds a polygonisation by scanning points vertically (primary sort by X,
	 * secondary by Y) and then removing edge crossings via a 2-opt (segment
	 * reversal) uncrossing routine.
	 * <p>
	 * The produced polygon has "vertical scanline" characteristics.
	 *
	 * @param points the input point set (may be {@code null}). If {@code null} an
	 *               empty {@link processing.core.PShape PShape} is returned. If the
	 *               set contains fewer than three distinct points a degenerate
	 *               {@link processing.core.PShape PShape} containing the input
	 *               points is returned.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon constructed by vertical scan and uncrossing.
	 * @since 2.2
	 */
	public static PShape vertical(Collection<PVector> points) {
		return scanAndResolve(points, false);
	}

	/**
	 * Produces a polygonisation by ordering points according to a Hilbert curve
	 * (space-filling curve) ordering, then applying a local uncrossing (2-opt) pass
	 * to remove any segment intersections.
	 * <p>
	 * Hilbert ordering tends to preserve locality and thus often produces visually
	 * compact, low-crossing initial orderings which the uncrossing step refines
	 * into a simple polygon.
	 *
	 * @param points the input point set (must not be {@code null}). If the set
	 *               contains fewer than three distinct points an appropriate
	 *               degenerate {@link processing.core.PShape PShape} containing the
	 *               input points will be returned.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon obtained from Hilbert ordering + uncrossing.
	 * @since 2.2
	 */
	public static PShape hilbert(Collection<PVector> points) {
		var seq = PGS_PointSet.hilbertSort(new ArrayList<PVector>(points));
		Uncrossing2Opt.uncross(seq);
		return toPolygon(seq);
	}

	/**
	 * Builds a polygonisation by grouping points into concentric "rings" around the
	 * centroid, ordering points within each ring by polar angle, and stitching the
	 * rings together into a single sequence.
	 * <p>
	 * This is a heuristic polygonisation: it favors "circular" or banded structures
	 * (concentric/clustered layouts) and often produces visually compact,
	 * low-crossing initial orders that the uncrossing step refines into a simple
	 * polygon.
	 *
	 * @param points the input point set (may be {@code null}). If {@code null} an
	 *               empty {@link processing.core.PShape PShape} is returned. If the
	 *               set contains fewer than three distinct points a degenerate
	 *               {@link processing.core.PShape PShape} containing the input
	 *               points is returned.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon constructed by concentric ring (circular) ordering and
	 *         subsequent uncrossing.
	 * @since 2.2
	 */
	public static PShape circular(Collection<PVector> points) {
		if (points == null) {
			return new PShape();
		}
		final int n = points.size();
		if (n < 3)
			return PGS_Conversion.fromPVector(new ArrayList<>(points));

		// center = centroid
		double cx = 0, cy = 0;
		for (PVector p : points) {
			cx += p.x;
			cy += p.y;
		}
		cx /= n;
		cy /= n;

		List<Info> info = new ArrayList<>(n);
		for (PVector p : points) {
			double dx = p.x - cx, dy = p.y - cy;
			info.add(new Info(p, Math.sqrt(dx * dx + dy * dy), FastMath.atan2(dy, dx)));
		}

		// sort by radius and split into rings (equal-size quantiles)
		Collections.sort(info, (a, b) -> Double.compare(a.r, b.r));
		int numRings = Math.max(1, (int) Math.round(Math.sqrt(n))); // heuristic
		List<List<Info>> rings = new ArrayList<>(numRings);
		for (int i = 0; i < numRings; i++)
			rings.add(new ArrayList<>());

		for (int i = 0; i < n; i++) {
			int bucket = (int) ((long) i * numRings / n); // maps 0..n-1 into 0..numRings-1
			rings.get(bucket).add(info.get(i));
		}

		// sort each ring by angle
		for (List<Info> ring : rings) {
			Collections.sort(ring, (a, b) -> Double.compare(a.theta, b.theta));
		}

		// concatenate rings, aligning each ring to the nearest start point and
		// alternating direction
		List<PVector> seq = new ArrayList<>(n);
		boolean forward = true;
		for (List<Info> ring : rings) {
			if (ring.isEmpty())
				continue;
			if (seq.isEmpty()) {
				// first ring: optionally start at smallest theta, and maybe reverse for parity
				if (!forward)
					Collections.reverse(ring);
				for (Info it : ring)
					seq.add(it.p);
			} else {
				// find index in ring nearest to last appended point
				PVector last = seq.get(seq.size() - 1);
				int start = 0;
				double best = Double.POSITIVE_INFINITY;
				for (int k = 0; k < ring.size(); k++) {
					double dx = last.x - ring.get(k).p.x;
					double dy = last.y - ring.get(k).p.y;
					double d2 = dx * dx + dy * dy;
					if (d2 < best) {
						best = d2;
						start = k;
					}
				}
				// append ring starting at start, in forward or reverse direction
				if (forward) {
					for (int k = 0; k < ring.size(); k++) {
						seq.add(ring.get((start + k) % ring.size()).p);
					}
				} else {
					for (int k = 0; k < ring.size(); k++) {
						int idx = (start - k) % ring.size();
						if (idx < 0)
							idx += ring.size();
						seq.add(ring.get(idx).p);
					}
				}
			}
			forward = !forward;
		}

		Uncrossing2Opt.uncross(seq);

		return toPolygon(seq);
	}

	/**
	 * Generates a polygonisation by angular (radial) sorting: points are sorted by
	 * angle around the centroid, tie-broken by distance from the centroid, and then
	 * a 2-opt uncrossing pass is applied.
	 * <p>
	 * The angular sort usually gives a star-shaped output.
	 *
	 * @param points the input point set (may be {@code null}). If {@code null} an
	 *               empty {@link processing.core.PShape PShape} is returned. If the
	 *               set contains fewer than three distinct points a degenerate
	 *               {@link processing.core.PShape PShape} containing the input
	 *               points is returned.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon constructed by radial sorting and uncrossing.
	 * @since 2.2
	 */
	public static PShape angular(Collection<PVector> points) {
		if (points == null) {
			return new PShape();
		}
		final int n = points.size();
		if (n == 0) {
			return PGS_Conversion.fromPVector(points);
		}
		if (n < 3) {
			return PGS_Conversion.fromPVector(new ArrayList<>(points));
		}

		// compute centroid as center for angular sort
		double cx = 0.0, cy = 0.0;
		for (PVector p : points) {
			cx += p.x;
			cy += p.y;
		}
		cx /= n;
		cy /= n;

		List<Info> info = new ArrayList<>(n);
		for (PVector p : points) {
			double dx = p.x - cx;
			double dy = p.y - cy;
			double theta = FastMath.atan2(dy, dx);
			double r = Math.sqrt(dx * dx + dy * dy);
			info.add(new Info(p, r, theta));
		}

		// sort by angle, tie-break by radius (closer first)
		Collections.sort(info, (a, b) -> {
			int c = Double.compare(a.theta, b.theta);
			if (c != 0)
				return c;
			return Double.compare(a.r, b.r);
		});

		List<PVector> seq = new ArrayList<>(n);
		for (Info it : info) {
			seq.add(it.p);
		}

		Uncrossing2Opt.uncross(seq);

		return toPolygon(seq);
	}

	/**
	 * Constructs a polygonisation using the "onion" (convex-layers) strategy:
	 * repeatedly peel convex hull layers (outermost first), stitch hull layers into
	 * a single cyclic order (alternating directions for continuity), insert any
	 * leftover points with a cheapest-insertion heuristic, and finally apply a
	 * 2-opt uncrossing pass.
	 * <p>
	 * This approach tends to respect global convex structure and produces
	 * spiral-like polygonisations that use all points as vertices.
	 *
	 * @param points the input point set (may be {@code null}). If {@code null} an
	 *               empty {@link processing.core.PShape PShape} is returned. If the
	 *               set contains fewer than three distinct points a degenerate
	 *               {@link processing.core.PShape PShape} containing the input
	 *               points is returned.
	 * @return a new {@link processing.core.PShape PShape} representing a simple
	 *         polygon constructed by convex-layer peeling, stitching and
	 *         uncrossing.
	 * @since 2.2
	 */
	public static PShape onion(Collection<PVector> points) {
		if (points == null)
			return new PShape();
		int n0 = points.size();
		if (n0 < 3)
			return PGS_Conversion.fromPVector(new ArrayList<>(points));

		// mutable working set
		List<PVector> remaining = new ArrayList<>(points);

		// peel convex hull layers
		List<List<PVector>> layers = new ArrayList<>();
		while (remaining.size() >= 3) {
			List<PVector> hull = convexHullMonotoneChain(remaining);
			if (hull.size() < 3)
				break; // degenerate (collinear etc.)
			layers.add(hull);

			// remove hull points (identity-based)
			var onHull = Collections.newSetFromMap(new IdentityHashMap<PVector, Boolean>());
			onHull.addAll(hull);
			remaining.removeIf(onHull::contains);
		}

		// stitch layers into one cyclic order (spiral-ish)
		List<PVector> seq = new ArrayList<>(points.size());
		boolean forward = true;

		for (List<PVector> layer : layers) {
			if (layer.isEmpty())
				continue;

			if (seq.isEmpty()) {
				if (!forward)
					java.util.Collections.reverse(layer);
				seq.addAll(layer);
			} else {
				PVector last = seq.get(seq.size() - 1);
				List<PVector> rotated = rotateToNearest(layer, last);
				if (!forward)
					Collections.reverse(rotated);
				seq.addAll(rotated);
			}
			forward = !forward;
		}

		// if anything left (0,1,2 points or collinear residue), insert cheaply
		for (PVector p : remaining) {
			insertCheapest(seq, p);
		}

		Uncrossing2Opt.uncross(seq);

		return toPolygon(seq);
	}

	/**
	 * Generic scan-based polygonisation. If primaryIsY is true, points are sorted
	 * primarily by Y then X (horizontal scanlines). Otherwise sorted primarily by X
	 * then Y (vertical scanlines). After sorting, a 2-opt style crossing removal is
	 * applied by iteratively reversing segments that cause segment intersections.
	 */
	private static PShape scanAndResolve(Collection<PVector> points, boolean primaryIsY) {
		// defensive handling
		if (points == null) {
			return new PShape();
		}

		final int n = points.size();
		if (n == 0) {
			return PGS_Conversion.fromPVector(points);
		}
		if (n < 3) {
			// trivial: nothing to polygonise
			return PGS_Conversion.fromPVector(new ArrayList<>(points));
		}

		// make a mutable copy
		List<PVector> seq = new ArrayList<>(points);

		// comparator depending on primary axis
		Comparator<PVector> cmp = primaryIsY ? Comparator.comparingDouble((PVector p) -> p.y).thenComparingDouble(p -> p.x)
				: Comparator.comparingDouble((PVector p) -> p.x).thenComparingDouble(p -> p.y);

		Collections.sort(seq, cmp);
		Uncrossing2Opt.uncross(seq);

		return toPolygon(seq);
	}

	private static List<PVector> convexHullMonotoneChain(List<PVector> pts) {
		// returns CCW hull without repeating the first point
		int n = pts.size();
		if (n < 3)
			return new ArrayList<>();

		// sort by x then y
		List<PVector> p = new ArrayList<>(pts);
		p.sort((a, b) -> {
			int cx = Float.compare(a.x, b.x);
			if (cx != 0)
				return cx;
			return Float.compare(a.y, b.y);
		});

		List<PVector> lower = new ArrayList<>();
		for (PVector v : p) {
			while (lower.size() >= 2 && cross(lower.get(lower.size() - 2), lower.get(lower.size() - 1), v) <= 0) {
				lower.remove(lower.size() - 1);
			}
			lower.add(v);
		}

		List<PVector> upper = new ArrayList<>();
		for (int i = p.size() - 1; i >= 0; i--) {
			PVector v = p.get(i);
			while (upper.size() >= 2 && cross(upper.get(upper.size() - 2), upper.get(upper.size() - 1), v) <= 0) {
				upper.remove(upper.size() - 1);
			}
			upper.add(v);
		}

		// remove last of each (it's the start of the other list)
		lower.remove(lower.size() - 1);
		upper.remove(upper.size() - 1);

		List<PVector> hull = new ArrayList<>(lower.size() + upper.size());
		hull.addAll(lower);
		hull.addAll(upper);
		return hull;
	}

	private static float cross(PVector o, PVector a, PVector b) {
		return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
	}

	private static List<PVector> rotateToNearest(List<PVector> ring, PVector target) {
		int m = ring.size();
		if (m == 0)
			return new ArrayList<>();

		int best = 0;
		double bestD2 = Double.POSITIVE_INFINITY;
		for (int i = 0; i < m; i++) {
			PVector p = ring.get(i);
			double dx = p.x - target.x;
			double dy = p.y - target.y;
			double d2 = dx * dx + dy * dy;
			if (d2 < bestD2) {
				bestD2 = d2;
				best = i;
			}
		}

		List<PVector> out = new ArrayList<>(m);
		for (int k = 0; k < m; k++)
			out.add(ring.get((best + k) % m));
		return out;
	}

	private static void insertCheapest(List<PVector> cycle, PVector p) {
		int n = cycle.size();
		if (n == 0) {
			cycle.add(p);
			return;
		}
		if (n == 1) {
			cycle.add(p);
			return;
		}

		int bestIdx = 0;
		double bestDelta = Double.POSITIVE_INFINITY;

		for (int i = 0; i < n; i++) {
			PVector a = cycle.get(i);
			PVector b = cycle.get((i + 1) % n);
			double delta = dist(a, p) + dist(p, b) - dist(a, b);
			if (delta < bestDelta) {
				bestDelta = delta;
				bestIdx = i + 1;
			}
		}
		cycle.add(bestIdx, p);
	}

	private static double dist(PVector a, PVector b) {
		double dx = a.x - b.x, dy = a.y - b.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	private static record Info(PVector p, double r, double theta) {
	}

	private static PShape toPolygon(List<PVector> points) {
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			points.add(points.get(0)); // close
		}
		return PGS_Conversion.fromPVector(points);
	}

}
