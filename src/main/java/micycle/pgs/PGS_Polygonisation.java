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
import micycle.pgs.commons.PEdge;
import micycle.pgs.commons.Uncrossing2Opt;
import net.jafama.FastMath;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Generates simple polygonisations of point sets.
 * <p>
 * A polygonisation is a simple polygon whose vertex set is exactly the given
 * point set, i.e. a non-self-intersecting Hamiltonian cycle through all points.
 * Different algorithms may produce different polygonizations of the same point
 * set.
 * <p>
 * This class includes polygonization variants that optimize geometric
 * objectives such as minimizing or maximizing the enclosed area, as well as
 * heuristic constructions based on ordering, proximity, or tour improvement.
 * <p>
 * Polygonizations are distinct from geometric hulls: hulls may select a
 * <b>subset</b> of extreme points to form an enclosing boundary, whereas
 * polygonizations use all points as vertices.
 *
 * @author Michael Carleton
 * @since 2.2
 */
public class PGS_Polygonisation {

	// https://discourse.processing.org/t/shape-generator-help-armin-hofmanns-rubber-band-shape-generator/33190/40

	// methods to create polygons from point sets
	// min/max area (heuristic)
	// tsp (min perimeter)
	// angular sort
	// horizontal sort (2-opt)
	// circular sort (2-opt)
	// Hamiltonian cycle?

	public static PShape minArea(Collection<PVector> points) {
		var coords = points.stream().map(p -> PGS.coordFromPVector(p)).toList();
		var g = AreaOptimalPolygonizer.polygonize(coords, AreaObjective.MINIMIZE);
		return toPShape(g);
	}

	public static PShape maxArea(Collection<PVector> points) {
		var coords = points.stream().map(p -> PGS.coordFromPVector(p)).toList();
		var g = AreaOptimalPolygonizer.polygonize(coords, AreaObjective.MAXIMIZE);
		return toPShape(g);
	}

	public static PShape minPerimeter(Collection<PVector> points) {
		return PGS_PointSet.findShortestTour(points);
	}

	/**
	 * Creates a polygonisation by scanning horizontally (i.e. sort primarily by Y,
	 * then X) then removing crossings with a 2-opt (segment reversal) routine.
	 */
	public static PShape horizontal(Collection<PVector> points) {
		// horizontal scanlines
		return scanAndResolve(points, true);
	}

	/**
	 * Creates a polygonisation by scanning vertically (i.e. sort primarily by X,
	 * then Y) then removing crossings with a 2-opt (segment reversal) routine.
	 */
	public static PShape vertical(Collection<PVector> points) {
//		vertical scanlines
		return scanAndResolve(points, false);
	}

	public static PShape hilbert(Collection<PVector> points) {
		var seq = PGS_PointSet.hilbertSort(new ArrayList<PVector>(points));
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
		Comparator<PVector> cmp;
		if (primaryIsY) {
			cmp = new Comparator<PVector>() {
				@Override
				public int compare(PVector a, PVector b) {
					if (a.y < b.y)
						return -1;
					if (a.y > b.y)
						return 1;
					if (a.x < b.x)
						return -1;
					if (a.x > b.x)
						return 1;
					return 0;
				}
			};
		} else {
			cmp = new Comparator<PVector>() {
				@Override
				public int compare(PVector a, PVector b) {
					if (a.x < b.x)
						return -1;
					if (a.x > b.x)
						return 1;
					if (a.y < b.y)
						return -1;
					if (a.y > b.y)
						return 1;
					return 0;
				}
			};
		}

		Collections.sort(seq, cmp);
		Uncrossing2Opt.uncross(seq);

		return toPolygon(seq);
	}

	/**
	 * instead of sorting purely by angle (angular/radial sort = sort by θ around a
	 * chosen center), "circular sorting" usually means grouping points into
	 * concentric rings (or convex layers) and then ordering points inside each ring
	 * by angle and stitching the rings together. The result is loopy/circular bands
	 * rather than one long star-shaped sweep.
	 * 
	 * @param points
	 * @return
	 */
	public static PShape circular(Collection<PVector> points) {
		// (concentric rings with simple stitching + 2‑opt)
		// other options: onion convex layers; spiral variant
//	    if (points == null) return PGS_Conversion.fromPVector(null);
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
	 * Angular / radial sort: sort points by angle around the centroid (atan2),
	 * tie-break by radius. Then run a 2-opt crossing removal.
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
			double r = FastMath.hypot(dx, dy);
			info.add(new Info(p, r, theta));
		}

		// sort by angle, tie-break by radius (closer first)
		Collections.sort(info, new Comparator<Info>() {
			@Override
			public int compare(Info a, Info b) {
				int c = Double.compare(a.theta, b.theta);
				if (c != 0)
					return c;
				return Double.compare(a.r, b.r);
			}
		});

		List<PVector> seq = new ArrayList<>(n);
		for (Info it : info) {
			seq.add(it.p);
		}

		Uncrossing2Opt.uncross(seq);

		return toPolygon(seq);
	}

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

		// uncross
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

	// Reverse elements in list from index i to j inclusive
	private static void reverse(List<PVector> list, int i, int j) {
		while (i < j) {
			PVector tmp = list.get(i);
			list.set(i, list.get(j));
			list.set(j, tmp);
			i++;
			j--;
		}
	}

	private static void uncross2Opt(List<PVector> seq) {
	    int n = seq.size();
	    if (n < 4) return;

	    PVector[] arr = seq.toArray(new PVector[n]);
	    
	    // Cache coordinates in primitive arrays (better CPU cache)
	    double[] x = new double[n];
	    double[] y = new double[n];
	    for (int k = 0; k < n; k++) {
	        x[k] = arr[k].x;
	        y[k] = arr[k].y;
	    }
	    
	    boolean improved = true;
	    while (improved) {
	        improved = false;
	        
	        for (int i = 0; i < n - 2; i++) {
	            double ax = x[i], ay = y[i];
	            double bx = x[i+1], by = y[i+1];
	            double abx = bx - ax;
	            double aby = by - ay;
	            
	            int jMax = (i == 0) ? n - 1 : n;
	            
	            for (int j = i + 2; j < jMax; j++) {
	                int j1 = (j + 1) % n;
	                
	                double cx = x[j], cy = y[j];
	                double dx = x[j1], dy = y[j1];
	                
	                double o1 = abx * (cy - ay) - aby * (cx - ax);
	                double o2 = abx * (dy - ay) - aby * (dx - ax);
	                
	                if ((o1 > 0 && o2 > 0) || (o1 < 0 && o2 < 0)) continue;
	                
	                double cdx = dx - cx;
	                double cdy = dy - cy;
	                double o3 = cdx * (ay - cy) - cdy * (ax - cx);
	                double o4 = cdx * (by - cy) - cdy * (bx - cx);
	                
	                if ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0)) {
	                    // Reverse coordinates
	                    reverseDoubles(x, i + 1, j);
	                    reverseDoubles(y, i + 1, j);
	                    improved = true;
	                    j = i + 1;
	                }
	            }
	        }
	    }
	    
	    // Rebuild PVectors
	    for (int i = 0; i < n; i++) {
	        arr[i].x = (float)x[i];
	        arr[i].y = (float)y[i];
	    }
	}

	private static void reverseDoubles(double[] arr, int start, int end) {
	    while (start < end) {
	        double tmp = arr[start];
	        arr[start] = arr[end];
	        arr[end] = tmp;
	        start++;
	        end--;
	    }
	}

	private static void reverseArray(PVector[] arr, int start, int end) {
	    while (start < end) {
	        PVector tmp = arr[start];
	        arr[start] = arr[end];
	        arr[end] = tmp;
	        start++;
	        end--;
	    }
	}

	private static boolean segmentsIntersect(final PVector a, final PVector b, final PVector c, final PVector d) {
		// Compute deltas once
		double abx = b.x - a.x;
		double aby = b.y - a.y;
		double acx = c.x - a.x;
		double acy = c.y - a.y;
		double adx = d.x - a.x;
		double ady = d.y - a.y;

		// Orient(a, b, c) and orient(a, b, d) share the ab vector
		double o1 = abx * acy - aby * acx;
		double o2 = abx * ady - aby * adx;

		// Early exit if same side
		if (o1 * o2 > 0)
			return false;

		double cdx = d.x - c.x;
		double cdy = d.y - c.y;
		double cax = a.x - c.x;
		double cay = a.y - c.y;
		double cbx = b.x - c.x;
		double cby = b.y - c.y;

		double o3 = cdx * cay - cdy * cax;
		double o4 = cdx * cby - cdy * cbx;

		// Check opposite sides
		return o3 * o4 < 0;
	}

	private static boolean segmentsIntersect2(PVector a, PVector b, PVector c, PVector d) {
		double o1 = orient(a, b, c);
		double o2 = orient(a, b, d);
		double o3 = orient(c, d, a);
		double o4 = orient(c, d, b);

		// proper intersection
		if (o1 * o2 < 0 && o3 * o4 < 0) {
			return true;
		}
//		final double EPS = 1e-9;
		// handle collinear / endpoint cases
//		if (Math.abs(o1) < EPS && onSegment(a, b, c))
//			return true;
//		if (Math.abs(o2) < EPS && onSegment(a, b, d))
//			return true;
//		if (Math.abs(o3) < EPS && onSegment(c, d, a))
//			return true;
//		if (Math.abs(o4) < EPS && onSegment(c, d, b))
//			return true;

		return false;
	}

	private static double orient(PVector a, PVector b, PVector c) {
		return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	}

	private static boolean onSegment(PVector a, PVector b, PVector p) {
		final double EPS = 1e-9;
		return p.x >= Math.min(a.x, b.x) - EPS && p.x <= Math.max(a.x, b.x) + EPS && p.y >= Math.min(a.y, b.y) - EPS && p.y <= Math.max(a.y, b.y) + EPS;
	}

	private static PShape toPolygon(List<PVector> points) {
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			points.add(points.get(0)); // close
		}
		return PGS_Conversion.fromPVector(points);
	}

}
