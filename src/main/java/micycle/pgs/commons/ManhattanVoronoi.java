package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;

/**
 * Computes Manhattan (L1) Voronoi diagrams clipped to an axis-aligned
 * rectangular bounding box.
 * <p>
 * This implementation is a Java port of the original JavaScript reference
 * implementation ({@code JDragovich/manhattan-voronoi}) and is intended to
 * preserve its behavior and structure for correctness and parity.
 * <p>
 * The underlying approach follows the divide-and-conquer scheme described by:
 * <blockquote> Lee, Der-Tsai, and C. K. Wong. <i>"Voronoui Diagrams in L_1(L_∞)
 * Metrics with 2-Dimensional Storage Applications."</i> SIAM Journal on
 * Computing 9, no. 1 (1980): 200–211. </blockquote>
 *
 * <h2>High-level algorithm (graph construction)</h2>
 * <ol>
 * <li><b>Preprocess</b> (optional): nudge input sites slightly to avoid the
 * degenerate "square bisector" case (|dx| = |dy|).</li>
 * <li><b>Sort</b> sites by x-coordinate (breaking ties by y).</li>
 * <li><b>Recursive split</b>: recursively divide the sorted sites into left and
 * right subsets; for base cases, compute and attach the L1 bisector for a pair
 * of sites.</li>
 * <li><b>Merge step</b>:
 * <ol>
 * <li>Choose an initial cross-subset bisector between a candidate site from the
 * left subset and its nearest neighbor in the right subset.</li>
 * <li><b>Walk the merge chain</b> upward and downward: repeatedly intersect the
 * current merge bisector with existing bisectors, trim bisectors at
 * intersection points, and "hop" to adjacent sites to continue the merge.</li>
 * <li><b>Orphan handling</b>: detect and remove bisectors that become trapped
 * by the newly formed merge boundary; add replacement merge bisectors when
 * needed.</li>
 * <li>Attach all merge bisectors to their incident sites and remove bisectors
 * invalidated by trapping.</li>
 * </ol>
 * </li>
 * <li><b>Post-process per site</b>: derive an ordered list of boundary vertices
 * for the clipped cell (and neighbor sites) from the site’s incident
 * bisectors.</li>
 * </ol>
 *
 * <p>
 * <b>Note:</b> This class produces a topological graph (sites with incident
 * bisectors and neighbor relations) and a set of boundary vertices for each
 * cell.
 * 
 * @author Original javascript implementation by Joe Dragovich
 * @author Java port by Michael Carleton
 */
public final class ManhattanVoronoi {

	static final LineIntersector li = new RobustLineIntersector();

	private final double minX;
	private final double minY;
	private final double maxX;
	private final double maxY;

	private ManhattanVoronoi(Envelope bounds) {
		if (bounds == null) {
			throw new IllegalArgumentException("Bounds must not be null.");
		}
		this.minX = bounds.getMinX();
		this.minY = bounds.getMinY();
		this.maxX = bounds.getMaxX();
		this.maxY = bounds.getMaxY();
	}

	public static final class Site {
		public final Coordinate site; // the actual point
		public List<Bisector> bisectors = new ArrayList<>();

		// outputs filled by generate
		public List<Coordinate> polygonPoints = new ArrayList<>();
		public List<Site> neighbors = new ArrayList<>();

		public Site(Coordinate site) {
			this.site = site;
		}

		/**
		 * Builds a JTS Polygon from this site's polygonPoints.
		 *
		 * @param gf geometry factory
		 */
		public Polygon toPolygon(GeometryFactory gf) {
			if (polygonPoints == null || polygonPoints.isEmpty()) {
				return gf.createPolygon(); // empty
			}

			List<Coordinate> pts = new ArrayList<>(polygonPoints.size());
			for (int i = 0; i < polygonPoints.size(); i++) {
				Coordinate c = polygonPoints.get(i);
				pts.add(c);
			}

			// Need at least 3 vertices
			if (pts.size() < 3) {
				return gf.createPolygon();
			}

			if (pts.size() < 3) {
				return gf.createPolygon();
			}

			// Close ring
			pts.add(new Coordinate(pts.get(0)));

			return gf.createPolygon(pts.toArray(Coordinate[]::new));
		}

		@Override
		public String toString() {
			return "Site{" + "site=" + fmt(site) + '}';
		}
	}

	public static final class Bisector {
		public final Site[] sites = new Site[2];
		public boolean up; // same meaning as JS
		public List<Coordinate> points = new ArrayList<>();
		public List<Coordinate> intersections = new ArrayList<>();
		public boolean compound = false;
		public int mergeLine = 0;

		public Bisector(Site a, Site b) {
			sites[0] = a;
			sites[1] = b;
		}

		public boolean sitesContainBoth(Site a, Site b) {
			return (sites[0] == a || sites[1] == a) && (sites[0] == b || sites[1] == b);
		}

		@Override
		public String toString() {
			return "Bisector{" + "sites=" + sites[0] + "," + sites[1] + ", up=" + up + ", points="
					+ points.stream().map(ManhattanVoronoi::fmt).collect(Collectors.toList()) + '}';
		}
	}

	public static List<Site> generate(Collection<Coordinate> sitePoints, double width, double height) {
		return generate(sitePoints, new Envelope(0, width, 0, height));
	}

	public static List<Site> generate(Collection<Coordinate> sitePoints, double width, double height, boolean nudgeData) {
		return generate(sitePoints, new Envelope(0, width, 0, height), nudgeData);
	}

	public static List<Site> generate(Collection<Coordinate> sitePoints, Envelope bounds) {
		return generate(sitePoints, bounds, true);
	}

	public static List<Site> generate(Collection<Coordinate> sitePoints, Envelope bounds, boolean nudgeData) {
		return new ManhattanVoronoi(bounds).generate(sitePoints, nudgeData);
	}

	private List<Site> generate(Collection<Coordinate> sitePoints, boolean nudgeData) {
		final int n = sitePoints.size();

		List<Coordinate> points = sitePoints.stream().map(Coordinate::copy).collect(Collectors.toList());

		if (nudgeData) {
			cleanData(points);
		}

		// Sort points by x then y
		points.sort((a, b) -> {
			int cx = Double.compare(a.x, b.x);
			return (cx != 0) ? cx : Double.compare(a.y, b.y);
		});

		// Create sites into an array (faster for in-place range recursion)
		Site[] sites = new Site[n];
		for (int i = 0; i < n; i++) {
			sites[i] = new Site(points.get(i));
		}

		// Build graph in-place (no BiFunction, no list slicing/copies)
		recursiveSplit(sites, 0, n);

		List<Site> graph = new ArrayList<>(n);
		Collections.addAll(graph, sites);

		postProcessSites(graph);
		return graph;
	}

	private void postProcessSites(List<Site> graph) {
		// Pre-create corners once
		final Coordinate[] corners = new Coordinate[] { new Coordinate(minX, minY), new Coordinate(maxX, minY), new Coordinate(maxX, maxY),
				new Coordinate(minX, maxY) };

		graph.parallelStream().forEach(site -> {
			buildPolygonPointsByChaining(site);
			injectCornersIfNeeded(site, corners);

			if (!site.polygonPoints.isEmpty()) {
				// Sort around site
				site.polygonPoints.sort((p1, p2) -> Double.compare(angle(site.site, p1), angle(site.site, p2)));
			}
//			computeNeighbors(site); // NOTE
		});
	}

	private void recursiveSplit(Site[] sites, int from, int to) {
		final int size = to - from;

		if (size > 2) {
			final int half = (size - (size % 2)) / 2;
			final int splitPoint = from + half;

			// recurse
			recursiveSplit(sites, from, splitPoint);
			recursiveSplit(sites, splitPoint, to);

			// working sites
			final Site lLast = sites[splitPoint - 1];

			// Find nearest neighbor in the right half WITHOUT sorting the entire R
			Site nearest = sites[splitPoint];
			double best = distance(lLast.site, nearest.site);
			for (int i = splitPoint + 1; i < to; i++) {
				Site s = sites[i];
				double d = distance(lLast.site, s.site);
				if (d < best) {
					best = d;
					nearest = s;
				}
			}

			StartingInfo startingInfo = determineStartingBisector(lLast, nearest, null);

			Bisector initialBisector = startingInfo.startingBisector;
			Site initialR = startingInfo.nearestNeighbor;
			Site initialL = startingInfo.w;

			// Single merge list; walkMergeLine appends as needed
			List<Bisector> mergeArray = new ArrayList<>();
			mergeArray.add(initialBisector);

			walkMergeLine(initialR, initialL, initialBisector, new Coordinate(maxX, maxY), true, null, mergeArray);
			walkMergeLine(initialR, initialL, initialBisector, new Coordinate(minX, minY), false, null, mergeArray);

			// attach merge bisectors
			for (int i = 0; i < mergeArray.size(); i++) {
				Bisector bisector = mergeArray.get(i);

				bisector.mergeLine = size;

				bisector.sites[0].bisectors = clearOutOrphans(bisector.sites[0], bisector.sites[1]);
				bisector.sites[1].bisectors = clearOutOrphans(bisector.sites[1], bisector.sites[0]);

				bisector.sites[0].bisectors.add(bisector);
				bisector.sites[1].bisectors.add(bisector);
			}

		} else if (size == 2) {
			Bisector bisector = findL1Bisector(sites[from], sites[from + 1]);
			sites[from].bisectors.add(bisector);
			sites[from + 1].bisectors.add(bisector);
		} else {
			// size 0/1: nothing to do
		}
	}

	private void walkMergeLine(Site currentR, Site currentL, Bisector currentBisector, Coordinate currentCropPoint, boolean goUp, Bisector crossedBorder,
			List<Bisector> mergeArray) {
		while (true) {

			// ensure bisector matches current sites; if not, create and trim
			if (!currentBisector.sitesContainBoth(currentR, currentL)) {
				currentBisector = findL1Bisector(currentR, currentL);
				trimBisector(currentBisector, crossedBorder, currentCropPoint);
				mergeArray.add(currentBisector);
			}

			List<CropCandidate> cropLArray = buildCropCandidatesForSide(currentL, currentBisector, currentR, currentCropPoint, goUp, crossedBorder, true);

			List<CropCandidate> cropRArray = buildCropCandidatesForSide(currentR, currentBisector, currentL, currentCropPoint, goUp, crossedBorder, false);

			CropCandidate cropL = (!cropLArray.isEmpty() && cropLArray.get(0).bisector != currentBisector) ? cropLArray.get(0)
					: new CropCandidate(null, goUp ? new Coordinate(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)
							: new Coordinate(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY));

			CropCandidate cropR = (!cropRArray.isEmpty() && cropRArray.get(0).bisector != currentBisector) ? cropRArray.get(0)
					: new CropCandidate(null, goUp ? new Coordinate(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)
							: new Coordinate(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY));

			// done?
			if (cropL.bisector == null && cropR.bisector == null) {

				Bisector leftOrphan = checkForOphans(currentR, currentL, goUp);
				Bisector rightOrphan = checkForOphans(currentL, currentR, goUp);

				if (leftOrphan != null) {
					// Remove trapped bisector
					for (Site s : leftOrphan.sites) {
						s.bisectors.removeIf(b -> b == leftOrphan);
					}

					Site hopTo = findHopTo(leftOrphan, currentL);
					currentR = findCorrectW(currentR, hopTo);

					Bisector newMergeBisector = findL1Bisector(hopTo, currentR);
					mergeArray.add(newMergeBisector);

					// continue with updated state
					currentL = hopTo;
					currentBisector = newMergeBisector;
					continue;

				} else if (rightOrphan != null) {
					for (Site s : rightOrphan.sites) {
						s.bisectors.removeIf(b -> b == rightOrphan);
					}

					Site hopTo = findHopTo(rightOrphan, currentR);
					currentL = findCorrectW(currentL, hopTo);

					Bisector newMergeBisector = findL1Bisector(hopTo, currentL);
					mergeArray.add(newMergeBisector);

					// continue with updated state
					currentR = hopTo;
					currentBisector = newMergeBisector;
					continue;
				}

				return; // finished
			}

			FirstBorderCross first = determineFirstBorderCross(cropR, cropL, currentCropPoint);

			if (first == FirstBorderCross.RIGHT) {
				trimBisector(cropR.bisector, currentBisector, cropR.point);
				trimBisector(currentBisector, cropR.bisector, cropR.point);
				currentBisector.intersections.add(cropR.point);

				crossedBorder = cropR.bisector;
				currentR = findOtherSite(cropR.bisector, currentR);
				currentCropPoint = cropR.point;

			} else if (first == FirstBorderCross.LEFT) {
				trimBisector(cropL.bisector, currentBisector, cropL.point);
				trimBisector(currentBisector, cropL.bisector, cropL.point);
				currentBisector.intersections.add(cropL.point);

				crossedBorder = cropL.bisector;
				currentL = findOtherSite(cropL.bisector, currentL);
				currentCropPoint = cropL.point;

			} else { // BOTH
				trimBisector(cropR.bisector, currentBisector, cropR.point);
				trimBisector(currentBisector, cropR.bisector, cropR.point);
				currentBisector.intersections.add(cropR.point);

				crossedBorder = cropR.bisector;
				currentR = findOtherSite(cropR.bisector, currentR);
				currentCropPoint = cropR.point;

				trimBisector(cropL.bisector, currentBisector, cropL.point);
				trimBisector(currentBisector, cropL.bisector, cropL.point);
				currentBisector.intersections.add(cropL.point);

				crossedBorder = cropL.bisector;
				currentL = findOtherSite(cropL.bisector, currentL);
				currentCropPoint = cropL.point;
			}
		}
	}

	private StartingInfo determineStartingBisector(Site w, Site nearestNeighbor, Coordinate lastIntersect) {
		while (true) {
			if (lastIntersect == null) {
				lastIntersect = w.site;
			}

			// horizontal ray to the right boundary
			final Coordinate z = new Coordinate(maxX, w.site.y);

			IntersectionHit hit = null;
			for (Bisector b : nearestNeighbor.bisectors) {
				Coordinate p = segmentBisectorIntersection(w.site, z, b);
				if (p != null) {
					hit = new IntersectionHit(p, b);
					break;
				}
			}

			if (hit != null && distance(w.site, hit.point) > distance(nearestNeighbor.site, hit.point)) {
				Bisector startingBisector = findL1Bisector(w, nearestNeighbor);
				return new StartingInfo(startingBisector, w, nearestNeighbor, hit.point);

			} else if (hit != null && distance(w.site, hit.point) < distance(nearestNeighbor.site, hit.point) && hit.point.x > lastIntersect.x) {

				nearestNeighbor = findOtherSite(hit.bisector, nearestNeighbor);
				lastIntersect = hit.point;
				continue;

			} else {
				w = findCorrectW(w, nearestNeighbor);
				Bisector startingBisector = findL1Bisector(w, nearestNeighbor);
				return new StartingInfo(startingBisector, w, nearestNeighbor, hit != null ? hit.point : w.site);
			}
		}
	}

	private static Coordinate segmentBisectorIntersection(Coordinate s0, Coordinate s1, Bisector b) {
		List<Coordinate> pts = b.points;
		for (int i = 0; i < pts.size() - 1; i++) {
			Coordinate p0 = pts.get(i);
			Coordinate p1 = pts.get(i + 1);
			Coordinate ip = segmentIntersection(s0, s1, p0, p1);
			if (ip != null) {
				return ip;
			}
		}
		return null;
	}

	private Site findCorrectW(Site w, Site nearestNeighbor) {
		while (true) {
			Bisector startingBisector = findL1Bisector(w, nearestNeighbor);

			Site bestHop = null;
			double bestDist = Double.POSITIVE_INFINITY;

			// find closest hopTo that traps the starting bisector
			List<Bisector> wb = w.bisectors;
			for (int i = 0; i < wb.size(); i++) {
				Bisector b = wb.get(i);
				Site hopTo = findHopTo(b, w);

				if (isBisectorTrapped(hopTo, startingBisector)) {
					double d = distance(hopTo.site, nearestNeighbor.site);
					if (d < bestDist) {
						bestDist = d;
						bestHop = hopTo;
					}
				}
			}

			if (bestHop == null) {
				return w;
			}

			w = bestHop;
		}
	}

	private Bisector checkForOphans(Site trapper, Site trapped, boolean goUp) {
		Bisector bestBisector = null;
		double bestExtreme = goUp ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;

		List<Bisector> tb = trapped.bisectors;
		for (int i = 0; i < tb.size(); i++) {
			Bisector b = tb.get(i);

			Site hopTo = findHopTo(b, trapped);
			boolean directionOk = (goUp == (hopTo.site.y < trapped.site.y));
			if (!directionOk) {
				continue;
			}

			if (!isBisectorTrapped(trapper, b)) {
				continue;
			}

			Bisector mergeLine = findL1Bisector(hopTo, trapper);
			double extreme = getExtremePoint(mergeLine, goUp);

			if (goUp) {
				if (extreme > bestExtreme) {
					bestExtreme = extreme;
					bestBisector = b;
				}
			} else {
				if (extreme < bestExtreme) {
					bestExtreme = extreme;
					bestBisector = b;
				}
			}
		}

		return bestBisector;
	}

	private List<CropCandidate> buildCropCandidatesForSide(Site sideSite, Bisector currentBisector, Site otherMergeSite, Coordinate currentCropPoint,
			boolean goUp, Bisector crossedBorder, boolean isLeftSide) {
		// Build candidate list + precompute hopTo for each candidate
		List<CropCandidate> candidates = new ArrayList<>();
		List<Site> hopTos = new ArrayList<>();

		List<Bisector> bis = sideSite.bisectors;
		for (int i = 0; i < bis.size(); i++) {
			Bisector b = bis.get(i);

			Coordinate p = bisectorIntersection(currentBisector, b);
			if (p == null) {
				continue;
			}

			Site hopTo = findHopTo(b, sideSite);
			boolean upward = isNewBisectorUpward(hopTo, sideSite, otherMergeSite, goUp);

			boolean sameCropAndSameBorder = samePoint(p, currentCropPoint) && b == crossedBorder;
			if ((goUp == upward) && !sameCropAndSameBorder) {
				candidates.add(new CropCandidate(b, p));
				hopTos.add(hopTo);
			}
		}

		// Sort stage (matches original JS directionality)
		candidates.sort((a, b) -> {
			if (isLeftSide) {
				Site hopToB = findHopTo(b.bisector, sideSite);
				Site hopToA = findHopTo(a.bisector, sideSite);
				return Double.compare(angle(sideSite.site, hopToB.site), angle(sideSite.site, hopToA.site));
			} else {
				Site hopToA = findHopTo(a.bisector, sideSite);
				Site hopToB = findHopTo(b.bisector, sideSite);
				return Double.compare(angle(sideSite.site, hopToA.site), angle(sideSite.site, hopToB.site));
			}
		});

		// Rebuild hopTos in the same order as candidates (since we sorted candidates)
		hopTos.clear();
		for (int i = 0; i < candidates.size(); i++) {
			hopTos.add(findHopTo(candidates.get(i).bisector, sideSite));
		}

		// JS-like "every(...)" filter (still O(n^2), but fewer repeated calls)
		List<CropCandidate> filtered = new ArrayList<>();
		for (int i = 0; i < candidates.size(); i++) {
			CropCandidate e = candidates.get(i);
			Site hopTo = hopTos.get(i);

			Bisector newMergeLine = findL1Bisector(otherMergeSite, hopTo);
			trimBisector(newMergeLine, e.bisector, e.point);

			boolean ok = true;
			for (int j = 0; j < candidates.size(); j++) {
				Site hopToD = hopTos.get(j);
				if (hopToD == hopTo) {
					continue;
				}
				if (isBisectorTrapped(hopToD, newMergeLine)) {
					ok = false;
					break;
				}
			}

			if (ok) {
				filtered.add(e);
			}
		}

		return filtered;
	}

	private void buildPolygonPointsByChaining(Site site) {

		if (site.bisectors.isEmpty()) {
			site.polygonPoints = new ArrayList<>();
			return;
		}

		Bisector startBisector = findStartBisectorOnEdge(site);
		if (startBisector == null) {
			startBisector = site.bisectors.get(0);
		}

		// used set: identity semantics
		Set<Bisector> used = Collections.newSetFromMap(new IdentityHashMap<>());
		used.add(startBisector);

		List<Coordinate> polygon = new ArrayList<>(startBisector.points.size());
		polygon.addAll(startBisector.points);

		// reverse if last point is on edge (matches JS)
		if (!polygon.isEmpty() && isPointOnEdge(polygon.get(polygon.size() - 1))) {
			Collections.reverse(polygon);
		}

		// chain remaining bisectors
		while (used.size() < site.bisectors.size()) {
			Coordinate last = polygon.get(polygon.size() - 1);

			Bisector next = null;
			double bestDist = Double.POSITIVE_INFINITY;

			for (Bisector candidate : site.bisectors) {
				if (used.contains(candidate)) {
					continue;
				}

				Coordinate a = candidate.points.get(0);
				Coordinate b = candidate.points.get(candidate.points.size() - 1);
				double candDist = Math.min(distance(last, a), distance(last, b));

				if (candDist < bestDist) {
					bestDist = candDist;
					next = candidate;
				}
			}

			if (next == null) {
				break; // JS assumes exists; keep safe
			}
			used.add(next);

			// append points, reversing if needed
			List<Coordinate> nextPts = new ArrayList<>(next.points);
			if (!nextPts.isEmpty() && samePoint(nextPts.get(nextPts.size() - 1), last)) {
				Collections.reverse(nextPts);
			}
			polygon.addAll(nextPts);
		}

		site.polygonPoints = polygon;
	}

	private Bisector findStartBisectorOnEdge(Site site) {
		for (Bisector b : site.bisectors) {
			for (Coordinate element : b.points) {
				if (isPointOnEdge(element)) {
					return b;
				}
			}
		}
		return null;
	}

	private void injectCornersIfNeeded(Site site, Coordinate[] corners) {
		if (site.polygonPoints.isEmpty()) {
			return;
		}

		Coordinate first = site.polygonPoints.get(0);
		Coordinate last = site.polygonPoints.get(site.polygonPoints.size() - 1);

		if (!(isPointOnEdge(first) && isPointOnEdge(last) && !arePointsOnSameEdge(first, last))) {
			return;
		}

		List<Coordinate> filteredCorners = new ArrayList<>(4);

		for (Coordinate corner : corners) {
			boolean ok = true;
			for (Bisector b : site.bisectors) {
				// you already have this helper
				if (segmentIntersectsBisector(corner, site.site, b)) {
					ok = false;
					break;
				}
			}

			if (ok) {
				filteredCorners.add(new Coordinate(corner));
			}
		}

		if (!filteredCorners.isEmpty()) {
			site.polygonPoints.addAll(filteredCorners);
		}
	}

	private static void computeNeighbors(Site site) {
		if (site.bisectors.isEmpty()) {
			site.neighbors = new ArrayList<>();
			return;
		}

		List<Site> neigh = new ArrayList<>(site.bisectors.size());
		for (int bi = 0; bi < site.bisectors.size(); bi++) {
			neigh.add(findHopTo(site.bisectors.get(bi), site));
		}
		site.neighbors = neigh;
	}

	private static boolean segmentIntersectsBisector(Coordinate s0, Coordinate s1, Bisector b) {
		// Intersect the segment (s0->s1) with every segment of the bisector polyline
		List<Coordinate> pts = b.points;
		for (int i = 0; i < pts.size() - 1; i++) {
			Coordinate p0 = pts.get(i);
			Coordinate p1 = pts.get(i + 1);
			if (segmentIntersection(s0, s1, p0, p1) != null) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Nudge points to hopefully eliminate square bisectors. Mutates the given
	 * coordinates, matching JS behavior.
	 */
	public static List<Coordinate> cleanData(List<Coordinate> data) {
		final double epsOffset = 1e-10;
		for (int i = 0; i < data.size(); i++) {
			Coordinate e = data.get(i);
			for (int j = 0; j < data.size(); j++) {
				Coordinate d = data.get(j);
				if (i != j && Math.abs(d.x - e.x) == Math.abs(d.y - e.y)) {
					d.x += epsOffset;
					d.y -= epsOffset;
				}
			}
		}
		return data;
	}

	private static Site findOtherSite(Bisector bisector, Site current) {
		return bisector.sites[0] == current ? bisector.sites[1] : bisector.sites[0];
	}

	private static double angle(Coordinate p1, Coordinate p2) {
		double ang = FastAtan2.atan2(p2.y - p1.y, p2.x - p1.x);
		if (ang < 0) {
			ang = Math.PI + Math.PI + ang;
		}
		return ang;
	}

	private enum FirstBorderCross {
		RIGHT, LEFT, BOTH
	}

	private static FirstBorderCross determineFirstBorderCross(CropCandidate cropR, CropCandidate cropL, Coordinate currentCropPoint) {
		double dr = Math.abs(cropR.point.y - currentCropPoint.y);
		double dl = Math.abs(cropL.point.y - currentCropPoint.y);

		if (dr == dl) {
			return FirstBorderCross.BOTH;
		}
		return (dr < dl) ? FirstBorderCross.RIGHT : FirstBorderCross.LEFT;
	}

	private record StartingInfo(Bisector startingBisector, Site w, Site nearestNeighbor, Coordinate startingIntersection) {
	}

	private record IntersectionHit(Coordinate point, Bisector bisector) {
	}

	/**
	 * Hyperoptimized findL1Bisector - eliminates allocations and redundant
	 * operations.
	 */
	private Bisector findL1Bisector(Site P1, Site P2) {
		final double p1x = P1.site.x;
		final double p1y = P1.site.y;
		final double p2x = P2.site.x;
		final double p2y = P2.site.y;

		final double xDistance = p1x - p2x;
		final double yDistance = p1y - p2y;

		if (p1x == p2x && p1y == p2y) {
			throw new IllegalArgumentException("Duplicate point: Points " + P1 + " and " + P2 + " are duplicates.");
		}

		final double absX = Math.abs(xDistance);
		final double absY = Math.abs(yDistance);

		// Square bisector check
		if (absX == absY) {
			throw new IllegalArgumentException("Square bisector: Points " + P1 + " and " + P2
					+ " are points on a square (their vertical distance equals horizontal distance). Consider nudgeData.");
		}

		Bisector bisector = new Bisector(P1, P2);

		// Fast path: vertical line (xDistance == 0)
		if (absX == 0) {
			final double midY = (p1y + p2y) * 0.5;
			bisector.up = false;
			bisector.points = List.of(new Coordinate(minX, midY), new Coordinate(maxX, midY));
			return bisector;
		}

		// Fast path: horizontal line (yDistance == 0)
		if (absY == 0) {
			final double midX = (p1x + p2x) * 0.5;
			bisector.up = true;
			bisector.points = List.of(new Coordinate(midX, minY), new Coordinate(midX, maxY));
			return bisector;
		}

		// Pre-compute midpoint (used in both branches)
		final double midX = (p1x + p2x) * 0.5;
		final double midY = (p1y + p2y) * 0.5;

		// Determine slope
		final double slope = (yDistance * xDistance > 0) ? -1.0 : 1.0;
		final double intercept = midY - midX * slope;

		// Branch based on dominant direction
		if (absX >= absY) {
			// Horizontal-dominant (up = true)
			bisector.up = true;

			// Compute vertices directly in sorted order
			final double v1x = (p1y - intercept) * slope; // slope is ±1, so division becomes multiplication
			final double v2x = (p2y - intercept) * slope;

			// Determine order based on y-coordinates
			if (p1y < p2y) {
				// p1y is lower, so v1 comes first
				bisector.points = List.of(new Coordinate(v1x, minY), new Coordinate(v1x, p1y), new Coordinate(v2x, p2y), new Coordinate(v2x, maxY));
			} else {
				// p2y is lower, so v2 comes first
				bisector.points = List.of(new Coordinate(v2x, minY), new Coordinate(v2x, p2y), new Coordinate(v1x, p1y), new Coordinate(v1x, maxY));
			}

		} else {
			// Vertical-dominant (up = false)
			bisector.up = false;

			// Compute vertices directly in sorted order
			final double v1y = p1x * slope + intercept;
			final double v2y = p2x * slope + intercept;

			// Determine order based on x-coordinates
			if (p1x < p2x) {
				// p1x is leftmost, so v1 comes first
				bisector.points = List.of(new Coordinate(minX, v1y), new Coordinate(p1x, v1y), new Coordinate(p2x, v2y), new Coordinate(maxX, v2y));
			} else {
				// p2x is leftmost, so v2 comes first
				bisector.points = List.of(new Coordinate(minX, v2y), new Coordinate(p2x, v2y), new Coordinate(p1x, v1y), new Coordinate(maxX, v1y));
			}
		}

		return bisector;
	}

	private static List<Bisector> clearOutOrphans(Site orphanage, Site trapPoint) {
		orphanage.bisectors.removeIf(b -> isBisectorTrapped(trapPoint, b));
		return orphanage.bisectors;
	}

	private static Site findHopTo(Bisector bisector, Site hopFrom) {
		return bisector.sites[0] == hopFrom ? bisector.sites[1] : bisector.sites[0];
	}

	/** Manhattan (L1) distance. */
	public static double distance(final Coordinate p1, final Coordinate p2) {
		return Math.abs(p1.x - p2.x) + Math.abs(p1.y - p2.y);
	}

	private static boolean isBisectorTrapped(Site trapPoint, Bisector bisector) {
		for (Coordinate point : bisector.points) {
			double dTrap = distance(trapPoint.site, point);
			if (!(dTrap <= distance(bisector.sites[0].site, point) && dTrap <= distance(bisector.sites[1].site, point))) {
				return false;
			}
		}
		return true;
	}

	private static double getExtremePoint(Bisector bisector, boolean goUp) {
		double acc = goUp ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
		for (Coordinate e : bisector.points) {
			acc = goUp ? Math.max(e.y, acc) : Math.min(e.y, acc);
		}
		return acc;
	}

	private static void trimBisector(Bisector target, Bisector intersector, Coordinate intersection) {
		if (intersector == null) {
			return;
		}

		// Find the "polygon site" (the intersector site not in target.sites)
		Site polygonSite = null;
		Site t0 = target.sites[0];
		Site t1 = target.sites[1];

		Site i0 = intersector.sites[0];
		Site i1 = intersector.sites[1];

		if (i0 != t0 && i0 != t1) {
			polygonSite = i0;
		} else if (i1 != t0 && i1 != t1) {
			polygonSite = i1;
		} else {
			return;
		}

		final Coordinate poly = polygonSite.site;

		List<Coordinate> src = target.points;
		List<Coordinate> newPoints = new ArrayList<>(src.size() + 1);

		final Coordinate s0 = t0.site;
		final Coordinate s1 = t1.site;

		for (int i = 0; i < src.size(); i++) {
			Coordinate p = src.get(i);

			// keep point if it's closer to both target sites than to polygonSite
			if (distance(p, s0) < distance(p, poly) && distance(p, s1) < distance(p, poly)) {
				newPoints.add(p);
			}
		}

		newPoints.add(new Coordinate(intersection.x, intersection.y));

		if (target.up) {
			newPoints.sort(Comparator.comparingDouble(c -> c.y));
		} else {
			newPoints.sort(Comparator.comparingDouble(c -> c.x));
		}

		target.points = newPoints;
	}

	private record CropCandidate(Bisector bisector, Coordinate point) {
	}

	private static boolean isNewBisectorUpward(Site hopTo, Site hopFrom, Site site, boolean goUpUnused) {
		double slope = (hopTo.site.y - site.site.y) / (hopTo.site.x - site.site.x);
		double intercept = hopTo.site.y - (slope * hopTo.site.x);

		if (Double.isInfinite(slope)) {
			return site.site.y > hopTo.site.y;
		}

		return hopFrom.site.y > (slope * hopFrom.site.x) + intercept;
	}

	private static Coordinate bisectorIntersection(Bisector b1, Bisector b2) {
		if (b1 == b2) {
			return null;
		}

		final List<Coordinate> pts1 = b1.points;
		final List<Coordinate> pts2 = b2.points;
		final int size1 = pts1.size();
		final int size2 = pts2.size();

		// Early exit for empty bisectors
		if (size1 < 2 || size2 < 2) {
			return null;
		}

		// Bounding box check
		final Coordinate p1_0 = pts1.get(0);
		final Coordinate p1_last = pts1.get(size1 - 1);
		final Coordinate p2_0 = pts2.get(0);
		final Coordinate p2_last = pts2.get(size2 - 1);

		final double b1_minX = Math.min(p1_0.x, p1_last.x);
		final double b1_maxX = Math.max(p1_0.x, p1_last.x);
		final double b1_minY = Math.min(p1_0.y, p1_last.y);
		final double b1_maxY = Math.max(p1_0.y, p1_last.y);

		final double b2_minX = Math.min(p2_0.x, p2_last.x);
		final double b2_maxX = Math.max(p2_0.x, p2_last.x);
		final double b2_minY = Math.min(p2_0.y, p2_last.y);
		final double b2_maxY = Math.max(p2_0.y, p2_last.y);

		// No overlap = no intersection
		if (b1_maxX < b2_minX || b2_maxX < b1_minX || b1_maxY < b2_minY || b2_maxY < b1_minY) {
			return null;
		}

		// Cache segment count (avoid repeated size() - 1)
		final int seg1Count = size1 - 1;
		final int seg2Count = size2 - 1;

		// Cache-friendly iteration with segment bounding box tests
		for (int i = 0; i < seg1Count; i++) {
			final Coordinate a0 = pts1.get(i);
			final Coordinate a1 = pts1.get(i + 1);

			// Segment 1 bounds (extract to locals for cache efficiency)
			final double a_minX = Math.min(a0.x, a1.x);
			final double a_maxX = Math.max(a0.x, a1.x);
			final double a_minY = Math.min(a0.y, a1.y);
			final double a_maxY = Math.max(a0.y, a1.y);

			for (int j = 0; j < seg2Count; j++) {
				final Coordinate b0 = pts2.get(j);
				final Coordinate bb1 = pts2.get(j + 1);

				// Cheap segment bounding box test (eliminates 90% of intersection checks)
				final double b_minX = Math.min(b0.x, bb1.x);
				final double b_maxX = Math.max(b0.x, bb1.x);

				if (a_maxX < b_minX || b_maxX < a_minX) {
					continue;
				}

				final double b_minY = Math.min(b0.y, bb1.y);
				final double b_maxY = Math.max(b0.y, bb1.y);

				if (a_maxY < b_minY || b_maxY < a_minY) {
					continue;
				}

				// Bounding boxes overlap - now do the expensive intersection test
				Coordinate intersect = segmentIntersection(a0, a1, b0, bb1);
				if (intersect != null) {
					return intersect; // Early termination on first intersection
				}
			}
		}

		return null;
	}

	/**
	 * Segment intersection ported from JS segementIntersection. - denom == 0 =>
	 * null (JS returned null for parallel/collinear) - if no intersection within
	 * segment bounds => null (JS returned false)
	 */
	private static Coordinate segmentIntersection(final Coordinate l10, final Coordinate l11, final Coordinate l20, final Coordinate l21) {
		final double denom = (l21.y - l20.y) * (l11.x - l10.x) - (l21.x - l20.x) * (l11.y - l10.y);
		if (denom == 0) {
			return null;
		}

		final double ua = ((l21.x - l20.x) * (l10.y - l20.y) - (l21.y - l20.y) * (l10.x - l20.x)) / denom;
		final double ub = ((l11.x - l10.x) * (l10.y - l20.y) - (l11.y - l10.y) * (l10.x - l20.x)) / denom;

		if (!(ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1)) {
			return null;
		}

		return new Coordinate(l10.x + ua * (l11.x - l10.x), l10.y + ua * (l11.y - l10.y));
	}

	private static Coordinate segmentIntersectionRobust(Coordinate a0, Coordinate a1, Coordinate b0, Coordinate b1) {
		li.computeIntersection(a0, a1, b0, b1);

		if (!li.hasIntersection()) {
			return null;
		}
		if (li.getIntersectionNum() != 1) {
			return null;
		}

		Coordinate ip = li.getIntersection(0);
		return new Coordinate(ip.x, ip.y);
	}

	private static boolean samePoint(final Coordinate p1, final Coordinate p2) {
		return p1.x == p2.x && p1.y == p2.y;
	}

	private boolean isPointOnEdge(final Coordinate p) {
		return p.x == minX || p.x == maxX || p.y == minY || p.y == maxY;
	}

	private boolean arePointsOnSameEdge(final Coordinate p1, final Coordinate p2) {
		return (p1.x == p2.x && (p1.x == minX || p1.x == maxX)) || (p1.y == p2.y && (p1.y == minY || p1.y == maxY));
	}

	private static String fmt(final Coordinate c) {
		return "[" + c.x + "," + c.y + "]";
	}
}