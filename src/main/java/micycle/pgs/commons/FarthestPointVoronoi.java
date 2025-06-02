package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Triangle;
import org.locationtech.jts.operation.polygonize.Polygonizer;

import micycle.pgs.commons.FarthestPointVoronoi.DCELVertex;

/**
 * Farthest-Point Voronoi Diagram
 * <p>
 * The region for a site ( p ) in a farthest-point Voronoi diagram (FPVD) is the
 * set of all points in the plane for which ( p ) is the farthest site among all
 * sites. In contrast, in a regular Voronoi diagram, regions hug their site—they
 * surround and contain the generating point, making it visually obvious which
 * region belongs to which site.
 * <p>
 * In the FPVD, regions do not hug their site; the generator site is far away
 * and not contained in its region. Only vertices of the convex hull have
 * regions in the FPVD, as only they can be farthest from some location in the
 * plane.
 * <p>
 * A useful interpretation is that all points within a given FPVD region
 * mutually share the same farthest site—the region’s generator. However, the
 * site for a region is not visually apparent, since it is not located within or
 * necessarily even near its region.
 *
 * @author Michael Carleton
 */
public class FarthestPointVoronoi {

	private Envelope clipEnv = null; // user‐supplied envelope (optional)
	private List<Coordinate> sites = null;
	private GeometryFactory gf = new GeometryFactory();
	private DCEL dcel;
	private Envelope env;

	/**
	 * If set, every edge in the final diagram will be clipped to the larger of
	 * (this envelope) and (an envelope surrounding the input sites).
	 */
	public void setClipEnvelope(Envelope clipEnv) {
		this.clipEnv = new Envelope(clipEnv);
	}

	/**
	 * Supply your sites as JTS Coordinates.
	 */
	public void setSites(Collection<Coordinate> coords) {
		this.sites = new ArrayList<>(coords);
	}

	/**
	 * Supply your sites as a JTS Geometry (e.g. a MultiPoint, or a
	 * GeometryCollection of Points).
	 */
	public void setSites(Geometry geom) {
		this.sites = new ArrayList<>();
		for (int i = 0; i < geom.getNumGeometries(); i++) {
			Geometry g = geom.getGeometryN(i);
			if (g instanceof Point) {
				sites.add(((Point) g).getCoordinate());
			}
		}
	}

	/**
	 * Compute the farthest‐site Voronoi diagram of the sites, as regions clipped
	 * according to the user-defined envelope (if present).
	 */
	public Geometry getDiagram() {
		getDCEL();

		var clipPoly = ((Polygon) gf.toGeometry(env)).getExteriorRing();
		var segs = dcel.getEdges().stream().map(e -> gf.createLineString(new Coordinate[] { e.origVertex, e.destVertex }))
				.collect(Collectors.toList());
		var segsGeom = gf.createMultiLineString(segs.toArray(new LineString[0]));
		var geom = segsGeom.union(clipPoly); // node

		var polygonizer = new Polygonizer(false);
		polygonizer.add(geom);
		return gf.createMultiPolygon(GeometryFactory.toPolygonArray(polygonizer.getPolygons()));
	}

	public DCEL getDCEL() {
		if (dcel == null) { // lazily
			// 1) compute the sites' bounding‐box
			env = new Envelope();
			for (Coordinate c : sites) {
				env.expandToInclude(c);
			}
			if (clipEnv != null) {
				env.expandToInclude(clipEnv);
			}

			double far = env.getDiameter() * 2;

			// 3) compute convex hull of the sites
			MultiPoint mp = gf.createMultiPointFromCoords(sites.toArray(new Coordinate[0]));
			Geometry ch = new ConvexHull(mp).getConvexHull();
			List<Coordinate> hullCoords = Arrays.asList(ch.getCoordinates());

			// 4) build the FPVD
			dcel = computeFPVD(hullCoords, far);
		}
		return dcel;
	}

	private static DCEL computeFPVD(List<Coordinate> hull, final double far) {
		// routine from
		// https://dccg.upc.edu/people/vera/Applet/smallest_enclosing_circle.html
		int n = hull.size();
		DCEL dcel = new DCEL(n);

		// build a closed ring array
		Coordinate[] ring = new Coordinate[n + 1];
		for (int i = 0; i < n; i++) {
			ring[i] = hull.get(i);
		}
		ring[n] = hull.get(0);

		// ensure CW ordering for correct bisector direction
		if (Orientation.isCCW(ring)) {
			Collections.reverse(hull);
		}

		// build initial site list S and the “infinite” bisector starts
		n -= 1;
		List<Coordinate> S = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			Coordinate c0 = hull.get(i);
			Coordinate c1 = hull.get((i + 1) % n);
			S.add(new Coordinate(c0.x, c0.y, i));

			double dx = c1.x - c0.x, dy = c1.y - c0.y;
			double len = Math.hypot(dx, dy);
			// outward normal for a CCW hull
			double nx = dy / len;
			double ny = -dx / len;

			double mx = 0.5 * (c0.x + c1.x) + far * nx;
			double my = 0.5 * (c0.y + c1.y) + far * ny;
			Coordinate infinitePt = new Coordinate(mx, my);

			dcel.vertices.add(new DCELVertex(infinitePt));
			dcel.points[i] = i;
		}

		// iteratively remove sites with the largest circumcircle
		// simple n^2 loop (fine since convex hull usually has very few points)
		while (S.size() > 2) {
			int p = maximize(S);
			int q = (p - 1 + S.size()) % S.size();
			int r = (p + 1) % S.size();

			Coordinate Pq = S.get(q), Pp = S.get(p), Pr = S.get(r);
			Coordinate center = Triangle.circumcentre(Pq, Pp, Pr);

			int newVid = dcel.vertices.size();
			int vidPp = dcel.points[(int) Pp.z];
			int vidPq = dcel.points[(int) Pq.z];

			dcel.edges.add(new DCELEdge(center, dcel.vertices.get(vidPp).point, newVid, vidPp));
			dcel.edges.add(new DCELEdge(center, dcel.vertices.get(vidPq).point, newVid, vidPq));

			int e2 = dcel.edges.size() - 1;
			int e1 = e2 - 1;

			// create new real vertex
			DCELVertex v = new DCELVertex(center, e2, e1, true);
			v.setCircumcircle(Pq, Pp, Pr);
			dcel.vertices.add(v);
			dcel.points[(int) Pp.z] = newVid;

			// update Pq‐vertex to point to these half‐edges
			int vq = dcel.points[(int) Pq.z];
			dcel.vertices.set(vq, v);
			dcel.points[(int) Pq.z] = vq;

			S.remove(p);
		}

		// if two sites remain, close them
		if (S.size() == 2) {
			Coordinate A = S.get(0), B = S.get(1);
			int va = dcel.points[(int) A.z];
			int vb = dcel.points[(int) B.z];
			dcel.edges.add(new DCELEdge(dcel.vertices.get(va).point, dcel.vertices.get(vb).point, va, vb));
		}

		return dcel;
	}

	/**
	 * Finds index i of the vertex triplet (i-1,i,i+1) that form the circumcircle
	 * having the possible largest radius.
	 */
	private static int maximize(List<Coordinate> S) {
		double bestR = -1;
		int bestI = -1, N = S.size();
		for (int i = 0; i < N; i++) {
			Coordinate prev = S.get((i - 1 + N) % N);
			Coordinate cur = S.get(i);
			Coordinate next = S.get((i + 1) % N);
			double r = circumradiusMetric(prev, cur, next);
			if (r > bestR) {
				bestR = r;
				bestI = i;
			}
		}
		return bestI;
	}

	/**
	 * Computes a value proportional to the square of the circumradius of the
	 * triangle defined by three coordinates. The actual circumradius is not
	 * computed; instead, the returned value can be used for circumradius
	 * comparisons between triangles (lower values correspond to smaller
	 * circumradii).
	 * <p>
	 * This method is optimized for speed and does not compute any square roots or
	 * trigonometric functions.
	 */
	public static double circumradiusMetric(Coordinate a, Coordinate b, Coordinate c) {
		// Edge vectors
		double abx = b.x - a.x, aby = b.y - a.y;
		double acx = c.x - a.x, acy = c.y - a.y;
		double bcx = c.x - b.x, bcy = c.y - b.y;

		// Squared edge lengths
		double ab2 = abx * abx + aby * aby;
		double bc2 = bcx * bcx + bcy * bcy;
		double ca2 = acx * acx + acy * acy; // |CA|² = |AC|²

		// Squared twice-area (safer than uu*vv - uv*uv)
		double cross = abx * acy - aby * acx;
		double cross2 = cross * cross;
		if (cross2 == 0.0) // collinear or degenerate
			return Double.POSITIVE_INFINITY;

		// metric ∝ R² (the constant 4 is dropped)
		return ab2 * bc2 * ca2 / cross2;
	}

	/**
	 * A directed DCEL edge from origVertex → destVertex.
	 * <p>
	 * However, {@link #hashCode()} and {@link #equals(Object)} are implemented such
	 * that this edge is considered equal (and hashes equally) to the edge from
	 * destVertex → origVertex for undirected comparisons.
	 */
	public static class DCELEdge {
		public final Coordinate origVertex;
		public final Coordinate destVertex;
		public final int origIndex;
		public final int destIndex;

		DCELEdge(Coordinate origVertex, Coordinate destVertex, int origIndex, int destIndex) {
			this.origVertex = origVertex;
			this.destVertex = destVertex;
			this.origIndex = origIndex;
			this.destIndex = destIndex;
		}

		@Override
		public int hashCode() {
			// For vertices, order agnostic
			int verticesHash = origVertex.hashCode() ^ destVertex.hashCode();
			// For indices, order agnostic
			int indicesHash = Integer.hashCode(origIndex) ^ Integer.hashCode(destIndex);
			// Combine
			return 31 * verticesHash + indicesHash;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null || getClass() != obj.getClass())
				return false;
			DCELEdge other = (DCELEdge) obj;

			// Direction-agnostic comparison
			boolean case1 = origVertex.equals(other.origVertex) && destVertex.equals(other.destVertex) && origIndex == other.origIndex
					&& destIndex == other.destIndex;
			boolean case2 = origVertex.equals(other.destVertex) && destVertex.equals(other.origVertex) && origIndex == other.destIndex
					&& destIndex == other.origIndex;

			return case1 || case2;
		}
	}

	/**
	 * A DCEL vertex, storing a geometric point, the indices of the outgoing (next)
	 * and incoming (prev) half‐edges, and a flag real.
	 */
	public static class DCELVertex {
		public final Coordinate point;
		public Coordinate[] circumcircle;
		/** whether this is a “real” Voronoi vertex vs. the dummy ray‐start */
		public final boolean real;

		DCELVertex(Coordinate p, int next, int prev, boolean real) {
			this.point = p;
			this.real = real;
		}

		DCELVertex(Coordinate p) {
			this(p, -1, -1, false);
		}

		void setCircumcircle(Coordinate a, Coordinate b, Coordinate c) {
			this.circumcircle = new Coordinate[3];
			circumcircle[0] = a;
			circumcircle[1] = b;
			circumcircle[2] = c;
		}
	}

	/**
	 * A simple container for the DCEL. points[i] = index in vertices of the
	 * DCEL‐vertex corresponding to original hull‐vertex i.
	 */
	public static class DCEL {
		public final List<DCELEdge> edges;
		/**
		 * Non-infinite vertices comprising the FPVD.
		 */
		private final List<DCELVertex> vertices;
		private final int[] points;

		DCEL(int nHullVertices) {
			this.edges = new ArrayList<>();
			this.vertices = new ArrayList<>();
			this.points = new int[nHullVertices];
		}

		/**
		 * Returns all (inner/non-infinite) vertices comprising the FPVD edges.
		 */
		public List<DCELVertex> getVertices() {
			return new ArrayList<>(vertices.stream().filter(v -> v.real).collect(Collectors.toSet()));
		}

		/**
		 * Returns all unique edges comprising the FPVD.
		 */
		public List<DCELEdge> getEdges() {
			return new ArrayList<>(new HashSet<>(edges));
		}
	}

}