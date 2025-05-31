package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Triangle;
import org.locationtech.jts.operation.polygonize.Polygonizer;

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
		var segs = dcel.edges.stream().map(e -> (Geometry) gf.createLineString(new Coordinate[] { e.origVertex, e.destVertex })).collect(Collectors.toList());
		segs.add(clipPoly);
		var geom = segs.stream().reduce((a, b) -> a.union(b)).get();

		var polygonizer = new Polygonizer();
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

			// 3) compute convex hull of the sites (so hull is in CCW order)
			MultiPoint mp = gf.createMultiPointFromCoords(sites.toArray(new Coordinate[0]));
			Geometry ch = new ConvexHull(mp).getConvexHull();
			List<Coordinate> hullCoords = Arrays.asList(ch.getCoordinates());

			// 4) build the FPVD
			dcel = computeFPVD(hullCoords, far);
		}
		return dcel;
	}

	private static DCEL computeFPVD(List<Coordinate> hull, double far) {
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
			dcel.vertices.add(new DCELVertex(center, e2, e1, true));
			dcel.points[(int) Pp.z] = newVid;

			// update Pq‐vertex to point to these half‐edges
			int vq = dcel.points[(int) Pq.z];
			dcel.vertices.set(vq, new DCELVertex(center, e2, e1, true));
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

	private static int maximize(List<Coordinate> S) {
		double bestR = -1;
		int bestI = -1, N = S.size();
		for (int i = 0; i < N; i++) {
			Coordinate prev = S.get((i - 1 + N) % N);
			Coordinate cur = S.get(i);
			Coordinate next = S.get((i + 1) % N);
			double r = Triangle.circumradius(prev, cur, next);
			if (r > bestR) {
				bestR = r;
				bestI = i;
			}
		}
		return bestI;
	}

	/**
	 * A directed DCEL edge from origVertex → destVertex.
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
	}

	/**
	 * A DCEL vertex, storing a geometric point, the indices of the outgoing (next)
	 * and incoming (prev) half‐edges, and a flag real.
	 */
	public static class DCELVertex {
		public final Coordinate point;
		/** whether this is a “real” Voronoi vertex vs. the dummy ray‐start */
		private final boolean real;

		DCELVertex(Coordinate p, int next, int prev, boolean real) {
			this.point = p;
			this.real = real;
		}

		DCELVertex(Coordinate p) {
			this(p, -1, -1, false);
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
		public final List<DCELVertex> vertices;
		private final int[] points;

		DCEL(int nHullVertices) {
			this.edges = new ArrayList<>();
			this.vertices = new ArrayList<>();
			this.points = new int[nHullVertices];
		}
	}

}