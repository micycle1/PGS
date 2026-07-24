package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.ItemVisitor;
import org.locationtech.jts.index.hprtree.HPRtree;

import net.jafama.FastMath;

/**
 * Computes 2D additively-weighted Voronoi cells (the <i>Apollonius diagram</i>)
 * for input sites encoded as {@link Coordinate}s: (x,y) is the site centre and
 * z is the additive weight (radius). The distance to site <i>i</i> is
 * {@code d_i(p) = |p - p_i| - w_i}; cell <i>i</i> is the locus where no other
 * site is closer.
 *
 * <h2>Formulation</h2>
 *
 * Around site <i>i</i>, a competitor <i>j</i> bounds the cell in direction
 * {@code u(t) = (cos t, sin t)} at radius
 *
 * <pre>
 * t_j(t) = N / (2 (dw + &lt;u, c&gt;)), c = p_j - p_i, dw = w_j - w_i, N = |c|&sup2; - dw&sup2;
 * </pre>
 *
 * and the cell boundary is the lower envelope {@code min_j t_j}. Rather than
 * sampling that envelope, this implementation works with the reciprocal radius:
 *
 * <pre>
 * rho_j(t) = 1 / t_j(t) = (dw + c_x cos t + c_y sin t) / (N/2)
 * </pre>
 *
 * which is affine in (cos t, sin t), hence a sinusoid. The lower envelope of
 * conics becomes an upper envelope of sinusoids,
 * {@code t(t) = 1 / max_j rho_j(t)}.
 *
 * <ul>
 * <li>Two sinusoids cross at most twice, so every envelope breakpoint is
 * available in closed form and nothing is discovered by sampling.</li>
 * <li>A clip half-plane contributes the same kind of reciprocal-radius
 * constraint.</li>
 * <li>For a non-dominated pair N is positive, and non-positive reciprocal
 * radius is naturally ignored by the upper envelope.</li>
 * </ul>
 *
 * <h2>Structure</h2>
 *
 * <ol>
 * <li>Coincident sites are deduplicated.</li>
 * <li>An R-tree over disk bounding boxes supplies a certified superset of the
 * competitors that can bind each cell.</li>
 * <li>The exact sinusoidal envelope determines the symbolic boundary.</li>
 * <li>Vertices are keyed by cyclically ordered triples of constraints.</li>
 * <li>Edges are keyed by their curve and symbolic endpoints, globally interned,
 * and flattened once.</li>
 * <li>Incident cells consume the same flattened edge in opposite directions,
 * making shared polygon boundaries watertight by construction.</li>
 * </ol>
 *
 * <p>
 * JTS stores the resulting polygons. If sites lie outside a supplied clipping
 * envelope, cells are first built in an enclosing polar box and JTS performs
 * the final exceptional trim against the requested envelope.
 *
 * @author Michael Carleton
 */
public class AdditivelyWeightedVoronoi {

	private static final int CLIP_XMAX = -1, CLIP_XMIN = -2, CLIP_YMAX = -3, CLIP_YMIN = -4;
	private static final double TWO_PI = 2 * Math.PI;
	private static final double ANG_EPS = 1e-12;
	private static final double MAX_PATCH_SPAN = Math.PI / 3;
	private static final double BOUNDS_SCALE = 1.5;

	private final GeometryFactory gf;
	private final double tol;

	/**
	 * @param gf             geometry factory used to create output polygons
	 * @param errorTolerance maximum sagitta of each emitted conic chord
	 */
	public AdditivelyWeightedVoronoi(GeometryFactory gf, double errorTolerance) {
		if (gf == null) {
			throw new IllegalArgumentException("GeometryFactory must not be null");
		}
		this.gf = gf;
		tol = Double.isFinite(errorTolerance) && errorTolerance > 0 ? errorTolerance : 1e-3;
	}

	/**
	 * Computes clipped polygonal approximations of the weighted Voronoi cells.
	 *
	 * @param sites  input sites; (x,y) is the centre and z the additive weight
	 * @param bounds optional clipping envelope
	 * @return non-empty cells in ascending site order, carrying the site index as
	 *         user data
	 */
	public List<Polygon> computeCells(List<Coordinate> sites, Envelope bounds) {
		if (sites == null || sites.isEmpty()) {
			return List.of();
		}
		return new Ctx(sites, bounds).compute();
	}

	// ------------------------------------------------------------------
	// per-invocation state
	// ------------------------------------------------------------------

	private final class Ctx {

		final int n;
		final double[] x, y, w;
		final boolean[] dead;
		final Envelope clip;
		final double bMinX, bMaxX, bMinY, bMaxY;
		final boolean trim;
		final double distEps, seedRadius, boxSpan;
		final DiskIndex index;

		Ctx(List<Coordinate> sites, Envelope bounds) {
			n = sites.size();
			x = new double[n];
			y = new double[n];
			w = new double[n];

			Envelope dataEnv = new Envelope();
			double maxAbsW = 0, maxAbsC = 0;

			for (int i = 0; i < n; i++) {
				Coordinate c = sites.get(i);
				if (c == null || !Double.isFinite(c.x) || !Double.isFinite(c.y)) {
					throw new IllegalArgumentException("Site " + i + " has non-finite coordinates");
				}
				x[i] = c.x;
				y[i] = c.y;
				w[i] = Double.isFinite(c.getZ()) ? c.getZ() : 0;
				maxAbsW = Math.max(maxAbsW, Math.abs(w[i]));
				maxAbsC = Math.max(maxAbsC, Math.max(Math.abs(x[i]), Math.abs(y[i])));
				dataEnv.expandToInclude(x[i], y[i]);
			}

			clip = bounds == null ? defaultClip(dataEnv, maxAbsW) : new Envelope(bounds);
			if (clip.isNull() || !(clip.getWidth() > 0) || !(clip.getHeight() > 0) || !Double.isFinite(clip.getMinX()) || !Double.isFinite(clip.getMaxX())
					|| !Double.isFinite(clip.getMinY()) || !Double.isFinite(clip.getMaxY())) {
				throw new IllegalArgumentException("Bounds must be finite and non-empty");
			}

			double scale = Math.max(Math.max(clip.getWidth(), clip.getHeight()), maxAbsC);
			if (!(scale > 0) || !Double.isFinite(scale)) {
				scale = 1;
			}
			distEps = 1e-12 * scale;

			/*
			 * The polar form requires every site strictly inside its box. When the
			 * requested clip does not provide that, use an enclosing box and trim the
			 * finished cells afterward.
			 */
			double margin = 1e-6 * scale + tol;
			boolean outside = false;
			for (int i = 0; i < n && !outside; i++) {
				outside = x[i] <= clip.getMinX() + margin || x[i] >= clip.getMaxX() - margin || y[i] <= clip.getMinY() + margin
						|| y[i] >= clip.getMaxY() - margin;
			}

			Envelope box = new Envelope(clip);
			if (outside) {
				box.expandToInclude(dataEnv);
				box.expandBy(margin);
			}

			trim = outside;
			bMinX = box.getMinX();
			bMaxX = box.getMaxX();
			bMinY = box.getMinY();
			bMaxY = box.getMaxY();
			boxSpan = box.getWidth() + box.getHeight();
			seedRadius = 1.5 * Math.sqrt(box.getWidth() * box.getHeight() / Math.max(1, n));
			index = new DiskIndex(x, y, w);

			/*
			 * Keep the heaviest coincident site, using the lowest index on ties.
			 */
			dead = new boolean[n];
			IntStream.range(0, n).parallel().forEach(i -> {
				IntList near = new IntList(8);
				index.query(x[i], y[i], 0, near);
				for (int t = 0; t < near.size; t++) {
					int j = near.get(t);
					if (j != i && x[j] == x[i] && y[j] == y[i] && (w[j] > w[i] || w[j] == w[i] && j < i)) {
						dead[i] = true;
						return;
					}
				}
			});
		}

		/**
		 * Builds symbolic cells in parallel, interns their shared geometry, flattens
		 * each distinct edge once, then assembles the polygons.
		 */
		List<Polygon> compute() {
			SCell[] cells = IntStream.range(0, n).parallel().mapToObj(this::symbolicCell).toArray(SCell[]::new);
			Graph graph = new Graph();
			for (SCell cell : cells) {
				if (cell != null) {
					graph.add(cell);
				}
			}
			graph.validate();
			graph.flatten();

			return IntStream.range(0, n).parallel().mapToObj(i -> polygons(cells[i], graph)).flatMap(List::stream).toList();
		}

		private final class Graph {

			final Map<VKey, double[]> vertices = new HashMap<>();
			final Map<EKey, Edge> edges = new HashMap<>();

			void add(SCell cell) {
				for (Use use : cell.boundary) {
					EKey key = use.key;
					double[] p0 = vertices.computeIfAbsent(key.v0, Ctx.this::vertex);
					double[] p1 = vertices.computeIfAbsent(key.v1, Ctx.this::vertex);
					Edge edge = edges.computeIfAbsent(key, k -> new Edge(k, p0, p1));
					edge.uses++;
					if (use.forward) {
						edge.forwardUses++;
					}
				}
			}

			void validate() {
				for (Edge edge : edges.values()) {
					if (edge.key.isClip()) {
						if (edge.uses != 1) {
							throw failure("Clip edge has " + edge.uses + " uses: " + edge.key);
						}
					} else if (edge.uses != 2 || edge.forwardUses != 1) {
						throw failure("Inconsistent shared edge: " + edge.key);
					}
				}
			}

			void flatten() {
				edges.values().parallelStream().forEach(e -> e.line = flattenEdge(e));
			}
		}

		// --------------------------------------------------------------
		// certified candidates and sinusoidal envelope
		// --------------------------------------------------------------

		private Cell buildCell(int i) {
			IntList ids = new IntList(16);
			seed(i, ids);

			for (int round = 0; round < 64; round++) {
				Constraint[] constraints = constraints(i, ids);
				if (constraints == null) {
					return null;
				}

				Arrays.sort(constraints);
				Cell cell = envelope(constraints);
				double radius = cell.maxRadius();

				IntList query = new IntList(8);
				query(i, radius, query);
				boolean grew = false;

				for (int t = 0; t < query.size; t++) {
					int j = query.get(t);
					if (!ids.contains(j)) {
						ids.add(j);
						grew = true;
					}
				}
				if (!grew) {
					return cell;
				}
			}
			throw failure("Candidate certification did not converge for site " + i);
		}

		private Constraint[] constraints(int i, IntList ids) {
			List<Constraint> out = new ArrayList<>(ids.size + 4);
			out.add(clipConstraint(i, CLIP_XMAX));
			out.add(clipConstraint(i, CLIP_XMIN));
			out.add(clipConstraint(i, CLIP_YMAX));
			out.add(clipConstraint(i, CLIP_YMIN));

			for (int t = 0; t < ids.size; t++) {
				int j = ids.get(t);
				if (j == i || dead[j]) {
					continue;
				}

				double dx = x[j] - x[i], dy = y[j] - y[i], dw = w[j] - w[i];
				double d = FastMath.hypot(dx, dy);
				if (dw >= d - distEps) {
					return null;
				}

				double den = 0.5 * (d * d - dw * dw);
				if (den > 0) {
					out.add(new Constraint(j, dw, dx, dy, den));
				}
			}
			return out.toArray(Constraint[]::new);
		}

		/**
		 * Opens the certified loop with a cheap radius guess. Nothing here needs to be
		 * complete; buildCell certifies the final candidate set.
		 */
		private void seed(int i, IntList out) {
			for (double radius = seedRadius;; radius *= 2) {
				out.clear();
				index.query(x[i], y[i], radius, out);

				int usable = 0;
				for (int t = 0; t < out.size; t++) {
					int j = out.get(t);
					if (j != i && !dead[j]) {
						usable++;
					}
				}
				if (usable >= 4) {
					return;
				}
				if (!Double.isFinite(radius) || radius >= boxSpan) {
					out.clear();
					index.all(out);
					return;
				}
			}
		}

		/**
		 * All j that can bind within radius r: {@code (d - w_j + w_i)/2 < r}.
		 */
		private void query(int i, double radius, IntList out) {
			IntList raw = new IntList(16);
			double expand = Double.isFinite(radius) ? Math.max(0, 2 * radius - w[i]) : Double.POSITIVE_INFINITY;

			if (Double.isFinite(expand)) {
				index.query(x[i], y[i], expand, raw);
			} else {
				index.all(raw);
			}

			for (int t = 0; t < raw.size; t++) {
				int j = raw.get(t);
				if (j == i || dead[j]) {
					continue;
				}
				double d = FastMath.hypot(x[j] - x[i], y[j] - y[i]);
				if (0.5 * (d - w[j] + w[i]) < radius) {
					out.add(j);
				}
			}
		}

		// --------------------------------------------------------------
		// symbolic topology
		// --------------------------------------------------------------

		private SCell symbolicCell(int i) {
			if (dead[i]) {
				return null;
			}

			Cell cell = buildCell(i);
			if (cell == null) {
				return null;
			}
			if (cell.k < 2) {
				throw failure("Degenerate envelope for site " + i);
			}

			List<Use> boundary = new ArrayList<>(cell.k);
			for (int s = 0; s < cell.k; s++) {
				int prev = cell.cs[cell.owner[(s - 1 + cell.k) % cell.k]].id;
				int owner = cell.cs[cell.owner[s]].id;
				int next = cell.cs[cell.owner[(s + 1) % cell.k]].id;

				VKey from = VKey.of(i, prev, owner);
				VKey to = VKey.of(i, owner, next);
				if (from.equals(to)) {
					continue;
				}

				EKey key = EKey.of(i, owner, from, to);
				boundary.add(new Use(key, key.v0.equals(from)));
			}

			if (boundary.size() < 2) {
				throw failure("Degenerate symbolic boundary for site " + i);
			}
			return new SCell(i, List.copyOf(boundary));
		}

		/**
		 * Solves one cyclically identified vertex in the frame of its lowest numbered
		 * incident site.
		 */
		private double[] vertex(VKey key) {
			int site = key.a, previous = key.b, current = key.c;
			Constraint p = constraintOf(site, previous), q = constraintOf(site, current);
			if (p == null || q == null) {
				throw failure("Invalid vertex " + key);
			}

			double dA = q.a * p.den - p.a * q.den;
			double dB = q.b * p.den - p.b * q.den;
			double dC = q.c * p.den - p.c * q.den;
			double magnitude = FastMath.hypot(dB, dC);

			if (!(magnitude > 0) || !Double.isFinite(magnitude)) {
				throw failure("Unresolved vertex " + key);
			}

			double ratio = Math.max(-1, Math.min(1, -dA / magnitude));
			double theta = FastMath.atan2(dC, dB) - FastMath.acos(ratio);
			double cos = FastMath.cos(theta), sin = FastMath.sin(theta);
			double num = p.num(cos, sin);

			if (!(num > 0)) {
				throw failure("Non-positive vertex radius at " + key);
			}

			double radius = p.den / num;
			double px = x[site] + radius * cos, py = y[site] + radius * sin;
			if (!(radius > 0) || !Double.isFinite(px) || !Double.isFinite(py)) {
				throw failure("Non-finite vertex " + key);
			}
			return new double[] { px, py };
		}

		/**
		 * Flattens one globally interned edge in canonical endpoint order.
		 */
		private List<Coordinate> flattenEdge(Edge edge) {
			if (edge.key.isClip() || w[edge.key.a] == w[edge.key.b]) {
				return List.of(coord(edge.p0), coord(edge.p1));
			}

			List<double[]> arc = canonicalArc(edge.key.a, edge.key.b, edge.p0, edge.p1);
			if (arc == null || arc.size() < 2) {
				throw failure("Could not construct edge " + edge.key);
			}

			List<Coordinate> line = new ArrayList<>(arc.size());
			for (double[] p : arc) {
				if (!finite(p)) {
					throw failure("Non-finite point on edge " + edge.key);
				}
				line.add(coord(p));
			}

			/*
			 * Force the globally interned endpoint values, even if an equivalent
			 * parameterisation path differs in its last bit.
			 */
			line.set(0, coord(edge.p0));
			line.set(line.size() - 1, coord(edge.p1));
			return List.copyOf(line);
		}

		// --------------------------------------------------------------
		// polygon assembly
		// --------------------------------------------------------------

		private List<Polygon> polygons(SCell cell, Graph graph) {
			if (cell == null) {
				return List.of();
			}

			List<Coordinate> ring = new ArrayList<>(4 * cell.boundary.size() + 1);
			for (Use use : cell.boundary) {
				Edge edge = graph.edges.get(use.key);
				if (edge == null || edge.line == null) {
					throw failure("Missing flattened edge " + use.key);
				}
				append(ring, edge.line, use.forward);
			}

			if (ring.size() < 4 || !ring.get(0).equals2D(ring.get(ring.size() - 1))) {
				throw failure("Open or degenerate ring for site " + cell.site);
			}

			Polygon polygon = gf.createPolygon(ring.toArray(Coordinate[]::new));
			if (!trim) {
				polygon.setUserData(cell.site);
				return List.of(polygon);
			}

			/*
			 * Compatibility path for sites outside the requested clip. Shared diagram edges
			 * were still flattened only once; JTS performs only this final exceptional
			 * overlay.
			 */
			Geometry clipped = polygon.intersection(gf.toGeometry(clip));
			List<Polygon> out = new ArrayList<>();
			for (int i = 0; i < clipped.getNumGeometries(); i++) {
				Geometry geometry = clipped.getGeometryN(i);
				if (geometry instanceof Polygon p && !p.isEmpty()) {
					p.setUserData(cell.site);
					out.add(p);
				}
			}
			return out;
		}

		private void append(List<Coordinate> ring, List<Coordinate> edge, boolean forward) {
			if (forward) {
				int from = ring.isEmpty() ? 0 : 1;
				if (!ring.isEmpty() && !ring.get(ring.size() - 1).equals2D(edge.get(0))) {
					throw failure("Disconnected edge");
				}
				for (int i = from; i < edge.size(); i++) {
					ring.add(new Coordinate(edge.get(i)));
				}
			} else {
				int from = ring.isEmpty() ? edge.size() - 1 : edge.size() - 2;
				if (!ring.isEmpty() && !ring.get(ring.size() - 1).equals2D(edge.get(edge.size() - 1))) {
					throw failure("Disconnected reversed edge");
				}
				for (int i = from; i >= 0; i--) {
					ring.add(new Coordinate(edge.get(i)));
				}
			}
		}

		// --------------------------------------------------------------
		// conic construction
		// --------------------------------------------------------------

		/**
		 * Flattens the hyperbolic a|b arc from p0 to p2 as exact rational quadratic
		 * patches.
		 */
		private List<double[]> canonicalArc(int a, int b, double[] p0, double[] p2) {
			Constraint constraint = siteConstraint(a, b);
			if (constraint == null) {
				return null;
			}

			double xa = x[a], ya = y[a];
			double t0 = FastMath.atan2(p0[1] - ya, p0[0] - xa);
			double t2 = FastMath.atan2(p2[1] - ya, p2[0] - xa);

			/*
			 * rho is least at gap and non-positive there, so the correct branch is the
			 * sweep that excludes it.
			 */
			double gap = FastMath.atan2(constraint.c, constraint.b) + Math.PI;
			double ccw = wrap(t2 - t0);
			boolean fromP0 = !(wrap(gap - t0) < ccw);
			double start = fromP0 ? t0 : t2;
			double span = fromP0 ? ccw : TWO_PI - ccw;

			if (!(span > 1e-13) || !Double.isFinite(span)) {
				return null;
			}

			int parts = Math.max(1, (int) Math.ceil(span / MAX_PATCH_SPAN));
			double[][] node = new double[parts + 1][];
			node[0] = fromP0 ? p0 : p2;
			node[parts] = fromP0 ? p2 : p0;

			for (int i = 1; i < parts; i++) {
				node[i] = onCurve(a, constraint, start + span * i / parts);
				if (node[i] == null) {
					return null;
				}
			}

			List<double[]> out = new ArrayList<>(4 * parts + 4);
			out.add(node[0]);

			for (int i = 0; i < parts; i++) {
				double[] shoulder = onCurve(a, constraint, start + span * (i + 0.5) / parts);
				if (shoulder == null) {
					return null;
				}
				patch(a, b, node[i], node[i + 1], shoulder, out);
			}

			if (!fromP0) {
				Collections.reverse(out);
			}
			return out;
		}

		/**
		 * Constructs one rational quadratic from its endpoints, end tangents and one
		 * interior point, then flattens it.
		 */
		private void patch(int a, int b, double[] q0, double[] q2, double[] s, List<double[]> out) {
			double[] t0 = tangent(a, b, q0), t2 = tangent(a, b, q2);
			double[] p1 = t0 == null || t2 == null ? null : meet(q0, t0, q2, t2);
			if (p1 == null) {
				throw failure("Could not construct conic control point for " + a + "|" + b);
			}

			double alpha = cross(p1[0] - s[0], p1[1] - s[1], q2[0] - s[0], q2[1] - s[1]);
			double beta = cross(q2[0] - s[0], q2[1] - s[1], q0[0] - s[0], q0[1] - s[1]);
			double gamma = cross(q0[0] - s[0], q0[1] - s[1], p1[0] - s[0], p1[1] - s[1]);

			boolean ag = alpha > 0 && gamma > 0 || alpha < 0 && gamma < 0;
			boolean ab = alpha > 0 && beta > 0 || alpha < 0 && beta < 0;
			double scale = Math.max(Math.abs(alpha), Math.max(Math.abs(beta), Math.abs(gamma)));

			if (!ag || !ab || !(scale > 0) || !Double.isFinite(scale)) {
				throw failure("Invalid conic barycentrics for " + a + "|" + b);
			}

			alpha /= scale;
			beta /= scale;
			gamma /= scale;
			double weight = Math.abs(beta) / (2 * Math.sqrt(alpha * gamma));

			if (!(weight > 0) || !Double.isFinite(weight)) {
				throw failure("Invalid conic weight for " + a + "|" + b);
			}
			flatten(q0, p1, q2, weight, tol, out, 0);
		}

		private double[] tangent(int a, int b, double[] p) {
			double ax = p[0] - x[a], ay = p[1] - y[a];
			double bx = p[0] - x[b], by = p[1] - y[b];
			double la = FastMath.hypot(ax, ay), lb = FastMath.hypot(bx, by);

			if (!(la > 0) || !(lb > 0)) {
				return null;
			}

			double tx = ax / la + bx / lb, ty = ay / la + by / lb;
			return tx * tx + ty * ty > 0 && Double.isFinite(tx) && Double.isFinite(ty) ? new double[] { tx, ty } : null;
		}

		private double[] onCurve(int site, Constraint constraint, double theta) {
			double cos = FastMath.cos(theta), sin = FastMath.sin(theta);
			double num = constraint.num(cos, sin);
			if (!(num > 0)) {
				return null;
			}

			double radius = constraint.den / num;
			double px = x[site] + radius * cos, py = y[site] + radius * sin;
			return radius > 0 && Double.isFinite(px) && Double.isFinite(py) ? new double[] { px, py } : null;
		}

		private Constraint constraintOf(int frame, int id) {
			return id >= 0 ? siteConstraint(frame, id) : clipConstraint(frame, id);
		}

		private Constraint siteConstraint(int i, int j) {
			double dx = x[j] - x[i], dy = y[j] - y[i], dw = w[j] - w[i];
			double den = 0.5 * (dx * dx + dy * dy - dw * dw);
			return den > 0 ? new Constraint(j, dw, dx, dy, den) : null;
		}

		private Constraint clipConstraint(int i, int id) {
			return switch (id) {
				case CLIP_XMAX -> new Constraint(id, 0, 1, 0, bMaxX - x[i]);
				case CLIP_XMIN -> new Constraint(id, 0, -1, 0, x[i] - bMinX);
				case CLIP_YMAX -> new Constraint(id, 0, 0, 1, bMaxY - y[i]);
				case CLIP_YMIN -> new Constraint(id, 0, 0, -1, y[i] - bMinY);
				default -> throw new IllegalArgumentException("Unknown clip id " + id);
			};
		}
	}

	// ------------------------------------------------------------------
	// symbolic keys and shared edges
	// ------------------------------------------------------------------

	/**
	 * Cyclically ordered constraints meeting at a vertex. Rotation is
	 * canonicalised; reversal is not, since it distinguishes the two possible
	 * Apollonius vertices of the same triple.
	 */
	private record VKey(int a, int b, int c) implements Comparable<VKey> {

		static VKey of(int a, int b, int c) {
			int[] v = { a, b, c };
			int first = 0;
			for (int i = 1; i < 3; i++) {
				if (v[i] >= 0 && (v[first] < 0 || v[i] < v[first])) {
					first = i;
				}
			}
			return new VKey(v[first], v[(first + 1) % 3], v[(first + 2) % 3]);
		}

		@Override
		public int compareTo(VKey o) {
			int c0 = Integer.compare(a, o.a);
			if (c0 == 0) {
				c0 = Integer.compare(b, o.b);
			}
			return c0 == 0 ? Integer.compare(c, o.c) : c0;
		}
	}

	/**
	 * A curve and its symbolically ordered endpoints. For a site bisector a and b
	 * are the ordered site ids; for a clip edge both equal the negative clip id.
	 */
	private record EKey(int a, int b, VKey v0, VKey v1) {

		static EKey of(int site, int owner, VKey from, VKey to) {
			int a = owner < 0 ? owner : Math.min(site, owner);
			int b = owner < 0 ? owner : Math.max(site, owner);
			return from.compareTo(to) <= 0 ? new EKey(a, b, from, to) : new EKey(a, b, to, from);
		}

		boolean isClip() {
			return a < 0;
		}
	}

	private record Use(EKey key, boolean forward) {
	}

	private record SCell(int site, List<Use> boundary) {
	}

	private static final class Edge {

		final EKey key;
		final double[] p0, p1;
		int uses, forwardUses;
		List<Coordinate> line;

		Edge(EKey key, double[] p0, double[] p1) {
			this.key = key;
			this.p0 = p0;
			this.p1 = p1;
		}
	}

	// ------------------------------------------------------------------
	// upper envelope of sinusoids
	// ------------------------------------------------------------------

	private static Cell envelope(Constraint[] constraints) {
		Cell envelope = new Cell(constraints);
		envelope.k = 1;
		envelope.owner[0] = 0;
		envelope.start[0] = 0;
		double radius = envelope.maxRadius();

		IntList owners = new IntList(16);
		DoubleList starts = new DoubleList(16);
		double[] roots = new double[2], split = new double[2];

		for (int f = 1; f < constraints.length; f++) {
			if (constraints[f].minT >= radius) {
				break;
			}
			envelope.insert(f, owners, starts, roots, split);
			radius = envelope.maxRadius();
		}
		return envelope;
	}

	private static int crossings(Constraint f, Constraint g, double[] out) {
		double dA = f.a * g.den - g.a * f.den;
		double dB = f.b * g.den - g.b * f.den;
		double dC = f.c * g.den - g.c * f.den;
		double radius = FastMath.hypot(dB, dC);

		if (!(radius > Math.abs(dA))) {
			return 0;
		}

		double phi = FastMath.atan2(dC, dB);
		double halfWidth = FastMath.acos(-dA / radius);
		out[0] = phi - halfWidth;
		out[1] = phi + halfWidth;
		return 2;
	}

	private static final class Cell {

		final Constraint[] cs;
		int[] owner;
		double[] start;
		int k;

		Cell(Constraint[] constraints) {
			cs = constraints;
			owner = new int[4 * constraints.length + 8];
			start = new double[owner.length];
		}

		double end(int s) {
			return s + 1 < k ? start[s + 1] : start[0] + TWO_PI;
		}

		void insert(int f, IntList newOwners, DoubleList newStarts, double[] roots, double[] split) {
			newOwners.clear();
			newStarts.clear();
			Constraint candidate = cs[f];

			for (int s = 0; s < k; s++) {
				int incumbentId = owner[s];
				Constraint incumbent = cs[incumbentId];
				double a0 = start[s], a1 = end(s);
				int rootCount = crossings(candidate, incumbent, roots);
				int splitCount = 0;

				for (int t = 0; t < rootCount; t++) {
					double root = a0 + wrap(roots[t] - a0);
					if (root > a0 + ANG_EPS && root < a1 - ANG_EPS) {
						split[splitCount++] = root;
					}
				}
				if (splitCount == 2 && split[0] > split[1]) {
					double tmp = split[0];
					split[0] = split[1];
					split[1] = tmp;
				}

				double low = a0;
				for (int t = 0; t <= splitCount; t++) {
					double high = t < splitCount ? split[t] : a1;
					double middle = 0.5 * (low + high);
					double cos = FastMath.cos(middle), sin = FastMath.sin(middle);
					int winner = candidate.num(cos, sin) * incumbent.den > incumbent.num(cos, sin) * candidate.den ? f : incumbentId;

					if (newOwners.size == 0 || newOwners.last() != winner) {
						newOwners.add(winner);
						newStarts.add(low);
					}
					low = high;
				}
			}

			int size = newOwners.size;
			if (size == 0) {
				return;
			}

			int from = size > 1 && newOwners.get(0) == newOwners.get(size - 1) ? 1 : 0;
			k = size - from;

			if (k > owner.length) {
				owner = Arrays.copyOf(owner, 2 * k);
				start = Arrays.copyOf(start, 2 * k);
			}
			for (int i = 0; i < k; i++) {
				owner[i] = newOwners.get(i + from);
				start[i] = newStarts.get(i + from);
			}
			if (k == 1) {
				start[0] = 0;
			}
		}

		double maxRadius() {
			double minimum = Double.POSITIVE_INFINITY;

			for (int s = 0; s < k; s++) {
				Constraint constraint = cs[owner[s]];
				double a0 = start[s], a1 = end(s);

				minimum = Math.min(minimum, constraint.rho(FastMath.cos(a0), FastMath.sin(a0)));
				minimum = Math.min(minimum, constraint.rho(FastMath.cos(a1), FastMath.sin(a1)));

				double low = FastMath.atan2(constraint.c, constraint.b) + Math.PI;
				if (wrap(low - a0) < a1 - a0) {
					minimum = Math.min(minimum, constraint.rho(FastMath.cos(low), FastMath.sin(low)));
				}
			}
			return minimum > 0 ? 1 / minimum : Double.POSITIVE_INFINITY;
		}
	}

	private static final class Constraint implements Comparable<Constraint> {

		final int id;
		final double a, b, c, den, minT;

		Constraint(int id, double a, double b, double c, double den) {
			this.id = id;
			this.a = a;
			this.b = b;
			this.c = c;
			this.den = den;
			minT = den / (a + FastMath.hypot(b, c));
		}

		double num(double cos, double sin) {
			return a + b * cos + c * sin;
		}

		double rho(double cos, double sin) {
			return num(cos, sin) / den;
		}

		@Override
		public int compareTo(Constraint other) {
			return Double.compare(minT, other.minT);
		}
	}

	// ------------------------------------------------------------------
	// rational quadratic flattening
	// ------------------------------------------------------------------

	/**
	 * Flattens a rational quadratic by shoulder subdivision. A segment is accepted
	 * only when its analytic sagitta does not exceed tolerance.
	 */
	private static void flatten(double[] p0, double[] p1, double[] p2, double weight, double tolerance, List<double[]> out, int depth) {

		double ex = p2[0] - p0[0], ey = p2[1] - p0[1];
		double length = FastMath.hypot(ex, ey);
		if (!(length > 0) || !(weight > 0) || !Double.isFinite(weight)) {
			throw failure("Degenerate rational quadratic");
		}

		double controlDistance = Math.abs(ex / length * (p1[1] - p0[1]) - ey / length * (p1[0] - p0[0]));
		double sagitta = controlDistance * weight / (1 + weight);

		if (Double.isFinite(sagitta) && sagitta <= tolerance) {
			out.add(p2);
			return;
		}
		if (depth >= 64) {
			throw failure("Rational quadratic did not converge");
		}

		double alpha = weight / (1 + weight);
		double[] q1 = { p0[0] + alpha * (p1[0] - p0[0]), p0[1] + alpha * (p1[1] - p0[1]) };
		double[] q2 = { p2[0] + alpha * (p1[0] - p2[0]), p2[1] + alpha * (p1[1] - p2[1]) };
		double[] shoulder = { q1[0] + 0.5 * (q2[0] - q1[0]), q1[1] + 0.5 * (q2[1] - q1[1]) };
		double childWeight = Math.sqrt(0.5 * (1 + weight));

		if (!finite(q1) || !finite(q2) || !finite(shoulder) || !Double.isFinite(childWeight)) {
			throw failure("Non-finite rational quadratic");
		}

		flatten(p0, q1, shoulder, childWeight, tolerance, out, depth + 1);
		flatten(shoulder, q2, p2, childWeight, tolerance, out, depth + 1);
	}

	// ------------------------------------------------------------------
	// helpers
	// ------------------------------------------------------------------

	private static Envelope defaultClip(Envelope dataEnv, double maxAbsWeight) {
		double diagonal = dataEnv.getDiameter();
		if (!(diagonal > 0)) {
			diagonal = 1;
		}
		double half = BOUNDS_SCALE * (diagonal + 2 * maxAbsWeight + 1);
		Coordinate centre = dataEnv.centre();
		return new Envelope(centre.x - half, centre.x + half, centre.y - half, centre.y + half);
	}

	private static double[] meet(double[] p, double[] u, double[] q, double[] v) {
		double denominator = cross(u[0], u[1], v[0], v[1]);
		if (denominator == 0 || !Double.isFinite(denominator)) {
			return null;
		}
		double t = cross(q[0] - p[0], q[1] - p[1], v[0], v[1]) / denominator;
		double px = p[0] + t * u[0], py = p[1] + t * u[1];
		return Double.isFinite(px) && Double.isFinite(py) ? new double[] { px, py } : null;
	}

	private static double cross(double ax, double ay, double bx, double by) {
		return ax * by - ay * bx;
	}

	private static double wrap(double angle) {
		double value = angle % TWO_PI;
		return value < 0 ? value + TWO_PI : value;
	}

	private static boolean finite(double[] point) {
		return Double.isFinite(point[0]) && Double.isFinite(point[1]);
	}

	private static Coordinate coord(double[] point) {
		return new Coordinate(point[0], point[1]);
	}

	private static IllegalStateException failure(String message) {
		return new IllegalStateException(message);
	}

	// ------------------------------------------------------------------
	// disk R-tree
	// ------------------------------------------------------------------

	private static final class DiskIndex {

		private final HPRtree tree = new HPRtree();
		private final Integer[] boxed;

		DiskIndex(double[] x, double[] y, double[] w) {
			boxed = new Integer[x.length];
			for (int i = 0; i < x.length; i++) {
				boxed[i] = i;
				double radius = Math.max(w[i], 0);
				tree.insert(new Envelope(x[i] - radius, x[i] + radius, y[i] - radius, y[i] + radius), boxed[i]);
			}
			tree.build();
		}

		void query(double x, double y, double expand, IntList out) {
			if (!Double.isFinite(expand)) {
				all(out);
				return;
			}
			Envelope envelope = new Envelope(x, x, y, y);
			if (expand > 0) {
				envelope.expandBy(expand);
			}
			tree.query(envelope, new Collector(out));
		}

		void all(IntList out) {
			for (int i = 0; i < boxed.length; i++) {
				out.add(i);
			}
		}

		private record Collector(IntList out) implements ItemVisitor {
			@Override
			public void visitItem(Object item) {
				out.add((Integer) item);
			}
		}
	}

	// ------------------------------------------------------------------
	// primitive lists
	// ------------------------------------------------------------------

	private static final class IntList {

		int[] values;
		int size;

		IntList(int capacity) {
			values = new int[Math.max(4, capacity)];
		}

		void add(int value) {
			if (size == values.length) {
				values = Arrays.copyOf(values, size * 2);
			}
			values[size++] = value;
		}

		int get(int index) {
			return values[index];
		}

		int last() {
			return values[size - 1];
		}

		void clear() {
			size = 0;
		}

		boolean contains(int value) {
			for (int i = 0; i < size; i++) {
				if (values[i] == value) {
					return true;
				}
			}
			return false;
		}
	}

	private static final class DoubleList {

		double[] values;
		int size;

		DoubleList(int capacity) {
			values = new double[Math.max(4, capacity)];
		}

		void add(double value) {
			if (size == values.length) {
				values = Arrays.copyOf(values, size * 2);
			}
			values[size++] = value;
		}

		double get(int index) {
			return values[index];
		}

		void clear() {
			size = 0;
		}
	}
}