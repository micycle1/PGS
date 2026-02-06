package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import org.locationtech.jts.algorithm.CGAlgorithmsDD;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;

/**
 * Contour regularization for JTS geometries.
 *
 * <h2>Regularities reinforced</h2>
 * <ul>
 * <li><b>Parallelism</b> (by snapping edge orientations via the direction
 * model)</li>
 * <li><b>Orthogonality</b> (optional insertion step)</li>
 * <li><b>Collinearity</b> (merge consecutive collinear edges)</li>
 * </ul>
 *
 * <h2>Notes / limitations</h2>
 * <ul>
 * <li>No guarantee of topological validity (self-intersections can occur).</li>
 * <li>Best suited to man-made / rectilinear-ish shapes.</li>
 * </ul>
 *
 * <p>
 * Conceptually inspired by CGAL Shape Regularization package.
 * </p>
 * 
 * @author Michael Carleton
 */
public final class ContourRegularization {

	private ContourRegularization() {
	}

	/**
	 * Parameters controlling detection/merge thresholds and optional steps.
	 *
	 * <p>
	 * Defaults are chosen to be similar in spirit to CGAL defaults, but are not
	 * numerically identical.
	 * </p>
	 */
	public static final class RegParameters {
		/**
		 * Angle threshold (degrees) used when testing "near parallel" between
		 * consecutive edges. Typical values: 3..10.
		 */
		public final double parallelAngleThresholdDeg;

		/**
		 * Maximum orthogonal distance between two consecutive parallel edges to
		 * consider them collinear (mergeable). Units are coordinate units.
		 */
		public final double maximumOffset;

		/**
		 * Edges shorter than this length are dropped before optimization.
		 */
		public final double minEdgeLength;

		/**
		 * If true, inserts an orthogonal edge when two consecutive edges are parallel.
		 */
		public final boolean insertOrthogonalWhenParallel;

		/**
		 * If true, for open contours, output keeps the original first and last vertex
		 * (XY). If false, endpoints can move due to rotation / reconnection.
		 */
		public final boolean preserveOpenEndpoints;

		public final ContourDirections directions;

		private RegParameters(Builder b) {
			this.parallelAngleThresholdDeg = b.parallelAngleThresholdDeg;
			this.maximumOffset = b.maximumOffset;
			this.minEdgeLength = b.minEdgeLength;
			this.insertOrthogonalWhenParallel = b.insertOrthogonalWhenParallel;
			this.preserveOpenEndpoints = b.preserveOpenEndpoints;
			this.directions = b.directions;
		}

		public static RegParameters defaults() {
			return builder().build();
		}

		public static Builder builder() {
			return new Builder();
		}

		public static final class Builder {
			private double parallelAngleThresholdDeg = 5.0;
			private double maximumOffset = 0.5;
			private double minEdgeLength = 1e-9;
			private boolean insertOrthogonalWhenParallel = true;
			private boolean preserveOpenEndpoints = true;
			private ContourDirections directions = new LongestEdgeDirections();

			public Builder parallelAngleThresholdDeg(double deg) {
				this.parallelAngleThresholdDeg = deg;
				return this;
			}

			public Builder maximumOffset(double v) {
				this.maximumOffset = v;
				return this;
			}

			public Builder minEdgeLength(double v) {
				this.minEdgeLength = v;
				return this;
			}

			public Builder insertOrthogonalWhenParallel(boolean v) {
				this.insertOrthogonalWhenParallel = v;
				return this;
			}

			public Builder preserveOpenEndpoints(boolean v) {
				this.preserveOpenEndpoints = v;
				return this;
			}

			public RegParameters build() {
				return new RegParameters(this);
			}

			public Builder directions(ContourDirections directions) {
				this.directions = directions;
				return this;
			}
		}
	}

	public interface ContourDirections {

		/**
		 * Creates an initialized direction model for a specific contour.
		 *
		 * <p>
		 * Recommended contract: return a NEW instance (do not mutate and return
		 * {@code this}), so that a single {@link RegParameters} instance can be reused
		 * safely across calls/threads.
		 * </p>
		 *
		 * @param coordinates ordered coordinates (for closed rings may include
		 *                    duplicate last==first)
		 * @param closed      whether to treat the contour as closed
		 * @return initialized direction model for this contour
		 */
		ContourDirections init(Coordinate[] coordinates, boolean closed);

		/**
		 * Orients an edge (segment) in-place toward the principal direction assigned to
		 * that edge.
		 *
		 * @param edgeIndex edge index in [0..numEdges-1]
		 * @param segment   segment to mutate
		 */
		void orient(int edgeIndex, LineSegment segment);
	}

	/**
	 * Default direction model: uses the orientation of the longest edge as the
	 * principal direction.
	 */
	public static final class LongestEdgeDirections implements ContourDirections {

		private final double refOrientationDeg; // only meaningful after init()
		private final boolean initialized;

		/** Prototype constructor (no contour yet). */
		public LongestEdgeDirections() {
			this.refOrientationDeg = 0.0;
			this.initialized = false;
		}

		private LongestEdgeDirections(double refOrientationDeg) {
			this.refOrientationDeg = refOrientationDeg;
			this.initialized = true;
		}

		@Override
		public ContourDirections init(Coordinate[] coordinates, boolean closed) {
			Coordinate[] pts = closed ? sanitizeClosed(coordinates) : copy(coordinates);
			double ref = computeLongestEdgeOrientationDeg(pts, closed);
			return new LongestEdgeDirections(ref);
		}

		@Override
		public void orient(int edgeIndex, LineSegment segment) {
			if (!initialized)
				throw new IllegalStateException("Not initialized; call init() first");
			double segOri = Geometry2D.orientationDeg(segment);
			double rot = Geometry2D.mod90AngleDifferenceDeg(segOri, refOrientationDeg);
			Geometry2D.rotateSegmentCCWAroundMidpointInPlace(segment, rot);
		}

		// same helper as before (use your existing)
		private static double computeLongestEdgeOrientationDeg(Coordinate[] pts, boolean closed) {
			double bestLen2 = -1.0, bestOri = 0.0;
			int n = pts.length, limit = closed ? n : (n - 1);
			for (int i = 0; i < limit; i++) {
				Coordinate a = pts[i];
				Coordinate b = pts[closed ? (i + 1) % n : (i + 1)];
				double dx = b.x - a.x, dy = b.y - a.y;
				double len2 = dx * dx + dy * dy;
				if (len2 > bestLen2) {
					bestLen2 = len2;
					bestOri = Geometry2D.orientationDeg(new LineSegment(a, b));
				}
			}
			return bestOri;
		}
	}

	/**
	 * Principal-direction model that automatically infers one or more dominant
	 * directions from the contour.
	 *
	 * <h2>Purpose</h2>
	 * <p>
	 * This strategy estimates a small set of dominant axes (principal directions)
	 * from the contour itself, then assigns each edge to one of these axes,
	 * enabling subsequent snapping (rotation) of edges to the inferred structure.
	 * </p>
	 *
	 * <h2>Algorithm (CGAL-inspired, simplified)</h2>
	 * <ol>
	 * <li>Compute all edge orientations and lengths.</li>
	 * <li>Mark edges shorter than {@link Options#minimumLength} as invalid for
	 * seeding axes.</li>
	 * <li>Sort edges by length descending.</li>
	 * <li>Iteratively pick the next longest unused valid edge as a new axis
	 * (seed).</li>
	 * <li>Assign other unused valid edges to that axis if they are near-parallel or
	 * near-orthogonal (within {@link Options#maximumAngleDeg}).</li>
	 * <li>Edges not assigned during axis discovery are filled in by propagating the
	 * nearest assigned axis along the contour ("unify along contour"), ensuring
	 * every edge has an axis index.</li>
	 * <li>If {@link Options#adjustDirections} is enabled, each discovered axis is
	 * slightly rotated by the average residual (mod-90) error of the edges assigned
	 * to it, yielding a better fit.</li>
	 * </ol>
	 *
	 * <h2>How this differs from {@link UserDefinedDirections}</h2>
	 * <ul>
	 * <li>{@code MultipleDirections} discovers axes from the contour data.</li>
	 * <li>{@code UserDefinedDirections} uses axes supplied by the user.</li>
	 * </ul>
	 *
	 * <h2>Behaviour with many edges</h2>
	 * <p>
	 * If the contour has 100 edges, this strategy will typically infer a small
	 * number of axes (often 1–4, depending on shape and thresholds). Each of the
	 * 100 edges is assigned to one inferred axis and then snapped to it during
	 * {@link #orient(int, LineSegment)}.
	 * </p>
	 *
	 * <h2>Fallback</h2>
	 * <p>
	 * If fewer than two meaningful axes are discovered (e.g., the contour is mostly
	 * one-directional), the strategy falls back to a single axis based on the
	 * longest edge and assigns all edges to it.
	 * </p>
	 */
	public static final class MultipleDirections implements ContourDirections {

		/**
		 * Configuration for {@link MultipleDirections}.
		 */
		public static final class Options {
			/**
			 * Maximum angle deviation (degrees) used during axis discovery.
			 *
			 * <p>
			 * An edge is considered compatible with an axis if it is either:
			 * <ul>
			 * <li>near-parallel: angle ≤ maximumAngleDeg</li>
			 * <li>near-orthogonal: angle ≥ 90 - maximumAngleDeg</li>
			 * </ul>
			 * </p>
			 */
			public final double maximumAngleDeg;
			/**
			 * Minimum edge length required for an edge to be used as an axis seed and to
			 * participate in axis estimation. Shorter edges are treated as weak/noisy
			 * evidence.
			 */
			public final double minimumLength;
			/**
			 * If true, each inferred axis is readjusted by the average residual angle of
			 * its assigned edges, improving fit to the input geometry.
			 */
			public final boolean adjustDirections;

			private Options(Builder b) {
				this.maximumAngleDeg = b.maximumAngleDeg;
				this.minimumLength = b.minimumLength;
				this.adjustDirections = b.adjustDirections;
			}

			public static Options defaults() {
				return builder().build();
			}

			public static Builder builder() {
				return new Builder();
			}

			public static final class Builder {
				private double maximumAngleDeg = 10.0;
				private double minimumLength = 3.0;
				private boolean adjustDirections = true;

				public Builder maximumAngleDeg(double v) {
					this.maximumAngleDeg = v;
					return this;
				}

				public Builder minimumLength(double v) {
					this.minimumLength = v;
					return this;
				}

				public Builder adjustDirections(boolean v) {
					this.adjustDirections = v;
					return this;
				}

				public Options build() {
					return new Options(this);
				}
			}
		}

		private final Options opt; // prototype options

		// initialized state:
		private final double[] axesDeg;
		private final int[] assigned;
		private final boolean initialized;

		/**
		 * Creates a prototype multiple-direction model using the provided options.
		 *
		 * <p>
		 * The actual axes and edge assignments are computed during
		 * {@link #init(Coordinate[], boolean)}.
		 * </p>
		 *
		 * @param opt axis discovery and adjustment options
		 */
		public MultipleDirections(Options opt) {
			this.opt = opt;
			this.axesDeg = null;
			this.assigned = null;
			this.initialized = false;
		}

		private MultipleDirections(Options opt, double[] axesDeg, int[] assigned) {
			this.opt = opt;
			this.axesDeg = axesDeg;
			this.assigned = assigned;
			this.initialized = true;
		}

		@Override
		public ContourDirections init(Coordinate[] coordinates, boolean closed) {
			Coordinate[] pts = closed ? sanitizeClosed(coordinates) : copy(coordinates);
			int edgeCount = closed ? pts.length : Math.max(0, pts.length - 1);

			double[] edgeOri = new double[edgeCount];
			double[] edgeLen = new double[edgeCount];
			boolean[] valid = new boolean[edgeCount];

			for (int e = 0; e < edgeCount; e++) {
				LineSegment s = edgeSegment(pts, closed, e);
				edgeOri[e] = Geometry2D.orientationDeg(s);
				edgeLen[e] = s.getLength();
				valid[e] = edgeLen[e] >= opt.minimumLength;
			}

			Integer[] idx = new Integer[edgeCount];
			for (int i = 0; i < edgeCount; i++)
				idx[i] = i;
			Arrays.sort(idx, (a, b) -> Double.compare(edgeLen[b], edgeLen[a]));

			int[] assigned = new int[edgeCount];
			Arrays.fill(assigned, -1);
			boolean[] used = new boolean[edgeCount];

			ArrayList<Double> axes = new ArrayList<>();
			int groupIndex = 0;

			while (true) {
				int seed = -1;
				for (int e : idx) {
					if (!used[e] && valid[e]) {
						seed = e;
						break;
					}
				}
				if (seed == -1)
					break;

				double axis = edgeOri[seed];
				axes.add(axis);
				assigned[seed] = groupIndex;
				used[seed] = true;

				for (int e = 0; e < edgeCount; e++) {
					if (e == seed || used[e] || !valid[e])
						continue;
					if (satisfiesAxisCondition(edgeOri[seed], edgeOri[e], opt.maximumAngleDeg)) {
						assigned[e] = groupIndex;
						used[e] = true;
					}
				}
				groupIndex++;
			}

			if (axes.size() <= 1) {
				// fallback to longest
				int longest = 0;
				for (int e = 1; e < edgeCount; e++)
					if (edgeLen[e] > edgeLen[longest])
						longest = e;
				double[] ax = new double[] { normalize180(edgeOri[longest]) };
				int[] as = new int[edgeCount];
				Arrays.fill(as, 0);
				return new MultipleDirections(opt, ax, as);
			}

			unifyAndCorrectAssignments(assigned, closed);

			double[] axesDeg = new double[axes.size()];
			for (int i = 0; i < axesDeg.length; i++)
				axesDeg[i] = normalize180(axes.get(i));

			if (opt.adjustDirections) {
				double[] sum = new double[axesDeg.length];
				double[] cnt = new double[axesDeg.length];

				for (int e = 0; e < edgeCount; e++) {
					if (!valid[e])
						continue;
					int a = assigned[e];
					double resid = Geometry2D.mod90AngleDifferenceDeg(edgeOri[e], axesDeg[a]);
					sum[a] += resid;
					cnt[a] += 1.0;
				}
				for (int a = 0; a < axesDeg.length; a++) {
					if (cnt[a] == 0.0)
						continue;
					axesDeg[a] = normalize180(axesDeg[a] + sum[a] / cnt[a]);
				}
			}

			return new MultipleDirections(opt, axesDeg, assigned);
		}

		@Override
		public void orient(int edgeIndex, LineSegment segment) {
			if (!initialized)
				throw new IllegalStateException("Not initialized; call init() first");
			int a = assigned[edgeIndex];
			double segOri = Geometry2D.orientationDeg(segment);
			double rot = Geometry2D.mod90AngleDifferenceDeg(segOri, axesDeg[a]);
			Geometry2D.rotateSegmentCCWAroundMidpointInPlace(segment, rot);
		}

		private static boolean satisfiesAxisCondition(double refOriDeg, double segOriDeg, double maxAngleDeg) {
			double a = angle0to90(refOriDeg, segOriDeg);
			return a <= maxAngleDeg || a >= (90.0 - maxAngleDeg);
		}
	}

	/**
	 * Principal-direction model where the user supplies one or more desired axis
	 * directions.
	 *
	 * <h2>How axes are interpreted</h2>
	 * <p>
	 * The {@code axesDeg} varargs defines a set of <em>global principal axes</em>
	 * (orientations) in degrees. Each value represents an axis direction normalized
	 * into the range {@code [0,180)}. For each axis, the model also implicitly
	 * accepts its orthogonal direction (axis + 90°), matching the CGAL concept of
	 * "parallel or orthogonal" fits.
	 * </p>
	 *
	 * <h2>Edge assignment (important)</h2>
	 * <p>
	 * Let the contour have {@code E} edges (segments), and you provide {@code M}
	 * axes in {@code axesDeg}. This does <b>not</b> mean you must provide {@code E}
	 * axes. Instead, each edge chooses an axis by classification:
	 * </p>
	 *
	 * <ol>
	 * <li>Compute the edge orientation in {@code [0,180)}.</li>
	 * <li>Scan axes in the order provided and select the <em>first</em> axis for
	 * which the edge is either:
	 * <ul>
	 * <li>near-parallel to the axis (angle ≤ {@code maxSnapAngleDeg}), or</li>
	 * <li>near-orthogonal to the axis (angle ≥ {@code 90 - maxSnapAngleDeg}).</li>
	 * </ul>
	 * </li>
	 * <li>If no axis matches, the edge is initially unassigned; unassigned edges
	 * are then filled-in by propagating the nearest assigned axis along the contour
	 * (a "unify along contour" pass), so that every edge ends up assigned.</li>
	 * </ol>
	 *
	 * <p>
	 * <b>Example:</b> If you pass two axes {@code axesDeg = {0, 45}} and the
	 * contour has 100 edges: each of the 100 edges will be assigned to axis 0° or
	 * 45° according to the test above. No per-edge axis array is required.
	 * </p>
	 *
	 * <h2>Orientation step</h2>
	 * <p>
	 * When {@link #orient(int, LineSegment)} is called for an edge, the edge is
	 * rotated around its midpoint by a signed "mod-90" difference so it becomes
	 * exactly aligned with the chosen axis (or its orthogonal), i.e., it snaps to
	 * the nearest of {axis, axis+90}.
	 * </p>
	 *
	 * <h2>Notes</h2>
	 * <ul>
	 * <li>If {@code axesDeg} is empty, this strategy cannot classify edges
	 * meaningfully. The implementation should either behave as a no-op or fall back
	 * to a default axis (implementation-defined).</li>
	 * <li>Because axes are checked in order, earlier axes have priority when an
	 * edge could match multiple axes.</li>
	 * </ul>
	 */
	public static final class UserDefinedDirections implements ContourDirections {

		private final double maxSnapAngleDeg;
		private final double[] axesDeg; // prototype options

		// initialized state:
		private final int[] assigned; // per edge
		private final boolean initialized;

		/**
		 * Creates a prototype user-defined axis model.
		 *
		 * @param maxSnapAngleDeg maximum angular deviation (degrees) for an edge to be
		 *                        considered parallel or orthogonal to an axis. Typical
		 *                        values: 3..10.
		 * @param axesDeg         one or more axis orientations in degrees. Each edge of
		 *                        the contour is classified against these axes (see
		 *                        class Javadoc). You do <b>not</b> need to supply one
		 *                        axis per edge.
		 */
		public UserDefinedDirections(double maxSnapAngleDeg, double... axesDeg) {
			this.maxSnapAngleDeg = maxSnapAngleDeg;
			this.axesDeg = new double[axesDeg.length];
			for (int i = 0; i < axesDeg.length; i++)
				this.axesDeg[i] = normalize180(axesDeg[i]);

			this.assigned = null;
			this.initialized = false;
		}

		private UserDefinedDirections(double maxSnapAngleDeg, double[] axesDeg, int[] assigned) {
			this.maxSnapAngleDeg = maxSnapAngleDeg;
			this.axesDeg = axesDeg;
			this.assigned = assigned;
			this.initialized = true;
		}

		@Override
		public ContourDirections init(Coordinate[] coordinates, boolean closed) {
			Coordinate[] pts = closed ? sanitizeClosed(coordinates) : copy(coordinates);
			int edgeCount = closed ? pts.length : Math.max(0, pts.length - 1);

			int[] assigned = new int[edgeCount];
			Arrays.fill(assigned, -1);

			for (int e = 0; e < edgeCount; e++) {
				LineSegment seg = edgeSegment(pts, closed, e);
				double segOri = Geometry2D.orientationDeg(seg);

				for (int d = 0; d < axesDeg.length; d++) {
					if (satisfiesAxisCondition(segOri, axesDeg[d], maxSnapAngleDeg)) {
						assigned[e] = d;
						break;
					}
				}
			}

			if (axesDeg.length == 0 || allUnassigned(assigned)) {
				Arrays.fill(assigned, 0); // fallback
			} else {
				unifyAndCorrectAssignments(assigned, closed);
			}

			// axesDeg is already normalized; safe to share (immutable)
			return new UserDefinedDirections(maxSnapAngleDeg, axesDeg, assigned);
		}

		@Override
		public void orient(int edgeIndex, LineSegment segment) {
			if (!initialized)
				throw new IllegalStateException("Not initialized; call init() first");
			int d = assigned[edgeIndex];
			double segOri = Geometry2D.orientationDeg(segment);
			double rot = Geometry2D.mod90AngleDifferenceDeg(segOri, axesDeg[d]);
			Geometry2D.rotateSegmentCCWAroundMidpointInPlace(segment, rot);
		}

		private static boolean satisfiesAxisCondition(double segOriDeg, double axisOriDeg, double maxAngleDeg) {
			double a = angle0to90(segOriDeg, axisOriDeg);
			return a <= maxAngleDeg || a >= (90.0 - maxAngleDeg);
		}

		private static boolean allUnassigned(int[] a) {
			for (int v : a)
				if (v != -1)
					return false;
			return true;
		}
	}

	/**
	 * Regularizes a geometry using {@link RegParameters#defaults()} and
	 * {@link LongestEdgeDirections}.
	 *
	 * @param input {@link LineString}, {@link LinearRing}, or {@link Polygon}
	 * @return a new geometry of the same runtime type
	 */
	public static Geometry regularize(Geometry input) {
		return regularize(input, RegParameters.defaults());
	}

	/**
	 * Regularizes a geometry.
	 *
	 * <p>
	 * Supported types:
	 * <ul>
	 * <li>{@link LineString} (open or closed)</li>
	 * <li>{@link LinearRing}</li>
	 * <li>{@link Polygon} (shell + holes each regularized independently)</li>
	 * </ul>
	 * </p>
	 *
	 * @param input  geometry to regularize
	 * @param params configuration
	 * @return a new geometry
	 */
	public static Geometry regularize(Geometry input, RegParameters params) {
		Objects.requireNonNull(input, "input");
		Objects.requireNonNull(params, "params");

		if (input instanceof Polygon p) {
			return regularize(p, params);
		}
		if (input instanceof LinearRing r) {
			return regularize(r, params);
		}
		if (input instanceof LineString ls) {
			return regularize(ls, params);
		}
		throw new IllegalArgumentException("Unsupported geometry type: " + input.getGeometryType());
	}

	/**
	 * Regularize a {@link LineString}. If the LineString is closed, it is treated
	 * as a ring.
	 */
	public static LineString regularize(LineString line, RegParameters params) {
		Objects.requireNonNull(line, "line");
		GeometryFactory gf = line.getFactory();

		boolean closed = line.isClosed();
		Coordinate[] coords = line.getCoordinates();
		Coordinate[] out = regularizeCoordinates(coords, closed, params);

		LineString result = gf.createLineString(out);
		result.setUserData(line.getUserData());
		return result;
	}

	/**
	 * Regularize a {@link LinearRing}. Always treated as closed.
	 */
	public static LinearRing regularize(LinearRing ring, RegParameters params) {
		Objects.requireNonNull(ring, "ring");
		GeometryFactory gf = ring.getFactory();

		Coordinate[] coords = ring.getCoordinates();
		Coordinate[] out = regularizeCoordinates(coords, true, params);

		LinearRing result = gf.createLinearRing(out);
		result.setUserData(ring.getUserData());
		return result;
	}

	/**
	 * Regularize a {@link Polygon}. Shell and each hole ring are processed
	 * independently.
	 */
	public static Polygon regularize(Polygon poly, RegParameters params) {
		Objects.requireNonNull(poly, "poly");
		GeometryFactory gf = poly.getFactory();

		LinearRing shell = poly.getExteriorRing();
		LinearRing shellR = regularize(shell, params);

		LinearRing[] holesR = null;
		int nh = poly.getNumInteriorRing();
		if (nh > 0) {
			holesR = new LinearRing[nh];
			for (int i = 0; i < nh; i++) {
				holesR[i] = regularize(poly.getInteriorRingN(i), params);
			}
		}

		Polygon result = gf.createPolygon(shellR, holesR);
		result.setUserData(poly.getUserData());
		return result;
	}

	/**
	 * Lowest-level entry point: regularize an ordered coordinate sequence.
	 *
	 * @param coordinates input vertices; for closed input may be with or without
	 *                    duplicate last==first
	 * @param closed      treat as ring if true; as open polyline if false
	 * @param directions  direction model (principal directions snapping)
	 * @param params      configuration
	 * @return regularized coordinates; if closed, guaranteed last==first
	 */
	public static Coordinate[] regularizeCoordinates(Coordinate[] coordinates, boolean closed, RegParameters params) {
		Objects.requireNonNull(coordinates, "coordinates");
		Objects.requireNonNull(params, "params");

		if (coordinates.length < (closed ? 4 : 2)) {
			// allow "already has closing point" case; sanitize below
			return copy(coordinates);
		}

		Coordinate[] pts = closed ? sanitizeClosed(coordinates) : copy(coordinates);
		if (pts.length < (closed ? 3 : 2)) {
			return closed ? ensureClosed(copy(pts)) : copy(pts);
		}

		// Keep originals if requested (open endpoints).
		final Coordinate openStart = (!closed && params.preserveOpenEndpoints) ? new Coordinate(pts[0]) : null;
		final Coordinate openEnd = (!closed && params.preserveOpenEndpoints) ? new Coordinate(pts[pts.length - 1]) : null;

		// 1) Build segments
		List<LineSegment> segs = buildSegments(pts, closed);

		// 2) Rotate toward principal directions
		ContourDirections dirs = params.directions.init(coordinates, closed);
		for (int i = 0; i < segs.size(); i++) {
			dirs.orient(i, segs.get(i));
		}

		// 3) Remove tiny edges
		segs = removeShort(segs, params.minEdgeLength);

		// CGAL-like early exits
		if (closed && segs.size() < 4) {
			return ensureClosed(copy(pts));
		}
		if (!closed && segs.size() < 1) {
			return copy(pts);
		}

		// 4) Merge consecutive collinear
		segs = mergeConsecutiveCollinear(segs, closed, params.parallelAngleThresholdDeg, params.maximumOffset);

		if (closed && segs.size() < 4) {
			return ensureClosed(copy(pts));
		}
		if (!closed && segs.size() < 1) {
			return copy(pts);
		}

		// 5) Insert orth edges if needed
		if (params.insertOrthogonalWhenParallel) {
			segs = insertOrthogonalEdgesWhenParallel(segs, closed, params.parallelAngleThresholdDeg);
		}

		if (closed && segs.size() < 4) {
			return ensureClosed(copy(pts));
		}
		if (!closed && segs.size() < 1) {
			return copy(pts);
		}

		// 6) Reconnect by intersections
		Coordinate[] out = reconnectToCoordinates(segs, closed);

		if (!closed && params.preserveOpenEndpoints) {
			out[0] = openStart;
			out[out.length - 1] = openEnd;
		}
		return out;
	}

	private static List<LineSegment> buildSegments(Coordinate[] pts, boolean closed) {
		int n = pts.length;
		int count = closed ? n : (n - 1);
		ArrayList<LineSegment> segs = new ArrayList<>(count);
		for (int i = 0; i < count; i++) {
			Coordinate a = pts[i];
			Coordinate b = pts[closed ? (i + 1) % n : (i + 1)];
			segs.add(new LineSegment(new Coordinate(a), new Coordinate(b)));
		}
		return segs;
	}

	private static List<LineSegment> removeShort(List<LineSegment> segs, double minLen) {
		double min2 = minLen * minLen;
		ArrayList<LineSegment> out = new ArrayList<>(segs.size());
		for (LineSegment s : segs) {
			double dx = s.p1.x - s.p0.x;
			double dy = s.p1.y - s.p0.y;
			if (dx * dx + dy * dy > min2) {
				out.add(s);
			}
		}
		return out;
	}

	private static List<LineSegment> mergeConsecutiveCollinear(List<LineSegment> segs, boolean closed, double angleThresholdDeg, double maxOffset) {
		if (segs.isEmpty()) {
			return segs;
		}

		List<List<LineSegment>> groups = new ArrayList<>();
		List<LineSegment> curr = new ArrayList<>();
		curr.add(segs.get(0));

		for (int i = 1; i < segs.size(); i++) {
			LineSegment ref = curr.get(0);
			LineSegment next = segs.get(i);

			if (Geometry2D.areParallel(ref, next, angleThresholdDeg) && Geometry2D.isCollinearEnough(ref, next, maxOffset)) {
				curr.add(next);
			} else {
				groups.add(curr);
				curr = new ArrayList<>();
				curr.add(next);
			}
		}
		groups.add(curr);

		// closed wrap merge
		if (closed && groups.size() >= 2) {
			List<LineSegment> first = groups.get(0);
			List<LineSegment> last = groups.get(groups.size() - 1);
			LineSegment a = last.get(0);
			LineSegment b = first.get(0);
			if (Geometry2D.areParallel(a, b, angleThresholdDeg) && Geometry2D.isCollinearEnough(a, b, maxOffset)) {
				ArrayList<LineSegment> merged = new ArrayList<>(last.size() + first.size());
				merged.addAll(last);
				merged.addAll(first);
				groups.set(0, merged);
				groups.remove(groups.size() - 1);
			}
		}
		List<LineSegment> out = new ArrayList<>(groups.size());
		for (List<LineSegment> g : groups) {
			out.add(Geometry2D.mergeCollinearGroup(g));
		}
		return out;
	}

	private static List<LineSegment> insertOrthogonalEdgesWhenParallel(List<LineSegment> segs, boolean closed, double angleThresholdDeg) {
		int n = segs.size();
		ArrayList<LineSegment> out = new ArrayList<>(n * 2);

		for (int i = 0; i < n; i++) {
			LineSegment si = segs.get(i);
			out.add(si);

			int j = closed ? (i + 1) % n : i + 1;
			if (!closed && j >= n) {
				break;
			}

			LineSegment sj = segs.get(j);
			if (Geometry2D.areParallel(si, sj, angleThresholdDeg)) {
				out.add(Geometry2D.createAverageOrth(si, sj));
			}
		}
		return out;
	}

	private static Coordinate[] reconnectToCoordinates(List<LineSegment> segs, boolean closed) {
		int n = segs.size();
		if (closed) {
			Coordinate[] out = new Coordinate[n + 1];
			for (int i = 0; i < n; i++) {
				int im = (i + n - 1) % n;
				out[i] = Geometry2D.infiniteLineIntersectionOrFallback(segs.get(im), segs.get(i), segs.get(i).p0);
			}
			out[n] = new Coordinate(out[0]);
			return out;
		} else {
			Coordinate[] out = new Coordinate[n + 1];
			out[0] = new Coordinate(segs.get(0).p0);
			for (int i = 1; i < n; i++) {
				out[i] = Geometry2D.infiniteLineIntersectionOrFallback(segs.get(i - 1), segs.get(i), segs.get(i).p0);
			}
			out[n] = new Coordinate(segs.get(n - 1).p1);
			return out;
		}
	}

	private static Coordinate[] sanitizeClosed(Coordinate[] input) {
		Coordinate[] c = copy(input);
		if (c.length >= 2 && c[0].equals2D(c[c.length - 1])) {
			return Arrays.copyOf(c, c.length - 1);
		}
		return c;
	}

	private static Coordinate[] ensureClosed(Coordinate[] ptsNoDupLast) {
		if (ptsNoDupLast.length == 0) {
			return ptsNoDupLast;
		}
		Coordinate[] out = Arrays.copyOf(ptsNoDupLast, ptsNoDupLast.length + 1);
		out[out.length - 1] = new Coordinate(out[0]);
		return out;
	}

	private static Coordinate[] copy(Coordinate[] in) {
		Coordinate[] out = new Coordinate[in.length];
		for (int i = 0; i < in.length; i++) {
			out[i] = new Coordinate(in[i]);
		}
		return out;
	}

	private static double normalize180(double ang) {
		ang %= 180.0;
		if (ang < 0)
			ang += 180.0;
		return ang;
	}

	/** Absolute angle between two orientations mapped into [0,90]. */
	private static double angle0to90(double aDeg, double bDeg) {
		double diff = Math.abs(aDeg - bDeg);
		diff = Math.min(diff, 180.0 - diff); // now in [0,90]
		return diff;
	}

	private static LineSegment edgeSegment(Coordinate[] pts, boolean closed, int edgeIndex) {
		int n = pts.length;
		int i = edgeIndex;
		int j = closed ? (i + 1) % n : (i + 1);
		return new LineSegment(pts[i], pts[j]);
	}

	/**
	 * CGAL-like "unify along contour" then "correct directions" pass. Operates
	 * in-place on {@code assigned}, where -1 indicates unassigned.
	 */
	private static void unifyAndCorrectAssignments(int[] assigned, boolean closed) {
		if (assigned.length == 0)
			return;

		if (closed) {
			unifyClosed(assigned);
			correctClosed(assigned);
		} else {
			unifyOpen(assigned);
			correctOpen(assigned);
		}
	}

	private static void unifyClosed(int[] assigned) {
		int n = assigned.length;
		for (int i = 0; i < n; i++) {
			if (assigned[i] != -1)
				continue;

			int im = (i + n - 1) % n;
			int ip = (i + 1) % n;

			boolean stop = false;
			int steps = 0;
			while (!stop && steps < n) {
				if (assigned[im] != -1) {
					assigned[i] = assigned[im];
					break;
				}
				if (assigned[ip] != -1) {
					assigned[i] = assigned[ip];
					break;
				}

				im = (im + n - 1) % n;
				ip = (ip + 1) % n;
				if (im == i || ip == i)
					stop = true;
				steps++;
			}
			if (assigned[i] == -1)
				assigned[i] = 0;
		}
	}

	private static void correctClosed(int[] assigned) {
		int n = assigned.length;
		int[] clean = new int[n];
		for (int i = 0; i < n; i++) {
			int im = (i + n - 1) % n;
			int ip = (i + 1) % n;
			int dm = assigned[im];
			int di = assigned[i];
			int dp = assigned[ip];
			clean[i] = (dm != -1 && dm == dp && di != dm) ? dm : di;
		}
		System.arraycopy(clean, 0, assigned, 0, n);
	}

	private static void unifyOpen(int[] assigned) {
		int n = assigned.length;
		for (int i = 0; i < n; i++) {
			if (assigned[i] != -1)
				continue;

			int im = (i > 0) ? i - 1 : -1;
			int ip = (i < n - 1) ? i + 1 : -1;

			boolean stop = false;
			int steps = 0;
			while (steps < n) {
				if (im != -1 && assigned[im] != -1) {
					assigned[i] = assigned[im];
					break;
				}
				if (ip != -1 && assigned[ip] != -1) {
					assigned[i] = assigned[ip];
					break;
				}

				if (stop)
					break;
				if (im > 0)
					im--;
				if (ip != -1 && ip < n - 1)
					ip++;

				if (im == 0 || ip == n - 1)
					stop = true;
				steps++;
			}
			if (assigned[i] == -1)
				assigned[i] = 0;
		}
	}

	private static void correctOpen(int[] assigned) {
		int n = assigned.length;
		if (n == 1)
			return;

		int[] clean = new int[n];

		// first
		clean[0] = (assigned[0] != assigned[1]) ? assigned[1] : assigned[0];

		// middle
		for (int i = 1; i < n - 1; i++) {
			int dm = assigned[i - 1];
			int di = assigned[i];
			int dp = assigned[i + 1];
			clean[i] = (dm != -1 && dm == dp && di != dm) ? dm : di;
		}

		// last
		clean[n - 1] = (assigned[n - 1] != assigned[n - 2]) ? assigned[n - 2] : assigned[n - 1];

		System.arraycopy(clean, 0, assigned, 0, n);
	}

	static final class Geometry2D {
		private Geometry2D() {
		}

		/**
		 * Orientation in degrees, normalized to [0,180) with CGAL-like sign convention.
		 */
		static double orientationDeg(LineSegment s) {
			double dx = s.p1.x - s.p0.x;
			double dy = s.p1.y - s.p0.y;
			if (dy < 0 || (dy == 0 && dx < 0)) {
				dx = -dx;
				dy = -dy;
			}
			double ang = Math.toDegrees(Math.atan2(dy, dx));
			if (ang < 0) {
				ang += 180.0;
			}
			return ang;
		}

		/** Signed mod-90 difference (CGAL internal::mod90_angle_difference_2). */
		static double mod90AngleDifferenceDeg(double angleI, double angleJ) {
			double diff = angleI - angleJ;
			int diff90 = (int) Math.floor(diff / 90.0);
			double toLower = 90.0 * (diff90 + 0.0) - diff;
			double toUpper = 90.0 * (diff90 + 1.0) - diff;
			return (Math.abs(toLower) < Math.abs(toUpper)) ? toLower : toUpper;
		}

		static void rotateSegmentCCWAroundMidpointInPlace(LineSegment s, double angleDeg) {
			double rad = Math.toRadians(angleDeg);
			double sin = Math.sin(rad);
			double cos = Math.cos(rad);

			double mx = 0.5 * (s.p0.x + s.p1.x);
			double my = 0.5 * (s.p0.y + s.p1.y);

			rotatePointCCWInPlace(s.p0, mx, my, cos, sin);
			rotatePointCCWInPlace(s.p1, mx, my, cos, sin);
		}

		private static void rotatePointCCWInPlace(Coordinate p, double cx, double cy, double cos, double sin) {
			double x = p.x - cx;
			double y = p.y - cy;
			double rx = x * cos - y * sin;
			double ry = y * cos + x * sin;
			p.x = cx + rx;
			p.y = cy + ry;
		}

		static boolean areParallel(LineSegment a, LineSegment b, double angleThresholdDeg) {
			double oa = orientationDeg(a);
			double ob = orientationDeg(b);
			double diff = Math.abs(oa - ob);
			diff = Math.min(diff, 180.0 - diff);
			diff = Math.min(diff, 90.0 - Math.abs(90.0 - diff)); // map to [0,90]
			return diff <= angleThresholdDeg;
		}

		/**
		 * CGAL-like test: distance from midpoint of s to projection on ref supporting
		 * line.
		 */
		static boolean isCollinearEnough(LineSegment ref, LineSegment s, double maxOffset) {
			double mx = 0.5 * (s.p0.x + s.p1.x);
			double my = 0.5 * (s.p0.y + s.p1.y);
			Coordinate proj = projectPointToLine(mx, my, ref.p0, ref.p1);
			double dx = mx - proj.x;
			double dy = my - proj.y;
			return (dx * dx + dy * dy) <= maxOffset * maxOffset;
		}

		static Coordinate projectPointToLine(double px, double py, Coordinate a, Coordinate b) {
			double ux = b.x - a.x;
			double uy = b.y - a.y;
			double denom = ux * ux + uy * uy;
			if (denom == 0.0) {
				return new Coordinate(a);
			}
			double t = ((px - a.x) * ux + (py - a.y) * uy) / denom;
			return new Coordinate(a.x + t * ux, a.y + t * uy);
		}

		/**
		 * Intersection of infinite supporting lines using JTS DD arithmetic. Falls back
		 * to {@code fallback} if parallel/collinear/undefined.
		 */
		static Coordinate infiniteLineIntersectionOrFallback(LineSegment s1, LineSegment s2, Coordinate fallback) {
			Coordinate p = CGAlgorithmsDD.intersection(s1.p0, s1.p1, s2.p0, s2.p1);
			return (p != null) ? p : new Coordinate(fallback);
		}

		/** Port of CGAL Contour_base_2::create_average_orth. */
		static LineSegment createAverageOrth(LineSegment segmentI, LineSegment segmentJ) {
			Coordinate p = projectPointToLine(segmentJ.p0.x, segmentJ.p0.y, segmentI.p0, segmentI.p1);
			Coordinate source = midpoint(p, segmentI.p1);

			Coordinate q = projectPointToLine(segmentI.p1.x, segmentI.p1.y, segmentJ.p0, segmentJ.p1);
			Coordinate target = midpoint(q, segmentJ.p0);

			return new LineSegment(source, target);
		}

		static Coordinate midpoint(Coordinate a, Coordinate b) {
			return new Coordinate(0.5 * (a.x + b.x), 0.5 * (a.y + b.y));
		}

		/**
		 * Merge a consecutive collinear group into one best-fit segment: - choose
		 * longest as reference direction - compute weighted (squared length) average
		 * normal offset of midpoints - project all endpoints onto shifted line and take
		 * min/max along reference axis
		 */
		static LineSegment mergeCollinearGroup(List<LineSegment> group) {
			if (group.size() == 1) {
				LineSegment s = group.get(0);
				return new LineSegment(new Coordinate(s.p0), new Coordinate(s.p1));
			}

			LineSegment ref = group.get(0);
			double best = -1;
			for (LineSegment s : group) {
				double len = s.getLength();
				if (len > best) {
					best = len;
					ref = s;
				}
			}

			double ux = ref.p1.x - ref.p0.x;
			double uy = ref.p1.y - ref.p0.y;
			double un = ref.p1.distance(ref.p0);
			if (un == 0.0) {
				return new LineSegment(new Coordinate(ref.p0), new Coordinate(ref.p1));
			}
			ux /= un;
			uy /= un;

			double nx = -uy;
			double ny = ux;

			double r0x = 0.5 * (ref.p0.x + ref.p1.x);
			double r0y = 0.5 * (ref.p0.y + ref.p1.y);

			double wsum = 0.0;
			double dsum = 0.0;
			for (LineSegment s : group) {
				double w = s.getLength();
				w = w * w; // squared length weights
				double mx = 0.5 * (s.p0.x + s.p1.x);
				double my = 0.5 * (s.p0.y + s.p1.y);
				double dx = mx - r0x;
				double dy = my - r0y;
				double signed = dx * nx + dy * ny;
				wsum += w;
				dsum += w * signed;
			}
			double d = (wsum == 0.0) ? 0.0 : (dsum / wsum);

			double rx = r0x + d * nx;
			double ry = r0y + d * ny;

			double tMin = Double.POSITIVE_INFINITY;
			double tMax = Double.NEGATIVE_INFINITY;
			for (LineSegment s : group) {
				tMin = Math.min(tMin, projT(rx, ry, ux, uy, s.p0));
				tMax = Math.max(tMax, projT(rx, ry, ux, uy, s.p0));
				tMin = Math.min(tMin, projT(rx, ry, ux, uy, s.p1));
				tMax = Math.max(tMax, projT(rx, ry, ux, uy, s.p1));
			}

			Coordinate a = new Coordinate(rx + tMin * ux, ry + tMin * uy);
			Coordinate b = new Coordinate(rx + tMax * ux, ry + tMax * uy);
			return new LineSegment(a, b);
		}

		private static double projT(double rx, double ry, double ux, double uy, Coordinate p) {
			return (p.x - rx) * ux + (p.y - ry) * uy;
		}
	}
}