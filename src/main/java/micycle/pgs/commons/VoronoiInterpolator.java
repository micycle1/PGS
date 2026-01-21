package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.locationtech.jts.algorithm.Distance;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.geom.util.GeometryCombiner;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;

import com.github.micycle1.geoblitz.HilbertParallelPolygonUnion;

import net.jafama.FastMath;

/**
 * Implements the Voronoi morph described in <em>Abstract morphing using the
 * Hausdorff distance and Voronoi diagrams</em> (de Kogel, van Kreveld,
 * Vermeulen).
 * <p>
 * Given two planar shapes {@code A} and {@code B} (as JTS {@link Geometry}),
 * the morph is evaluated for {@code alpha ∈ [0,1]} by:
 * <ul>
 * <li>computing the overlap {@code A ∩ B} (kept fixed throughout the
 * morph),</li>
 * <li>partitioning the non-overlapping parts {@code A \ B} and {@code B \ A} by
 * Voronoi cells of sampled boundary sites of the opposite shape,</li>
 * <li>moving each partition piece toward its associated closest-site on the
 * other shape, using a piecewise affine transform:
 * <ul>
 * <li>for vertex-sites: uniform scaling toward the vertex,</li>
 * <li>for edge-sites: scaling perpendicular to the edge’s supporting line,</li>
 * </ul>
 * </li>
 * <li>unioning (optionally) the transformed pieces with the fixed overlap.</li>
 * </ul>
 * <p>
 * This class separates an expensive, one-time preprocessing step
 * ({@link #prepareVoronoiPartition(Geometry, Geometry, double, double)
 * prepareVoronoiPartition()}) from the evaluation step
 * ({@link #interpolateVoronoi(VoronoiPartition, double, boolean)
 * interpolateVoronoi()}), so the same partition can be reused to render many
 * frames for different {@code alpha}.
 * 
 * @author Michael Carleton
 */
public final class VoronoiInterpolator {

	private VoronoiInterpolator() {
	}

	/**
	 * Cached preprocessing result for Voronoi morph evaluation.
	 * <p>
	 * Contains:
	 * <ul>
	 * <li>{@link #overlap}: the fixed overlap {@code A ∩ B} (unchanged throughout
	 * the morph),</li>
	 * <li>{@link #aPieces}: partitioned pieces of {@code A \ B} tagged with
	 * closest-site information on {@code B},</li>
	 * <li>{@link #bPieces}: partitioned pieces of {@code B \ A} tagged with
	 * closest-site information on {@code A}.</li>
	 * </ul>
	 * The expensive work is the geometric partitioning/overlay done during
	 * preparation. The interpolation step only applies affine transforms to these
	 * cached pieces and optionally unions them.
	 */
	public static final class VoronoiPartition {
		/** Geometry factory used to create output geometries. */
		public final GeometryFactory gf;

		/** Fixed overlap region {@code A ∩ B} (unchanged throughout the morph). */
		public final Geometry overlap;

		/**
		 * Partition pieces from {@code A \ B}, each tagged with closest-site info on
		 * {@code B}.
		 */
		public final List<Piece> aPieces;

		/**
		 * Partition pieces from {@code B \ A}, each tagged with closest-site info on
		 * {@code A}.
		 */
		public final List<Piece> bPieces;

		/**
		 * Creates a partition container. Intended to be produced by
		 * {@link VoronoiInterpolator#prepareVoronoiPartition(Geometry, Geometry, double, double)
		 * prepareVoronoiPartition()}.
		 *
		 * @param gf      geometry factory for outputs
		 * @param overlap fixed overlap {@code A ∩ B}
		 * @param aPieces partition pieces from {@code A \ B}
		 * @param bPieces partition pieces from {@code B \ A}
		 */
		private VoronoiPartition(GeometryFactory gf, Geometry overlap, List<Piece> aPieces, List<Piece> bPieces) {
			this.gf = gf;
			this.overlap = overlap;
			this.aPieces = aPieces;
			this.bPieces = bPieces;
		}
	}

	/**
	 * A single partition piece together with the closest-site information that
	 * determines how it moves during the morph.
	 * <p>
	 * {@link #geom} is the static (alpha-independent) geometry of the piece;
	 * {@link #site} defines the affine transform applied at evaluation time.
	 */
	private static final class Piece {
		/** Partitioned polygonal piece (independent of {@code alpha}). */
		public final Geometry geom;

		/**
		 * Closest-site information (vertex/edge) used to derive the affine transform.
		 */
		public final SiteInfo site;

		/**
		 * Creates a piece record.
		 *
		 * @param geom piece geometry (typically a {@link Polygon})
		 * @param site closest-site info describing the motion model for this piece
		 */
		private Piece(Geometry geom, SiteInfo site) {
			this.geom = geom;
			this.site = site;
		}
	}

	/**
	 * Performs the expensive, one-time preprocessing step: fixes inputs, computes
	 * the fixed overlap {@code A ∩ B}, and partitions the non-overlapping parts
	 * using Voronoi cells of sampled boundary sites.
	 * <p>
	 * The returned {@link VoronoiPartition} is intended to be reused to evaluate
	 * many intermediate shapes for different {@code alpha} values via
	 * {@link #interpolateVoronoi(VoronoiPartition, double, boolean)
	 * interpolateVoronoi()}.
	 *
	 * @param a                input shape {@code A}
	 * @param b                input shape {@code B}
	 * @param maxSegmentLength maximum segment length used to densify the boundary
	 *                         when sampling Voronoi sites; {@code <= 0} disables
	 *                         densification
	 * @param clipExpand       non-negative expansion applied to the combined
	 *                         envelope of {@code A} and {@code B} to form the
	 *                         Voronoi clip envelope
	 * @return a reusable partition containing the fixed overlap and tagged
	 *         partition pieces
	 * @throws NullPointerException if {@code a} or {@code b} is null
	 */
	public static VoronoiPartition prepareVoronoiPartition(Geometry a, Geometry b, double maxSegmentLength, double clipExpand) {
		Objects.requireNonNull(a, "a");
		Objects.requireNonNull(b, "b");

		Geometry aFix = GeometryFixer.fix(a);
		Geometry bFix = GeometryFixer.fix(b);
		GeometryFactory gf = aFix.getFactory();

		Geometry overlap = aFix.intersection(bFix);
		Geometry aOutsideB = aFix.difference(bFix);
		Geometry bOutsideA = bFix.difference(aFix);

		Envelope clip = new Envelope(aFix.getEnvelopeInternal());
		clip.expandToInclude(bFix.getEnvelopeInternal());
		clip.expandBy(Math.max(clipExpand, 0.0));

		// Precompute partition pieces once:
		List<Piece> aPieces = partitionPieces(aOutsideB, bFix, maxSegmentLength, clip, gf);
		List<Piece> bPieces = partitionPieces(bOutsideA, aFix, maxSegmentLength, clip, gf);

		return new VoronoiPartition(gf, overlap, aPieces, bPieces);
	}

	/**
	 * Fast evaluation step: computes the intermediate shape for the given
	 * {@code alpha} by applying affine transforms to the cached partition pieces
	 * and combining them with the fixed overlap.
	 * <p>
	 * Pieces originating from {@code A \ B} are transformed with fraction
	 * {@code alpha}; pieces originating from {@code B \ A} are transformed with
	 * fraction {@code 1 - alpha}.
	 *
	 * @param p       cached preprocessing result returned by
	 *                {@link #prepareVoronoiPartition(Geometry, Geometry, double, double)}
	 * @param alpha   morph parameter in {@code [0, 1]}
	 * @param doUnion if {@code true}, unions all parts into a clean area geometry
	 *                (slowest); if {@code false}, combines without union (fast) and
	 *                may leave overlaps/seams
	 * @return the interpolated geometry for {@code alpha}
	 * @throws NullPointerException     if {@code p} is null
	 * @throws IllegalArgumentException if {@code alpha} is outside {@code [0,1]}
	 */
	public static Geometry interpolateVoronoi(VoronoiPartition p, double alpha, boolean doUnion) {
		Objects.requireNonNull(p, "partition");
		if (alpha < 0.0 || alpha > 1.0) {
			throw new IllegalArgumentException("alpha must be in [0,1]");
		}

		if (alpha == 0.0) {
			// Reconstruct something close to A: overlap + aPieces at fraction 0 + bPieces
			// at fraction 1
			// If you need exact A, just keep A separately. This is meant for animation.
		}
		if (alpha == 1.0) {
			// same note as alpha==0
		}

		double fracA = alpha;
		double fracB = 1.0 - alpha;

		// Transform pieces in parallel (cheap compared to overlay ops)
		List<Geometry> movedA = p.aPieces.isEmpty() ? List.of()
				: p.aPieces.parallelStream().map(pc -> buildTransformForSite(pc.site, fracA).transform(pc.geom)).filter(g -> !g.isEmpty()).toList();

		List<Geometry> movedB = p.bPieces.isEmpty() ? List.of()
				: p.bPieces.parallelStream().map(pc -> buildTransformForSite(pc.site, fracB).transform(pc.geom)).filter(g -> !g.isEmpty()).toList();

		ArrayList<Geometry> out = new ArrayList<>(1 + movedA.size() + movedB.size());
		if (p.overlap != null && !p.overlap.isEmpty()) {
			out.add(p.overlap);
		}
		out.addAll(movedA);
		out.addAll(movedB);

		if (out.isEmpty()) {
			return p.gf.createGeometryCollection();
		}

		if (!doUnion) {
			return GeometryCombiner.combine(out);
		}
		Geometry u = HilbertParallelPolygonUnion.union(out);
		return u;
	}

	private enum SiteType {
		VERTEX, EDGE
	}

	private static final class SiteInfo {
		final Coordinate site;
		final SiteType type;

		// EDGE only
		final Coordinate edgeOrigin;
		final double edgeAngle;

		private SiteInfo(Coordinate site, SiteType type, Coordinate edgeOrigin, double edgeAngle) {
			this.site = site;
			this.type = type;
			this.edgeOrigin = edgeOrigin;
			this.edgeAngle = edgeAngle;
		}

		static SiteInfo vertex(Coordinate c) {
			return new SiteInfo(new Coordinate(c), SiteType.VERTEX, null, 0.0);
		}

		static SiteInfo edge(Coordinate c, Coordinate origin, double angle) {
			return new SiteInfo(new Coordinate(c), SiteType.EDGE, new Coordinate(origin), angle);
		}
	}

	/**
	 * One-time: partition `moving` by Voronoi cells of `other`, returning
	 * (piece,site) pairs. This is where the expensive `moving.intersection(cell)`
	 * happens.
	 */
	private static List<Piece> partitionPieces(Geometry moving, Geometry other, double maxSegmentLength, Envelope clip, GeometryFactory gf) {
		if (moving.isEmpty()) {
			return List.of();
		}

		List<SiteInfo> sites = sampleBoundarySites(other, maxSegmentLength);
		if (sites.isEmpty()) {
			return List.of();
		}

		Geometry cells = buildVoronoiCells(sites, clip, gf);

		Map<String, SiteInfo> siteByKey = new HashMap<>(sites.size() * 2);
		for (SiteInfo si : sites) {
			siteByKey.put(key(si.site), si);
		}

		int n = cells.getNumGeometries();
		Envelope movEnv = moving.getEnvelopeInternal();

		return IntStream.range(0, n).parallel().boxed().flatMap(i -> {
			Geometry cell = cells.getGeometryN(i);
			if (!(cell instanceof Polygon) || cell.isEmpty()) {
				return Stream.empty();
			}
			if (!movEnv.intersects(cell.getEnvelopeInternal())) {
				return Stream.empty();
			}

			Object ud = cell.getUserData();
			if (!(ud instanceof Coordinate)) {
				return Stream.empty();
			}

			SiteInfo si = siteByKey.get(key((Coordinate) ud));
			if (si == null) {
				return Stream.empty();
			}

			Geometry piece = moving.intersection(cell);
			if (piece.isEmpty()) {
				return Stream.empty();
			}

			// Fast path: already a single polygon
			if (piece instanceof Polygon p) {
				return p.isEmpty() ? Stream.empty() : Stream.of(new Piece(p, si));
			}

			// Slower path only when needed: explode collections/multipolygons to polygons
			@SuppressWarnings("unchecked")
			List<Polygon> polys = PolygonExtracter.getPolygons(piece);
			if (polys.isEmpty()) {
				return Stream.empty(); // line/point/etc
			}

			return polys.stream().filter(p -> !p.isEmpty()).map(p -> new Piece(p, si));
		}).toList();
	}

	/** Same as before */
	private static AffineTransformation buildTransformForSite(SiteInfo si, double fraction) {
		double s = 1.0 - fraction;

		if (si.type == SiteType.VERTEX) {
			return AffineTransformation.scaleInstance(s, s, si.site.x, si.site.y);
		}

		double angle = si.edgeAngle;
		Coordinate o = si.edgeOrigin;

		AffineTransformation t = new AffineTransformation();
		t.translate(-o.x, -o.y);
		t.rotate(-angle);
		t.scale(1.0, s);
		t.rotate(angle);
		t.translate(o.x, o.y);
		return t;
	}

	private static Geometry buildVoronoiCells(List<SiteInfo> sites, Envelope clip, GeometryFactory gf) {
		Coordinate[] coords = new Coordinate[sites.size()];
		for (int i = 0; i < sites.size(); i++) {
			coords[i] = sites.get(i).site;
		}

		Geometry mp = gf.createMultiPointFromCoords(coords);

		VoronoiDiagramBuilder vdb = new VoronoiDiagramBuilder();
		vdb.setSites(mp);
		vdb.setClipEnvelope(clip);

		// Returns a GeometryCollection of polygons; each polygon has
		// userData=Coordinate(site)
		return vdb.getDiagram(gf);
	}

	private static List<SiteInfo> sampleBoundarySites(Geometry g, double maxSegmentLength) {
		Geometry boundary = g.getBoundary();
		if (boundary.isEmpty()) {
			return List.of();
		}

		Geometry boundaryToSample = boundary;
		if (maxSegmentLength > 0.0) {
			boundaryToSample = Densifier.densify(boundary, maxSegmentLength);
		}

		// Extract segments from the ORIGINAL boundary (for EDGE direction
		// classification).
		List<LineSegment> segs = extractSegments(boundary);

		// Extract vertices of the ORIGINAL boundary.
		Set<String> vertexKeys = new HashSet<>();
		for (Coordinate c : boundary.getCoordinates()) {
			vertexKeys.add(key(c));
		}

		// Create sites from sampled coordinates.
		// If near a vertex -> VERTEX site, else -> EDGE site using closest boundary
		// segment direction.
		Set<String> seen = new HashSet<>();
		List<SiteInfo> sites = new ArrayList<>();

		for (Coordinate c : boundaryToSample.getCoordinates()) {
			String k = key(c);
			if (!seen.add(k)) {
				continue;
			}

			if (vertexKeys.contains(k)) {
				sites.add(SiteInfo.vertex(c));
			} else {
				SegmentMatch m = closestSegmentMatch(c, segs);
				if (m != null) {
					double ang = FastMath.atan2(m.seg.p1.y - m.seg.p0.y, m.seg.p1.x - m.seg.p0.x);
					sites.add(SiteInfo.edge(c, m.seg.p0, ang));
				} else {
					// Fallback: treat as vertex site
					sites.add(SiteInfo.vertex(c));
				}
			}
		}

		return sites;
	}

	private static final class SegmentMatch {
		final LineSegment seg;
		final double dist2;

		SegmentMatch(LineSegment seg, double dist2) {
			this.seg = seg;
			this.dist2 = dist2;
		}
	}

	private static SegmentMatch closestSegmentMatch(Coordinate p, List<LineSegment> segs) {
		LineSegment best = null;
		double bestD2 = Double.POSITIVE_INFINITY;

		for (LineSegment s : segs) {
			double d2 = Distance.pointToSegmentSq(p, s.p0, s.p1);
			if (d2 < bestD2) {
				bestD2 = d2;
				best = s;
			}
		}
		return best == null ? null : new SegmentMatch(best, bestD2);
	}

	private static List<LineSegment> extractSegments(Geometry boundary) {
		List<LineSegment> segs = new ArrayList<>();
		for (int i = 0; i < boundary.getNumGeometries(); i++) {
			Geometry gi = boundary.getGeometryN(i);
			if (!(gi instanceof LineString)) {
				continue;
			}
			Coordinate[] cs = ((LineString) gi).getCoordinates();
			for (int j = 0; j + 1 < cs.length; j++) {
				if (!cs[j].equals2D(cs[j + 1])) {
					segs.add(new LineSegment(cs[j], cs[j + 1]));
				}
			}
		}
		return segs;
	}

	// coordinate key with rounding to make Coordinate usable as map key
	private static String key(Coordinate c) {
		double eps = 1e-9;
		long ix = Math.round(c.x / eps);
		long iy = Math.round(c.y / eps);
		return ix + ":" + iy;
	}
}