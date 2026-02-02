package micycle.pgs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SplittableRandom;
import java.util.stream.Collectors;

import org.jgrapht.alg.interfaces.MatchingAlgorithm;
import org.jgrapht.alg.matching.blossom.v5.KolmogorovWeightedMatching;
import org.jgrapht.alg.matching.blossom.v5.KolmogorovWeightedPerfectMatching;
import org.jgrapht.alg.matching.blossom.v5.ObjectiveSense;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.util.LinearComponentExtracter;
import org.locationtech.jts.noding.MCIndexSegmentSetMutualIntersector;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.SegmentIntersector;
import org.locationtech.jts.noding.SegmentSetMutualIntersector;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.SegmentStringUtil;
import org.tinfour.common.IIncrementalTin;

import com.github.micycle1.geoblitz.IndexedLengthIndexedLine;

import micycle.pgs.color.Colors;
import micycle.pgs.commons.FastAtan2;
import micycle.pgs.commons.Nullable;
import micycle.pgs.commons.PEdge;
import net.jafama.FastMath;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Generation of random sets of <b>non-intersecting</b> line segments (and
 * associated functions).
 * <p>
 * Methods in this class output segments as collections of
 * {@link micycle.pgs.commons.PEdge PEdges}; such collections can be converted
 * into LINES PShapes with {@link #toPShape(Collection) toPShape(edges)}.
 * 
 * @author Michael Carleton
 * @since 1.3.0
 */
public class PGS_SegmentSet {

	private PGS_SegmentSet() {
	}

	/**
	 * Generates N non-intersecting segments via intersection and noding.
	 * <p>
	 * The segments are generated within a bounding box anchored at (0, 0) having
	 * the width and height specified.
	 * <p>
	 * Algorithm:
	 * <ul>
	 * <li>Generate a set of N random segments (will have many intersections)</li>
	 * <li>Node the random segment set and split segments at intersections</li>
	 * <li>Shrink split segments by ~30% (to increase distance between each
	 * other)</li>
	 * <li>The split segments will be very numerous; sort by length and return the
	 * longest N segments</li>
	 * </ul>
	 * 
	 * @param width  width of the bounds in which to generate segments; segment x
	 *               coordinates will not be greater than this value
	 * @param height height of the bounds in which to generate segments; segment y
	 *               coordinates will not be greater than this value
	 * @param n      number of segments to generate
	 * @param seed   number used to initialize the underlying pseudorandom number
	 *               generator
	 * @return set of N random non-intersecting line segments
	 */
	public static List<PEdge> nodedSegments(double width, double height, int n, long seed) {
		double maxLength = 80; // suitable value
		int effectiveN = n > 500 ? 2 * n / 3 : n;
		final SplittableRandom random = new SplittableRandom(seed);
		final List<SegmentString> segments = new ArrayList<>(effectiveN);
		for (int i = 0; i < effectiveN; i++) {
			final double length = maxLength;
			final Coordinate a = new Coordinate(random.nextDouble() * (width - length * 2) + length, random.nextDouble() * (height - length * 2) + length);
			final double theta = random.nextDouble(Math.PI * 2);
			final double x = FastMath.cosQuick(theta) * length;
			final double y = FastMath.sinQuick(theta) * length;
			final Coordinate b = new Coordinate(a.x + x, a.y + y);

			segments.add(new NodedSegmentString(new Coordinate[] { a, b }, null));
		}

		final List<LineSegment> ls = new ArrayList<>(segments.size());

		PGS.nodeSegmentStrings(segments).forEach(s -> {
			LineSegment lineSegment = new LineSegment(s.getCoordinate(0), s.getCoordinate(1));
			ls.add(lineSegment);
		});

		final double fract = 0.15;
		ls.sort((a, b) -> Double.compare(b.getLength(), a.getLength()));

		List<PEdge> edges = new ArrayList<>(segments.size());
		ls.subList(0, Math.min(n, ls.size())).forEach(lineSegment -> {
			Coordinate a = lineSegment.pointAlong(fract);
			Coordinate b = lineSegment.pointAlong(1 - fract);
			PEdge e = new PEdge(a.x, a.y, b.x, b.y);
			edges.add(e);
		});

		return edges;
	}

	/**
	 * Generates N non-intersecting segments via a <i>Perfect matching</i> algorithm
	 * applied to a triangulation populated with random points.
	 * <p>
	 * The segments are generated within a bounding box anchored at (0, 0) having
	 * the width and height specified.
	 * <p>
	 * The <code>graphMatchedSegments</code> methods are arguably the best
	 * approaches for random segment set generation.
	 * 
	 * @param width  width of the bounds in which to generate segments; segment x
	 *               coordinates will not be greater than this value
	 * @param height height of the bounds in which to generate segments; segment y
	 *               coordinates will not be greater than this value
	 * @param n      number of segments to generate
	 * @param seed   number used to initialize the underlying pseudorandom number
	 *               generator
	 * @return set of N random non-intersecting line segments
	 */
	public static List<PEdge> graphMatchedSegments(double width, double height, int n, long seed) {
		// pjLDS chosen as a comprimise between random and poisson
		n *= 2; // since #segments = #vertices/2
		return graphMatchedSegments(PGS_PointSet.plasticJitteredLDS(0, 0, width, height, n, seed));
	}

	/**
	 * Generates non-intersecting segments via a <i>Perfect matching</i> algorithm
	 * applied to a triangulation populated with the given points.
	 * <p>
	 * Generates N/2 segments, where N is the number of points in the input. If the
	 * number of points is odd, the last point is discarded.
	 * <p>
	 * The <code>graphMatchedSegments</code> methods are arguably the best
	 * approaches for random segment set generation.
	 * 
	 * @param points point set from which to compute segments
	 * @return set of non-intersecting line segments
	 */
	public static List<PEdge> graphMatchedSegments(List<PVector> points) {
		if ((points.size() & 1) == 1) { // coerce size to be even
			points = points.subList(0, points.size() - 1);
		}
		return graphMatchedSegments(PGS_Triangulation.delaunayTriangulationMesh(points));
	}

	/**
	 * Generates non-intersecting segments via a <i>Perfect matching</i> algorithm
	 * applied to the given triangulation.
	 * <p>
	 * Generates N/2 segments, where N is the number of vertices in the
	 * triangulation. If the number of points is odd, the last point is discarded.
	 * <p>
	 * The <code>graphMatchedSegments</code> methods are arguably the best
	 * approaches for random segment set generation.
	 * 
	 * @param triangulation triangulation object to use for vertices and graph
	 *                      connectivity
	 * @return set of non-intersecting line segments
	 */
	public static List<PEdge> graphMatchedSegments(IIncrementalTin triangulation) {
		// explained here https://stackoverflow.com/a/72565245/
		final var g = PGS_Triangulation.toGraph(triangulation);
		MatchingAlgorithm<PVector, PEdge> m;
		try {
			m = new KolmogorovWeightedPerfectMatching<>(g, ObjectiveSense.MAXIMIZE);
			m.getMatching();
		} catch (Exception e2) {
			m = new KolmogorovWeightedMatching<>(g, ObjectiveSense.MAXIMIZE);
		}
		return new ArrayList<>(m.getMatching().getEdges());
	}

	/**
	 * Generates a set of N random non-intersecting line segments via brute-forcing.
	 * Plentifully fast enough for many applications.
	 * <p>
	 * The segments are generated within a bounding box anchored at (0, 0) having
	 * the width and height specified.
	 * 
	 * @param width  width of the bounds in which to generate segments; segment x
	 *               coordinates will not be greater than this value
	 * @param height height of the bounds in which to generate segments; segment y
	 *               coordinates will not be greater than this value
	 * @param n      number of segments to generate
	 * @return set of N random non-intersecting line segments
	 * @see #stochasticSegments(double, double, int, double, double, long)
	 */
	public static List<PEdge> stochasticSegments(double width, double height, int n) {
		return stochasticSegments(width, height, n, 1, Math.min(width, height), System.nanoTime());
	}

	/**
	 * Generates a set of N random non-intersecting line segments of the given
	 * length via brute-forcing. Plentifully fast enough for many applications.
	 * <p>
	 * The segments are generated within a bounding box anchored at (0, 0) having
	 * the width and height specified.
	 * 
	 * @param width  width of the bounds in which to generate segments; segment x
	 *               coordinates will not be greater than this value
	 * @param height height of the bounds in which to generate segments; segment y
	 *               coordinates will not be greater than this value
	 * @param n      number of segments to generate
	 * @param length segment length (for all segments)
	 * @return set of N random non-intersecting line segments
	 * @see #stochasticSegments(double, double, int, double, double, long)
	 */
	public static List<PEdge> stochasticSegments(double width, double height, int n, double length) {
		return stochasticSegments(width, height, n, length, length, System.nanoTime());
	}

	/**
	 * Generates a set of N random non-intersecting line segments via brute-forcing.
	 * Plentifully fast enough for many applications.
	 * <p>
	 * The segments are generated within a bounding box anchored at (0, 0) having
	 * the width and height specified.
	 * 
	 * @param width     width of the bounds in which to generate segments; segment x
	 *                  coordinates will not be greater than this value
	 * @param height    height of the bounds in which to generate segments; segment
	 *                  y coordinates will not be greater than this value
	 * @param n         number of segments to generate
	 * @param minLength minimum segment length (inclusive)
	 * @param maxLength maximum segment length (exclusive)
	 * @param seed      number used to initialize the underlying pseudorandom number
	 *                  generator
	 * @return set of N random non-intersecting line segments
	 */
	public static List<PEdge> stochasticSegments(final double width, final double height, final int n, final double minLength, final double maxLength,
			long seed) {
		boolean monoLength = minLength == maxLength;
		final SplittableRandom random = new SplittableRandom(seed);
		List<LineSegment> segments = new ArrayList<>(n);

		find: while (segments.size() < n) {
			final double length = monoLength ? minLength : random.nextDouble(minLength, maxLength);
			final Coordinate a = new Coordinate(random.nextDouble() * (width - length * 2) + length, random.nextDouble() * (height - length * 2) + length);
			final double theta = random.nextDouble(Math.PI * 2);
			final double x = FastMath.cosQuick(theta) * length;
			final double y = FastMath.sinQuick(theta) * length;
			final Coordinate b = new Coordinate(a.x + x, a.y + y);

			for (LineSegment s : segments) {
				if (intersect(a, b, s.p0, s.p1)) {
					continue find;
				}
			}
			LineSegment candidate = new LineSegment(a, b);
			segments.add(candidate);

			/*
			 * Sort by longest segment first. Could increase speed but sorting cancels out
			 * any iteration improvements.
			 */
//			if (!uniLength) {
//				segments.sort((q,w) -> Double.compare(w.getLength(), q.getLength()));
//			}
		}

		return toPEdges(segments);
	}

	/**
	 * Generates a set of N straight parallel segments, centered on a given point.
	 * 
	 * @param centerX the x coordinate of the center of the segments arrangment
	 * @param centerY the y coordinate of the center of the segments arrangment
	 * @param length  length of each segment
	 * @param spacing distance between successive segments
	 * @param angle   line angle in radians, where 0 is parallel to x-axis
	 *                (horizontal)
	 * @param n       number of segments to generate. if odd then the middle segment
	 *                lies on the center point; if even, then the first two segments
	 *                are spaced evenly from the center point
	 * @return set of N parallel line segments
	 */
	public static List<PEdge> parallelSegments(double centerX, double centerY, double length, double spacing, double angle, int n) {
		List<PEdge> edges = new ArrayList<>(n);
		if (n < 1) {
			return edges;
		}
		PVector center = new PVector((float) centerX, (float) centerY);

		float dx = (float) (Math.cos(angle + PConstants.HALF_PI) * spacing);
		float dy = (float) (Math.sin(angle + PConstants.HALF_PI) * spacing);
		float cos = (float) Math.cos(angle);
		float sin = (float) Math.sin(angle);
		float l = (float) length;

		int i;

		float offX = 0, offY = 0;
		if (n % 2 == 1) { // odd number
			PVector a = center.copy().add(new PVector(cos * l, sin * l));
			PVector b = center.copy().add(new PVector(cos * -l, sin * -l));
			edges.add(new PEdge(a, b));
			i = 1;
		} else {
			// half the distance for first two if no middle
			dx /= 2;
			dy /= 2;
			i = 1;
			PVector a = PVector.add(center, new PVector(dx * i, dy * i));
			PVector b = PVector.add(a, new PVector(cos * l, sin * l));
			a.add(cos * -l, sin * -l);
			edges.add(new PEdge(a, b));
			a = PVector.add(center, new PVector(dx * -i, dy * -i));
			b = PVector.add(a, new PVector(cos * l, sin * l));
			a.add(cos * -l, sin * -l);
			edges.add(new PEdge(a, b));
			i = 2;
			dx *= 2;
			dy *= 2;
			offX = dx / 2;
			offY = dy / 2;
		}

		for (; i < 1 + (n / 2); i++) {
			PVector a = PVector.add(center, new PVector(dx * i - offX, dy * i - offY));
			PVector b = PVector.add(a, new PVector(cos * l, sin * l));
			a.add(cos * -l, sin * -l);
			edges.add(new PEdge(a, b));

			a = PVector.add(center, new PVector(dx * -i + offX, dy * -i + offY));
			b = PVector.add(a, new PVector(cos * l, sin * l));
			a.add(cos * -l, sin * -l);
			edges.add(new PEdge(a, b));
		}
		return edges;
	}

	/**
	 * Function that supplies the perpendicular segment length at each sampled
	 * position along a shape component (optionally varying by location, phase, or
	 * normal angle).
	 */
	@FunctionalInterface
	public interface SegmentLengthFn {
		/**
		 * @param x        sampled boundary point x
		 * @param y        sampled boundary point y
		 * @param posFrac  fractional position along the current component in [0,1)
		 *                 (includes {@code startOffset} phase)
		 * @param angleRad angle of the (outward) unit normal in radians at the sample.
		 *                 (Tangent angle is {@code angleRad - PI/2}).
		 * @return desired segment length L. If {@code <= 0}, the segment is skipped.
		 */
		double length(double x, double y, double posFrac, double angleRad);
	}

	/**
	 * Extracts perpendicular segments along each linear component of {@code shape},
	 * with each segment centered on the path/outline.
	 *
	 * <p>
	 * This method samples positions along each linear component (each path or
	 * polygon boundary) at approximately {@code interSegmentDistance} spacing. At
	 * every sampled position it builds a length {@code L} segment that is
	 * perpendicular to the local direction and centered on the sampled point (i.e.
	 * it extends {@code L/2} to each side).
	 *
	 * <p>
	 * {@code startOffset} is wrapped into {@code [0,1)} and acts like a phase along
	 * each component. For closed boundaries, increasing {@code startOffset}
	 * advances sampling counterclockwise and decreasing it advances clockwise. For
	 * open paths, it advances forward/backward along the path direction.
	 *
	 * @param shape                the input {@link PShape} containing rings or line
	 *                             strings
	 * @param interSegmentDistance spacing between successive segments along each
	 *                             component (arc-length units)
	 * @param L                    length of each perpendicular segment (must be
	 *                             &gt; 0)
	 * @param startOffset          fractional phase along each component (0..1);
	 *                             values outside this range are wrapped
	 * @return a list of {@link PEdge} segments from every linear component; empty
	 *         if none produce segments
	 * @since 2.2
	 */
	public static List<PEdge> perpendicularPathSegments(PShape shape, double interSegmentDistance, double L, double startOffset) {
		return perpendicularPathSegments(shape, interSegmentDistance, (x, y, t, a) -> L, startOffset);
	}

	/**
	 * Extracts perpendicular segments along each linear component of {@code shape},
	 * with each segment centered on the path/outline, using a user-supplied
	 * function to vary segment length.
	 *
	 * <p>
	 * This method samples positions along each linear component (each path or
	 * polygon boundary) at approximately {@code interSegmentDistance} spacing. At
	 * every sampled position it builds a length {@code L} segment that is
	 * perpendicular to the local direction and centered on the sampled point (i.e.
	 * it extends {@code L/2} to each side).
	 *
	 * <p>
	 * {@code startOffset} is wrapped into {@code [0,1)} and acts like a phase along
	 * each component. For closed boundaries, increasing {@code startOffset}
	 * advances sampling counterclockwise and decreasing it advances clockwise. For
	 * open paths, it advances forward/backward along the path direction.
	 *
	 * @param shape                the input {@link PShape} containing rings or line
	 *                             strings
	 * @param interSegmentDistance spacing between successive segments along each
	 *                             component (arc-length units)
	 * @param lengthFn             function that returns a per-sample segment length
	 * @param startOffset          fractional phase along each component (0..1);
	 *                             values outside this range are wrapped
	 * @return a list of {@link PEdge} segments; empty if none produce segments
	 * @since 2.2
	 */
	public static List<PEdge> perpendicularPathSegments(PShape shape, double interSegmentDistance, SegmentLengthFn lengthFn, double startOffset) {
		List<PEdge> edges = new ArrayList<>();
		if (interSegmentDistance <= 0 || lengthFn == null) {
			return edges;
		}

		final double startNorm = ((startOffset % 1.0) + 1.0) % 1.0;

		PGS.applyToLinealGeometries(shape, line -> {

			final boolean closed = line.isClosed();

			// Only meaningful for closed components; keeps "outward" consistent.
			if (closed && Orientation.isCCW(line.getCoordinates())) {
				line = line.reverse();
			}

			final IndexedLengthIndexedLine l = new IndexedLengthIndexedLine(line);
			final double end = l.getEndIndex();
			if (end <= 0) {
				return line;
			}

			final int count = Math.max(1, (int) Math.round(end / interSegmentDistance));
			final double increment = 1.0 / count;

			// ε chosen automatically; inline clamp
			final double eps = Math.max(end * 1e-4, Math.min(interSegmentDistance * 0.35, end * 0.01));

			for (int i = 0; i < count; i++) {
				final double posFrac = (startNorm + i * increment) % 1.0;
				final double idx = posFrac * end;

				final double idxM = closed ? wrapIndex(idx - eps, end) : Math.max(0.0, Math.min(end, idx - eps));

				final double idxP = closed ? wrapIndex(idx + eps, end) : Math.max(0.0, Math.min(end, idx + eps));

				final Coordinate pm = l.extractPoint(idxM);
				final Coordinate pp = l.extractPoint(idxP);
				final Coordinate pc = l.extractPoint(idx);

				final double tx = pp.x - pm.x;
				final double ty = pp.y - pm.y;
				final double tlen = Math.sqrt(tx * tx + ty * ty);
				if (tlen == 0) {
					continue;
				}

				// Left-hand unit normal (-ty, tx). For normalized CW rings, this points
				// outward.
				final double nx = -ty / tlen;
				final double ny = tx / tlen;
				final double normalAngle = FastAtan2.atan2(ny, nx);

				final double L = lengthFn.length(pc.x, pc.y, posFrac, normalAngle);
				if (!(L > 0)) {
					continue;
				}

				final double half = L * 0.5;

				PVector a = new PVector((float) (pc.x - nx * half), (float) (pc.y - ny * half));
				PVector b = new PVector((float) (pc.x + nx * half), (float) (pc.y + ny * half));
				edges.add(new PEdge(a, b));
			}

			return line;
		});

		return edges;
	}

	/**
	 * Creates a fabric-like layout of horizontal and vertical segments on a regular
	 * cell grid, controlled by the A (horizontal run), B (vertical run), and C (row
	 * shift) weave parameters.
	 *
	 * <p>
	 * The pattern is built on a rectangular grid of cells (size {@code cellSize}).
	 * Each cell is assigned one of two states (“horizontal on top” or “vertical on
	 * top”) so the grid reads like a simple over/under weaving diagram. Contiguous
	 * same-state cells in a row become one horizontal segment; contiguous
	 * same-state cells in a column become one vertical segment.
	 * </p>
	 *
	 * <ul>
	 * <li><b>A</b> - how many cells in a row the horizontal strand stays on top.
	 * Increasing A lengthens horizontal runs (longer horizontal elements).</li>
	 * <li><b>B</b> - how many cells in a row the vertical strand stays on top.
	 * Increasing B lengthens vertical runs (longer vertical elements).</li>
	 * <li><b>C</b> - how far each successive row is shifted (a phase offset).
	 * Changing C slides the pattern row-by-row and can change where segments meet
	 * or form longer/shorter junctions. (C is applied modulo A + B.)</li>
	 * </ul>
	 *
	 * <p>
	 * Many traditional weaves are expressible with A–B–C (e.g., plain weave 1–1–1,
	 * twill 2–2–1). Patterns where the shift and period are “coprime” (gcd(A+B,
	 * C)=1) tend to produce a single connected repeating motif (“hang together”);
	 * if not, the repeat unit can be larger or the motif can repeat in bands.
	 * </p>
	 *
	 * @param width    domain width
	 * @param height   domain height
	 * @param cellSize size of a grid cell (world units)
	 * @param A        consecutive weft-visible (horizontal on-top) cells per
	 *                 period; controls horizontal run length (A >= 1)
	 * @param B        consecutive warp-visible (vertical on-top) cells per period;
	 *                 controls vertical run length (B >= 1)
	 * @param C        row-to-row offset (phase shift) in cells (C >= 0)
	 * @return list of {@link PEdge} segments representing the weaves
	 * @throws IllegalArgumentException if cellSize <= 0, A <= 0, or B <= 0
	 * @since 2.2
	 */
	public static List<PEdge> weaveSegments(final double width, final double height, final double cellSize, final int A, final int B, final int C) {
		return weaveSegments(width, height, cellSize, A, B, C, false, 1, true);
	}

	/**
	 * Creates a fabric-like layout of horizontal and vertical segments on a regular
	 * cell grid, controlled by the A (horizontal run), B (vertical run), and C (row
	 * shift) weave parameters.
	 *
	 * <p>
	 * The pattern is built on a rectangular grid of cells (size {@code cellSize}).
	 * Each cell is assigned one of two states (“horizontal on top” or “vertical on
	 * top”) so the grid reads like a simple over/under weaving diagram. Contiguous
	 * same-state cells in a row become one horizontal segment; contiguous
	 * same-state cells in a column become one vertical segment.
	 * </p>
	 *
	 * <ul>
	 * <li><b>A</b> — how many cells in a row the horizontal strand stays on top.
	 * Increasing A lengthens horizontal runs (longer horizontal elements).</li>
	 * <li><b>B</b> — how many cells in a row the vertical strand stays on top.
	 * Increasing B lengthens vertical runs (longer vertical elements).</li>
	 * <li><b>C</b> — how far each successive row is shifted (a phase offset).
	 * Changing C slides the pattern row-by-row and can change where segments meet
	 * or form longer/shorter junctions. (C is applied modulo A + B.)</li>
	 * </ul>
	 *
	 * <p>
	 * Many traditional weaves are expressible with A–B–C (e.g., plain weave 1–1–1,
	 * twill 2–2–1). Patterns where the shift and period are “coprime” (gcd(A+B,
	 * C)=1) tend to produce a single connected repeating motif (“hang together”);
	 * if not, the repeat unit can be larger or the motif can repeat in bands.
	 * </p>
	 * 
	 * </ul>
	 * <h3>Endpoint placement and edge behavior</h3>
	 * <p>
	 * The {@code cellFraction} value controls how segment endpoints sit inside the
	 * terminal cell of each run: 0.5 places endpoints at cell centers (short
	 * segments), values closer to 1 move endpoints toward cell edges (longer
	 * segments). If {@code extendSingletonsToEdge} is true, single-cell runs that
	 * touch the outer domain boundary are extended to the grid edge so they are not
	 * rendered as tiny isolated dots at the border.
	 * </p>
	 *
	 * @param width                  domain width
	 * @param height                 domain height
	 * @param cellSize               size of a grid cell (world units)
	 * @param A                      consecutive weft-visible (horizontal on-top)
	 *                               cells per period; controls horizontal run
	 *                               length (A >= 1)
	 * @param B                      consecutive warp-visible (vertical on-top)
	 *                               cells per period; controls vertical run length
	 *                               (B >= 1)
	 * @param C                      row-to-row offset (phase shift) in cells (C >=
	 *                               0)
	 * @param swapColors             if true, swap roles of weft/warp (horizontal ↔
	 *                               vertical)
	 * @param cellFraction           endpoint position inside the run end cells
	 *                               (0.5..1.0 typical)
	 * @param extendSingletonsToEdge if true, extend single-cell runs on the domain
	 *                               boundary to the edge
	 * @return list of {@link PEdge} segments representing the weaves
	 * @since 2.2
	 */
	private static List<PEdge> weaveSegments(final double width, final double height, final double cellSize, final int A, final int B, final int C,
			final boolean swapColors, final double cellFraction, final boolean extendSingletonsToEdge) {
		/*
		 * Implements 'ABC-Auxetics: An Implicit Design Approach for Negative Poisson’s
		 * Ratio Materials'
		 */
		if (cellSize <= 0) {
			throw new IllegalArgumentException("cellSize must be > 0");
		}
		if (A <= 0 || B <= 0) {
			throw new IllegalArgumentException("A and B must be > 0");
		}
		if (!(cellFraction > 0.0)) {
			throw new IllegalArgumentException("cellFraction must be > 0");
		}

		final int P = A + B;

		final double f = Math.max(0.0, Math.min(1.0, cellFraction));

		// Fit an integer grid inside the domain; center it.
		final int cols = Math.max(1, (int) Math.floor(width / cellSize));
		final int rows = Math.max(1, (int) Math.floor(height / cellSize));
		final double gridW = cols * cellSize;
		final double gridH = rows * cellSize;
		final double dx = (width - gridW) * 0.5;
		final double dy = (height - gridH) * 0.5;

		// World-space bounds of the actual grid (may be inset if dx/dy != 0)
		final double gridLeft = dx;
		final double gridRight = dx + gridW;
		final double gridBottom = dy;
		final double gridTop = dy + gridH;

		// Build 2-color matrix: true = weft (horizontal), false = warp (vertical)
		final boolean[][] weft = new boolean[rows][cols];
		for (int r = 0; r < rows; r++) {
			final int shift = Math.floorMod(r * C, P);
			for (int c = 0; c < cols; c++) {
				final int idx = Math.floorMod(c + shift, P);
				boolean isWeft = idx < A;
				if (swapColors) {
					isWeft = !isWeft;
				}
				weft[r][c] = isWeft;
			}
		}

		final List<PEdge> segs = new ArrayList<>();

		// Horizontal segments from contiguous WEFT runs per row
		for (int r = 0; r < rows; r++) {
			int c = 0;
			while (c < cols) {
				if (!weft[r][c]) {
					c++;
					continue;
				}

				final int start = c;
				while (c + 1 < cols && weft[r][c + 1]) {
					c++;
				}
				final int end = c;

				final int runLen = end - start + 1;

				final double y = dy + (r + 0.5) * cellSize;

				// Endpoints controlled by cellFraction:
				// left endpoint is inside first cell at (1-f), right endpoint inside last cell
				// at f.
				double x0 = dx + (start + (1.0 - f)) * cellSize;
				double x1 = dx + (end + f) * cellSize;

				// Special case: single-cell run on boundary becomes half-cell from center to
				// boundary.
				if (extendSingletonsToEdge && runLen == 1) {
					final double cx = dx + (start + 0.5) * cellSize;
					if (start == 0) {
						x0 = gridLeft;
						x1 = cx;
					} else if (end == cols - 1) {
						x0 = cx;
						x1 = gridRight;
					}
				}

				// Clamp to grid bounds (keeps segments from bleeding into margins due to f)
				x0 = Math.max(gridLeft, Math.min(gridRight, x0));
				x1 = Math.max(gridLeft, Math.min(gridRight, x1));

				segs.add(new PEdge(x0, y, x1, y));
				c++;
			}
		}

		// Vertical segments from contiguous WARP runs per column
		for (int c = 0; c < cols; c++) {
			int r = 0;
			while (r < rows) {
				if (weft[r][c]) {
					r++;
					continue;
				} // warp = !weft

				final int start = r;
				while (r + 1 < rows && !weft[r + 1][c]) {
					r++;
				}
				final int end = r;

				final int runLen = end - start + 1;

				final double x = dx + (c + 0.5) * cellSize;

				double y0 = dy + (start + (1.0 - f)) * cellSize;
				double y1 = dy + (end + f) * cellSize;

				if (extendSingletonsToEdge && runLen == 1) {
					final double cy = dy + (start + 0.5) * cellSize;
					if (start == 0) {
						y0 = gridBottom;
						y1 = cy;
					} else if (end == rows - 1) {
						y0 = cy;
						y1 = gridTop;
					}
				}

				y0 = Math.max(gridBottom, Math.min(gridTop, y0));
				y1 = Math.max(gridBottom, Math.min(gridTop, y1));

				segs.add(new PEdge(x, y0, x, y1));
				r++;
			}
		}

		return segs;
	}

	/**
	 * Converts a collection of {@link micycle.pgs.commons.PEdge PEdges} into a
	 * <code>LINES</code> shape.
	 * 
	 * @param segments collection of segments
	 * @return <code>LINES</code> shape representing segments
	 */
	public static PShape toPShape(Collection<PEdge> segments) {
		return toPShape(segments, null, null, 4);
	}

	/**
	 * Converts a collection of {@link micycle.pgs.commons.PEdge PEdges} into a
	 * <code>LINES</code> shape, having the (optional) styling provided.
	 * 
	 * @param segments     collection of PEdge segments
	 * @param strokeColor  nullable/optional (default = {@link Colors#PINK})
	 * @param strokeCap    nullable/optional (default = <code>ROUND</code>)
	 * @param strokeWeight nullable/optional (default = <code>2</code>)
	 * @return shape representing segments
	 */
	public static PShape toPShape(Collection<PEdge> segments, @Nullable Integer strokeColor, @Nullable Integer strokeCap, @Nullable Integer strokeWeight) {
		PShape lines = PGS.prepareLinesPShape(strokeColor, strokeCap, strokeWeight);
		segments.forEach(s -> {
			lines.vertex(s.a.x, s.a.y);
			lines.vertex(s.b.x, s.b.y);
		});
		lines.endShape();

		return lines;
	}

	public static Map<PEdge, Integer> toBag(List<PEdge> edges) {
		Map<PEdge, Integer> edgeBag = new HashMap<>(edges.size());
		for (PEdge edge : edges) {
			edgeBag.put(edge, edgeBag.getOrDefault(edge, 0) + 1);
		}
		return edgeBag;
	}

	/**
	 * Dissolves the edges from a collection of {@link micycle.pgs.commons.PEdge
	 * PEdges} into a set of maximal-length LineStrings in which each unique segment
	 * appears only once. This method works by fusing segments that share endpoints
	 * into longer linear strings.
	 * <p>
	 * This method may be preferred to {@link #toPShape(Collection) toPShape()} when
	 * the input segments form a linear string(s).
	 * 
	 * @param segments Collection of PEdge objects to dissolve into maximal-length
	 *                 LineStrings
	 * @return A PShape object representing the dissolved LineStrings
	 * @since 1.4.0
	 */
	public static PShape dissolve(Collection<PEdge> segments) {
		Geometry g = SegmentStringUtil.toGeometry(fromPEdges(segments), PGS.GEOM_FACTORY);
		if (g.isEmpty()) {
			return new PShape();
		}
		Geometry dissolved = LineDissolver.dissolve(g);
		return PGS_Conversion.toPShape(dissolved);
	}

	/**
	 * 
	 * Computes all intersection points among the supplied edges.
	 * <p>
	 * Each PEdge in {@code edges} is treated as a line segment and intersections
	 * are computed pairwise between segment interiors. Endpoint-endpoint "touches"
	 * are not included.
	 * 
	 * @param edges collection of PEdge objects to test for intersections
	 * @return a List<PVector> containing intersection points; empty if none are
	 *         found
	 * @since 2.2
	 */
	public static List<PVector> intersections(Collection<PEdge> edges) {
		return PGS_Processing.intersections(fromPEdges(edges), false);
	}

	/**
	 * Computes all intersection points between two collections of line segments.
	 * <p>
	 * Given two collections of edges, this method finds and returns all points
	 * where a segment from the first collection intersects with a segment from the
	 * second collection, including intersection endpoints if applicable.
	 * </p>
	 *
	 * @param edgesA the first collection of line segments to compare
	 * @param edgesB the second collection of line segments to compare
	 * @return a list of {@link PVector} representing all intersection points found
	 *         between the two collections of segments
	 * @since 2.1
	 */
	public static List<PVector> intersections(Collection<PEdge> edgesA, Collection<PEdge> edgesB) {
		var segsA = fromPEdges(edgesA);
		var segsB = fromPEdges(edgesB);
		return PGS_Processing.intersections(segsA, segsB);
	}

	/**
	 * Extracts a list of unique PEdge segments representing the given shape.
	 * <p>
	 * This method iterates through all the child shapes of the input shape,
	 * creating PEdge segments for each pair of consecutive vertices.
	 *
	 * @param shape The shape from which to extract the edges. Supports holes and
	 *              GROUP shapes.
	 * @return A list of unique PEdge segments representing the edges of the input
	 *         shape and its child shapes.
	 * @since 1.4.0
	 */
	public static List<PEdge> fromPShape(PShape shape) {
		List<PEdge> edges = new ArrayList<>(shape.getFamily() != PConstants.GROUP ? shape.getVertexCount() : shape.getChildCount() * 4);
		@SuppressWarnings("unchecked")
		List<LineString> strings = LinearComponentExtracter.getLines(PGS_Conversion.fromPShape(shape));
		strings.forEach(s -> {
			Coordinate[] coords = s.getCoordinates();
			boolean closed = coords[0].equals2D(coords[coords.length - 1]);
			for (int i = 0; i < coords.length - (closed ? 0 : 1); i++) {
				Coordinate a = coords[i];
				Coordinate b = coords[(i + 1) % coords.length];
				if (a.equals(b)) {
					continue;
				}
				final PEdge e = new PEdge(a.x, a.y, b.x, b.y);
				edges.add(e);
			}
		});
		return edges;
	}

	/**
	 * Stretches each PEdge segment in the provided list by a specified factor. The
	 * stretch is applied by scaling the distance between the edge's vertices, while
	 * keeping the midpoint of the edge constant.
	 *
	 * @param segments The list of PEdges to be stretched.
	 * @param factor   The factor by which to stretch each PEdge. A value greater
	 *                 than 1 will stretch the edges, while a value between 0 and 1
	 *                 will shrink them.
	 * @return A new List of PEdges representing the stretched edges.
	 * @since 1.4.0
	 */
	public static List<PEdge> stretch(List<PEdge> segments, double factor) {
		List<PEdge> stretchedEdges = new ArrayList<>(segments.size());

		for (PEdge edge : segments) {
			PVector midpoint = PVector.add(edge.a, edge.b).mult(0.5f);
			PVector newA = PVector.add(midpoint, PVector.sub(edge.a, midpoint).mult((float) factor));
			PVector newB = PVector.add(midpoint, PVector.sub(edge.b, midpoint).mult((float) factor));
			stretchedEdges.add(new PEdge(newA, newB));
		}

		return stretchedEdges;
	}

	/**
	 * Removes segments having a length less than the given length from a collection
	 * of segmensts.
	 * 
	 * @param segments  list of segments to filter
	 * @param minLength the minimum segment length to keep
	 * @return a filtered copy of input collection
	 */
	public static List<PEdge> filterByMinLength(List<PEdge> segments, double minLength) {
		return segments.stream().filter(s -> s.length() >= minLength).collect(Collectors.toList());
	}

	/**
	 * Removes segments having a length either less than some fraction or more than
	 * <code>1/fraction</code> of the mean segment length from a collection of
	 * segments.
	 * <p>
	 * This method can be used to homogenise a segment set.
	 * 
	 * @param segments list of segments to filter
	 * @param fraction fraction of mean length to keep segments
	 * @return a filtered copy of input collection
	 */
	public static List<PEdge> filterByAverageLength(List<PEdge> segments, double fraction) {
		final double lenAvg = segments.stream().mapToDouble(PEdge::length).average().orElse(0);
		return segments.stream().filter(e -> e.length() > fraction * lenAvg && e.length() < 1 / fraction * lenAvg).collect(Collectors.toList());
	}

	/**
	 * Removes axis-aligned (horizontal and vertical) segments (within a given angle
	 * tolerance) from a collection of segments.
	 * 
	 * @param segments   list of segments to filter
	 * @param angleDelta angle tolerance, in radians
	 * @return a filtered copy of the input where axis-aligned segments have been
	 *         removed
	 */
	public static List<PEdge> filterAxisAligned(List<PEdge> segments, double angleDelta) {
		List<PEdge> filtered = new ArrayList<>(segments);
		filtered.removeIf(s -> {
			double angle = Math.abs(FastAtan2.atan2(s.b.y - s.a.y, s.b.x - s.a.x)) % (Math.PI / 2);
			return (angle + angleDelta > Math.PI / 2 || angle - angleDelta < 0);
		});
		return filtered;
	}

	/**
	 * Removes segments that are near others. More specifically, removes any segment
	 * that lies closer than a given threshold to a previously retained segment.
	 * <p>
	 * This method uses a simple greedy, O(n²) algorithm: it iterates through the
	 * input {@code segments} in order, keeps the first segment unconditionally, and
	 * for each subsequent segment rejects it if its distance to any segment already
	 * in the output list is less than {@code distance}.
	 * </p>
	 *
	 * @param segments the list of segments to filter.
	 * @param distance the minimum allowed distance between any two segments in the
	 *                 returned list; must be non-negative.
	 * @return a new {@link List} containing those segments from the input (in
	 *         original order) such that each pair of segments in the returned list
	 *         is at least {@code distance} apart.
	 * @since 2.1
	 */
	public static List<PEdge> filterNear(List<PEdge> segments, double distance) {
		List<PEdge> filtered = new ArrayList<>();
		for (PEdge seg : segments) {
			boolean tooClose = false;
			for (PEdge kept : filtered) {
				if (seg.distance(kept) < distance) {
					tooClose = true;
					break;
				}
			}
			if (!tooClose) {
				filtered.add(seg);
			}
		}
		return filtered;
	}

	/**
	 * Retains line segments from a set of line segments that are wholly contained
	 * within a given shape.
	 *
	 * @param segments a list of line segments to check for containment within the
	 *                 shape
	 * @param shape    the polygonal shape to check for interior segments
	 * @return a list of interior segments contained within the shape
	 * @since 1.4.0
	 */
	public static List<PEdge> getPolygonInteriorSegments(List<PEdge> segments, PShape shape) {
		Geometry g = PGS_Conversion.fromPShape(shape);
		final Collection<?> segmentStrings = SegmentStringUtil.extractBasicSegmentStrings(g);
		final SegmentSetMutualIntersector mci = new MCIndexSegmentSetMutualIntersector(segmentStrings);

		Map<SegmentString, PEdge> map = new HashMap<>(segments.size());
		segments.forEach(e -> map.put(PGS.createSegmentString(e.a, e.b), e));

		List<PEdge> interiorSegments = new ArrayList<>();
		IndexedPointInAreaLocator locator = new IndexedPointInAreaLocator(g);
		RobustLineIntersector lineIntersector = new RobustLineIntersector();
		Set<SegmentString> segSet = new HashSet<>(fromPEdges(segments));
		Set<SegmentString> segSet2 = new HashSet<>();

		mci.process(fromPEdges(segments), new SegmentIntersector() {
			@Override
			public void processIntersections(SegmentString e0, int segIndex0, SegmentString e1, int segIndex1) {
				lineIntersector.computeIntersection(e0.getCoordinate(0), e0.getCoordinate(1), e1.getCoordinate(0), e1.getCoordinate(1));

				if (lineIntersector.getIntersectionNum() > 0) { // no intersection -- either inside or outside
					if (locator.locate(e0.getCoordinate(0)) != Location.EXTERIOR) { // one point is inside
						interiorSegments.add(map.get(e1));
					}
				} else {
					segSet2.add(e0);
				}
			}

			@Override
			public boolean isDone() {
				return false;
			}
		});

		segSet.removeAll(segSet2);
		segSet.removeIf(s -> locator.locate(s.getCoordinate(1)) == Location.EXTERIOR || locator.locate(s.getCoordinate(0)) == Location.EXTERIOR);

		return fromSegmentString(segSet);
	}

	private static List<PEdge> toPEdges(Collection<LineSegment> segments) {
		List<PEdge> edges = new ArrayList<>(segments.size());
		segments.forEach(s -> {
			PEdge e = new PEdge(s.p0.x, s.p0.y, s.p1.x, s.p1.y);
			edges.add(e);
		});

		return edges;
	}

	private static List<PEdge> fromSegmentString(Collection<SegmentString> segments) {
		List<PEdge> edges = new ArrayList<>(segments.size());
		segments.forEach(s -> {
			PEdge e = new PEdge(s.getCoordinate(0).x, s.getCoordinate(0).y, s.getCoordinate(1).x, s.getCoordinate(1).y);
			edges.add(e);
		});

		return edges;
	}

	private static List<SegmentString> fromPEdges(Collection<PEdge> edges) {
		List<SegmentString> segments = new ArrayList<>(edges.size());
		edges.forEach(e -> {
			SegmentString s = PGS.createSegmentString(e.a, e.b);
			segments.add(s);
		});

		return segments;
	}

	private static boolean ccw(Coordinate A, Coordinate B, Coordinate C) {
		return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x);
	}

	private static boolean intersect(Coordinate A, Coordinate B, Coordinate C, Coordinate D) {
		return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);
	}

	private static double wrapIndex(double idx, double end) {
		double r = idx % end;
		if (r < 0) {
			r += end;
		}
		return r;
	}
}
