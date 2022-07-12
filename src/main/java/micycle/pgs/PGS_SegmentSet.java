package micycle.pgs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.SplittableRandom;

import org.jgrapht.alg.interfaces.MatchingAlgorithm;
import org.jgrapht.alg.matching.blossom.v5.KolmogorovWeightedPerfectMatching;
import org.jgrapht.alg.matching.blossom.v5.ObjectiveSense;
import org.jgrapht.graph.SimpleGraph;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.SegmentString;
import org.tinfour.common.IIncrementalTin;

import micycle.pgs.color.RGB;
import micycle.pgs.commons.Nullable;
import micycle.pgs.commons.PEdge;
import net.jafama.FastMath;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Generation of random sets of <b>non-intersecting</b> line segments (and
 * associated functions).
 * <p>
 * Methods in this class output segments as collections of
 * {@link micycle.pgs.commons.PEdge PEdges}; such collections can be converted
 * into LINES PShapes with {@link #toPShape(Collection)
 * toPShape(Collection<PEdge>)}.
 * 
 * @author Michael Carleton
 * @since 1.2.1
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
			final Coordinate a = new Coordinate(random.nextDouble() * (width - length * 2) + length,
					random.nextDouble() * (height - length * 2) + length);
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
	 * approaches for random segment set generation. Graph matched / graph perfect
	 * matching. In this method, the input point set is triangulated and a matching
	 * ran on that.
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
	 * Number of segments = #vertices/2
	 * 
	 * Let P be a set of n points, not all on the same line. Let k be the number of
	 * points on the boundary of the convex hull of P. Any triangulation of P has 1)
	 * 2n − 2 − k triangles, and 2) 3n − 3 − k edges
	 * 
	 * @param triangulation
	 * @return
	 */

	/**
	 * Generates non-intersecting segments via a <i>Perfect matching</i> algorithm
	 * applied to the given triangulation.
	 * <p>
	 * Generates N/2 segments, where N is the number of vertices in the
	 * triangulation. If the number of points is odd, the last point is discarded.
	 * <p>
	 * The <code>graphMatchedSegments</code> methods are arguably the best
	 * approaches for random segment set generation. Graph matched / graph perfect
	 * matching. In this method, the input point set is triangulated and a matching
	 * ran on that.
	 * 
	 * @param points point set from which to compute segments
	 * @return set of non-intersecting line segments
	 */
	public static List<PEdge> graphMatchedSegments(IIncrementalTin triangulation) {
		// explained here https://stackoverflow.com/a/72565245/
		final SimpleGraph<PVector, PEdge> g = PGS_Triangulation.toGraph(triangulation);
		MatchingAlgorithm<PVector, PEdge> m;
		try {
			m = new KolmogorovWeightedPerfectMatching<>(g, KolmogorovWeightedPerfectMatching.DEFAULT_OPTIONS, ObjectiveSense.MAXIMIZE);
			return new ArrayList<>(m.getMatching().getEdges());
		} catch (IllegalArgumentException e) {
			// catch exception if no perfect matching is possible
			return new ArrayList<>();
		}
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
		return stochasticSegments(width, height, n, 1, Math.min(width, height), System.currentTimeMillis());
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
		return stochasticSegments(width, height, n, length, length, System.currentTimeMillis());
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
	public static List<PEdge> stochasticSegments(final double width, final double height, final int n, final double minLength,
			final double maxLength, long seed) {
		boolean monoLength = minLength == maxLength;
		final SplittableRandom random = new SplittableRandom(seed);
		List<LineSegment> segments = new ArrayList<>(n);

		find: while (segments.size() < n) {
			final double length = monoLength ? minLength : random.nextDouble(minLength, maxLength);
			final Coordinate a = new Coordinate(random.nextDouble() * (width - length * 2) + length,
					random.nextDouble() * (height - length * 2) + length);
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
	 * Converts a collection of {@link micycle.pgs.commons.PEdge PEdges} into a
	 * <code>LINES</code> shape.
	 * 
	 * @param segments collection of segments
	 * @return shape representing segments
	 */
	public static PShape toPShape(Collection<PEdge> segments) {
		return toPShape(segments, null, null, 4);
	}

	/**
	 * Converts a collection of {@link micycle.pgs.commons.PEdge PEdges} into a
	 * <code>LINES</code> shape, having the (optional) styling provided.
	 * 
	 * @param segments     collection of segments
	 * @param strokeColor  nullable/optional (default = {@link RGB#PINK})
	 * @param strokeCap    nullable/optional (default = <code>ROUND</code>)
	 * @param strokeWeight nullable/optional (default = <code>2</code>)
	 * @return shape representing segments
	 */
	public static PShape toPShape(Collection<PEdge> segments, @Nullable Integer strokeColor, @Nullable Integer strokeCap,
			@Nullable Integer strokeWeight) {
		PShape lines = PGS.prepareLinesPShape(strokeColor, strokeCap, strokeWeight);
		segments.forEach(s -> {
			lines.vertex(s.a.x, s.a.y);
			lines.vertex(s.b.x, s.b.y);
		});
		lines.endShape();

		return lines;
	}

	private static List<PEdge> toPEdges(List<LineSegment> segments) {
		List<PEdge> edges = new ArrayList<>(segments.size());
		segments.forEach(s -> {
			PEdge e = new PEdge(s.p0.x, s.p0.y, s.p1.x, s.p1.y);
			edges.add(e);
		});

		return edges;
	}
	

	private static boolean ccw(Coordinate A, Coordinate B, Coordinate C) {
		return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x);
	}

	private static boolean intersect(Coordinate A, Coordinate B, Coordinate C, Coordinate D) {
		return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);
	}
}
