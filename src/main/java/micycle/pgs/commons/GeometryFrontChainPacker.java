package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;
import java.util.SplittableRandom;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.quadtree.Quadtree;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;

import com.github.micycle1.geoblitz.CircleIndex;
import com.github.micycle1.geoblitz.PointDistanceIndex;
import com.github.micycle1.geoblitz.YStripesPointInAreaLocator;

/**
 * Circle packing of an arbitrary polygonal boundary using the front-chain
 * packing algorithm from 'Visualization of Large Hierarchical Data by Circle
 * Packing' (adapted from Mike Bostock's rectangle implementation).
 * <p>
 * Circles are seeded at the centre of the polygon's envelope and grown outwards
 * in a spiral. Unlike the rectangle variant, this version respects the true
 * shape of the input {@link Polygon} (including holes and concave regions):
 * <ul>
 * <li>A candidate circle whose disc would cross the polygon boundary is
 * iteratively <em>shrunk</em> (down to, at most, {@code radiusMin}) and
 * re-seated tangentially against its two parent circles so that it fits inside
 * the boundary.</li>
 * <li>Candidates that still cannot fit (centre outside the polygon, or the
 * available clearance is smaller than {@code radiusMin}) are kept in the
 * front-chain — so the spiral can continue past concavities and holes — but are
 * excluded from the output.</li>
 * </ul>
 * Boundary queries are accelerated with {@link IndexedFacetDistance} (STR-tree
 * over the boundary linework) and {@link IndexedPointInAreaLocator}; overlap
 * guarding of boundary-shrunk circles uses an incremental {@link Quadtree}.
 *
 * @author Mike Bostock (original front-chain algorithm)
 * @author JTS/arbitrary-geometry adaptation
 */
public class GeometryFrontChainPacker {

	private static final double EPSILON = 1e-9;
	/** Slack used for tangency/overlap tests (matches intersects()). */
	private static final double OVERLAP_EPSILON = 1e-6;
	/** Max tangent re-seating iterations when shrinking a circle to fit. */
	private static final int MAX_SHRINK_ITERATIONS = 16;

	private final double radiusMin, radiusMax;
	private final SplittableRandom rand;

	private final double centerX, centerY; // spiral origin (envelope centre)
	private final double envWidth, envHeight;

	/**
	 * Square of the max distance between a circle centre (of radiusMax) and the
	 * spiral origin. The packing terminates when the next candidate centre falls
	 * outside this distance (i.e. the spiral has covered the whole envelope).
	 */
	private final double maxDistSq;

	private final PointOnGeometryLocator interiorLocator;
	private final PointDistanceIndex boundaryDistance;

	/**
	 * Spatial index over every *accepted* circle. Used to guard boundary-shrunk
	 * circles against overlapping already-placed circles: shrinking pulls a circle
	 * back towards its two parents, where it may collide with interior (non-front)
	 * circles that the front-chain intersection scan cannot see.
	 */
	private final CircleIndex<Coordinate> circleIndex = new CircleIndex<>();

	private final List<Coordinate> circles;

	/**
	 * Creates a packer and computes the packing immediately.
	 *
	 * @param polygon   polygonal region to pack (may be concave / contain holes)
	 * @param radiusMin minimum circle radius (must be &gt; 0)
	 * @param radiusMax maximum circle radius
	 * @param seed      RNG seed (packings are deterministic per seed)
	 */
	public GeometryFrontChainPacker(Polygon polygon, double radiusMin, double radiusMax, long seed) {
		this.radiusMin = Math.max(EPSILON, Math.min(radiusMin, radiusMax));
		this.radiusMax = Math.max(EPSILON, Math.max(radiusMin, radiusMax));
		this.rand = new SplittableRandom(seed);

		final Envelope e = polygon.getEnvelopeInternal();
		this.centerX = (e.getMinX() + e.getMaxX()) / 2d;
		this.centerY = (e.getMinY() + e.getMaxY()) / 2d;
		this.envWidth = e.getWidth();
		this.envHeight = e.getHeight();
		this.maxDistSq = envWidth * envHeight / 2 + this.radiusMax * this.radiusMax;

		this.interiorLocator = new YStripesPointInAreaLocator(polygon);
		this.boundaryDistance = new PointDistanceIndex(polygon.getBoundary());

		this.circles = pack();
	}

	/**
	 * @return packed circles; x/y are the circle centre, z is its radius
	 */
	public List<Coordinate> getCircles() {
		return circles;
	}

	private List<Coordinate> pack() {
		final List<Coordinate> chain = new ArrayList<>();
		final Set<Coordinate> discarded = Collections.newSetFromMap(new IdentityHashMap<>());

		// init first chain of 3, seeded at the envelope centre
		chain.add(new Coordinate(0, 0, randomRadius()));
		chain.add(new Coordinate(0, 0, randomRadius()));
		chain.add(new Coordinate(0, 0, randomRadius()));

		final Coordinate A = chain.get(0);
		final Coordinate B = chain.get(1);
		Coordinate C = chain.get(2);
		A.x = centerX - B.z;
		A.y = centerY;
		B.x = centerX + A.z;
		B.y = centerY;
		tangentPosition(B, A, C);

		// the first three circles may themselves straddle the boundary
		accept(A, fitInPlace(A), discarded);
		accept(B, fitInPlace(B), discarded);
		accept(C, fit(B, A, C, false), discarded);

		Node a = new Node(A);
		Node b = new Node(B);
		Node c = new Node(C);
		a.next = c.previous = b;
		b.next = a.previous = c;
		c.next = b.previous = a;

		chain.add(new Coordinate(0, 0, randomRadius()));

		boolean fitted;
		int iter = 0;
		double sampledRadius = chain.get(chain.size() - 1).z;
		pack: while (true) {
			C = chain.get(chain.size() - 1);
			// restore the sampled radius: a retry after a front splice gets new
			// parents, which may allow the full radius again
			C.z = sampledRadius;
			tangentPosition(a.c, b.c, C);
			if (!withinBounds(C)) {
				break; // spiral has covered the whole envelope
			}
			fitted = fit(a.c, b.c, C, true);
			c = new Node(C);

			// Find the closest intersecting circle on the front-chain, if any.
			// "Closeness" is determined by linear distance along the front-chain.
			Node j = b.next;
			Node k = a.previous;
			double sj = b.c.z;
			double sk = a.c.z;
			do {
				// The scan may legitimately walk a large part of the front, which
				// grows with the packing; the bound below only trips on a corrupted
				// chain. Crucially, when it trips the candidate is *discarded* from
				// the output rather than inserted as a visibly overlapping circle.
				if (++iter > chain.size() + 16) {
					fitted = false;
					break;
				}
				if (sj <= sk) {
					if (intersects(j.c, c.c)) {
						b = j;
						a.next = b;
						b.previous = a;
						continue pack;
					}
					sj += j.c.z;
					j = j.next;
				} else {
					if (intersects(k.c, c.c)) {
						a = k;
						a.next = b;
						b.previous = a;
						continue pack;
					}
					sk += k.c.z;
					k = k.previous;
				}
			} while (j != k.next);
			iter = 0;

			// Success! Insert the new circle c between a and b.
			accept(C, fitted, discarded);
			c.previous = a;
			c.next = b;
			a.next = b.previous = b = c;

			// Compute the new closest circle pair to the spiral origin.
			double aa = score(a);
			double ca;
			while ((c = c.next) != b) {
				if ((ca = score(c)) < aa) {
					a = c;
					aa = ca;
				}
			}
			b = a.next;
			chain.add(new Coordinate(0, 0, randomRadius())); // last circle fit within bounds; add another
			sampledRadius = chain.get(chain.size() - 1).z;
		}

		chain.remove(chain.size() - 1); // trailing candidate that fell outside the spiral bounds

		final List<Coordinate> out = new ArrayList<>(chain.size());
		for (Coordinate circle : chain) {
			if (!discarded.contains(circle)) {
				out.add(circle);
			}
		}
		return out;
	}

	/**
	 * Records a circle's fate: accepted circles enter the spatial index (so later
	 * boundary-shrunk circles can avoid them); rejected circles are excluded from
	 * the output but remain on the front-chain.
	 */
	private void accept(Coordinate c, boolean fitted, Set<Coordinate> discarded) {
		if (fitted) {
			circleIndex.insert(c.x, c.y, c.z, c);
		} else {
			discarded.add(c);
		}
	}

	/**
	 * Positions circle c tangentially against circles a and b (using c's current
	 * radius), writing the result into c.
	 */
	private static void tangentPosition(Coordinate b, Coordinate a, Coordinate c) {
		final double dx = b.x - a.x;
		final double dy = b.y - a.y;
		final double d2 = dx * dx + dy * dy;
		if (d2 > EPSILON) {
			final double a2 = (a.z + c.z) * (a.z + c.z);
			final double b2 = (b.z + c.z) * (b.z + c.z);
			if (a2 > b2) {
				final double x = (d2 + b2 - a2) / (2 * d2);
				final double y = Math.sqrt(Math.max(0d, b2 / d2 - x * x));
				c.x = b.x - x * dx - y * dy;
				c.y = b.y - x * dy + y * dx;
			} else {
				final double x = (d2 + a2 - b2) / (2 * d2);
				final double y = Math.sqrt(Math.max(0d, a2 / d2 - x * x));
				c.x = a.x + x * dx - y * dy;
				c.y = a.y + x * dy + y * dx;
			}
		} else {
			c.x = a.x + c.z;
			c.y = a.y;
		}
	}

	/**
	 * Makes circle c respect the polygon boundary. If c's disc would cross the
	 * boundary, its radius is shrunk (bounded below by radiusMin) and c is
	 * re-seated tangentially against its parents a and b, since a smaller radius
	 * yields a different tangent position.
	 * <p>
	 * Because shrinking pulls the centre back towards the parents — into space that
	 * may already be occupied by interior circles the front-chain scan cannot
	 * detect — any circle that has been shrunk is additionally checked against the
	 * spatial index of accepted circles, and shrunk further if it overlaps one.
	 * This is what prevents "double layer" artefacts near the boundary. Virgin
	 * (unshrunk) placements are left untouched so the normal front-chain
	 * intersection/splicing logic keeps driving the packing.
	 *
	 * @param guardOverlaps whether to apply the accepted-circle overlap guard once
	 *                      shrinking has occurred
	 * @return true if c (possibly shrunk) lies fully inside the polygon and
	 *         overlaps no accepted circle
	 */
	private boolean fit(Coordinate b, Coordinate a, Coordinate c, boolean guardOverlaps) {
		boolean shrunk = false;
		for (int i = 0; i < MAX_SHRINK_ITERATIONS; i++) {
			if (isExterior(c)) {
				// A large radius can push the tangent position outside the polygon.
				// A smaller radius pulls the centre back towards the parents, so
				// shrink and re-seat instead of failing outright. Use the distance
				// to the boundary as a hint for how much to cut.
				if (c.z <= radiusMin + EPSILON) {
					return false; // already minimal and still outside
				}
				final double overshoot = distanceToBoundary(c);
				c.z = Math.max(radiusMin, Math.min(c.z * 0.5, c.z - overshoot));
				shrunk = true;
				tangentPosition(b, a, c);
				continue;
			}
			final double allowed = allowedRadius(c, guardOverlaps && shrunk);
			if (allowed >= c.z - OVERLAP_EPSILON) {
				return true; // fits
			}
			final double target = Math.max(radiusMin, allowed);
			if (target >= c.z - EPSILON) {
				return false; // at radiusMin and still overlapping
			}
			c.z = target;
			shrunk = true;
			tangentPosition(b, a, c);
		}
		return !isExterior(c) && allowedRadius(c, guardOverlaps && shrunk) >= c.z - OVERLAP_EPSILON;
	}

	/**
	 * Largest radius circle c could have at its current centre: the clearance to
	 * the polygon boundary, optionally further constrained by the distance to
	 * already-accepted circles.
	 */
	private double allowedRadius(Coordinate c, boolean guardOverlaps) {
		double allowed = distanceToBoundary(c);
		if (guardOverlaps && allowed > 0 && circleIndex.size() > 0) {
			final var n = circleIndex.nearest(c.x, c.y);
			if (n != null && n.clearance < allowed) {
				allowed = n.clearance;
			}
		}
		return allowed;
	}

	/**
	 * Boundary-fits a circle whose position is fixed (used for the seed circles):
	 * shrinks the radius only, without tangential re-seating.
	 */
	private boolean fitInPlace(Coordinate c) {
		if (isExterior(c)) {
			return false;
		}
		final double clearance = distanceToBoundary(c);
		if (clearance >= c.z - EPSILON) {
			return true;
		}
		if (clearance < radiusMin) {
			return false;
		}
		c.z = clearance;
		return true;
	}

	private boolean isExterior(Coordinate c) {
		return interiorLocator.locate(c) == Location.EXTERIOR;
	}

	private double distanceToBoundary(Coordinate c) {
		return boundaryDistance.distance(c);
	}

	/**
	 * Determines whether the candidate circle centre is still within the spiral's
	 * working region (the polygon's envelope, roughly). Used as the packing
	 * termination condition.
	 */
	private boolean withinBounds(Coordinate v) {
		final double dx = v.x - centerX;
		final double dy = v.y - centerY;
		return dx * dx + dy * dy < maxDistSq;
	}

	private static boolean intersects(Coordinate a, Coordinate b) {
		final double dr = a.z + b.z - OVERLAP_EPSILON;
		final double dx = b.x - a.x;
		final double dy = b.y - a.y;
		return dr > 0 && dr * dr > dx * dx + dy * dy;
	}

	private double score(Node node) {
		final Coordinate a = node.c;
		final Coordinate b = node.next.c;
		final double ab = a.z + b.z;
		final double cx = (a.x * b.z + b.x * a.z) / ab - centerX;
		final double cy = (a.y * b.z + b.y * a.z) / ab - centerY;
		return Math.max(Math.abs(cx * envHeight), Math.abs(cy * envWidth));
	}

	private double randomRadius() {
		return radiusMin == radiusMax ? radiusMin : rand.nextDouble(radiusMin, radiusMax);
	}

	private static class Node {

		final Coordinate c;
		Node next, previous;

		Node(Coordinate circle) {
			this.c = circle;
		}
	}

}