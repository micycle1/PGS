package micycle.pgs.commons;

import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedPolygon;
import org.locationtech.jts.index.strtree.STRtree;

import net.metaopt.swarm.FitnessFunction;
import net.metaopt.swarm.pso.Particle;
import net.metaopt.swarm.pso.Swarm;

/**
 * Finds an approximate largest area triangle of arbitrary orientation in a
 * concave polygon via particle swarm optimisation.
 * 
 * @author Michael Carleton
 *
 */
public class MaximumInscribedTriangle {

	private static final int SWARM_SIZE = 3000;
	private static final int MAX_GENERATIONS = 1000;
	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory();
	private static final Random RANDOM = ThreadLocalRandom.current();

	private final Swarm swarm;
	private final TriangleFitness fitnessFunction;

	double minArea = 0;

	public MaximumInscribedTriangle(Polygon polygon) {
		fitnessFunction = new TriangleFitness(polygon);
		// For triangle, our candidate vector is of dimension 6: [x1,y1,x2,y2,x3,y3].
		swarm = new ParallelSwarm(SWARM_SIZE, new TriangleCandidate(), fitnessFunction);

		// Use the envelope of the polygon as the search bounds.
		Envelope e = polygon.getEnvelopeInternal();
		double[] minPosition = new double[] { e.getMinX(), e.getMinY(), e.getMinX(), e.getMinY(), e.getMinX(), e.getMinY() };
		double[] maxPosition = new double[] { e.getMaxX(), e.getMaxY(), e.getMaxX(), e.getMaxY(), e.getMaxX(), e.getMaxY() };

		var r = MaximumInscribedCircle.getRadiusLine(polygon, 1).getLength();
		/*
		 * Area of equilateral triangle inscribed in MIC. Sets a minimum bound to the
		 * MIT area.
		 */
		minArea = (3 * Math.sqrt(3) / 4) * Math.pow(r, 2);

		swarm.setMinPosition(minPosition);
		swarm.setMaxPosition(maxPosition);
		swarm.init();
		initializeSwarm(swarm, e);
	}

	/**
	 * Initialize the swarm particles. Several groups use different heuristics: -
	 * Some choose points on polygon envelope edges. - Others choose random points
	 * inside the envelope.
	 */
	private void initializeSwarm(Swarm swarm, Envelope e) {
		Particle[] particles = swarm.getParticles();
		int particlesPerGroup = SWARM_SIZE / 4;
		double minX = e.getMinX();
		double minY = e.getMinY();
		double width = e.getWidth();
		double height = e.getHeight();

		for (int i = 0; i < SWARM_SIZE; i++) {
			double[] position = new double[6];

			if (i < particlesPerGroup) {
				// Group 1: Pick vertices on the envelope boundary.
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = RANDOM.nextBoolean() ? e.getMinY() : e.getMaxY();
				position[2] = minX + RANDOM.nextDouble() * width;
				position[3] = RANDOM.nextBoolean() ? e.getMinY() : e.getMaxY();
				position[4] = minX + RANDOM.nextDouble() * width;
				position[5] = RANDOM.nextBoolean() ? e.getMinY() : e.getMaxY();
			} else if (i < 2 * particlesPerGroup) {
				// Group 2: Random points in the envelope.
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = minY + RANDOM.nextDouble() * height;
				position[2] = minX + RANDOM.nextDouble() * width;
				position[3] = minY + RANDOM.nextDouble() * height;
				position[4] = minX + RANDOM.nextDouble() * width;
				position[5] = minY + RANDOM.nextDouble() * height;
			} else if (i < 3 * particlesPerGroup) {
				// Group 3: Two points from one edge and one from the opposite edge.
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = e.getMinY();
				position[2] = minX + RANDOM.nextDouble() * width;
				position[3] = e.getMinY();
				position[4] = minX + RANDOM.nextDouble() * width;
				position[5] = e.getMaxY();
			} else {
				// Group 4: Random exploration.
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = minY + RANDOM.nextDouble() * height;
				position[2] = minX + RANDOM.nextDouble() * width;
				position[3] = minY + RANDOM.nextDouble() * height;
				position[4] = minX + RANDOM.nextDouble() * width;
				position[5] = minY + RANDOM.nextDouble() * height;
			}

			particles[i].setPosition(position);
			particles[i].setVelocity(new double[] { 0, 0, 0, 0, 0, 0 });
			particles[i].setBestPosition(position.clone());
		}
	}

	/**
	 * The main entry method which runs the swarm evolution, then refines the best
	 * candidate using Apache Commons Math's Nelder–Mead optimizer.
	 */
	public Polygon computeMIT() {
		int gen = 0, stableCount = 0;
		double lastFitness = Double.MIN_VALUE;

		// Global search via swarm evolution.
		while (gen++ < MAX_GENERATIONS) {
			swarm.evolve();

			// Occasionally reinitialize worst performing particles.
			if (gen % 100 == 0) {
				reinitializeWorstParticles(swarm);
			}
			double bestFitness = swarm.getBestFitness();
			if (bestFitness == lastFitness) {
				if (stableCount++ > 75) {
					break;
				}
			} else {
				stableCount = 0;
			}
			lastFitness = bestFitness;
		}

		double[] bestCandidate = swarm.getBestPosition().clone();
		bestCandidate = refineCandidate(bestCandidate);
		return triangleFromCoords(bestCandidate);
	}

	/**
	 * Reinitialize the bottom 10% of swarm particles around the best candidate.
	 */
	private void reinitializeWorstParticles(Swarm swarm) {
		Particle[] particles = swarm.getParticles();
		double[] bestPos = swarm.getBestPosition();
		int reinitCount = SWARM_SIZE / 10;
		double[] maxPos = swarm.getMaxPosition();
		double[] minPos = swarm.getMinPosition();
		for (int i = 0; i < reinitCount; i++) {
			double[] newPos = new double[6];
			for (int j = 0; j < 6; j++) {
				double range = (maxPos[j] - minPos[j]) * 0.1;
				newPos[j] = bestPos[j] + (RANDOM.nextDouble() - 0.5) * range;
				newPos[j] = Math.max(minPos[j], Math.min(maxPos[j], newPos[j]));
			}
			particles[i].setPosition(newPos);
			particles[i].setVelocity(new double[] { 0, 0, 0, 0, 0, 0 });
		}
	}

	/**
	 * Uses Nelder–Mead simplex optimization (via Apache Commons Math) to refine the
	 * best candidate. The objective function is defined as the negative fitness.
	 */
	private double[] refineCandidate(double[] candidate) {
		MultivariateFunction objective = point -> -fitnessFunction.evaluate(point);

		SimplexOptimizer optimizer = new SimplexOptimizer(1e-8, 1e-8);
		NelderMeadSimplex simplex = new NelderMeadSimplex(candidate.length);
		try {
			PointValuePair result = optimizer.optimize(new MaxEval(1000), new ObjectiveFunction(objective), GoalType.MINIMIZE, new InitialGuess(candidate),
					simplex);
			return result.getPoint();
		} catch (Exception e) {
			return candidate;
		}
	}

	/**
	 * Constructs a triangle polygon from a candidate vector. Candidate parameters:
	 * [x1, y1, x2, y2, x3, y3]
	 */
	private static Polygon triangleFromCoords(double[] position) {
		// Create three coordinates for the triangle.
		Coordinate p0 = new Coordinate(position[0], position[1]);
		Coordinate p1 = new Coordinate(position[2], position[3]);
		Coordinate p2 = new Coordinate(position[4], position[5]);
		// Ensure the ring is closed.
		Coordinate[] coords = new Coordinate[] { p0, p1, p2, p0 };
		return GEOM_FACTORY.createPolygon(coords);
	}

	/**
	 * The fitness function for candidate triangles in a concave polygon. Instead of
	 * calling PreparedGeometry.containsProperly(…), we use a custom “fast–path”
	 * test: 1. All triangle vertices must be inside the polygon (using
	 * ray–casting). 2. None of the triangle edges may intersect any polygon
	 * boundary edge.
	 *
	 * If both tests pass, return the triangle’s absolute area.
	 */
	private class TriangleFitness extends FitnessFunction {
		// Cached polygon vertices (exterior ring) for fast point testing.
		private double[] polyX;
		private double[] polyY;
		// Spatial index over the polygon’s edges.
		private STRtree polygonEdgeIndex;
		private final boolean hasHoles;
		private PreparedPolygon cache;

		TriangleFitness(Polygon polygon) {
			hasHoles = polygon.getNumInteriorRing() > 0;
			if (hasHoles) {
				cache = new PreparedPolygon(polygon);
			} else { // faster approach (but doesn't detect holes)
				Coordinate[] coords = polygon.getExteriorRing().getCoordinates();
				int n = coords.length - 1; // last coordinate is a duplicate of the first
				polyX = new double[n];
				polyY = new double[n];
				for (int i = 0; i < n; i++) {
					polyX[i] = coords[i].x;
					polyY[i] = coords[i].y;
				}
				// Build an STRtree on polygon edges.
				polygonEdgeIndex = new STRtree();
				for (int i = 0; i < n; i++) {
					Coordinate p0 = coords[i];
					Coordinate p1 = coords[(i + 1) % n];
					LineSegment seg = new LineSegment(p0, p1);
					polygonEdgeIndex.insert(new Envelope(seg.p0, seg.p1), seg);
				}
				polygonEdgeIndex.build();
			}
		}

		@Override
		public double evaluate(double[] position) {
			double x0 = position[0], y0 = position[1];
			double x1 = position[2], y1 = position[3];
			double x2 = position[4], y2 = position[5];
			// Compute absolute area using the cross–product formula.
			double area = Math.abs(0.5 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1)));
			if (area < minArea) {
				return 0;
			}

			if (hasHoles) {
				Polygon tri = triangleFromCoords(position);
				if (cache.containsProperly(tri)) {
					return area;
				} else {
					return 0;
				}
			} else {
				// 1. Check that all three vertices are strictly inside the polygon.
				if (!pointInPolygon(x0, y0, polyX, polyY) || !pointInPolygon(x1, y1, polyX, polyY)
						|| !pointInPolygon(x2, y2, polyX, polyY)) {
					return 0;
				}

				// 2. Check for any intersection between triangle edges and polygon boundary.
				// Define the triangle’s three edges.
				if (edgeIntersectsPolygon(x0, y0, x1, y1) || edgeIntersectsPolygon(x1, y1, x2, y2) || edgeIntersectsPolygon(x2, y2, x0, y0)) {
					return 0;
				}

				return area;
			}
		}

		/**
		 * Returns true if the triangle edge from (ax,ay) to (bx,by) intersects a
		 * polygon edge.
		 */
		private boolean edgeIntersectsPolygon(double ax, double ay, double bx, double by) {
			Envelope edgeEnv = new Envelope(ax, bx, ay, by);
			List<?> candidates = polygonEdgeIndex.query(edgeEnv);
			LineSegment triangleEdge = new LineSegment(new Coordinate(ax, ay), new Coordinate(bx, by));
			for (Object obj : candidates) {
				LineSegment polyEdge = (LineSegment) obj;
				if (segmentsIntersect(triangleEdge, polyEdge)) {
					return true;
				}
			}
			return false;
		}

		/**
		 * Determines whether the two given non-collinear line segments intersect.
		 *
		 * <p>
		 * This fast geometric method assumes that the bounding boxes (envelopes) of the
		 * segments already intersect, so no additional envelope intersection test is
		 * performed. It also assumes that no three endpoints are collinear, ensuring
		 * that none of the computed cross products are zero.
		 * </p>
		 *
		 * <p>
		 * The algorithm proceeds by computing the cross products to determine the
		 * relative orientations of the endpoints of each segment with respect to the
		 * line defined by the other segment. Specifically, let segment 1 be defined by
		 * points A and B, and segment 2 by points C and D. The method computes:
		 * </p>
		 *
		 * <ul>
		 * <li>d1 = cross product of (D - C) and (A - C)</li>
		 * <li>d2 = cross product of (D - C) and (B - C)</li>
		 * <li>d3 = cross product of (B - A) and (C - A)</li>
		 * <li>d4 = cross product of (B - A) and (D - A)</li>
		 * </ul>
		 *
		 * <p>
		 * If d1 and d2 have the same sign, then both A and B lie on the same side of
		 * the line through C and D, meaning segment 1 does not cross segment 2.
		 * Similarly, if d3 and d4 have the same sign, segment 2 does not cross segment
		 * 1. Therefore, the segments intersect if and only if A and B lie on opposite
		 * sides of the line through C and D, and C and D lie on opposite sides of the
		 * line through A and B.
		 * </p>
		 *
		 * @param seg1 the first line segment
		 * @param seg2 the second line segment
		 * @return {@code true} if the segments intersect; {@code false} otherwise
		 */
		private boolean segmentsIntersect(LineSegment seg1, LineSegment seg2) {
			// Cache coordinates
			double ax = seg1.p0.x, ay = seg1.p0.y, bx = seg1.p1.x, by = seg1.p1.y;
			double cx = seg2.p0.x, cy = seg2.p0.y, dx = seg2.p1.x, dy = seg2.p1.y;

			// Compute differences for segment 2 (t)
			double rdx = dx - cx, rdy = dy - cy;
			// Compute cross products for seg1 endpoints relative to seg2's line
			double d1 = rdx * (ay - cy) - rdy * (ax - cx);
			double d2 = rdx * (by - cy) - rdy * (bx - cx);

			// If d1 and d2 are of the same sign, seg1 does not cross seg2.
			if ((d1 > 0) == (d2 > 0)) {
				return false;
			}

			// Compute differences for segment 1 (s)
			double sdx = bx - ax, sdy = by - ay;
			// Compute cross products for seg2 endpoints relative to seg1's line
			double d3 = sdx * (cy - ay) - sdy * (cx - ax);
			double d4 = sdx * (dy - ay) - sdy * (dx - ax);

			// The segments intersect if and only if seg2's endpoints lie on opposite
			// sides of seg1's line.
			return (d3 > 0) != (d4 > 0);
		}

		/**
		 * A standard ray-casting point-in-polygon test.
		 *
		 * @param x     the test point x coordinate
		 * @param y     the test point y coordinate
		 * @param polyX the array of polygon vertex x coordinates
		 * @param polyY the array of polygon vertex y coordinates
		 * @return true if the point is inside
		 */
		private boolean pointInPolygon(double x, double y, double[] polyX, double[] polyY) {
			boolean inside = false;
			int n = polyX.length;
			for (int i = 0, j = n - 1; i < n; j = i++) {
				// Check if point is between the y-interval of the edge.
				if (((polyY[i] > y) != (polyY[j] > y)) && (x < (polyX[j] - polyX[i]) * (y - polyY[i]) / (polyY[j] - polyY[i]) + polyX[i])) {
					inside = !inside;
				}
			}
			return inside;
		}
	}

	/**
	 * A Particle subclass for triangle candidates. The dimensionality is 6.
	 */
	private class TriangleCandidate extends Particle {
		public TriangleCandidate() {
			super(6);
		}

		@Override
		public Object selfFactory() {
			return new TriangleCandidate();
		}
	}
}
