package micycle.pgs.commons;

import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;

import net.jafama.FastMath;
import net.metaopt.swarm.FitnessFunction;
import net.metaopt.swarm.pso.Particle;
import net.metaopt.swarm.pso.Swarm;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.PointValuePair;

import java.util.Random;

/**
 * Finds an approximate largest area rectangle of arbitrary orientation in a
 * concave polygon via particle swarm optimisation.
 * 
 * @author Michael Carleton
 *
 */
public class MaximumInscribedRectangle {

	private static final int SWARM_SIZE = 2500; // swarm size is tunable
	private static final int GENERATIONS = 1000; // maximum PSO generations
	private static final double ASPECT_WEIGHT = 0.05; // bonus for non-square shapes
	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory();
	private static final Random RANDOM = new Random();

	private final Swarm swarm;
	private final RectangleFitness fitnessFunction; // used by both PSO and local optimization

	public MaximumInscribedRectangle(Polygon polygon) {
		// We compute a lower bound based on the maximum inscribed circle.
		MaximumInscribedCircle mic = new MaximumInscribedCircle(polygon, 2);
		double minSquare = mic.getRadiusLine().getLength() / Math.sqrt(2);

		fitnessFunction = new RectangleFitness(polygon);
		swarm = new ParallelSwarm(SWARM_SIZE, new RectangleCandidate(), fitnessFunction);

		Envelope e = polygon.getEnvelopeInternal();
		double[] maxPosition = new double[] { e.getMaxX(), e.getMaxY(), e.getWidth(), e.getHeight(), Math.PI };
		double[] minPosition = new double[] { e.getMinX(), e.getMinY(), minSquare, minSquare, -Math.PI };

		swarm.setMaxPosition(maxPosition);
		swarm.setMinPosition(minPosition);
		swarm.init();
		initializeSwarm(swarm, e, minSquare);
	}

	private void initializeSwarm(Swarm swarm, Envelope e, double minSquare) {
		Particle[] particles = swarm.getParticles();
		final int particlesPerRegion = SWARM_SIZE / 4; // Four groups
		double minX = e.getMinX();
		double minY = e.getMinY();
		double width = e.getWidth();
		double height = e.getHeight();

		for (int i = 0; i < SWARM_SIZE; i++) {
			double[] position = new double[5];

			if (i < particlesPerRegion) {
				// Group 1: Along envelope edges
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = RANDOM.nextDouble() < 0.5 ? e.getMinY() : e.getMaxY();
				position[2] = width * (0.3 + RANDOM.nextDouble() * 0.7);
				position[3] = height * (0.3 + RANDOM.nextDouble() * 0.7);
				position[4] = RANDOM.nextDouble() * Math.PI - Math.PI / 2;
			} else if (i < 2 * particlesPerRegion) {
				// Group 2: Vertical rectangles
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = minY + RANDOM.nextDouble() * height;
				position[2] = width * (0.1 + RANDOM.nextDouble() * 0.3);
				position[3] = height * (0.7 + RANDOM.nextDouble() * 0.3);
				position[4] = RANDOM.nextDouble() * (Math.PI / 6) - Math.PI / 12;
			} else if (i < 3 * particlesPerRegion) {
				// Group 3: Horizontal rectangles
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = minY + RANDOM.nextDouble() * height;
				position[2] = width * (0.7 + RANDOM.nextDouble() * 0.3);
				position[3] = height * (0.1 + RANDOM.nextDouble() * 0.3);
				position[4] = Math.PI / 2 + RANDOM.nextDouble() * (Math.PI / 6) - Math.PI / 12;
			} else {
				// Group 4: Random exploration
				position[0] = minX + RANDOM.nextDouble() * width;
				position[1] = minY + RANDOM.nextDouble() * height;
				position[2] = minSquare + RANDOM.nextDouble() * (width - minSquare);
				position[3] = minSquare + RANDOM.nextDouble() * (height - minSquare);
				position[4] = RANDOM.nextDouble() * Math.PI - Math.PI / 2;
			}
			particles[i].setPosition(position);
			particles[i].setVelocity(new double[] { 0, 0, 0, 0, 0 });
			particles[i].setBestPosition(position.clone());
		}
	}

	/**
	 * This method runs the swarm evolution and then refines the best candidate
	 * using Apache Commons Math Nelder–Mead simplex optimization.
	 */
	public Polygon computeMIR() {
		int gen = 0, stableCount = 0;
		double lastFitness = Double.MIN_VALUE;

		// Evolve the swarm until convergence criterion or maximum generation reached.
		while (gen++ < GENERATIONS) {
			swarm.evolve();
			if (gen % 50 == 0) {
				reinitializeWorstParticles(swarm);
			}
			double bestFitness = swarm.getBestFitness();
			if (bestFitness == lastFitness) {
				if (stableCount++ > 50) {
					break;
				}
			} else {
				stableCount = 0;
			}
			lastFitness = bestFitness;
		}

		// Retrieve the best candidate from the swarm
		double[] bestCandidate = swarm.getBestPosition().clone();
		// Refine using Apache Commons Math
		bestCandidate = refineCandidate(bestCandidate);
		return rectFromCoords(bestCandidate);
	}

	private void reinitializeWorstParticles(Swarm swarm) {
		Particle[] particles = swarm.getParticles();
		double[] bestPos = swarm.getBestPosition();
		int reinitCount = SWARM_SIZE / 10;
		double[] maxPos = swarm.getMaxPosition();
		double[] minPos = swarm.getMinPosition();
		for (int i = 0; i < reinitCount; i++) {
			double[] newPos = new double[5];
			for (int j = 0; j < 5; j++) {
				double range = (maxPos[j] - minPos[j]) * 0.1;
				newPos[j] = bestPos[j] + (RANDOM.nextDouble() - 0.5) * range;
				// Clamp within bounds.
				newPos[j] = Math.max(minPos[j], Math.min(maxPos[j], newPos[j]));
			}
			particles[i].setPosition(newPos);
			particles[i].setVelocity(new double[] { 0, 0, 0, 0, 0 });
		}
	}

	/**
	 * Uses Apache Commons Math optimizer (Nelder–Mead Simplex) to refine the best
	 * candidate. We define a MultivariateFunction that returns the negative fitness
	 * for minimization.
	 */
	private double[] refineCandidate(double[] candidate) {
		// We want to maximize the fitness, so we minimize negative fitness.
		MultivariateFunction objective = new MultivariateFunction() {
			public double value(double[] point) {
				// Return the negative fitness
				return -fitnessFunction.evaluate(point);
			}
		};

		// No additional bounds mapping is used here since our candidate is already
		// inside the feasible region.
		// However, you can wrap this with a MultivariateFunctionMappingAdapter if
		// needed.
		MultivariateOptimizer optimizer = new SimplexOptimizer(1e-8, 1e-8);
		NelderMeadSimplex simplex = new NelderMeadSimplex(candidate.length);
		try {
			PointValuePair result = optimizer.optimize(new MaxEval(1000), new ObjectiveFunction(objective), GoalType.MINIMIZE, new InitialGuess(candidate),
					simplex);
			return result.getPoint();
		} catch (Exception e) {
			// If the optimizer fails, return the original candidate.
			return candidate;
		}
	}

	/**
	 * Builds a rectangle Polygon based on a candidate parameter vector.
	 *
	 * The candidate parameters are: [x, y, width, height, angle] and the rectangle
	 * corners are computed from these values.
	 */
	private static Polygon rectFromCoords(double[] position) {
		double x = position[0];
		double y = position[1];
		double w = position[2];
		double h = position[3];
		double a = position[4];

		double cosA = FastMath.cos(a);
		double sinA = FastMath.sin(a);
		double dx = cosA * w;
		double dy = sinA * w;
		double dx2 = -sinA * h;
		double dy2 = cosA * h;

		Coordinate[] coords = new Coordinate[5];
		coords[0] = new Coordinate(x, y);
		coords[1] = new Coordinate(x + dx, y + dy);
		coords[2] = new Coordinate(x + dx + dx2, y + dy + dy2);
		coords[3] = new Coordinate(x + dx2, y + dy2);
		coords[4] = coords[0];

		return GEOM_FACTORY.createPolygon(coords);
	}

	/**
	 * Evaluates each candidate rectangle. A candidate rectangle (constructed from
	 * position) is given a fitness equal to its area, with a small bonus for
	 * non-square aspect ratios. If the rectangle is not properly contained in the
	 * geometry, it returns 0.
	 */
	private class RectangleFitness extends FitnessFunction {
		private final PreparedGeometry geometry;

		RectangleFitness(Geometry geom) {
			this.geometry = PreparedGeometryFactory.prepare(geom);
		}

		@Override
		public double evaluate(double[] position) {
			double w = position[2];
			double h = position[3];
			Polygon rect = rectFromCoords(position);
			if (!geometry.containsProperly(rect)) {
				return 0;
			}
			double aspectRatio = (w > h) ? w / h : h / w;
			return w * h * (1 + ASPECT_WEIGHT * (aspectRatio - 1));
		}
	}

	// Particle candidate for the swarm; note the dimension is 5.
	private class RectangleCandidate extends Particle {
		public RectangleCandidate() {
			super(5);
		}

		@Override
		public Object selfFactory() {
			return new RectangleCandidate();
		}
	}
}
