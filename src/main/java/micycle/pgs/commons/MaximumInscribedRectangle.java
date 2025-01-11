package micycle.pgs.commons;

import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;

import net.jafama.FastMath;
import net.metaopt.swarm.FitnessFunction;
import net.metaopt.swarm.pso.Particle;
import net.metaopt.swarm.pso.Swarm;
import processing.core.PVector;

/**
 * Finds an approximate largest area rectangle of arbitrary orientation in a
 * polygon via particle swarm optimisation.
 * 
 * @author Michael Carleton
 *
 */
public class MaximumInscribedRectangle {

	private static final int SWARM_SIZE = 2500; // Reduced from 2500
	private static final int GENERATIONS = 1000; // Reduced from 1000
	private static final double ASPECT_WEIGHT = 0.05; // Small bonus for non-square shapes
	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory();
	private final Swarm swarm;

	public MaximumInscribedRectangle(Polygon polygon) {
		final MaximumInscribedCircle mic = new MaximumInscribedCircle(polygon, 2);
		final double minSquare = mic.getRadiusLine().getLength() / Math.sqrt(2);

		// Create multiple sub-swarms with different initial conditions
		final FitnessFunction fitnessFunction = new RectangleFitness(polygon);
		swarm = new Swarm(SWARM_SIZE, new RectangleCandidate(), fitnessFunction);

		final Envelope e = polygon.getEnvelopeInternal();
		final double[] maxPosition = new double[] { e.getMaxX(), e.getMaxY(), e.getWidth(), e.getHeight(), Math.PI };
		final double[] minPosition = new double[] { e.getMinX(), e.getMinY(), minSquare, minSquare, -Math.PI };

		swarm.setMaxPosition(maxPosition);
		swarm.setMinPosition(minPosition);

		// Initialize particles in promising regions
		swarm.init();
		initializeSwarm(swarm, e, minSquare);
	}

	private void initializeSwarm(Swarm swarm, Envelope e, double minSquare) {
		Particle[] particles = swarm.getParticles();
		int particlesPerRegion = SWARM_SIZE / 4;

		for (int i = 0; i < SWARM_SIZE; i++) {
			double[] position = new double[5];

			if (i < particlesPerRegion) {
				// Group 1: Initialize along envelope edges
				position[0] = e.getMinX() + Math.random() * e.getWidth();
				position[1] = Math.random() < 0.5 ? e.getMinY() : e.getMaxY();
				position[2] = e.getWidth() * (0.3 + Math.random() * 0.7);
				position[3] = e.getHeight() * (0.3 + Math.random() * 0.7);
				position[4] = Math.random() * Math.PI - Math.PI / 2;
			} else if (i < 2 * particlesPerRegion) {
				// Group 2: Try vertical rectangles
				position[0] = e.getMinX() + Math.random() * e.getWidth();
				position[1] = e.getMinY() + Math.random() * e.getHeight();
				position[2] = e.getWidth() * (0.1 + Math.random() * 0.3);
				position[3] = e.getHeight() * (0.7 + Math.random() * 0.3);
				position[4] = Math.random() * Math.PI / 6 - Math.PI / 12;
			} else if (i < 3 * particlesPerRegion) {
				// Group 3: Try horizontal rectangles
				position[0] = e.getMinX() + Math.random() * e.getWidth();
				position[1] = e.getMinY() + Math.random() * e.getHeight();
				position[2] = e.getWidth() * (0.7 + Math.random() * 0.3);
				position[3] = e.getHeight() * (0.1 + Math.random() * 0.3);
				position[4] = Math.PI / 2 + Math.random() * Math.PI / 6 - Math.PI / 12;
			} else {
				// Group 4: Random exploration
				position[0] = e.getMinX() + Math.random() * e.getWidth();
				position[1] = e.getMinY() + Math.random() * e.getHeight();
				position[2] = minSquare + Math.random() * (e.getWidth() - minSquare);
				position[3] = minSquare + Math.random() * (e.getHeight() - minSquare);
				position[4] = Math.random() * Math.PI - Math.PI / 2;
			}

			particles[i].setPosition(position);
			particles[i].setVelocity(new double[] { 0, 0, 0, 0, 0 });
			particles[i].setBestPosition(position.clone());
		}
	}

	public Polygon computeMIR() {
		int i = 0;
		int same = 0;
		double lastFitness = Double.MIN_VALUE;

		while (i++ < GENERATIONS) {
			swarm.evolve();

			// Periodically reinitialize worst performing particles
			if (i % 50 == 0) {
				reinitializeWorstParticles(swarm);
			}

			if (swarm.getBestFitness() == lastFitness) {
				if (same++ > 50) {
					break;
				}
			} else {
				same = 0;
			}
			lastFitness = swarm.getBestFitness();
		}

		return getBestRectangleResult(swarm);
	}

	private Polygon getBestRectangleResult(Swarm swarm) {
		return rectFromCoords(swarm.getBestPosition());
	}

	private void reinitializeWorstParticles(Swarm swarm) {
		Particle[] particles = swarm.getParticles();
		double[] bestPos = swarm.getBestPosition();

		// Reinitialize bottom 10% of particles
		int reinitCount = SWARM_SIZE / 10;
		for (int i = 0; i < reinitCount; i++) {
			double[] newPos = new double[5];
			// Explore around current best solution
			for (int j = 0; j < 5; j++) {
				double range = (swarm.getMaxPosition()[j] - swarm.getMinPosition()[j]) * 0.1;
				newPos[j] = bestPos[j] + (Math.random() - 0.5) * range;
				// Ensure bounds
				newPos[j] = Math.max(swarm.getMinPosition()[j], Math.min(swarm.getMaxPosition()[j], newPos[j]));
			}
			particles[i].setPosition(newPos);
			particles[i].setVelocity(new double[] { 0, 0, 0, 0, 0 });
		}
	}

	private class RectangleFitness extends FitnessFunction {
		private PreparedGeometry geometry;

		RectangleFitness(Geometry geometry) {
			this.geometry = PreparedGeometryFactory.prepare(geometry);
		}

		@Override
		public double evaluate(double[] position) {
			final double w = position[2];
			final double h = position[3];
			if (!geometry.containsProperly(rectFromCoords(position))) {
				return 0;
			}
			// Add small bonus for non-square shapes to encourage elongated rectangles
			double aspectRatio = Math.max(w / h, h / w);
			return h * w * (1 + ASPECT_WEIGHT * (aspectRatio - 1));
		}
	}

	private class RectangleCandidate extends Particle {
		public RectangleCandidate() {
			super(5);
		}

		@Override
		public Object selfFactory() {
			return new RectangleCandidate();
		}
	}

	private static final Polygon rectFromCoords(double[] position) {
		final double x = position[0];
		final double y = position[1];
		final double w = position[2];
		final double h = position[3];
		final double a = position[4];

		Coordinate[] coords = new Coordinate[5];
		PVector base = new PVector((float) x, (float) y);
		PVector dir1 = new PVector((float) FastMath.cos(a), (float) FastMath.sin(a));
		PVector dir2 = dir1.copy().rotate((float) Math.PI * 0.5f);
		PVector base2 = base.copy().add(dir1.copy().mult((float) w));
		coords[0] = coordFromPVector(base);
		coords[1] = coordFromPVector(base2);
		dir2.mult((float) h);
		coords[2] = coordFromPVector(base2.add(dir2));
		coords[3] = coordFromPVector(base.add(dir2));
		coords[4] = coords[0];

		return GEOM_FACTORY.createPolygon(coords);
	}

	private static final Coordinate coordFromPVector(PVector p) {
		return new Coordinate(p.x, p.y);
	}

}
