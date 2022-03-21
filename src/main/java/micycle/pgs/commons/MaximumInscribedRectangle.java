package micycle.pgs.commons;

import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;

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

	private static final int SWARM_SIZE = 2500;
	private static final int GENERATIONS = 1000;
	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory();

	private final Swarm swarm;

	public MaximumInscribedRectangle(Polygon polygon) {
		/*
		 * Find MIC to give a best-guess minimum of the MIR's side lengths (we know the
		 * polygon can at least contain a square that is inscribed within the MIC).
		 */
		final MaximumInscribedCircle mic = new MaximumInscribedCircle(polygon, 2);
		final double minSquare = mic.getRadiusLine().getLength() / Math.sqrt(2);

		final FitnessFunction fitnessFunction = new RectangleFitness(PreparedGeometryFactory.prepare(polygon));
		swarm = new Swarm(SWARM_SIZE, new RectangleCandidate(), fitnessFunction);

		final Envelope e = polygon.getEnvelopeInternal();
		final double[] maxPosition = new double[] { e.getMaxX(), e.getMaxY(), e.getWidth(), e.getHeight(), Math.PI };
		swarm.setMaxPosition(maxPosition);
		final double[] minPosition = new double[] { e.getMinX(), e.getMinY(), minSquare, minSquare, -Math.PI };
		swarm.setMinPosition(minPosition);

	}

	/**
	 * Computes the Returns the
	 * 
	 * @return a rectangle polygon
	 */
	public Polygon computeMIR() {
		int i = 0;
		int same = 0;
		double lastFitness = Double.MIN_VALUE;
		while (i++ < GENERATIONS) {
			swarm.evolve();
			i++;
			if (swarm.getBestFitness() == lastFitness) {
				if (same++ > 75) { // evolution is somewhat stuck on local minimum -- exit early
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

	private class RectangleFitness extends FitnessFunction {

		private PreparedGeometry geometry;

		RectangleFitness(PreparedGeometry geometry) {
			this.geometry = geometry;
		}

		@Override
		public double evaluate(double[] position) {
			final double w = position[2];
			final double h = position[3];
			return geometry.containsProperly(rectFromCoords(position)) ? h * w : 0;
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
		PVector dir1 = new PVector((float) Math.cos(a), (float) Math.sin(a));
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
