package micycle.pgs.commons;

import java.util.stream.IntStream;

import net.metaopt.swarm.FitnessFunction;
import net.metaopt.swarm.pso.Particle;
import net.metaopt.swarm.pso.Swarm;

/**
 * A particle swarm that evaluates particle fitness in parallel for improved
 * performance.
 *
 * @author Michael Carleton
 */
public class ParallelSwarm extends Swarm {

	public ParallelSwarm(int numberOfParticles, Particle sampleParticle, FitnessFunction fitnessFunction) {
		super(numberOfParticles, sampleParticle, fitnessFunction);
	}

	@Override
	public void evaluate() {
		// Initialize bestFitness on first run.
		if (Double.isNaN(bestFitness)) {
			bestFitness = (fitnessFunction.isMaximize() ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY);
			bestParticleIndex = -1;
		}

		final int len = particles.length;
		final double[] fitnesses = new double[len];

		// Evaluate each particle's fitness in parallel.
		IntStream.range(0, len).parallel().forEach(i -> {
			double fit = fitnessFunction.evaluate(particles[i]);
			fitnesses[i] = fit;
		});

		// Now process the results sequentially to update the global best and
		// neighborhood.
		for (int i = 0; i < len; i++) {
			numEvaluations++;
			double fit = fitnesses[i];

			// Update 'best global' position if this particle is better.
			if (fitnessFunction.isBetterThan(bestFitness, fit)) {
				bestFitness = fit;
				bestParticleIndex = i;
				if (bestPosition == null) {
					bestPosition = new double[sampleParticle.getDimension()];
				}
				particles[i].copyPosition(bestPosition);
			}

			// Update 'best neighborhood' information.
			if (neighborhood != null) {
				neighborhood.update(this, particles[i]);
			}

		}
	}

}
