package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.SplittableRandom;

import org.jgrapht.Graph;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.util.CollectionUtil;

/**
 * Finds a solution to a graph coloring using a genetic algorithm.
 * <p>
 * This class implements the technique described in <i>Genetic Algorithm Applied
 * to the Graph Coloring Problem</i> by <i>Musa M. Hindi and Roman V.
 * Yampolskiy</i>.
 * <p>
 * The genetic algorithm process continues until it either finds a solution
 * (i.e. 0 conflicts between adjacent vertices) or the algorithm has been run
 * for the predefined number of generations.
 * <p>
 * The goal of the algorithm is to improve the fitness of the population (a
 * coloring) by mating its fittest individuals to produce superior offspring
 * that offer a better solution to the problem. This process continues until a
 * terminating condition is reached which could be simply that the total number
 * of generations has been run or any other parameter like non-improvement of
 * fitness over a certain number of generations or that a solution for the
 * problem has been found.
 * 
 * @author Soroush Javadi
 * @author Refactored for JGraphT by Michael Carleton
 *
 * @param <V> the graph vertex type
 * @param <E> the graph edge type
 */
public class GeneticColoring<V, E> implements VertexColoringAlgorithm<V> {

	private final int vertexCount;
	private final int maxGenerations;
	private final int populationSize;
	// fitness threshold for choosing a parent selection and mutation algorithm
	private final int fitnessThreshold;
	private SplittableRandom rand;
	private int colorsCount;

	final Map<V, Integer> colors;
	private final List<int[]> neighborCache;
	final Map<V, Integer> vertexIds;

	/**
	 * Creates with a population size of 50; "the value was chosen after testing a
	 * number of different population sizes. The value 50 was the least value that
	 * produced the desired results".
	 * 
	 * @param graph
	 */
	public GeneticColoring(Graph<V, E> graph) {
		this(graph, 100, 50, 4);
	}

	public GeneticColoring(Graph<V, E> graph, int maxGenerations, int populationSize, int fitnessThreshold) {
		if (graph == null || maxGenerations < 1 || populationSize < 2) {
			throw new IllegalArgumentException();
		}

		this.vertexCount = graph.vertexSet().size();
		this.maxGenerations = maxGenerations;
		this.populationSize = populationSize;
		this.fitnessThreshold = fitnessThreshold;
		this.rand = new SplittableRandom(); // NOTE unseeded

		this.colors = CollectionUtil.newHashMapWithExpectedSize(graph.vertexSet().size());

		this.vertexIds = new HashMap<>();
		int i = 0;
		for (V v : graph.vertexSet()) {
			vertexIds.put(v, i);
		}

		this.neighborCache = new ArrayList<>();
		final NeighborCache<V, E> neighborCacheJT = new NeighborCache<>(graph);
		for (V v : graph.vertexSet()) {
			neighborCache.add(neighborCacheJT.neighborsOf(v).stream().map(vertexIds::get) // map vertex objects to Integer ID
					.mapToInt(Integer::intValue) // Integer -> int
					.toArray());
		}
	}

	@Override
	public Coloring<V> getColoring() {
		/*
		 * For PGS, attempt to find a 4-color (optimal) solution within the
		 * maxGenerations threshold. If a solution is not found, find a 5-color solution
		 * instead (much more attainable/faster).
		 */
		if (!getSolution(4)) {
			getSolution(5);
		}
		return new ColoringImpl<>(colors, colorsCount);
	}

	private boolean getSolution(int colors) {
		this.colorsCount = colors;
		Population population = new Population();
		while (population.bestFitness() != 0 && population.generation() < maxGenerations) {
			population.nextGeneration();
		}
		if (population.bestFitness() == 0) {
			vertexIds.forEach((v, i) -> this.colors.put(v, population.bestIndividual()[i]));
			return true;
		}
		return false;
	}

	private int[] neighborsOf(int v) {
		return neighborCache.get(v);
	}

	private class Population {

		private List<Individual> population; // NOTE heap suitable?
		private int generation = 0;

		Population() {
			population = new ArrayList<>(populationSize);
			for (int i = 0; i < populationSize; i++) {
				population.add(new Individual());
			}
			sort();
		}

		/**
		 * With each generation the bottom performing half of the population is removed
		 * and new randomly generated chromosomes are added.
		 */
		public void nextGeneration() {
			final int halfSize = populationSize / 2;
			List<Individual> children = new ArrayList<>(halfSize);
			for (int i = 0; i < halfSize; i++) {
				Parents parents = selectParents();
				Individual child = new Individual(parents);
				child.mutate();
				children.add(child);
			}
			for (int i = 0; i < halfSize; i++) {
				population.set(populationSize - i - 1, children.get(i));
			}
			sort();
			generation++;
		}

		/**
		 * 
		 * @return the best/fittest color assignment
		 */
		public int[] bestIndividual() {
			return population.get(0).chromosome;
		}

		public int bestFitness() {
			return population.get(0).fitness;
		}

		public int generation() {
			return generation;
		}

		/**
		 * Choosing a parent selection and mutation method depends on the state of the
		 * population and how close it is to finding a solution.
		 * <p>
		 * If the best fitness is greater than {@code fitnessThreshold} then
		 * parentSelection1 and mutation1 are used. Otherwise, parentSelection2 and
		 * mutation2 are used. This alteration is the result of experimenting with the
		 * different data sets. It was observed that when the best fitness score is low
		 * (i.e. approaching an optimum) the usage of parent selection 2 (which copies
		 * the best chromosome as the new child) along with mutation2 (which randomly
		 * selects a color for the violating vertex) results in a solution more often
		 * and more quickly than using the other two respective methods.
		 */
		private Parents selectParents() {
			return bestFitness() > fitnessThreshold ? selectParents1() : selectParents2();
		}

		private Parents selectParents1() {
			Individual tempParent1, tempParent2, parent1, parent2;
			tempParent1 = population.get(rand.nextInt(populationSize));
			do {
				tempParent2 = population.get(rand.nextInt(populationSize));
			} while (tempParent1 == tempParent2);
			parent1 = (tempParent1.fitness > tempParent2.fitness ? tempParent2 : tempParent1);
			do {
				tempParent1 = population.get(rand.nextInt(populationSize));
				do {
					tempParent2 = population.get(rand.nextInt(populationSize));
				} while (tempParent1 == tempParent2);
				parent2 = (tempParent1.fitness > tempParent2.fitness ? tempParent2 : tempParent1);
			} while (parent1 == parent2);
			return new Parents(parent1, parent2);
		}

		private Parents selectParents2() {
			return new Parents(population.get(0), population.get(1));
		}

		private void sort() {
			population.sort(Comparator.comparingInt(m -> m.fitness));
		}

		/**
		 * A candidate graph coloring.
		 */
		private class Individual {
			/**
			 * each element of chromosome represents a color of a vertex
			 */
			private int[] chromosome;
			/**
			 * fitness is defined as the number of 'bad' edges, i.e., edges connecting two
			 * vertices with the same color
			 */
			private int fitness;

			/**
			 * Instantiate a random individual.
			 */
			Individual() {
				chromosome = new int[vertexCount];
				for (int i = 0; i < vertexCount; i++) {
					chromosome[i] = rand.nextInt(colorsCount);
				}
				scoreFitness();
			}

			// crossover
			Individual(Parents parents) {
				chromosome = new int[vertexCount];
				final int crosspoint = rand.nextInt(vertexCount);
				System.arraycopy(parents.parent1.chromosome, 0, chromosome, 0, crosspoint);
				System.arraycopy(parents.parent2.chromosome, crosspoint, chromosome, crosspoint, vertexCount - crosspoint);
				scoreFitness();
			}

			public void mutate() {
				if (bestFitness() > fitnessThreshold) {
					mutate1();
				} else {
					mutate2();
				}
			}

			private void mutate1() {
				for (int v = 0; v < vertexCount; v++) {
					for (int w : neighborsOf(v)) {
						if (chromosome[v] == chromosome[w]) {
							HashSet<Integer> validColors = new HashSet<>();
							for (int c = 0; c < colorsCount; c++) {
								validColors.add(c);
							}
							for (int u : neighborsOf(v)) {
								validColors.remove(chromosome[u]);
							}
							if (!validColors.isEmpty()) {
								chromosome[v] = (int) validColors.toArray()[rand.nextInt(validColors.size())];
							}
							break;
						}
					}
				}
				scoreFitness();
			}

			private void mutate2() {
				for (int v = 0; v < vertexCount; v++) {
					for (int w : neighborsOf(v)) {
						if (chromosome[v] == chromosome[w]) {
							chromosome[v] = rand.nextInt(colorsCount);
							break;
						}
					}
				}
				scoreFitness();
			}

			/**
			 * A bad edge is defined as an edge connecting two vertices that have the same
			 * color. The number of bad edges is the fitness score for the chromosome
			 * (higher number is worse fitness).
			 */
			private void scoreFitness() {
				int f = 0;
				for (int v = 0; v < vertexCount; v++) {
					for (int w : neighborsOf(v)) {
						if (chromosome[v] == chromosome[w]) {
							f++;
						}
					}
				}
				/*
				 * Divide by 2 to account for double counting. Doesn't really matter if we're
				 * simply using fitness score to sort individuals.
				 */
				fitness = f / 2;
			}
		}

		private class Parents {
			public final Individual parent1;
			public final Individual parent2;

			public Parents(Individual parent1, Individual parent2) {
				this.parent1 = parent1;
				this.parent2 = parent2;
			}
		}
	}
}