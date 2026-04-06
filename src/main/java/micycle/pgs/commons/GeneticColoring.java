package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.SplittableRandom;

import org.jgrapht.Graph;
import org.jgrapht.alg.interfaces.VertexColoringAlgorithm;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.util.CollectionUtil;

/**
 * Memetic (GA + local search) solver for the k-Coloring problem on a JGraphT
 * graph.
 * <p>
 * The algorithm searches for a proper coloring by minimizing the number of
 * conflicting edges (adjacent vertices with the same color). It targets a
 * 4-coloring first; if unsuccessful within the generation limit, it repairs the
 * best found assignment into a proper coloring using up to 5 colors.
 *
 * <p>
 * Approach (concise):
 * <ul>
 * <li>Representation: int chromosome of length |V|; gene = color in [0,
 * k-1].</li>
 * <li>Fitness: number of conflicting edges (counted once per edge).</li>
 * <li>Initialization: strong seeds via DSATUR and degree-ordered greedy, plus
 * random individuals.</li>
 * <li>Selection: tournament; elitism preserves the best individuals each
 * generation.</li>
 * <li>Crossover: uniform (favoring agreement) and occasional 2-point
 * crossover.</li>
 * <li>Mutation: conflict-directed; for conflicted vertices choose a valid color
 * using a boolean mask (micro-optimized); otherwise minimize conflicts.</li>
 * <li>Local search: steepest-descent conflict repair applied to each child
 * (memetic step).</li>
 * <li>Diversification: partial re-seeding on stagnation to escape local
 * minima.</li>
 * <li>Precomputation: integer vertex IDs and neighbor arrays for fast inner
 * loops.</li>
 * <li>Termination: success when fitness = 0 or after maxGenerations; fallback
 * repair guarantees a proper coloring.</li>
 * </ul>
 *
 * @param <V> graph vertex type
 * @param <E> graph edge type
 */
public class GeneticColoring<V, E> implements VertexColoringAlgorithm<V> {

	private final int vertexCount;
	private final int maxGenerations;
	private final int populationSize;
	private final int fitnessThreshold;
	private SplittableRandom rand;
	private int colorsCount;

	final java.util.Map<V, Integer> colors;
	private final List<int[]> neighborCache;
	final java.util.Map<V, Integer> vertexIds;

	// Precomputed degrees and degree order
	private final int[] degrees;
	private final int[] degreeOrder;

	public GeneticColoring(Graph<V, E> graph, long seed) {
		this(graph, 200, 100, 4, seed); // larger defaults to give memetic search room
	}

	public GeneticColoring(Graph<V, E> graph, int maxGenerations, int populationSize, int fitnessThreshold, long seed) {
		if (graph == null || maxGenerations < 1 || populationSize < 2) {
			throw new IllegalArgumentException();
		}

		this.vertexCount = graph.vertexSet().size();
		this.maxGenerations = maxGenerations;
		this.populationSize = populationSize;
		this.fitnessThreshold = fitnessThreshold;
		this.rand = new SplittableRandom(seed);

		this.colors = CollectionUtil.newHashMapWithExpectedSize(graph.vertexSet().size());

		// Vertex IDs 0..n-1
		this.vertexIds = new HashMap<>();
		int i = 0;
		for (V v : graph.vertexSet()) {
			vertexIds.put(v, i++);
		}

		// Neighbor cache indexed by ID
		this.neighborCache = new ArrayList<>(vertexCount);
		for (int k = 0; k < vertexCount; k++) {
			neighborCache.add(null);
		}
		final NeighborCache<V, E> nc = new NeighborCache<>(graph);
		for (V v : graph.vertexSet()) {
			int id = vertexIds.get(v);
			int[] neigh = nc.neighborsOf(v).stream().map(vertexIds::get).mapToInt(Integer::intValue).toArray();
			neighborCache.set(id, neigh);
		}

		// Degrees and degree-descending order for greedy seeding
		this.degrees = new int[vertexCount];
		for (int v = 0; v < vertexCount; v++) {
			degrees[v] = neighborCache.get(v).length;
		}
		this.degreeOrder = new int[vertexCount];
		{
			Integer[] idx = new Integer[vertexCount];
			for (int v = 0; v < vertexCount; v++) {
				idx[v] = v;
			}
			Arrays.sort(idx, (a, b) -> Integer.compare(degrees[b], degrees[a]));
			for (int k = 0; k < vertexCount; k++) {
				degreeOrder[k] = idx[k];
			}
		}
	}

	@Override
	public Coloring<V> getColoring() {
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

	private int[] greedyColoringByOrder(int[] order) {
		int[] chrom = new int[vertexCount];
		Arrays.fill(chrom, -1);
		boolean[] used = new boolean[colorsCount];
		int[] candidates = new int[colorsCount];
		for (int idx = 0; idx < vertexCount; idx++) {
			int v = order[idx];
			Arrays.fill(used, false);
			for (int u : neighborsOf(v)) {
				int cu = chrom[u];
				if (cu >= 0) {
					used[cu] = true;
				}
			}
			int n = 0;
			for (int c = 0; c < colorsCount; c++) {
				if (!used[c]) {
					candidates[n++] = c;
				}
			}
			chrom[v] = (n > 0) ? candidates[rand.nextInt(n)] : rand.nextInt(colorsCount);
		}
		return chrom;
	}

	// DSATUR seeding (bounded to colorsCount; if blocked, choose least-conflicting
	// color)
	private int[] dsaturColoring() {
		int[] color = new int[vertexCount];
		Arrays.fill(color, -1);
		int[] sat = new int[vertexCount]; // saturation degree
		int[] usedCounts = new int[vertexCount]; // temporary reuse per vertex
		boolean[][] neighborColors = new boolean[vertexCount][colorsCount];

		for (int colored = 0; colored < vertexCount; colored++) {
			int v = -1, bestSat = -1, bestDeg = -1;
			for (int u = 0; u < vertexCount; u++) {
				if (color[u] != -1) {
					continue;
				}
				int s = sat[u];
				int d = degrees[u];
				if (s > bestSat || (s == bestSat && d > bestDeg)) {
					bestSat = s;
					bestDeg = d;
					v = u;
				}
			}
			Arrays.fill(neighborColors[v], false);
			for (int w : neighborsOf(v)) {
				if (color[w] >= 0) {
					neighborColors[v][color[w]] = true;
				}
			}
			int pick = -1, bestConf = Integer.MAX_VALUE;
			for (int c = 0; c < colorsCount; c++) {
				int conf = 0;
				for (int w : neighborsOf(v)) {
					if (color[w] == c) {
						conf++;
					}
				}
				if (!neighborColors[v][c]) {
					pick = c;
					break;
				} // conflict-free available
				if (conf < bestConf) {
					bestConf = conf;
					pick = c;
				}
			}
			color[v] = pick;
			// update saturation of neighbors
			for (int w : neighborsOf(v)) {
				if (color[w] == -1) {
					if (!neighborColors[w][pick]) {
						neighborColors[w][pick] = true;
						sat[w]++;
					}
				}
			}
		}
		return color;
	}

	private int[] shuffledOrderFrom(int[] base) {
		int[] order = base.clone();
		for (int i = order.length - 1; i > 0; i--) {
			int j = rand.nextInt(i + 1);
			int t = order[i];
			order[i] = order[j];
			order[j] = t;
		}
		return order;
	}

	private class Population {

		private List<Individual> population;
		private int generation = 0;
		private int bestSoFar = Integer.MAX_VALUE;
		private int stagnantGens = 0;

		// memetic parameters
		private final int eliteCount = Math.max(2, populationSize / 5);
		private final int tournamentK = 3;
		private final int localSearchCap = Math.max(200, vertexCount * 10);

		Population() {
			population = new ArrayList<>(populationSize);

			// Strong seeds
			population.add(new Individual(dsaturColoring()));
			population.add(new Individual(greedyColoringByOrder(degreeOrder)));

			int greedySeeds = Math.max(2, populationSize / 10);
			for (int s = 2; s < greedySeeds; s++) {
				int[] order = shuffledOrderFrom(degreeOrder);
				population.add(new Individual(greedyColoringByOrder(order)));
			}

			while (population.size() < populationSize) {
				population.add(new Individual());
			}

			sort();
			bestSoFar = bestFitness();
		}

		public void nextGeneration() {
			// Elitism: keep top eliteCount
			List<Individual> next = new ArrayList<>(populationSize);
			for (int i = 0; i < eliteCount; i++) {
				next.add(population.get(i));
			}

			// Fill the rest with children
			while (next.size() < populationSize) {
				Individual p1 = tournamentSelect();
				Individual p2 = tournamentSelect();
				Individual child = rand.nextBoolean() ? new Individual(new Parents(p1, p2)) : new Individual(twoPointCrossover(p1, p2));
				child.mutate();
				child.localSearch(localSearchCap); // memetic step
				next.add(child);
			}

			population = next;
			sort();
			generation++;

			int best = bestFitness();
			if (best < bestSoFar) {
				bestSoFar = best;
				stagnantGens = 0;
			} else {
				stagnantGens++;
				if (stagnantGens >= 20 && best > 0) {
					diversify();
					stagnantGens = 0;
					bestSoFar = bestFitness();
				}
			}
		}

		private Individual tournamentSelect() {
			Individual best = null;
			for (int i = 0; i < tournamentK; i++) {
				Individual cand = population.get(rand.nextInt(populationSize));
				if (best == null || cand.fitness < best.fitness) {
					best = cand;
				}
			}
			return best;
		}

		private Parents twoPointCrossover(Individual a, Individual b) {
			int i = rand.nextInt(vertexCount);
			int j = rand.nextInt(vertexCount);
			if (i > j) {
				int t = i;
				i = j;
				j = t;
			}
			int[] child = new int[vertexCount];
			System.arraycopy(a.chromosome, 0, child, 0, i);
			System.arraycopy(b.chromosome, i, child, i, j - i);
			System.arraycopy(a.chromosome, j, child, j, vertexCount - j);
			Individual wrapped = new Individual(child);
			return new Parents(wrapped, wrapped); // reuse Individual(Parents) signature via a wrapper
		}

		private void diversify() {
			// re-seed the weakest third
			int start = (int) (populationSize * 0.67);
			for (int i = start; i < populationSize; i++) {
				if (rand.nextDouble() < 0.5) {
					population.set(i, new Individual());
				} else {
					if (rand.nextBoolean()) {
						population.set(i, new Individual(dsaturColoring()));
					} else {
						population.set(i, new Individual(greedyColoringByOrder(shuffledOrderFrom(degreeOrder))));
					}
				}
			}
			sort();
		}

		public int[] bestIndividual() {
			return population.get(0).chromosome;
		}

		public int bestFitness() {
			return population.get(0).fitness;
		}

		public int generation() {
			return generation;
		}

		private void sort() {
			population.sort(Comparator.comparingInt(m -> m.fitness));
		}

		private class Individual {
			private int[] chromosome;
			private int fitness;

			Individual() {
				chromosome = new int[vertexCount];
				for (int i = 0; i < vertexCount; i++) {
					chromosome[i] = rand.nextInt(colorsCount);
				}
				scoreFitness();
			}

			Individual(int[] chromosome) {
				this.chromosome = chromosome.clone();
				scoreFitness();
			}

			// Uniform crossover favoring agreement
			Individual(Parents parents) {
				chromosome = new int[vertexCount];
				final int[] p1 = parents.parent1.chromosome;
				final int[] p2 = parents.parent2.chromosome;
				for (int i = 0; i < vertexCount; i++) {
					int c1 = p1[i], c2 = p2[i];
					chromosome[i] = (c1 == c2) ? c1 : (rand.nextBoolean() ? c1 : c2);
				}
				scoreFitness();
			}

			public void mutate() {
				// Conflict-directed mutation
				boolean[] used = new boolean[colorsCount];
				int[] candidates = new int[colorsCount];

				// early phase vs late phase adapts intensity
				int passes = (bestFitness() > fitnessThreshold) ? 2 : 1;
				for (int pass = 0; pass < passes; pass++) {
					for (int v = 0; v < vertexCount; v++) {
						int cv = chromosome[v];
						boolean conflict = false;
						for (int w : neighborsOf(v)) {
							if (cv == chromosome[w]) {
								conflict = true;
								break;
							}
						}
						if (!conflict) {
							continue;
						}

						// try best color for v
						Arrays.fill(used, false);
						for (int u : neighborsOf(v)) {
							used[chromosome[u]] = true;
						}
						int n = 0;
						for (int c = 0; c < colorsCount; c++) {
							if (!used[c]) {
								candidates[n++] = c;
							}
						}
						if (n > 0) {
							chromosome[v] = candidates[rand.nextInt(n)];
						} else {
							// pick color minimizing conflicts
							int bestC = cv, bestConf = Integer.MAX_VALUE;
							for (int c = 0; c < colorsCount; c++) {
								int conf = 0;
								for (int u : neighborsOf(v)) {
									if (chromosome[u] == c) {
										conf++;
									}
								}
								if (conf < bestConf) {
									bestConf = conf;
									bestC = c;
								}
							}
							chromosome[v] = bestC;
						}
					}
				}
				scoreFitness();
			}

			// Local search: steepest-descent conflict repair with cap
			void localSearch(int cap) {
				int moves = 0;
				boolean improved = true;
				int[] colorCounts = new int[colorsCount];

				while (improved && moves < cap && fitness > 0) {
					improved = false;

					// Build a random permutation of vertices to reduce bias
					int[] order = shuffledOrderFrom(degreeOrder);
					for (int idx = 0; idx < vertexCount && moves < cap; idx++) {
						int v = order[idx];
						int cv = chromosome[v];

						int currentConf = 0;
						for (int w : neighborsOf(v)) {
							if (chromosome[w] == cv) {
								currentConf++;
							}
						}
						if (currentConf == 0) {
							continue;
						}

						Arrays.fill(colorCounts, 0);
						for (int w : neighborsOf(v)) {
							colorCounts[chromosome[w]]++;
						}

						int bestC = cv;
						int bestScore = currentConf;
						for (int c = 0; c < colorsCount; c++) {
							int score = colorCounts[c];
							if (score < bestScore || (score == bestScore && c < bestC)) {
								bestScore = score;
								bestC = c;
							}
						}
						if (bestC != cv) {
							chromosome[v] = bestC;
							moves++;
							// Recompute fitness lazily after a batch; here recompute fully occasionally
							if ((moves & 15) == 0) {
								scoreFitness();
							}
							improved = true;
						}
					}
					// finalize fitness recompute
					scoreFitness();
				}
			}

			private void scoreFitness() {
				int f = 0;
				for (int v = 0; v < vertexCount; v++) {
					int cv = chromosome[v];
					for (int w : neighborsOf(v)) {
						if (w > v && cv == chromosome[w]) {
							f++;
						}
					}
				}
				fitness = f;
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