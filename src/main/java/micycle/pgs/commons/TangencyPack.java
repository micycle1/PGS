package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.graph.SimpleGraph;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;

import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import micycle.pgs.PGS_Triangulation;
import net.jafama.FastMath;
import processing.core.PVector;

/**
 * Implements a circle packing algorithm described by Collins and Stephenson
 * (2003) to find an arrangement of circles which corresponds to a graph of
 * desired circle tangencies.
 * <p>
 * The algorithm takes a graph (in triangulation form) which specifies a desired
 * pattern of circle tangencies and searches for an arrangement of circle
 * positions and sizes which satisfy that pattern.
 * <p>
 * Given any set of radii, it is possible to compute the angles of the triangles
 * using the law of cosines. The final radii are those for which the angles at
 * any vertex sum to exactly 2π. Thus, the algorithm searches for the radii of
 * the disks by making small incremental updates to the radii, increasing the
 * radius if the angle sum is more than 2π and decreasing the radius of the
 * angle sum is less than 2π.
 * <p>
 * This implementation (specifically circle coordinate placement) is based on an
 * implementation in the <i>packcircles</i> R package.
 * 
 * @author Michael Carleton
 */

public class TangencyPack {

	/*-
	 * Good explanation of this algorithm here:
	 * 		http://www.ams.org/publicoutreach/feature-column/fc-2015-12
	 * Thorough (yet ugly) implementation here (by Stephenson):
	 * 		https://github.com/kensmath/CirclePack/blob/CP-develop/src/rePack/EuclPacker.java
	 */

	private static final double TWO_PI = Math.PI * 2;

	private final IIncrementalTin triangulation;

	private int[][] flowersIds; // Each row contains neighbor IDs for a flower
	private int[] flowerSizes; // Precomputed sizes of each flower
	private double[] delArray; // Precomputed sin(π / n) for each flower
	private double[] factorDelArray; // Precomputed (1 - del)/del for each flower

	private double[] radiiArray;
	private Complex[] placementsArray;
	private List<PVector> circles;
	private Vertex centralVertex;
	private double[] boundaryRadii;

	private List<Vertex> allVertices;
	private Map<Vertex, Integer> vertexToId;
	private List<Vertex> interiorVertices;
	private Map<Vertex, Integer> interiorVertexIndex;
	private int[] interiorVertexIds;

	public TangencyPack(IIncrementalTin triangulation, double boundaryRadii) {
		this.triangulation = triangulation;
		this.boundaryRadii = new double[] { boundaryRadii };
		init();
	}

	public TangencyPack(IIncrementalTin triangulation, List<Double> boundaryRadii) {
		this.triangulation = triangulation;
		this.boundaryRadii = boundaryRadii.stream().mapToDouble(Double::doubleValue).toArray();
		init();
	}

	public TangencyPack(IIncrementalTin triangulation, double[] boundaryRadii) {
		this.triangulation = triangulation;
		this.boundaryRadii = boundaryRadii;
		init();
	}

	private void init() {
		Set<Vertex> perimeterVertices = new HashSet<>();
		triangulation.getPerimeter().forEach(e -> {
			perimeterVertices.add(e.getA());
			perimeterVertices.add(e.getB());
		});

		allVertices = new ArrayList<>(triangulation.getVertices());
		vertexToId = new HashMap<>();
		for (int i = 0; i < allVertices.size(); i++) {
			vertexToId.put(allVertices.get(i), i);
		}

		radiiArray = new double[allVertices.size()];
		flowersIds = new int[allVertices.size()][];
		interiorVertices = new ArrayList<>();
		interiorVertexIndex = new HashMap<>();

		SimpleGraph<Vertex, IQuadEdge> graph = PGS_Triangulation.toTinfourGraph(triangulation);
		NeighborCache<Vertex, IQuadEdge> neighbors = new NeighborCache<>(graph);

		int boundaryIndex = 0;
		PVector meanVertexPos = new PVector();
		int vi = 0;
		for (Vertex v : allVertices) {
			if (perimeterVertices.contains(v)) {
				radiiArray[vertexToId.get(v)] = boundaryRadii[boundaryIndex++ % boundaryRadii.length];
			} else {
				List<Vertex> flower = neighbors.neighborListOf(v);
				flower.sort(new RadialComparator(v));
				int[] flowerIds = new int[flower.size()];
				for (int i = 0; i < flowerIds.length; i++) {
					flowerIds[i] = vertexToId.get(flower.get(i));
				}
				flowersIds[vi++] = flowerIds;
				interiorVertices.add(v);
				int vid = vertexToId.get(v);
				radiiArray[vid] = boundaryRadii[0] / 10;
				meanVertexPos.add((float) v.x, (float) v.y);
			}
		}

		interiorVertexIds = new int[interiorVertices.size()];
		for (int i = 0; i < interiorVertices.size(); i++) {
			Vertex v = interiorVertices.get(i);
			interiorVertexIds[i] = vertexToId.get(v);
			interiorVertexIndex.put(v, i);
		}

		meanVertexPos.div(interiorVertices.size());
		double minDist = Double.MAX_VALUE;
		for (Vertex v : interiorVertices) {
			double dist = v.getDistanceSq(meanVertexPos.x, meanVertexPos.y);
			if (dist < minDist) {
				minDist = dist;
				centralVertex = v;
			}
		}

		// Precompute flower sizes and trigonometric terms
		flowerSizes = new int[interiorVertices.size()];
		delArray = new double[interiorVertices.size()];
		factorDelArray = new double[interiorVertices.size()];

		for (int i = 0; i < interiorVertices.size(); i++) {
			int n = flowersIds[i].length;
			flowerSizes[i] = n;
			double del = FastMath.sin(Math.PI / n); // del = sin(π / n)
			delArray[i] = del;
			factorDelArray[i] = (1 - del) / del; // Precompute (1-del)/del
		}
	}

	public List<PVector> pack() {
		computeRadiiSuperStep();
		computeCenters();
		return circles;
	}

	private void computeRadii() {
		double ttoler;
		// adapt tolerance based on input size. seems sufficient
		if (interiorVertices.size() <= 100) {
			ttoler = 1e-3; // Base case
		} else {
			double exponent = 3.5 + (interiorVertices.size()) / 100.0;
			ttoler = Math.pow(10, -exponent);
		}
		int key = 1;
		double accumErr2 = Double.MAX_VALUE;
		int localPasses = 0;

		while (accumErr2 > ttoler && localPasses < 3 * interiorVertices.size()) {
			Object2DoubleMap<Vertex> R1 = new Object2DoubleOpenHashMap<>(interiorVertices.size());
			for (Vertex v : interiorVertices) {
				R1.put(v, radiiArray[vertexToId.get(v)]);
			}

			double c1 = computeAngleSums();
			c1 = Math.sqrt(c1);
			if (c1 < ttoler) {
				break;
			}

			double factor = c1 / accumErr2;
			if (factor >= 1.0) {
				accumErr2 = c1;
				key = 1;
				localPasses++;
				continue;
			}

			Object2DoubleMap<Vertex> R2 = new Object2DoubleOpenHashMap<>(interiorVertices.size());
			for (Vertex v : interiorVertices) {
				R2.put(v, radiiArray[vertexToId.get(v)]);
			}

			double lmax = 10000;
			for (Vertex v : R1.keySet()) {
				double r1 = R1.getDouble(v);
				double r2 = R2.getDouble(v);
				double rat = r2 - r1;
				if (rat < 0) {
					double tr = -r2 / rat;
					lmax = Math.min(lmax, tr);
				}
			}
			lmax /= 2;

			double lambda;
			if (key == 1) {
				lambda = Math.min(0.75 / (1 - factor), factor);
			} else {
				lambda = factor;
			}
			lambda = Math.min(lambda, lmax);

			for (Vertex v : interiorVertices) {
				int vid = vertexToId.get(v);
				double r1 = R1.getDouble(v);
				double r2 = R2.getDouble(v);
				radiiArray[vid] = r2 + lambda * (r2 - r1);
			}

			// NOTE probably faster to not call computeAngleSums() again. simply use sum
			// from before
//            accumErr2 = computeAngleSums();
//            accumErr2 = Math.sqrt(accumErr2);
			accumErr2 = c1;

			localPasses++;
		}
	}

	private void computeRadiiSuperStep() {
		// Precompute tolerance based on problem size
		final double ttoler = interiorVertices.size() <= 10 ? 1e-2 : Math.pow(10, -(3.5 + interiorVertices.size() / 100.0));

		// Preallocate working arrays once
		final double[] R1 = new double[radiiArray.length];
		final double[] R2 = new double[radiiArray.length];

		int key = 1;
		double accumErr2 = Double.MAX_VALUE;
		int localPasses = 1;
		final int maxPasses = 3 * interiorVertices.size();

		while (accumErr2 > ttoler && localPasses < maxPasses) {
			// Phase 1: Standard iteration
			System.arraycopy(radiiArray, 0, R1, 0, radiiArray.length);
			double c1;
			double factor;

			// Single-pass factor calculation
			do {
				c1 = Math.sqrt(computeAngleSums());
				factor = c1 / accumErr2;
				if (factor >= 1.0) {
					accumErr2 = c1;
					key = 1;
				}
			} while (factor >= 1.0);

			// Phase 2: Super-step preparation
			System.arraycopy(radiiArray, 0, R2, 0, radiiArray.length);

			// Lambda calculation with precomputed values
			final double lambda = calculateLambda(R1, R2, factor, key);

			// Vectorized radius update
			updateRadii(R1, R2, lambda);

			// Error calculation with early exit check
			accumErr2 = Math.sqrt(computeAngleSums());

			// Convergence monitoring
			if (!updateState(R1, R2, c1, accumErr2, factor, key, lambda)) {
				key = (key == 1) ? 2 : 1;
			}

			localPasses++;
		}
	}

	private double calculateLambda(double[] R1, double[] R2, double factor, int key) {
		double lmax = Double.MAX_VALUE;

		// Parallel safe iteration (if needed)
		for (int vid : interiorVertexIds) {
			final double r1 = R1[vid];
			final double r2 = R2[vid];
			final double rat = r2 - r1;
			if (rat < 0) {
				lmax = Math.min(lmax, -r2 / rat);
			}
		}
		lmax /= 2;

		if (key == 1) {
			return Math.min(0.75 / (1 - factor), factor);
		} else {
			return Math.min(factor / (1 - factor), lmax);
		}
	}

	private void updateRadii(double[] R1, double[] R2, double lambda) {
		for (int vid : interiorVertexIds) {
			final double delta = R2[vid] - R1[vid];
			radiiArray[vid] = R2[vid] + lambda * delta;
		}
	}

	private boolean updateState(double[] R1, double[] R2, double c1, double accumErr2, double factor, int key, double lambda) {
		final double pred = FastMath.exp(lambda * FastMath.log(factor));
		final double act = accumErr2 / c1;

		if (act >= 1) {
			System.arraycopy(R2, 0, radiiArray, 0, radiiArray.length);
			return false;
		}
		return act <= pred;
	}

	private double computeAngleSums() {
		double error = 0;
		for (int i = 0; i < interiorVertices.size(); i++) {
			int vId = interiorVertexIds[i];
			int[] flower = flowersIds[i];
			double ra = radiiArray[vId];

			// Compute angle sum for the flower
			double angleSum = 0;
			int n = flowerSizes[i];
			// NOTE inlined angleSum
			for (int j = 0; j < n; j++) {
				int currentId = flower[j];
				int nextId = (j + 1 < n) ? flower[j + 1] : flower[0];
				double b = radiiArray[currentId];
				double c = radiiArray[nextId];
				double bc = b * c;
				double denominator = ra * ra + ra * (b + c) + bc;
				double x = 1 - (2 * bc) / denominator;
				x = Math.max(-1.0, Math.min(1.0, x));
				angleSum += FastMath.acos(x);
			}

			// Update radius using precomputed values
			double factorDel = factorDelArray[i];
			double bet = FastMath.sin(angleSum / (2 * n));
			double r2 = ra * bet * factorDel / (1 - bet);

			radiiArray[vId] = r2;
			error += (angleSum - TWO_PI) * (angleSum - TWO_PI);
		}
		return error;
	}

	private double angleSum(final double rc, final List<Integer> flower) {
		final int n = flower.size();
		double sum = 0.0;
		for (int j = 0; j < n; j++) {
			int currentId = flower.get(j);
			int nextId = flower.get((j + 1) % n);
			sum += tangentAngle(rc, radiiArray[currentId], radiiArray[nextId]);
		}
		return sum;
	}

	private void computeCenters() {
		if (interiorVertices.isEmpty()) {
			circles = new ArrayList<>();
			return;
		}

		placementsArray = new Complex[allVertices.size()];
		int centralId = vertexToId.get(centralVertex);
		placementsArray[centralId] = new Complex(0, 0);

		int centralIndex = interiorVertexIndex.get(centralVertex);
		int[] centralFlower = flowersIds[centralIndex];
		int k2Id = centralFlower[0];
		placementsArray[k2Id] = new Complex(radiiArray[centralId] + radiiArray[k2Id], 0);

		place(centralVertex);
		place(allVertices.get(k2Id));

		circles = new ArrayList<>(allVertices.size());
		for (Vertex v : allVertices) {
			int id = vertexToId.get(v);
			Complex pos = placementsArray[id];
			if (pos != null) {
				circles.add(new PVector((float) pos.getReal(), (float) pos.getImaginary(), (float) radiiArray[id]));
			}
		}
	}

	private void place(Vertex centre) {
		int centreId = vertexToId.get(centre);
		if (!interiorVertexIndex.containsKey(centre)) {
			return;
		}

		int centreIndex = interiorVertexIndex.get(centre);
		int[] flower = flowersIds[centreIndex];
		int nc = flower.length;
		double rcentre = radiiArray[centreId];
		Complex minusI = new Complex(0, -1);

		for (int i = 0; i < nc; i++) {
			int sId = flower[i];
			int tId = flower[(i + 1) % nc];

			if (placementsArray[sId] != null && placementsArray[tId] == null) {
				double rs = radiiArray[sId];
				double rt = radiiArray[tId];
				double theta = tangentAngle(rcentre, rs, rt);

				Complex offset = placementsArray[sId].subtract(placementsArray[centreId]).divide(new Complex(rs + rcentre));
				offset = offset.multiply(minusI.multiply(theta).exp());
				placementsArray[tId] = placementsArray[centreId].add(offset.multiply(rt + rcentre));

				place(allVertices.get(tId));
			}
		}
	}

	private static double tangentAngle(double a, double b, double c) {
		double q = b * c;
		double o = 1 - 2 * q / (a * a + a * (b + c) + q);
		return FastMath.acos(o);
	}

	private static class RadialComparator implements Comparator<Vertex> {
		private Vertex origin;

		public RadialComparator(Vertex origin) {
			this.origin = origin;
		}

		@Override
		public int compare(Vertex o1, Vertex o2) {
			return polarCompare(origin, o1, o2);
		}

		private static int polarCompare(Vertex o, Vertex p, Vertex q) {
			double dxp = p.x - o.x;
			double dyp = p.y - o.y;
			double dxq = q.x - o.x;
			double dyq = q.y - o.y;

			double alph = FastMath.atan2(dyp, dxp);
			double beta = FastMath.atan2(dyq, dxq);
			return Double.compare(alph, beta);
		}
	}
}
