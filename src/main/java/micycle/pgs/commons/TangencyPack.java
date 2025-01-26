package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
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
	private List<List<Integer>> flowersIds;
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
		flowersIds = new ArrayList<>();
		interiorVertices = new ArrayList<>();
		interiorVertexIndex = new HashMap<>();

		SimpleGraph<Vertex, IQuadEdge> graph = PGS_Triangulation.toTinfourGraph(triangulation);
		NeighborCache<Vertex, IQuadEdge> neighbors = new NeighborCache<>(graph);

		int boundaryIndex = 0;
		PVector meanVertexPos = new PVector();
		for (Vertex v : allVertices) {
			if (perimeterVertices.contains(v)) {
				radiiArray[vertexToId.get(v)] = boundaryRadii[boundaryIndex++ % boundaryRadii.length];
			} else {
				List<Vertex> flower = neighbors.neighborListOf(v);
				flower.sort(new RadialComparator(v));
				List<Integer> flowerIds = new ArrayList<>();
				for (Vertex neighbor : flower) {
					flowerIds.add(vertexToId.get(neighbor));
				}
				flowersIds.add(flowerIds);
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
		double ttoler;
	    if (interiorVertices.size() <= 10) {
	        ttoler = 1e-2; // Base case
	    } else {
	        double exponent = 3.5 + (interiorVertices.size()) / 100.0;
	        ttoler = Math.pow(10, -exponent);
	    }
	    int key = 1; // initial superstep type
	    double accumErr2 = Double.MAX_VALUE;
	    int localPasses = 1;
	    while ((accumErr2 > ttoler && localPasses < 3 * interiorVertices.size())) { // main loop
	        double[] R1 = Arrays.copyOf(radiiArray, radiiArray.length);
	        double c1;
	
	        double factor;
	
	        do { // Make sure factor < 1.0
	            c1 = computeAngleSums();
	            c1 = Math.sqrt(c1);
	
	            factor = c1 / accumErr2;
	            if (factor >= 1.0) {
	                accumErr2 = c1;
	                key = 1;
	            }
	        } while (factor >= 1.0);
	
	        // ================= superstep calculation ====================
	        double[] R2 = Arrays.copyOf(radiiArray, radiiArray.length);
	
	        // find maximum step one can safely take
	        double lmax = 10000;
	        for (int vid : interiorVertexIds) {
	            double r1 = R1[vid];
	            double r2 = R2[vid];
	            double rat = r2 - r1;
	            if (rat < 0) {
	                double tr = (-r2 / rat);
	                lmax = Math.min(lmax, tr); // to keep R>0
	            }
	        }
	        lmax = lmax / 2;
	
	        // do super step
	        double m = 1;
	        int sct = 1;
	        int fct = 2;
	        double lambda;
	        if (key == 1) { // type 1 SS
	            lambda = m * factor;
	            double mmax = 0.75 / (1 - factor); // upper limit on m
	            double mm = (1 + 0.8 / (sct + 1)) * m;
	            m = Math.min(mmax, mm);
	        } else { // type 2 SS
	            double fact0 = 0.0;
	            double ftol = 0.0;
	            if (sct > fct && Math.abs(factor - fact0) < ftol) {
	                lambda = factor / (1 - factor);
	                sct = -1;
	            } else {
	                lambda = factor;
	            }
	        }
	        lambda = Math.min(lambda, lmax);
	
	        // interpolate new radii labels
	        for (int vid : interiorVertexIds) {
	            double r1 = R1[vid];
	            double r2 = R2[vid];
	            radiiArray[vid] = r2 + lambda * (r2 - r1);
	        }
	
	        // check results
	        accumErr2 = computeAngleSums();
	        accumErr2 = Math.sqrt(accumErr2);
	
	        double pred = FastMath.exp(lambda * FastMath.log(factor));
	        double act = accumErr2 / c1;
	        if (act < 1) {
	            if (act > pred) { // not as good as expected: reset
	                if (key == 1) key = 2;
	            }
	        } else { // reset to before superstep
	            System.arraycopy(R2, 0, radiiArray, 0, radiiArray.length);
	            accumErr2 = c1;
	            if (key == 2) key = 1;
	        }
	
	        localPasses++;
	    }
	}

	private double computeAngleSums() {
		double error = 0;
		for (int i = 0; i < interiorVertices.size(); i++) {
			int vId = interiorVertexIds[i];
			List<Integer> flower = flowersIds.get(i);
			double ra = radiiArray[vId];
			double angleSum = angleSum(ra, flower);

			int N = 2 * flower.size();
			double del = FastMath.sin(TWO_PI / N);
			double bet = FastMath.sin(angleSum / N);
			double r2 = ra * bet * (1 - del) / (del * (1 - bet));

			radiiArray[vId] = r2;
			angleSum -= TWO_PI;
			error += angleSum * angleSum;
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
		List<Integer> centralFlower = flowersIds.get(centralIndex);
		int k2Id = centralFlower.get(0);
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
		List<Integer> flower = flowersIds.get(centreIndex);
		int nc = flower.size();
		double rcentre = radiiArray[centreId];
		Complex minusI = new Complex(0, -1);

		for (int i = 0; i < nc; i++) {
			int sId = flower.get(i);
			int tId = flower.get((i + 1) % nc);

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
