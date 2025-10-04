package micycle.pgs.commons;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;
import java.util.SplittableRandom;
import java.util.random.RandomGenerator;

import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.utils.TriangleCollector;

import micycle.pgs.PGS_Triangulation;
import processing.core.PShape;
import processing.core.PVector;

/**
 * @author Michael Carleton
 */
public final class ShapeRandomPointSampler {

	private final double[] ax, ay, bx, by, cx, cy; // triangle vertices
	private final double[] prob; // alias table probs
	private final int[] alias; // alias table indices
	private final int n; // number of triangles
	private RandomGenerator rng;

	public ShapeRandomPointSampler(final PShape shape, final long seed) {
		reseed(seed);

		// Build constrained Delaunay TIN
		final IIncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape);
		final boolean constrained = !tin.getConstraints().isEmpty();

		// Collect valid triangles and their areas
		final List<double[]> tris = new ArrayList<>();
		final List<Double> areas = new ArrayList<>();
		final double eps = 1e-15;

		TriangleCollector.visitSimpleTriangles(tin, (SimpleTriangle tri) -> {
			final IConstraint region = tri.getContainingRegion();
			final boolean inside = !constrained || (region != null && region.definesConstrainedRegion());
			if (!inside)
				return;

			final Vertex A = tri.getVertexA();
			final Vertex B = tri.getVertexB();
			final Vertex C = tri.getVertexC();

			// Robust, unsigned area = 0.5 * |(B-A) x (C-A)|
			final double area2 = (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y);
			final double area = 0.5 * Math.abs(area2);
			if (area > eps && Double.isFinite(area)) {
				tris.add(new double[] { A.x, A.y, B.x, B.y, C.x, C.y });
				areas.add(area);
			}
		});

		this.n = tris.size();
		if (n == 0) {
			throw new IllegalStateException("No valid triangles found for shape; cannot sample.");
		}

		this.ax = new double[n];
		this.ay = new double[n];
		this.bx = new double[n];
		this.by = new double[n];
		this.cx = new double[n];
		this.cy = new double[n];
		for (int i = 0; i < n; i++) {
			double[] t = tris.get(i);
			ax[i] = t[0];
			ay[i] = t[1];
			bx[i] = t[2];
			by[i] = t[3];
			cx[i] = t[4];
			cy[i] = t[5];
		}

		// Build alias table for area weights (Walker’s method)
		this.prob = new double[n];
		this.alias = new int[n];

		final double sumArea = areas.stream().mapToDouble(Double::doubleValue).sum();
		final double[] scaled = new double[n]; // weights scaled by n: w_i * n
		final Deque<Integer> small = new ArrayDeque<>();
		final Deque<Integer> large = new ArrayDeque<>();

		for (int i = 0; i < n; i++) {
			double w = areas.get(i) / sumArea;
			double p = w * n;
			scaled[i] = p;
			if (p < 1.0)
				small.add(i);
			else
				large.add(i);
		}

		while (!small.isEmpty() && !large.isEmpty()) {
			int s = small.removeLast();
			int l = large.removeLast();
			prob[s] = scaled[s];
			alias[s] = l;
			scaled[l] = (scaled[l] + scaled[s]) - 1.0;
			if (scaled[l] < 1.0)
				small.add(l);
			else
				large.add(l);
		}

		// Any leftover get prob=1
		while (!large.isEmpty())
			prob[large.removeLast()] = 1.0;
		while (!small.isEmpty())
			prob[small.removeLast()] = 1.0;
	}

	public void reseed(long seed) {
		this.rng = new SplittableRandom(seed);
	}

	public PVector getRandomPoint() {
		// Alias sampling: pick triangle
		final int i = rng.nextInt(n);
		final double u = rng.nextDouble();
		final int triIndex = (u < prob[i]) ? i : alias[i];

		// Uniform inside triangle using barycentric sqrt-trick
		final double r1 = rng.nextDouble();
		final double r2 = rng.nextDouble();
		final double sqrtR1 = Math.sqrt(r1);
		final double wA = 1.0 - sqrtR1;
		final double wB = sqrtR1 * (1.0 - r2);
		final double wC = sqrtR1 * r2;

		final float x = (float) (wA * ax[triIndex] + wB * bx[triIndex] + wC * cx[triIndex]);
		final float y = (float) (wA * ay[triIndex] + wB * by[triIndex] + wC * cy[triIndex]);
		return new PVector(x, y);
	}

	public List<PVector> getRandomPoints(int count) {
		List<PVector> pts = new ArrayList<>(count);
		for (int i = 0; i < count; i++) {
			pts.add(getRandomPoint());
		}
		return pts;
	}

	public int triangleCount() {
		return n;
	}
}