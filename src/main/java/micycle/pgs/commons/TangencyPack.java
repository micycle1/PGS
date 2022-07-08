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

import com.github.xmunkki.fixpoint.Fixed64;

import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import micycle.pgs.PGS_Triangulation;
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

	private static final double TOLERANCE = 1 + 1e-8;

	private final IIncrementalTin triangulation;
	/**
	 * Maps a vertex to a list of it neighbouring vertices; the neighbour list is
	 * ordered radially around the given vertex.
	 */
	private Map<Vertex, List<Vertex>> flowers;
	/**
	 * The radius of each circle (including boundary circles).
	 */
	private Object2DoubleMap<Vertex> radii;
	private double[] boundaryRadii;
	private Map<Vertex, Complex> placements = new HashMap<>();
	private List<PVector> circles;
	private Vertex centralVertex;

	/**
	 * Creates a circle packing using tangancies specified by a triangulation.
	 * 
	 * @param triangulation Pattern of tangencies; vertices connected by an edge in
	 *                      the triangulation represent tangent circles in the
	 *                      packing
	 * @param boundaryRadii Radii of the circles (same for every circle) associated
	 *                      with the boundary/perimeter vertices of the
	 *                      triangulation
	 */
	public TangencyPack(IIncrementalTin triangulation, double boundaryRadii) {
		this.triangulation = triangulation;
		this.boundaryRadii = new double[] { boundaryRadii };
		init();
	}

	/**
	 * Creates a circle packing using tangancies specified by a triangulation.
	 * 
	 * @param triangulation Pattern of tangencies; vertices connected by an edge in
	 *                      the triangulation represent tangent circles in the
	 *                      packing
	 * @param boundaryRadii List of radii values of the circles associated with the
	 *                      boundary/perimeter vertices of the triangulation. The
	 *                      list may have fewer radii than the number of boundary
	 *                      vertices; in this case, boundary radii will wrap around
	 *                      the list
	 */
	public TangencyPack(IIncrementalTin triangulation, List<Double> boundaryRadii) {
		this.triangulation = triangulation;
		this.boundaryRadii = boundaryRadii.stream().mapToDouble(Double::doubleValue).toArray();
		init();
	}

	/**
	 * Creates a circle packing using tangancies specified by a triangulation.
	 * 
	 * @param triangulation Pattern of tangencies; vertices connected by an edge in
	 *                      the triangulation represent tangent circles in the
	 *                      packing
	 * @param boundaryRadii Array of radii values of the circles associated with the
	 *                      boundary/perimeter vertices of the triangulation. The
	 *                      list may have fewer radii than the number of boundary
	 *                      vertices; in this case, boundary radii will wrap around
	 *                      the list
	 */
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

		flowers = new HashMap<>();
		radii = new Object2DoubleOpenHashMap<>(triangulation.getVertices().size());

		SimpleGraph<Vertex, IQuadEdge> graph = PGS_Triangulation.toGraph(triangulation, true);
		NeighborCache<Vertex, IQuadEdge> neighbors = new NeighborCache<>(graph);

		final PVector meanVertexPos = new PVector();
		int index = 0;
		for (Vertex v : graph.vertexSet()) {
			if (perimeterVertices.contains(v)) {
				radii.put(v, boundaryRadii[index++ % boundaryRadii.length]);
			} else {
				List<Vertex> flower = neighbors.neighborListOf(v);
				RadialComparator c = new RadialComparator(v);
				flower.sort(c);
				flowers.put(v, flower);
				radii.put(v, boundaryRadii[0] / 10);
				meanVertexPos.add((float) v.x, (float) v.y);
			}
		}
		
		// pick a rather central vertex, so output is same on identical input
		meanVertexPos.div(flowers.size());
		double maxDist = Double.MAX_VALUE;
		for (Vertex v : flowers.keySet()) {
			double dist = v.getDistanceSq(meanVertexPos.x, meanVertexPos.y);
			if 	(dist < maxDist) {
				maxDist = dist;
				centralVertex = v;
			}	
		}
	}

	/**
	 * Computes and returns a circle packing for the configuration of tangencies
	 * given by the triangulation.
	 * 
	 * @return a list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 */
	public List<PVector> pack() {
		computeRadii();
		computeCenters();
		return circles;
	}

	/**
	 * Find radii of circles using numerical relaxation. Circle radii converge
	 * rapidly to a unique fixed point for which all flower angles are are within a
	 * desired tolerance of 2π, at which point iteration stops and a packing is
	 * found.
	 */
	private void computeRadii() {
		double lastChange = TOLERANCE + 1;
		while (lastChange > TOLERANCE) {
			lastChange = 1.0;
			for (Vertex v : flowers.keySet()) {
				double theta = flower(v);
				lastChange = Math.max(lastChange, theta);
			}
		}
	}

	/**
	 * Determine the centers of the circles using radii of the interior circles.
	 * 
	 * @return
	 */
	private void computeCenters() {
		if (flowers.size() > 0) {
			placements = new HashMap<>();

			Vertex k1 = centralVertex; // pick one internal circle
			placements.put(k1, new Complex(0, 0)); // place it at the origin

			Vertex k2 = flowers.get(k1).get(0); // pick one of its neighbors
			placements.put(k2, new Complex(radii.getDouble(k1) + radii.getDouble(k2))); // place it on the real axis
			place(k1); // recursively place the rest
			place(k2);
		}

		circles = new ArrayList<>(radii.size());
		placements.forEach(
				(v, pos) -> circles.add(new PVector((float) pos.getReal(), (float) pos.getImaginary(), (float) radii.getDouble(v))));
	}

	/**
	 * Compute the angle sum for the petals surrounding the given vertex and update
	 * the radius of the vertex such that the angle sum would equal 2π.
	 * 
	 * @param center target vertex
	 * @return a measure of the error (difference between target angle sum (2π) and
	 *         the actual angle sum
	 */
	private double flower(final Vertex center) {
		List<Vertex> flower = flowers.get(center);
		final int n = flower.size();
		final double rc = radii.getDouble(center);
		double sum = 0.0;
		for (int i = 0; i < n; i++) {
			int j = i + 1 == n ? 0 : i + 1;
			sum += tangentAngle(rc, radii.getDouble(flower.get(i)), radii.getDouble(flower.get(j)));
		}

		double hat = rc / (1.0 / sin(sum / (2 * n)) - 1);
		double newrad = hat * (1.0 / sin(Math.PI / n) - 1);
		radii.put(center, newrad);

		return Math.max(newrad / rc, rc / newrad);
	}

	/**
	 * Recursively determine centers of all circles surrounding a given vertex.
	 */
	private void place(final Vertex centre) {

		if (!flowers.containsKey(centre)) {
			return; // boundary vertex
		}

		List<Vertex> flower = flowers.get(centre);
		final int nc = flower.size();
		final double rcentre = radii.getDouble(centre);

		final Complex minusI = new Complex(0.0, -1);

		for (int i = -nc; i < nc - 1; i++) {
			int ks = i < 0 ? nc + i : i;
			Vertex s = flower.get(ks);
			double rs = radii.getDouble(s);

			int kt = ks + 1 < nc ? ks + 1 : 0;
			Vertex t = flower.get(kt);
			double rt = radii.getDouble(t);

			if (placements.containsKey(s) && !placements.containsKey(t)) {
				double theta = tangentAngle(rcentre, rs, rt);
				Complex offset = (placements.get(s).subtract(placements.get(centre))).divide(new Complex(rs + rcentre));
				offset = offset.multiply(minusI.multiply(theta).exp());
				placements.put(t, placements.get(centre).add(offset.multiply(rt + rcentre)));

				place(t);
			}
		}
	}

	/**
	 * Computes the angle that circles y and z make with circle x (angle yxz). The
	 * circles are given by their radii and are mutually tangent.
	 * 
	 * @param rx radius of circle x, the circle of interest
	 * @param ry radius of circle y, a "petal" circle
	 * @param rz radius of circle z, a "petal" circle
	 * @return angle of yxz
	 */
	private static double tangentAngle(final double rx, final double ry, final double rz) {
		final double x = (ry * rz) / ((rx + ry) * (rx + rz));
//		return 2 * Fixed64.ToDouble(Fixed64.Asin(Fixed64.FromDouble(Math.sqrt(x)))); // slower
		return 2 * Fixed64.ToDouble(Fixed64.Atan(Fixed64.FromDouble(Math.sqrt(x) / (Math.sqrt(1 - x)))));
	}

	/**
	 * @deprecated slower than chosen method
	 */
	@Deprecated
	private static double tangentAngle2(double a, double b, double c) {
		final double q = b * c;
		final double o = 1 - 2 * q / (a * a + a * (b + c) + q);
		return Math.acos(o);
	}

	private static double sin(double x) {
		return Fixed64.ToDouble(Fixed64.Sin(Fixed64.FromDouble(x)));
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

		/**
		 * Given two points p and q compare them with respect to their radial ordering
		 * about point o. First checks radial ordering.
		 *
		 * @param o the origin
		 * @param p a point
		 * @param q another point
		 * @return -1, 0 or 1 depending on whether angle p is less than, equal to or
		 *         greater than angle q
		 */
		private static int polarCompare(Vertex o, Vertex p, Vertex q) {
			double dxp = p.x - o.x;
			double dyp = p.y - o.y;
			double dxq = q.x - o.x;
			double dyq = q.y - o.y;

			int result = 0;
			double alph = atan2Quick(dxp, dyp);
			double beta = atan2Quick(dxq, dyq);
			if (alph < beta) {
				result = -1;
			}
			if (alph > beta) {
				result = 1;
			}
			return result;
		}

		private static double atan2Quick(final double y, final double x) {
			final double THREE_QRTR_PI = Math.PI * 0.75;
			final double QRTR_PI = Math.PI * 0.25;

			double r, angle;
			final double abs_y = Math.abs(y) + 1e-10f; // kludge to prevent 0/0 condition

			if (x < 0.0f) {
				r = (x + abs_y) / (abs_y - x); // (3)
				angle = THREE_QRTR_PI; // (4)
			} else {
				r = (x - abs_y) / (x + abs_y); // (1)
				angle = QRTR_PI; // (2)
			}
			angle += (0.1963f * r * r - 0.9817f) * r; // (2 | 4)
			if (y < 0.0f) {
				return (-angle); // negate if in quad III or IV
			} else {
				return (angle);
			}
		}

	}

}
