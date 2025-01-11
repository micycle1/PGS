package micycle.pgs.commons;

import static micycle.pgs.PGS_Conversion.getChildren;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleWeightedGraph;

import micycle.pgs.PGS_Conversion;
import net.jafama.FastMath;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Models a GROUP PShape comprising faces of a 2D mesh as a topological entity,
 * ready for mesh optimisation via Laplacian smoothing, a form of diffusion
 * smoothing. In Laplacian smoothing, vertices are replaced with the average of
 * the positions of their adjacent vertices.
 * 
 * @author Michael Carleton
 * @since 1.4.0
 */
public class PMesh {

	/*-
	 * Implement HC-algorithm from "Improved Laplacian Smoothing of Noisy Surface Meshes"
	 */

	private final PShape mesh;
	private Map<PVector, PMeshVertex> meshVertices;
	private PMeshVertex[] meshVerticesList; // faster iteration than meshVertices.values()

	public PMesh(PShape mesh) {
		this.mesh = mesh;
		initMeshVertices();
	}

	private void initMeshVertices() {
		meshVertices = new HashMap<>(3 * mesh.getChildCount());
		Map<PEdge, Integer> edgeCounts = new HashMap<>();
		SimpleGraph<PVector, PEdge> g = new SimpleWeightedGraph<>(PEdge.class);
		for (PShape face : getChildren(mesh)) {
			for (int i = 0; i < face.getVertexCount(); i++) {
				final PVector a = face.getVertex(i);
				final PVector b = face.getVertex((i + 1) % face.getVertexCount());
				if (a.equals(b)) {
					continue;
				}
				final PEdge e = new PEdge(a, b);

				edgeCounts.merge(e, 1, Integer::sum);

				if (g.addVertex(a)) {
					meshVertices.put(a, new PMeshVertex(a));
				}
				if (g.addVertex(b)) {
					meshVertices.put(b, new PMeshVertex(b));
				}
				g.addEdge(a, b, e);
			}
		}

		// mark perimeter vertices
		List<PEdge> perimeterEdges = edgeCounts.entrySet().stream().filter(entry -> entry.getValue().intValue() == 1)
				.map(entry -> entry.getKey()).collect(Collectors.toList());
		perimeterEdges.forEach(e -> {
			meshVertices.get(e.a).onBoundary = true;
			meshVertices.get(e.b).onBoundary = true;
		});

		// populate vertex neighbors
		meshVertices.forEach((v, mv) -> {
			List<PMeshVertex> neighbors = new ArrayList<>(5);
			g.outgoingEdgesOf(v).forEach(e -> {
				if (e.a.equals(v)) {
					neighbors.add(meshVertices.get(e.b));
				} else {
					neighbors.add(meshVertices.get(e.a));
				}
			});
			mv.neighbors = neighbors;
		});

		meshVerticesList = new ArrayList<>(meshVertices.values()).toArray(new PMeshVertex[0]);
	}

	/**
	 * Performs one pass of simple laplacian smoothing on the mesh.
	 * <p>
	 * During simple laplacian smoothing, vertices are moved to the barycenter
	 * center of their neighbors (equally-weighted), generally resulting in more
	 * isotropic mesh.
	 * 
	 * @param excludeBoundaryVertices a boolean value indicating whether or not to
	 *                                exclude the boundary vertices from being
	 *                                smoothed. Generally this should be set to
	 *                                true, otherwise the mesh will shrink as it is
	 *                                smoothed.
	 * @return the average displacement distance of the smoothed vertices
	 */
	public float smooth(final boolean excludeBoundaryVertices) {
		float totalDisplacement = 0;
		int displacedVertices = 0;

		for (PMeshVertex mv : meshVerticesList) {
			if (excludeBoundaryVertices && mv.onBoundary) {
				continue;
			}

			PVector mean = new PVector(0, 0);
			/*
			 * The original point *could* be included (having a variable weighting) to bind
			 * the location of the smoothed vertex to its original position.
			 */
			for (PMeshVertex n : mv.neighbors) {
				PVector neighbor = n.smoothedVertex;
				mean.add(neighbor);
			}
			mean.div(mv.neighbors.size());
			totalDisplacement += mv.smoothedVertex.dist(mean);
			displacedVertices++;

			/*
			 * Modify the position immediately (this variant of laplacian smoothing is
			 * called "sequential", in constrast the the "simultaneous" version).
			 */
			mv.smoothedVertex.set(mean);
		}
		return totalDisplacement / displacedVertices;
	}

	/**
	 * Performs one pass of weighted laplacian smoothing on the mesh.
	 * <p>
	 * This slightly more complex approximation of the Laplacian uses weights
	 * proportional to the inverse distance between the vertices, such that a vertex
	 * is "pulled" more towards its farther-away neighbours.
	 * 
	 * @param excludeBoundaryVertices a boolean value indicating whether or not to
	 *                                exclude the boundary vertices from being
	 *                                smoothed. Generally this should be set to
	 *                                true, otherwise the mesh will shrink as it is
	 *                                smoothed.
	 * @return the displacement distance of the most displaced vertex
	 */
	public float smoothWeighted(final boolean excludeBoundaryVertices) {
		float totalDisplacement = 0;
		int displacedVertices = 0;
		float maxDisplacement = 0;

		for (PMeshVertex mv : meshVerticesList) {
			if (excludeBoundaryVertices && mv.onBoundary) {
				continue;
			}

			PVector mean = new PVector(0, 0);
			float totalWeight = 0;
			for (PMeshVertex n : mv.neighbors) {
				PVector neighbor = n.smoothedVertex;
				float w = neighbor.dist(mv.smoothedVertex);
				mean.add(PVector.mult(neighbor, w));
				totalWeight += w;
			}
			mean.div(totalWeight);
			final float displacement = mv.smoothedVertex.dist(mean);
			totalDisplacement += displacement;
			maxDisplacement = Math.max(maxDisplacement, displacement);
			displacedVertices++;

			mv.smoothedVertex.set(mean);
		}
//		return totalDisplacement / displacedVertices;
		return maxDisplacement;
	}

	/**
	 * Performs a single pass of Taubin smoothing. Each Taubin pass peforms two
	 * laplacian smoothing passes on the mesh (shrinking then inflating) in an
	 * effort to preserve the volume of the mesh.
	 * <p>
	 * In geometric terms, Taubin smoothing diffuses the mesh inwards and outwards
	 * to attenuate details while keeping the surface in roughly the same position.
	 * Although this approach is not guaranteed to preserve mesh volume, it does a
	 * good job if the parameters λ and µ are well chosen. On the down side, Taubin
	 * smoothing requires more iterations to achieve a level of smoothing comparable
	 * to other methods.
	 * <p>
	 * A good starting point for λ and µ is a µ negative value that is mildly larger
	 * than λ: i.e. <code>λ=0.2</code>, <code>µ=-0.201</code>.
	 *
	 * @param λ                       Controls the amount of inward diffusion /
	 *                                shrinkage. Should be positive.
	 * @param µ                       Controls the amount of outward diffusion /
	 *                                inflation. Should be negative, with an
	 *                                absolute value greater than lambda.
	 * @param excludeBoundaryVertices A boolean value indicating whether or not to
	 *                                exclude the boundary vertices from being
	 *                                smoothed, preserving their original location
	 *                                and the area of the mesh.
	 * @return the displacement distance of the most displaced vertex during the
	 *         shrinkage.
	 */
	public float smoothTaubin(double λ, double µ, boolean excludeBoundaryVertices) {
//		if (λ > Math.abs(µ) + 1e-5) {
//			µ = -(λ + 1e-5);
//		}
		float displacement = smoothScaled(λ, excludeBoundaryVertices);
		if (µ != 0) {
			smoothScaled(µ, excludeBoundaryVertices);
		}

		return displacement;
	}

	/**
	 * Performs one pass of scaled laplacian smoothing on the mesh.
	 * 
	 * <p>
	 * This smoothing affects vertices <i>simultaneously</i> (i.e. the displacement
	 * is first computed for every vertex against the original coordinates, and
	 * finally applied to every vertex at the end of the iteration).
	 * 
	 * @param scale                   A scalar parameter that controls the diffusion
	 *                                speed.
	 * @param excludeBoundaryVertices A boolean value indicating whether or not to
	 *                                exclude the boundary vertices from being
	 *                                smoothed, preserving their original location
	 *                                and the area of the mesh.
	 * @return the displacement distance of the most displaced vertex
	 */
	private float smoothScaled(final double scale, final boolean excludeBoundaryVertices) {
		// Use an array instead of ArrayList for displacements
		final PVector[] displacements = new PVector[meshVerticesList.length];

		float maxDisplacement = 0;
		final float s = (float) scale;

		for (int i = 0; i < meshVerticesList.length; i++) {
			final PMeshVertex mv = meshVerticesList[i];
			final PVector laplacian = new PVector();
			displacements[i] = laplacian;

			if (excludeBoundaryVertices && mv.onBoundary) {
				continue;
			}

			for (PMeshVertex n : mv.neighbors) {
				PVector neighbor = n.smoothedVertex;
				laplacian.add(neighbor);
			}
			laplacian.mult(s / mv.neighbors.size());
			laplacian.add(PVector.mult(mv.smoothedVertex, -s));

			float displacementSq = laplacian.mag();
			maxDisplacement = Math.max(displacementSq, maxDisplacement);
		}

		for (int j = 0; j < meshVerticesList.length; j++) {
			PMeshVertex mv = meshVerticesList[j];
			mv.smoothedVertex.add(displacements[j]);
		}

		return maxDisplacement;
	}

	private float smoothCotanWeighted(final boolean excludeBoundaryVertices) { // NOTE working?
		// https://rodolphe-vaillant.fr/entry/69/c-code-for-cotangent-weights-over-a-triangular-mesh
	
		final double eps = 1e-6f;
		final double cotan_max = FastMath.cos(eps) / FastMath.sin(eps);
	
		for (PMeshVertex mv : meshVertices.values()) {
			if (excludeBoundaryVertices && mv.onBoundary) {
				continue;
			}
			final PVector i = mv.smoothedVertex;
			final PVector mean = new PVector(0, 0);
			double totalWeight = 0;
	
			for (int k = 0; k < mv.neighbors.size(); k++) {
				PVector v_prev = mv.neighbors.get(k == 0 ? mv.neighbors.size() - 1 : k - 1).smoothedVertex;
				PVector v = mv.neighbors.get(k).smoothedVertex;
				PVector v_next = mv.neighbors.get((k + 1) % mv.neighbors.size()).smoothedVertex;
	
				// Calculate cotangent weights
				PVector v1 = PVector.sub(i, v_prev);
				PVector v2 = PVector.sub(v, v_prev);
				PVector v3 = PVector.sub(i, v_next);
				PVector v4 = PVector.sub(v, v_next);
	
				double cotan_alpha = cotan(v1, v2);
				double cotan_beta = cotan(v3, v4);
	
				double weight = (cotan_alpha + cotan_beta);
	
				if (Double.isNaN(weight)) {
					weight = 0;
				}
				/*
				 * Compute the cotangent value close to 0.0f and PI. As cotan approaches
				 * infinity close to those values, we clamp.
				 */
				weight = clamp(weight, -cotan_max, cotan_max);
				mean.add(PVector.mult(v, (float) weight));
				totalWeight += weight;
	
			}
			mean.div((float) totalWeight);
			mv.smoothedVertex.set(mean);
	
		}
	
		return 0;
	}

	private static float cotan(final PVector a, final PVector b) {
		return PVector.dot(a, b) / PVector.cross(a, b, null).normalize().mag();
	}

	private static double clamp(double value, double min, double max) {
		return Math.max(min, Math.min(value, max));
	}

	/**
	 * Smoothes the mesh until the maximum displacement between the original and
	 * smoothed vertices is less than or equal to the specified threshold value.
	 * 
	 * @param preservePerimeter a boolean value indicating whether or not to
	 *                          preserve the perimeter of the mesh during smoothing
	 * @param maxDisplacement   the maximum allowable average displacement between
	 *                          the original and smoothed vertices.
	 */
	int relax(final boolean preservePerimeter, double minDisplacement) {
		minDisplacement = Double.max(minDisplacement, 1e-9);
		int relaxations = 0;
		double displacement = Double.MAX_VALUE;

		do {
			displacement = smooth(preservePerimeter);
			relaxations++;
		} while (displacement > minDisplacement);

		return relaxations;
	}

	/**
	 * 
	 * Returns a PShape object representing the smoothed mesh, with vertices moved
	 * to their smoothed positions.
	 * 
	 * @return a PShape object representing the smoothed mesh
	 */
	public PShape getMesh() {
		PShape out = new PShape(PConstants.GROUP);
		PGS_Conversion.getChildren(mesh).forEach(f -> {
			PShape shape = PGS_Conversion.copy(f);
			for (int i = 0; i < f.getVertexCount(); i++) {
				shape.setVertex(i, meshVertices.get(shape.getVertex(i)).smoothedVertex);
			}
			out.addChild(shape);
		});
		return out;
	}

	private static class PMeshVertex {

		final PVector originalVertex;
		final PVector smoothedVertex;
		boolean onBoundary;
		List<PMeshVertex> neighbors;

		public PMeshVertex(PVector v) {
			originalVertex = v; // same reference
			smoothedVertex = v.copy();
		}

		@Override
		public int hashCode() {
			int result = 1;
			result = 31 * result + Float.floatToIntBits(originalVertex.x);
			result = 31 * result + Float.floatToIntBits(originalVertex.y);
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof PMeshVertex)) {
				return false;
			}
			final PMeshVertex v = (PMeshVertex) obj;
			return originalVertex.equals(v.originalVertex);
		}
	}
}
