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
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Models a GROUP PShape comprising faces of a 2D mesh as a topological entity,
 * ready for mesh optimisation via Laplacian smoothing. In Laplacian smoothing,
 * vertices are replaced with the average of the positions of their adjacent
 * vertices.
 * 
 * @author Michael Carleton
 * @since 1.4.0
 */
public class PMesh {

	/*-
	 * Implement HC-algorithm from "Improved Laplacian Smoothing of Noisy Surface Meshes"?
	 * Cotangent weights? http://rodolphe-vaillant.fr/entry/69/c-code-for-cotangent-weights-over-a-triangular-mesh
	 */

	private final PShape mesh;
	private Map<PVector, PMeshVertex> meshVertices;

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
	}

	/**
	 * Performs one pass of simple laplacian smoothing on the mesh.
	 * <p>
	 * During laplacian smoothing, vertices are moved to the geometric center of
	 * their incident vertices, generally resulting in more isotropic mesh.
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

		for (PMeshVertex mv : meshVertices.values()) {
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
	 * During laplacian smoothing, vertices are moved to the geometric center of
	 * their incident vertices, generally resulting in more isotropic mesh. In this
	 * weighted variant, the neighbours of a given vertex are weighted by their
	 * euclidean distance to the vertex (so the vertex is "pulled" more towards its
	 * farther-away neighbours).
	 * 
	 * @param excludeBoundaryVertices a boolean value indicating whether or not to
	 *                                exclude the boundary vertices from being
	 *                                smoothed. Generally this should be set to
	 *                                true, otherwise the mesh will shrink as it is
	 *                                smoothed.
	 * @return the displacement of the most displaced vertex
	 */
	public float smoothWeighted(final boolean excludeBoundaryVertices) {
		float totalDisplacement = 0;
		int displacedVertices = 0;
		float maxDisplacement = 0;

		for (PMeshVertex mv : meshVertices.values()) {
			if (excludeBoundaryVertices && mv.onBoundary) {
				continue;
			}

			PVector mean = new PVector(0, 0);
			float weight = 0;
			for (PMeshVertex n : mv.neighbors) {
				PVector neighbor = n.smoothedVertex;
				float w = neighbor.dist(mv.smoothedVertex);
				mean.add(PVector.mult(neighbor, w));
				weight += w;
			}
			mean.div(weight);
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
	 * Performs a single pass of Taubin smoothing. Each Taubin smoothing pass
	 * peforms two laplacian smoothing passes on the mesh (imwards then outwards) in
	 * an effort to prevent the mesh from shrinking during the smoothing operation.
	 *
	 * @param lamda first Taubin parameter. Should be positive
	 * @param mu    second Taubin parameter. Should be negative, with an absolute
	 *              value greater than lambda
	 */
	public void smoothTaubin(double lamda, double mu) {
		if (lamda > Math.abs(mu) + 1e-3) {
			mu = -(lamda + 1e-3);
		}
		smoothScaled(lamda);
		if (mu != 0) {
			smoothScaled(mu);
		}
	}

	/**
	 * Performs one pass of scaled laplacian smoothing on the mesh.
	 * 
	 * @param excludeBoundaryVertices
	 * @return the average displacement distance of the smoothed vertices
	 */
	private float smoothScaled(double scale) {
		List<PVector> L = new ArrayList<>();
		for (int i = 0; i < meshVertices.size(); i++) {
			L.add(new PVector());
		}

		final float s = (float) scale;
		float totalDisplacement = 0;
		int displacedVertices = 0;

		int i = 0;
		List<PMeshVertex> mvs = new ArrayList<>(meshVertices.values());
		for (PMeshVertex mv : mvs) {
			PVector lap = L.get(i++);
			for (PMeshVertex n : mv.neighbors) {
				PVector neighbor = n.smoothedVertex;
				lap.add(neighbor);
			}
			lap.mult(s / mv.neighbors.size());
			lap.add(PVector.mult(mv.smoothedVertex, -s));
			totalDisplacement += mv.smoothedVertex.dist(lap);
			displacedVertices++;
		}
		for (int j = 0; j < meshVertices.size(); j++) {
			mvs.get(j).smoothedVertex.add(L.get(j));
		}
		return totalDisplacement / displacedVertices;
	}

	/**
	 * 
	 * Smooths the mesh until the maximum displacement between the original and
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

	// see
	// http://rodolphe-vaillant.fr/entry/69/c-code-for-cotangent-weights-over-a-triangular-mesh
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
