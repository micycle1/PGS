package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Softer tesselations.
 *
 * A Java implementation of the algorithm described in the paper "A Generative
 * Approach to Smooth Tessellations" (PNAS Nexus, 2024).
 *
 * <p>
 * This class generates and renders smooth, curved tessellations using Bezier
 * curves to soften the edges of a base mesh. The algorithm supports various
 * tangent modes to control the direction and curvature of the edges, enabling a
 * wide range of artistic and geometric effects.
 * </p>
 *
 * <p>
 * The algorithm is highly customizable, supporting multiple tangent modes and
 * input meshes. It can be used for generative art, architectural design, and
 * scientific visualization.
 * </p>
 *
 * @see <a href=
 *      "https://academic.oup.com/pnasnexus/article/3/9/pgae311/7754698">Original
 *      Paper</a>
 * @author Michael Carleton
 * @author CLAUDIO ESPERANÇA
 */
public class SoftCells {

	// https://openprocessing.org/sketch/2419830

	/**
	 * Tangents are computed for each edge based on the selected tangent mode. These
	 * tangents control the curvature of the Bezier curves used to render the edges.
	 *
	 * Each tangent is scaled to half the length of the shortest edge incident to
	 * the source vertex. Unless otherwise noted, the direction is flipped per edge
	 * so the tangent generally points in the same hemisphere as the edge.
	 */
	public enum TangentMode {
		/** Fixed +X direction for all edges of the vertex; no per-edge flip. */
		HORIZONTAL,
		/** Fixed +Y direction for all edges of the vertex; no per-edge flip. */
		VERTICAL,
		/**
		 * 45° diagonal with slope +1 (vector (1,1)); flipped per edge to align via dot
		 * sign.
		 */
		DIAGONAL,
		/**
		 * 45° diagonal with slope −1 (vector (1,−1)); flipped per edge to align via dot
		 * sign.
		 */
		DIAGONAL2,
		/**
		 * Randomly picks one of the two 45° diagonals per vertex; flipped per edge to
		 * align.
		 */
		RANDOM_DIAGONAL,
		/**
		 * Picks one of three directions 60° apart per vertex; flipped per edge to
		 * align.
		 */
		RANDOM_60DEG,
		/**
		 * Uses the sum/average of incident edge directions; flipped per edge to align.
		 */
		ADAPTIVE,
		/** Uniformly random direction per vertex; flipped per edge to align. */
		RANDOM,
		/**
		 * Horizontal base; sign determined by row parity and the edge’s vertical sign.
		 * Edges that are nearly horizontal (abs(u.y) < 1e-3) align to their own
		 * direction.
		 */
		EVEN_ODD,
		/**
		 * Alternates between the two 45° diagonals by (col+row) parity; flipped per
		 * edge to align.
		 */
		ALT_DIAGONAL,
		/**
		 * Cycles through three 60° directions by (col+row)%3; flipped per edge to
		 * align.
		 */
		ALT_60DEG
	}

	private List<PVector> points = new ArrayList<>();
	private List<int[]> faceMap = new ArrayList<>();
	private Map<Integer, PVector> edgeTangentMap = new HashMap<>();
	private List<List<Integer>> vertexEdgeMap = new ArrayList<>();

	private float avgSize = 1f;
	private float minX = 0, minY = 0, maxX = 0, maxY = 0;

	// deterministic RNG
	private Random rng = new Random(0L);

	public SoftCells() {
		this(0L);
	}

	public SoftCells(long seed) {
		setSeed(seed);
	}

	/**
	 * Convenience method: load mesh, compute tangents, build shape in one call.
	 */
	public PShape generate(PShape mesh, TangentMode mode, float ratio) {
		loadMeshFromPShape(mesh);
		computeTangents(mode);
		return buildSoftFacesShape(ratio);
	}

	public void setSeed(long seed) {
		this.rng = new Random(seed);
	}

	/**
	 * Loads mesh data from a PShape and populates face/vertex data structures.
	 * Assumes the input mesh is conforming, meaning vertices are shared exactly
	 * between faces.
	 *
	 * @param mesh PShape containing the mesh (each child is a face polygon)
	 */
	public void loadMeshFromPShape(final PShape mesh) {
		// 1) reset
		resetDataStructures();

		// 2) collect unique vertices by value
		final Set<PVector> unique = new LinkedHashSet<>();
		for (int f = 0; f < mesh.getChildCount(); f++) {
			final PShape face = mesh.getChild(f);
			for (int v = 0; v < face.getVertexCount(); v++) {
				unique.add(face.getVertex(v));
			}
		}

		// 3) build and sort points (x, then y, then z)
		points = new ArrayList<>(unique);
		points.sort(Comparator.<PVector, Float>comparing(p -> p.x).thenComparing(p -> p.y));

		// 4) init adjacency and index map (value-based)
		vertexEdgeMap = new ArrayList<>(points.size());
		for (int i = 0; i < points.size(); i++) {
			vertexEdgeMap.add(new ArrayList<>());
		}
		final Map<PVector, Integer> index = new HashMap<>(points.size() * 2);
		for (int i = 0; i < points.size(); i++) {
			index.put(points.get(i), i);
		}

		// 5) faces + adjacency using sorted indices
		for (int f = 0; f < mesh.getChildCount(); f++) {
			final PShape face = mesh.getChild(f);
			final int n = face.getVertexCount();
			final int[] faceVerts = new int[n];

			for (int v = 0; v < n; v++) {
				final PVector p = face.getVertex(v);
				final Integer idx = index.get(p);
				if (idx == null) {
					throw new IllegalStateException("Vertex not in index map: " + p);
				}
				faceVerts[v] = idx;
			}
			faceMap.add(faceVerts);

			for (int i = 0; i < n; i++) {
				final int a = faceVerts[i];
				final int b = faceVerts[(i + 1) % n];
				addUnique(vertexEdgeMap.get(a), b);
				addUnique(vertexEdgeMap.get(b), a);
			}
		}

		// 6) finalise
		sortVertexNeighborsByAngle();
		computeBoundsAndAvgSize();
	}

	private void resetDataStructures() {
		points = new ArrayList<>();
		faceMap = new ArrayList<>();
		vertexEdgeMap = new ArrayList<>();
		edgeTangentMap = new HashMap<>();
		avgSize = 1f;
		minX = minY = maxX = maxY = 0;
	}

	private void addUnique(final List<Integer> list, final int value) {
		if (!list.contains(value)) {
			list.add(value);
		}
	}

	private void sortVertexNeighborsByAngle() {
		for (int i = 0; i < points.size(); i++) {
			final PVector center = points.get(i);
			final List<Integer> neighbors = vertexEdgeMap.get(i);

			neighbors.sort((a, b) -> {
				final PVector va = PVector.sub(points.get(a), center);
				final PVector vb = PVector.sub(points.get(b), center);
				return Float.compare(va.heading(), vb.heading());
			});
		}
	}

	private void computeBoundsAndAvgSize() {
		if (points.isEmpty()) {
			avgSize = 1f;
			return;
		}

		minX = maxX = points.get(0).x;
		minY = maxY = points.get(0).y;

		for (PVector p : points) {
			if (p.x < minX) {
				minX = p.x;
			}
			if (p.y < minY) {
				minY = p.y;
			}
			if (p.x > maxX) {
				maxX = p.x;
			}
			if (p.y > maxY) {
				maxY = p.y;
			}
		}

		// average unique edge length
		double sum = 0.0;
		int cnt = 0;
		for (int i = 0; i < points.size(); i++) {
			final PVector pi = points.get(i);
			for (int j : vertexEdgeMap.get(i)) {
				if (j > i) { // unique edge (i,j) with i<j
					final PVector pj = points.get(j);
					sum += PVector.dist(pi, pj);
					cnt++;
				}
			}
		}
		avgSize = (cnt > 0) ? (float) (sum / cnt) : Math.max(1f, (maxX - minX + maxY - minY) * 0.01f);
		if (avgSize <= 0) {
			avgSize = 1f;
		}
	}

	public void computeTangents(final TangentMode mode) {
		edgeTangentMap = new HashMap<>();
		for (int i = 0; i < points.size(); i++) {
			final TangentEstimatorFunc tangentFunc = getTangentEstimator(mode, i);
			final List<Integer> neighbors = vertexEdgeMap.get(i);
			if (neighbors != null) {
				for (int k = 0; k < neighbors.size(); k++) {
					final int j = neighbors.get(k);
					final PVector tangent = tangentFunc.estimateTangent(k, i, j);
					edgeTangentMap.put(getKey(i, j), tangent);
				}
			}
		}
	}

	private int getKey(final int src, final int dst) {
		return points.size() * src + dst;
	}

	/**
	 * Creates a tangent estimator function for a given vertex based on the
	 * specified tangent mode.
	 *
	 * <p>
	 * This method is a core part of the algorithm for generating smooth, curved
	 * tessellations. It determines the direction and magnitude of tangents for each
	 * edge connected to a vertex, which are later used to compute Bezier curves for
	 * rendering the tessellation. The tangent estimator function returned by this
	 * method is specific to a single vertex and is applied to all its neighboring
	 * edges.
	 * </p>
	 *
	 * <p>
	 * The tangent mode defines how tangents are calculated:
	 * <ul>
	 * <li><b>Fixed Directions:</b> Modes like {@link TangentMode#HORIZONTAL} and
	 * {@link TangentMode#DIAGONAL} use predefined directions (e.g., horizontal,
	 * vertical, or diagonal) to compute tangents.</li>
	 * <li><b>Adaptive Directions:</b> Modes like {@link TangentMode#ADAPTIVE}
	 * calculate tangents based on the average direction of neighboring edges,
	 * ensuring smooth transitions.</li>
	 * <li><b>Randomized Directions:</b> Modes like {@link TangentMode#RANDOM} and
	 * {@link TangentMode#RANDOM_60DEG} generate random directions for tangents,
	 * either uniformly or constrained to specific angles (e.g., 60°).</li>
	 * </ul>
	 * </p>
	 *
	 * <p>
	 * The tangent estimator function returned by this method takes the index of a
	 * neighboring edge and computes a tangent vector for that edge. The tangent
	 * vector is scaled to half the length of the shortest neighboring edge (to
	 * ensure smooth curvature) and is oriented based on the edge's direction and
	 * the chosen tangent mode.
	 * </p>
	 *
	 * @param mode The tangent mode to use for computing tangents. This determines
	 *             how tangent directions are calculated.
	 * @param i    The index of the vertex for which the tangent estimator is being
	 *             created.
	 * @return A {@link TangentEstimatorFunc} that computes tangent vectors for
	 *         edges connected to the vertex.
	 */
	private TangentEstimatorFunc getTangentEstimator(final TangentMode mode, final int i) {
		final List<Integer> neighbors = vertexEdgeMap.get(i);
		final List<PVector> neighborVectors = new ArrayList<>();
		final PVector srcPoint = points.get(i);

		// build the list of edge-vectors emanating from vertex i
		for (final int nbr : neighbors) {
			neighborVectors.add(PVector.sub(points.get(nbr), srcPoint));
		}

		// compute average length, minimum length and sum-direction
		float avgLen = 0;
		float minLen = Float.MAX_VALUE;
		final PVector avgDir = new PVector();
		if (!neighborVectors.isEmpty()) {
			for (final PVector v : neighborVectors) {
				final float m = v.mag();
				avgLen += m;
				minLen = Math.min(minLen, m);
				avgDir.add(v);
			}
			avgLen /= neighborVectors.size();
		} else {
			minLen = 0;
		}

		final float halfMin = minLen * 0.5f;
		// make any base directions we’ll need (already scaled to halfMin)
		final PVector horiz = new PVector(1, 0).mult(halfMin);
		final PVector vert = new PVector(0, 1).mult(halfMin);
		final PVector diag1 = new PVector(1, 1).normalize().mult(halfMin);
		final PVector diag2 = new PVector(1, -1).normalize().mult(halfMin);

		// PRE‐CHOOSING random direction once per vertex:
		// note: this is never mutated later
		final float TWO_PI = (float) (Math.PI * 2.0);
		final PVector randDir = PVector.fromAngle(rng.nextFloat() * TWO_PI).mult(halfMin);
		final PVector randDiagonal = (rng.nextFloat() < 0.5f ? diag1 : diag2);

		// adaptive means sum-direction
		final PVector adaptiveDir = avgDir.copy().setMag(halfMin);

		switch (mode) {
			case HORIZONTAL :
				return (k, src, dst) -> {
					// WORK ON A COPY:
					return horiz.copy();
				};

			case VERTICAL :
				return (k, src, dst) -> {
					return vert.copy();
				};

			case DIAGONAL :
				return (k, src, dst) -> {
					final PVector dir = diag1.copy();
					if (neighborVectors.get(k).dot(dir) < 0) {
						dir.mult(-1);
					}
					return dir;
				};

			case DIAGONAL2 :
				return (k, src, dst) -> {
					final PVector dir = diag2.copy();
					if (neighborVectors.get(k).dot(dir) < 0) {
						dir.mult(-1);
					}
					return dir;
				};

			case RANDOM_DIAGONAL :
				// randDiagonal was chosen once above, do NOT re‐roll per edge
				return (k, src, dst) -> {
					final PVector dir = randDiagonal.copy();
					if (neighborVectors.get(k).dot(dir) < 0) {
						dir.mult(-1);
					}
					return dir;
				};

			case RANDOM_60DEG :
				// pick one of three 60° directions once per vertex
				final PVector[] tris = { new PVector(0, 1), new PVector(0, 1).rotate((float) Math.toRadians(60)),
						new PVector(0, 1).rotate((float) Math.toRadians(-60)) };
				final PVector triDir = tris[rng.nextInt(tris.length)].normalize().mult(halfMin);
				return (k, src, dst) -> {
					final PVector dir = triDir.copy();
					if (neighborVectors.get(k).dot(dir) < 0) {
						dir.mult(-1);
					}
					return dir;
				};

			case ADAPTIVE :
				return (k, src, dst) -> {
					final PVector dir = adaptiveDir.copy();
					if (neighborVectors.get(k).dot(dir) < 0) {
						dir.mult(-1);
					}
					return dir;
				};

			case RANDOM :
				return (k, src, dst) -> {
					final PVector dir = randDir.copy();
					if (neighborVectors.get(k).dot(dir) < 0) {
						dir.mult(-1);
					}
					return dir;
				};

			case EVEN_ODD :
				return (k, src, dst) -> {
					final int col = (int) ((srcPoint.x - minX) / avgSize);
					final int row = (int) ((srcPoint.y - minY) / avgSize);
					final PVector base = horiz.copy();
					final PVector u = neighborVectors.get(k);
					// if perfectly horizontal edge, just align
					if (Math.abs(u.y) < 1e-3) {
						if (u.dot(base) < 0) {
							base.mult(-1);
						}
						return base;
					}
					if (row % 2 == 0) {
						if (u.y < 0) {
							base.mult(-1);
						}
					} else {
						if (u.y > 0) {
							base.mult(-1);
						}
					}
					return base;
				};

			case ALT_DIAGONAL :
				return (k, src, dst) -> {
					final int col = (int) ((srcPoint.x - minX) / avgSize);
					final int row = (int) ((srcPoint.y - minY) / avgSize);
					final PVector pick = ((col + row) % 2 == 0 ? diag1 : diag2).copy();
					if (neighborVectors.get(k).dot(pick) < 0) {
						pick.mult(-1);
					}
					return pick;
				};

			case ALT_60DEG :
				return (k, src, dst) -> {
					final int col = (int) ((srcPoint.x - minX) / avgSize);
					final int row = (int) ((srcPoint.y - minY) / avgSize);
					final PVector[] altTris = { new PVector(0, 1), new PVector(0, 1).rotate((float) Math.toRadians(60)),
							new PVector(0, 1).rotate((float) Math.toRadians(-60)) };
					final PVector pick = altTris[(Math.abs(col + row)) % 3].normalize().mult(halfMin);
					if (neighborVectors.get(k).dot(pick) < 0) {
						pick.mult(-1);
					}
					return pick;
				};

			default :
				return (k, src, dst) -> new PVector(0, 0);
		}
	}

	interface TangentEstimatorFunc {
		PVector estimateTangent(int k, int srcIndex, int dstIndex);
	}

	public PShape buildSoftFacesShape(final float ratio) {
		final PShape grp = new PShape(PConstants.GROUP);

		for (final int[] vtx : faceMap) {
			final PShape poly = new PShape(PShape.PATH);
			poly.setStroke(0);
			poly.setStroke(true);
			poly.setStrokeWeight(2);
			poly.setFill(255);
			poly.setFill(true);

			poly.beginShape();
			int i = vtx[vtx.length - 1];
			PVector p1 = points.get(i);
			poly.vertex(p1.x, p1.y);

			for (final int j : vtx) {
				final PVector p2 = points.get(j);
				final PVector t1 = edgeTangentMap.get(points.size() * i + j);
				final PVector t2 = edgeTangentMap.get(points.size() * j + i);

				if (t1 != null && t2 != null) {
					poly.bezierVertex(p1.x + t1.x * ratio, p1.y + t1.y * ratio, p2.x + t2.x * ratio, p2.y + t2.y * ratio, p2.x, p2.y);
				} else {
					poly.vertex(p2.x, p2.y);
				}
				i = j;
				p1 = p2;
			}

			poly.endShape(PConstants.CLOSE);
			grp.addChild(poly);
		}

		return grp;
	}
}