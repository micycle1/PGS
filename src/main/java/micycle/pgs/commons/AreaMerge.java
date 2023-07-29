package micycle.pgs.commons;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.util.NeighborCache;
import org.jgrapht.alg.util.Pair;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import micycle.pgs.PGS_Conversion;
import micycle.pgs.PGS_ShapeBoolean;
import micycle.pgs.PGS_ShapePredicates;
import processing.core.PShape;

/**
 * The AreaMerge class provides a method to recursive merge smaller faces of a
 * mesh into their adjacent faces. The procedure continues until there are no
 * resulting faces with an area smaller than the specified threshold.
 * 
 * @author Michael Carleton
 */
public class AreaMerge {

	private AreaMerge() {
	}

	/**
	 * Recursively merges smaller faces of a mesh into their adjacent faces. The
	 * procedure continues until there are no resulting faces with an area smaller
	 * than the specified threshold.
	 * 
	 * @param mesh
	 * @param areaThreshold
	 * @return
	 */
	public static PShape areaMerge(PShape mesh, double areaThreshold) {
		SimpleGraph<PShape, DefaultEdge> graph = PGS_Conversion.toDualGraph(mesh);
		NeighborCache<PShape, DefaultEdge> neighborCache = new NeighborCache<>(graph);

		Map<PShape, FaceGroup> initialFaceMap = new HashMap<>(graph.vertexSet().size());
		SimpleGraph<FaceGroup, DefaultEdge> groupsGraph = new SimpleGraph<>(DefaultEdge.class);
		TreeSet<FaceGroup> smallGroups = new TreeSet<>(); // groups having area < areaThreshold
		for (PShape face : graph.vertexSet()) {
			double area = PGS_ShapePredicates.area(face);
			FaceGroup f = new FaceGroup(face, area, neighborCache);
			initialFaceMap.put(face, f);
			groupsGraph.addVertex(f);

			if (area < areaThreshold) {
				smallGroups.add(f);
			}
		}

		graph.edgeSet().forEach(e -> {
			PShape a = graph.getEdgeSource(e);
			PShape b = graph.getEdgeTarget(e);
			/*
			 * Now add edges to the neighboring groups graph. Initially the groups have the
			 * same topology as the faces, since each group comprises one face.
			 */
			groupsGraph.addEdge(initialFaceMap.get(a), initialFaceMap.get(b));
		});

		while (!smallGroups.isEmpty()) {
			final FaceGroup toMerge = smallGroups.pollFirst();
			Collection<FaceGroup> neighbors = Graphs.neighborSetOf(groupsGraph, toMerge);

			// find smallest neighbor of the toMerge face
//			FaceGroup smallestNeighbor = neighbors.stream().min((a, b) -> Double.compare(a.area, b.area)).orElse(null);
			FaceGroup smallestNeighbor = Graphs.neighborListOf(groupsGraph, toMerge).get(0);
			if (smallestNeighbor == null) {
				break; // exit merging
			}

			smallestNeighbor.mergeWith(toMerge); // merge face groups
			mergeVertices(groupsGraph, smallestNeighbor, toMerge); // update topology

			if (smallestNeighbor.area > areaThreshold) {
				smallGroups.remove(smallestNeighbor);
			} else {
				// re-insert and order the new group
//				smallGroups.remove(smallestNeighbor);
//				smallGroups.add(smallestNeighbor);
			}
		}

		final boolean hasHoles = PGS_ShapePredicates.holes(mesh) > 0;
		return PGS_Conversion.flatten(groupsGraph.vertexSet().stream().map(g -> {
			if (hasHoles) {
				return PGS_ShapeBoolean.unionMesh(g.faces.keySet());
			} else {
				return PGS_ShapeBoolean.unionMeshWithoutHoles(g.faces.keySet());
			}
		}).collect(Collectors.toList()));
	}

	/**
	 * Removes a vertex <code>remove</code> from the graph, but points all of its
	 * incoming edges towards vertex <code>keep</code>, essentially merging the two
	 * vertices.
	 */
	private static <V, E> void mergeVertices(Graph<V, E> graph, V keep, V remove) {
		if (graph == null || keep == null || remove == null) {
			throw new IllegalArgumentException("Arguments cannot be null");
		}
		if (!graph.containsVertex(keep) || !graph.containsVertex(remove)) {
			throw new IllegalArgumentException("Both vertices should be in the graph");
		}

		// Move all edges to the 'keep' vertex and remove the 'remove' vertex
		for (E edge : graph.edgesOf(remove)) {
			Pair<V, V> endpoints = Pair.of(graph.getEdgeSource(edge), graph.getEdgeTarget(edge));

			// Determine the 'other' vertex
			V other = endpoints.getFirst().equals(remove) ? endpoints.getSecond() : endpoints.getFirst();

			// Only create a new edge if it doesn't already exist and isn't a self-loop
			if (!other.equals(keep) && !graph.containsEdge(keep, other)) {
				graph.addEdge(keep, other);
			}
		}

		// Finally, remove the 'remove' vertex
		graph.removeVertex(remove);
	}

	/**
	 * Models a merged group of faces.
	 */
	static class FaceGroup implements Comparable<FaceGroup> {

		private double area;
		/** A map of faces comprising this group and their respective areas. */
		private Map<PShape, Double> faces;
		private Set<PShape> neighborFaces; // union of neighbors of each face - faces themselves
		private final NeighborCache<PShape, DefaultEdge> neighborCache;

		FaceGroup(PShape initial, double area, NeighborCache<PShape, DefaultEdge> meshNeighborCache) {
			this.faces = new HashMap<>();
			this.neighborFaces = new HashSet<>();
			this.neighborCache = meshNeighborCache;
			addFace(initial, area);
		}

		/**
		 * @return true if the face is unique and new to the facegroup
		 */
		public boolean addFace(PShape face, double area) {
			if (face == null) {
//				System.err.println("Tried to add null face.");
				return false;
			}
//			if (!faces.isEmpty() && !neighborFaces.contains(face)) {
//				System.err.println("Tried to add a face that is not a neighboring face of the group.");
//				return false;
//			}
			if (!faces.containsKey(face)) {
				faces.put(face, area);
				this.area += area;
//				recomputeNeighbors(face);
				return true;
			}
			return false;
		}

		private void recomputeNeighbors(PShape face) {
			neighborFaces.addAll(neighborCache.neighborsOf(face));
			/*
			 * Now remove any of the group's own faces from the neighbours.
			 */
			neighborFaces.removeAll(faces.keySet());
		}

		private void recomputeNeighbors() {
			neighborFaces.clear();
			for (PShape face : faces.keySet()) {
				neighborFaces.addAll(neighborCache.neighborsOf(face));
			}
			neighborFaces.removeAll(faces.keySet());
		}

		/**
		 * Merges another facegroup into this one.
		 */
		public boolean mergeWith(FaceGroup other) {
//			if (!neighbors(other)) {
//				return false;
//			}
			boolean changed = false;
			for (Entry<PShape, Double> faceEntry : other.faces.entrySet()) {
				changed = changed | addFace(faceEntry.getKey(), faceEntry.getValue());
			}

			return changed;
		}

		/**
		 * @return true if this and the other facegroup have no faces in common
		 */
		public boolean disjoint(FaceGroup other) {
			Set<PShape> myFaces = new HashSet<>(faces.keySet());
			myFaces.retainAll(other.faces.keySet());
			return myFaces.isEmpty();
		}

		public boolean neighbors(PShape face) {
			return neighborFaces.contains(face);
		}

		/**
		 * @param facegroup
		 * @return true if the groups are disjoint and neighborFaces each other
		 */
		public boolean neighbors(FaceGroup facegroup) {
			if (!disjoint(facegroup)) {
				return false;
			}
			for (PShape face : facegroup.faces.keySet()) {
				if (neighborFaces.contains(face)) {
					return true;
				}
			}
			return false;
		}

		public boolean hasFace(PShape face) {
			return faces.keySet().contains(face);
		}

		@Override
		public int compareTo(FaceGroup o) {
			return Double.compare(this.area, o.area);
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null || getClass() != obj.getClass()) {
				return false;
			}
			FaceGroup other = (FaceGroup) obj;
			return Objects.equals(faces, other.faces);
		}
	}

}
