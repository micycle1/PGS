package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedPolygon;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.index.SpatialIndex;
import org.locationtech.jts.index.hprtree.HilbertEncoder;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.operation.overlayng.OverlayNG;

/**
 * Computes the combined area where two or more input shapes overlap. Optimised
 * for efficient calculation for many input shapes/intersections.
 * 
 * @author Michael Carleton
 */
public class FastOverlapRegions {

	private final SpatialIndex spatialIndex;
	private final List<IndexedGeom> indexedGeoms;
	private final GeometryFactory factory;

	@SuppressWarnings("unchecked")
	public FastOverlapRegions(Geometry geometry) {
		this(PolygonExtracter.getPolygons(geometry), geometry.getFactory());
	}

	public FastOverlapRegions(List<Geometry> geometries, GeometryFactory factory) {
		// build STRtree of valid input geometries
		spatialIndex = new STRtree();
		indexedGeoms = new ArrayList<>(geometries.size());
		this.factory = factory;

		int k = 0;
		for (int i = 0; i < geometries.size(); i++) {
			Geometry g = geometries.get(i);
			if (g == null || g.isEmpty() || g.getDimension() < 2 || !(g instanceof Polygonal)) {
				continue;
			}
			IndexedGeom ig = new IndexedGeom(g, k++);
			indexedGeoms.add(ig);
			spatialIndex.insert(g.getEnvelopeInternal(), ig);
		}

//		spatialIndex.build();
	}

	public Geometry get(boolean union) {
		// find overlapping regions
		var patches = findPairwiseOverlaps();

		if (!union) {
			return factory.buildGeometry(patches.stream().map(p -> p.geom).toList());
		}

		/*
		 * Before we would find and union connected component groups (groups of patches
		 * that are known to form an intersecting "blob") to reduce unneccessary union
		 * operations. However, this optimisation is negligible given the hilbert-sorted
		 * parallelStream approach.
		 */
//		var groups = findConnectedComponents(patches);

		var geoms = patches.stream().map(p -> p.geom).collect(Collectors.toList());
		return fastUnion(geoms);
	}

	/**
	 * Computes intersections between all geometries, with parallelism and efficient
	 * short-circuiting.
	 * 
	 * @return
	 */
	private List<IntersectPatch> findPairwiseOverlaps() {
		AtomicInteger patchId = new AtomicInteger(0);
		var allPairwiseIntersections = indexedGeoms.parallelStream().flatMap(igA -> {
			Stream.Builder<IntersectPatch> builder = Stream.builder();
			Geometry a = igA.geom;

			@SuppressWarnings("unchecked")
			List<IndexedGeom> candidates = (List<IndexedGeom>) spatialIndex.query(a.getEnvelopeInternal());

			for (IndexedGeom igB : candidates) {
				// skip self‐pair and duplicates
				if (igB.idx <= igA.idx) {
					continue;
				}

				// fast reject if they don’t actually meet
				Geometry b = igB.geom;
				if (!igA.preparedGeometry.intersects(b)) {
					continue;
				}

				// compute the real intersection
				Geometry inter = OverlayNG.overlay(a, b, OverlayNG.INTERSECTION); // polygonal
				if (!inter.isEmpty()) {
					builder.add(new IntersectPatch(patchId.getAndIncrement(), inter));
				}

			}
			return builder.build();
		}).toList();

		return allPairwiseIntersections;
	}

	/**
	 * Much faster CascadedPolygonUnion (similar idea, but faster sorting and
	 * parallel union).
	 */
	private Geometry fastUnion(List<Geometry> geoms) {
		HilbertEncoder.sort(geoms);
		return geoms.parallelStream().reduce((g1, g2) -> {
			var result = g1.union(g2);
			if (!(result instanceof Polygonal)) {
				var polygons = PolygonExtracter.getPolygons(result);
				if (polygons.size() == 1) {
					result = (Polygon) polygons.get(0);
				} else {
					result = factory.buildGeometry(polygons);
				}
			}
			return result;
		}).orElse(factory.createEmpty(2));
	}

	/**
	 * Finds groups of connected (intersecting) geometries from the input list.
	 * Returns a list of lists, where each inner list contains geometries of one
	 * connected component.
	 */
	private List<List<Geometry>> findConnectedComponents(List<IntersectPatch> patches) {
		STRtree patchIndex = new STRtree();
		for (IntersectPatch ip : patches) {
			patchIndex.insert(ip.env, ip);
		}
		patchIndex.build();

		UnionFind uf2 = new UnionFind(patches.size());

		// For every patch, look for neighbouring patches whose polygons actually
		// touch/overlap
		List<int[]> unionPairs = patches.parallelStream().flatMap(ip -> {
			@SuppressWarnings("unchecked")
			List<IntersectPatch> neigh = patchIndex.query(ip.env);
			return neigh.stream().filter(jp -> jp.idx > ip.idx).filter(jp -> ip.preparedGeometry.intersects(jp.geom)).map(jp -> new int[] { ip.idx, jp.idx });
		}).toList();

		// Do union in serial, to avoid concurrency issues
		for (int[] pair : unionPairs) {
			uf2.union(pair[0], pair[1]);
		}

		Map<Integer, List<Geometry>> islands = new HashMap<>();
		for (IntersectPatch ip : patches) {
			int root = uf2.find(ip.idx);
			islands.computeIfAbsent(root, k -> new ArrayList<>()).add(ip.geom);
		}
		return new ArrayList<>(islands.values());
	}

	private static record IndexedGeom(Geometry geom, int idx, PreparedGeometry preparedGeometry) {
		private IndexedGeom(Geometry geom, int idx) {
			this(geom, idx, new PreparedPolygon((Polygonal) geom));
		}
	}

	private static record IntersectPatch(int idx, // a dense 0..M−1 index
			Geometry geom, Envelope env, PreparedGeometry preparedGeometry) {
		public IntersectPatch(int idx, Geometry g) {
			this(idx, g, g.getEnvelopeInternal(), new PreparedPolygon((Polygonal) g));
		}
	}

	private static class UnionFind {
		private final int[] parent, size; // change rank to size

		public UnionFind(int n) {
			parent = new int[n];
			size = new int[n];
			for (int i = 0; i < n; i++) {
				parent[i] = i;
				size[i] = 1; // initialize to 1
			}
		}

		public int find(int x) {
			// iterative path‐halving
			while (parent[x] != x) {
				parent[x] = parent[parent[x]];
				x = parent[x];
			}
			return x;
		}

		public void union(int a, int b) {
			int ra = find(a), rb = find(b);
			if (ra == rb)
				return;
			if (size[ra] < size[rb]) {
				parent[ra] = rb;
				size[rb] += size[ra];
			} else {
				parent[rb] = ra;
				size[ra] += size[rb];
			}
		}
	}

}
