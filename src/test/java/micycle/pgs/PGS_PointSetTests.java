package micycle.pgs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.junit.jupiter.api.Test;

import processing.core.PVector;

class PGS_PointSetTests {

	@Test
	void testPrunePointsWithinDistance() {
		List<PVector> points = new ArrayList<>();
		points.add(new PVector(0, 0));
		points.add(new PVector(1, 0)); // within distance 5 of (0,0)
		points.add(new PVector(10, 0)); // farther than 5 from (0,0) and (1,0)
		points.add(new PVector(11, 0)); // within distance 5 of (10,0)

		List<PVector> pruned = PGS_PointSet.prunePointsWithinDistance(points, 5);
		assertEquals(2, pruned.size());
		assertEquals(new PVector(0, 0), pruned.get(0));
		assertEquals(new PVector(10, 0), pruned.get(1));
		List<PVector> unpruned = PGS_PointSet.prunePointsWithinDistance(List.of(new PVector(0, 0), new PVector(1, 0), new PVector(2, 0)), 0);
		assertEquals(3, unpruned.size());

		List<PVector> heavilyPruned = PGS_PointSet.prunePointsWithinDistance(List.of(new PVector(0, 0), new PVector(1, 0), new PVector(2, 0)), 1000);
		assertEquals(1, heavilyPruned.size());
	}

	@Test
	void testPruneRandomRemoveN() {
		List<PVector> points = new ArrayList<>();
		for (int i = 0; i < 20; i++) {
			points.add(new PVector(i, i));
		}

		List<PVector> removedNone = PGS_PointSet.pruneRandomRemoveN(points, 0, 42);
		assertEquals(points.size(), removedNone.size());

		List<PVector> reduced = PGS_PointSet.pruneRandomRemoveN(points, 5, 42);
		assertEquals(15, reduced.size());

		List<PVector> repeated = PGS_PointSet.pruneRandomRemoveN(points, 5, 42);
		assertEquals(reduced, repeated);

		Set<PVector> original = new HashSet<>(points);
		for (PVector p : reduced) {
			assertTrue(original.contains(p), "Result should only contain original points");
		}
	}

	@Test
	void testPruneRandomToN() {
		List<PVector> points = new ArrayList<>();
		for (int i = 0; i < 20; i++) {
			points.add(new PVector(i, i));
		}

		List<PVector> reduced = PGS_PointSet.pruneRandomToN(points, 8, 42);
		assertEquals(8, reduced.size());

		List<PVector> repeated = PGS_PointSet.pruneRandomToN(points, 8, 42);
		assertEquals(reduced, repeated);
	}

	@Test
	void testHilbertSortPreservesElements() {
		List<PVector> points = new ArrayList<>();
		for (int i = 0; i < 30; i++) {
			points.add(new PVector(i * 10, i * 5));
		}
		List<PVector> sorted = PGS_PointSet.hilbertSort(points);
		assertEquals(points.size(), sorted.size());
		// Same elements, different order
		assertTrue(new HashSet<>(sorted).containsAll(points));
		assertTrue(new HashSet<>(points).containsAll(sorted));
	}

	@Test
	void testHilbertSortSmallList() {
		// Lists < 24 points are returned unchanged
		List<PVector> points = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			points.add(new PVector(i, 0));
		}
		List<PVector> sorted = PGS_PointSet.hilbertSort(points);
		assertEquals(points, sorted); // unchanged
	}

	@Test
	void testHilbertSortSpatialLocality() {
		// After sorting, consecutive points should tend to be spatially close.
		// We verify the total path length of the Hilbert-sorted sequence is shorter
		// than a deliberately bad ordering (e.g. zigzag across the grid).
		List<PVector> grid = new ArrayList<>();
		for (int y = 0; y < 10; y++) {
			for (int x = 0; x < 10; x++) {
				grid.add(new PVector(x * 10, y * 10));
			}
		}
		// Zigzag: even rows left-to-right, odd rows right-to-left — but shuffled
		List<PVector> zigzag = new ArrayList<>(grid);
		Collections.shuffle(zigzag, new Random(999));

		List<PVector> sorted = PGS_PointSet.hilbertSort(grid);

		double sortedPathLen = pathLength(sorted);
		double zigzagPathLen = pathLength(zigzag);

		assertTrue(sortedPathLen < zigzagPathLen, "Hilbert-sorted path length (" + sortedPathLen + ") should be shorter than shuffled (" + zigzagPathLen + ")");
	}

	@Test
	void testClusterGroupCount() {
		List<PVector> points = new ArrayList<>();
		// 3 tight clusters
		for (int i = 0; i < 10; i++) {
			points.add(new PVector(i, 0)); // cluster 1
			points.add(new PVector(100 + i, 0)); // cluster 2
			points.add(new PVector(200 + i, 0)); // cluster 3
		}

		List<List<PVector>> clusters = PGS_PointSet.cluster(points, 3, 42);
		assertEquals(3, clusters.size());

		// All points accounted for
		int total = clusters.stream().mapToInt(List::size).sum();
		assertEquals(30, total);
	}

	@Test
	void testClusterDeterministic() {
		List<PVector> points = new ArrayList<>();
		for (int i = 0; i < 30; i++) {
			points.add(new PVector(i * 10, i * 5));
		}
		List<List<PVector>> a = PGS_PointSet.cluster(points, 3, 42);
		List<List<PVector>> b = PGS_PointSet.cluster(points, 3, 42);
		assertEquals(a.size(), b.size());
		for (int i = 0; i < a.size(); i++) {
			assertEquals(a.get(i).size(), b.get(i).size());
		}
	}

	private static double pathLength(List<PVector> pts) {
		double len = 0;
		for (int i = 1; i < pts.size(); i++) {
			len += pts.get(i).dist(pts.get(i - 1));
		}
		return len;
	}

}
