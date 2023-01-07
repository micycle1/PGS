package micycle.pgs.commons;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import processing.core.PConstants;
import processing.core.PShape;

/**
 * Generates random space-filling curves.
 * <p>
 * A space-filling curve is a continuous curve that traverses every point of an
 * plane exactly once.
 * <p>
 * This class implements the technique described in "Context-based space filling
 * curves" by Dafner et al, and is based on <a href=
 * "https://observablehq.com/@esperanc/random-space-filling-curves">this</a>
 * javascript implementation by Claudio Esperança.
 * 
 * @author Michael Carleton
 *
 */
public class RandomSpaceFillingCurve {

	private final Random random;

	private final int ncols, nrows;
	private boolean[] grid;
	private List<Edge> edges;
	private ArrayDeque<int[]> cellVisitQueue;
	private Connect[] connect;
	private List<Integer> path;

	private final int nrows2, ncols2;
	private List<Edge> edges2;

	/**
	 * Instantiates a random space-filling curve. The curve is created upon
	 * initialisation.
	 * 
	 * @param ncols number of cells in the x direction (columns)
	 * @param nrows number of cells in the y direction (rows)
	 */
	public RandomSpaceFillingCurve(int ncols, int nrows) {
		this(ncols, nrows, System.currentTimeMillis());
	}

	/**
	 * Instantiates a random space-filling curve, having a given random seed. The
	 * curve is created upon initialisation.
	 * 
	 * @param ncols number of cells in the x direction (columns)
	 * @param nrows number of cells in the y direction (rows)
	 * @param seed
	 */
	public RandomSpaceFillingCurve(int ncols, int nrows, long seed) {
		random = new XoRoShiRo128PlusRandom(seed);
		this.ncols = ncols;
		this.nrows = nrows;
		nrows2 = 2 * nrows;
		ncols2 = 2 * ncols;

		gridSpanningTree();
		hamiltonianFromSpanningTree();
	}

	public PShape getCurve(float cellWidth, float cellHeight) {
		cellWidth /= 2;
		cellHeight /= 2;
		PShape curve = new PShape(PShape.PATH);
		curve.setStrokeWeight(3);
		curve.setFill(true);
		curve.setStroke(true);
		curve.setFill(-1231234);
		curve.setStrokeJoin(PConstants.MITER);

		curve.beginShape();
		for (Integer i : path) {
			curve.vertex(((i % ncols2) + 0.5f) * cellWidth + 0.5f, ((i / ncols2) + 0.5f) * cellHeight + 0.5f);
		}
		curve.endShape(PConstants.CLOSE);

		return curve;
	}

	public PShape getSkeleton(float cellWidth, float cellHeight) {
		PShape skeleton = new PShape(PShape.GEOMETRY);
		skeleton.setStrokeWeight(2);
		skeleton.setStroke(true);
		skeleton.setStrokeJoin(PConstants.MITER);
		skeleton.setStrokeCap(PConstants.SQUARE);

		skeleton.beginShape(PConstants.LINES);
		for (Edge e : edges) {
			skeleton.vertex(((e.k % ncols) + 0.5f) * cellWidth + 0.5f, ((e.k / ncols) + 0.5f) * cellHeight + 0.5f);
			skeleton.vertex(((e.n % ncols) + 0.5f) * cellWidth + 0.5f, ((e.n / ncols) + 0.5f) * cellHeight + 0.5f);
		}
		skeleton.endShape();

		return skeleton;
	}

	private void gridSpanningTree() {
		grid = new boolean[ncols * nrows];
		edges = new ArrayList<>();

		cellVisitQueue = new ArrayDeque<>();
		cellVisitQueue.add(new int[] { random.nextInt(ncols * nrows), -1 });
		while (!cellVisitQueue.isEmpty()) {
			visit(cellVisitQueue.pollLast());
		}
		edges.remove(0);

		connect = Stream.iterate(0, x -> x + 1).limit(ncols * nrows).map(i -> new Connect()).toArray(Connect[]::new);

		for (Edge e : edges) {
			int i = e.k;
			int j = e.n;
			if (i > j) { // swap
				int tmp = i;
				i = j;
				j = tmp;
			}
			if (i / ncols == j / ncols) {
				connect[i].right = connect[j].left = true;
			} else {
				connect[i].down = connect[j].up = true;
			}
		}
	}

	private void visit(int[] in) {
		int k = in[0];
		int parentN = in[1];
		if (grid[k]) {
			return;
		}

		edges.add(new Edge(k, parentN));
		grid[k] = true;
		final int i = k % ncols;
		final int j = k / ncols;

		List<Integer> neighbors = new ArrayList<>(4);
		if (i > 0) {
			neighbors.add(k - 1);
		}
		if (j > 0) {
			neighbors.add(k - ncols);
		}
		if (i + 1 < ncols) {
			neighbors.add(k + 1);
		}
		if (j + 1 < nrows) {
			neighbors.add(k + ncols);
		}

		Collections.shuffle(neighbors, random);
		for (Integer n : neighbors) {
			if (!grid[n]) {
				cellVisitQueue.add(new int[] { n, k });
			}
		}
	}

	private void hamiltonianFromSpanningTree() {
		edges2 = new ArrayList<>();
		int i = 0;
		for (Connect cell : connect) {
			if (cell.right) {
				edges2.add(new Edge(index2(i, 1, 0), index2(i, 2, 0)));
				edges2.add(new Edge(index2(i, 1, 1), index2(i, 2, 1)));
			} else {
				edges2.add(new Edge(index2(i, 1, 0), index2(i, 1, 1)));
			}
			if (!cell.left) {
				edges2.add(new Edge(index2(i, 0, 0), index2(i, 0, 1)));
			}
			if (cell.down) {
				edges2.add(new Edge(index2(i, 0, 1), index2(i, 0, 2)));
				edges2.add(new Edge(index2(i, 1, 1), index2(i, 1, 2)));
			} else {
				edges2.add(new Edge(index2(i, 0, 1), index2(i, 1, 1)));
			}
			if (!cell.up) {
				edges2.add(new Edge(index2(i, 0, 0), index2(i, 1, 0)));
			}
			i++;
		}

		List<ArrayList<Integer>> link = Stream.iterate(0, x -> x + 1).limit(ncols2 * nrows2).map(q -> new ArrayList<Integer>())
				.collect(Collectors.toList());
		for (Edge e : edges2) {
			link.get(e.k).add(e.n);
			link.get(e.n).add(e.k);
		}

		int j = 0;
		path = new ArrayList<>();
		boolean[] visited = new boolean[edges2.size()];
		for (int k = edges2.size(); k > 0; k--) {
			path.add(j);
			visited[j] = true;
			j = visited[link.get(j).get(0)] ? link.get(j).get(1) : link.get(j).get(0);
		}
	}

	private int index2(int i, int dcol, int drow) {
		return ((i / ncols) * 2 + drow) * ncols2 + (i % ncols) * 2 + dcol;
	}

	private class Edge {
		int k, n;

		Edge(int k, int n) {
			this.k = k;
			this.n = n;
		}
	}

	private class Connect {
		boolean left, right, up, down;
	}

}
