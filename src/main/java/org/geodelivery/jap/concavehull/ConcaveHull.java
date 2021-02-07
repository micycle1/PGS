package org.geodelivery.jap.concavehull;

import java.util.ArrayList;
import java.util.Collections;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.planargraph.DirectedEdge;
import org.locationtech.jts.planargraph.DirectedEdgeStar;
import org.locationtech.jts.planargraph.Edge;
import org.locationtech.jts.planargraph.Node;
import org.locationtech.jts.planargraph.PlanarGraph;

/**
 * Concave hull algorithm by Pimin Konstantin Kefaloukos and Elias Lsfgren.
 * <p>
 * The algorithm works by extracting a point cloud from the input geometry, and
 * computing the Delaunay Triangulation of this point cloud. The triangulation
 * forms the basis of the output polygon, and the algorithm iteratively removes
 * triangles at the perimeter of the triangulation until a stop criterium is
 * met.
 * </p>
 * 
 * <p>
 * The algorithm takes a geometry and a tunable parameter &Alpha; &isin; {0,1}.
 * The parameter is used together with one of the three heuristics described
 * below.
 * </p>
 * 
 * <p>
 * The perimeter of the triangulation initially corresponds to the Convex Hull.
 * The algorithm circles the perimeter and targets edges that can be deleted,
 * and replaced by the two inner triangle edges. The algorithm applies three
 * rules to determine whether a perimeter edge should be deleted. The first two
 * rules ensure that the resulting polygon is valid. The third rule is a
 * threshold for the minimum length of edges that are deleted.
 * </p>
 * 
 * <p>
 * When computing the minimum edge length threshold, the algorithm uses one of
 * three heuristics:
 * <ul>
 * <li>Average of edge lengths on initial perimeter, times alpha</li>
 * <li>Edge lengths of initial perimeter, sorted. Alpha used as index into
 * sorted array. Alpha=0.5 corresponds to median value.</li>
 * <li>Compute the minimum spanning tree, and use longest MST edge time
 * alpha</li>
 * </ul>
 * </p>
 * 
 * @author Pimin Konstantin Kefaloukos
 */
public class ConcaveHull {

	private final static double DEFAULT_ALPHA = 0.5;
	private double _alpha;
	private ThresholdHeuristic _thresholdHeuristic;

	public ConcaveHull(double alpha) {
		this(ThresholdHeuristic.MED, alpha);
	}

	/**
	 * 
	 * @param thresholdHeuristic
	 * @param alpha 0...1
	 */
	public ConcaveHull(ThresholdHeuristic thresholdHeuristic, double alpha) {
		super();
		this._thresholdHeuristic = thresholdHeuristic;
		this._alpha = alpha;
	}

	public ConcaveHull() {
		this(ThresholdHeuristic.MED, DEFAULT_ALPHA);
	}

	/**
	 * @param compression
	 * @return
	 */
	public Geometry transform(Geometry geom) {

		// marked node means "exposed"
		// marked edge means "deleted"
		// System.out.println("Delaunay");
		DelaunayGraph delaunay = new DelaunayGraph();
		PlanarGraph graph = delaunay.transform(geom);

		// Find threshold
		// System.out.println("MST");
		// double threshold = getThreshold(graph);

		// Find perimeter
		// System.out.println("Perimeter");
		Perimeter perimeter = findPerimeter(graph);

		// Find starting point
		DirectedEdge successorEdge = perimeter.getStartEdge();
		int numExposed = perimeter.getNumExposed();
		// Note: threshold length for removing an edge... this is perhaps a bad
		// heuristic for threshold.
		// Have seen MST work well in other algorithm.
		Node start = successorEdge.getFromNode();
		Node from = start;
		Node to = successorEdge.getToNode();

		// System.out.println("Threshold fast");
		double threshold;
		switch (this._thresholdHeuristic) {
			case AVG:
				threshold = getThresholdAvg(successorEdge);
				break;
			case MED:
				threshold = getThresholdMed(successorEdge);
				break;
			// case MST:
			// threshold = getThresholdMst( graph );
			// break;
			default:
				threshold = getThresholdMed(successorEdge);
				break;
		}

		// do a number of laps
		// System.out.println("ConcaveHull");
		while (true) {

			boolean hasDeleted = false;

			// do a lap around the perimeter
			// delete edges that are already exposed, but not edges that become exposed
			do {
				from = successorEdge.getFromNode();
				to = successorEdge.getToNode();
				double length = from.getCoordinate().distance(to.getCoordinate());

				// compute triangle edges
				DirectedEdge triangleEdge1 = nextEdgeCW(successorEdge);
				Node triangleNode = triangleEdge1.getToNode();
				DirectedEdge inv = triangleEdge1.getEdge().getDirEdge(triangleNode);
				DirectedEdge triangleEdge2 = nextEdgeCW(inv);

				// can successor edge be deleted?
				// rule 1: from and to have degree > 2
				// rule 2: opposite node in CW triangle is not exposed
				// rule 3: edge is longer than or equal to "threshold"
				boolean rule1 = getDegree(successorEdge.getFromNode()) > 2 && getDegree(successorEdge.getToNode()) > 2;
				boolean rule2 = !triangleNode.isMarked();
				boolean rule3 = length >= threshold;

				if (rule1 && rule2 & rule3) {
					// delete edge and replace with two triangle edges
					from.setData(triangleEdge1);
					triangleNode.setData(triangleEdge2);
					markDeleted(successorEdge);
					triangleNode.setMarked(true);
					numExposed++;
					hasDeleted = true;
				}
				successorEdge = (DirectedEdge) to.getData();
				from = to;

			} while (from != start);

			if (!hasDeleted) {
				break;
			}
		}

		Geometry result = toPolygon(start, numExposed);
		return result;
	}

	/*
	 * private double getThresholdMst(PlanarGraph triangulationGraph) { // This
	 * method is slow, but pretty good. Uses MST PlanarGraph mst = new
	 * MinimumSpanningTree().computeGraph(triangulationGraph); double threshold =
	 * Double.MIN_VALUE; for(Object obj : mst.getEdges()) { LineMergeEdge edge =
	 * (LineMergeEdge) obj; threshold = Math.max(threshold,
	 * edge.getLine().getLength()); } return _alpha * threshold; }
	 */
	private double getThresholdAvg(DirectedEdge successorEdge) {
		// Faster threshold.
		// - Only considers edges on perimeter
		// - Uses average length of perimeter edges
		Node start = successorEdge.getFromNode();
		DirectedEdge dirEdge = successorEdge;
		double avg = 0;
		double count = 0;

		do {
			avg += dirEdge.getFromNode().getCoordinate().distance(dirEdge.getToNode().getCoordinate());
			count++;
			dirEdge = (DirectedEdge) dirEdge.getToNode().getData();
		} while (dirEdge.getFromNode() != start);

		return _alpha * avg / count;
	}

	/**
	 * Faster threshold. Only considers edges on perimeter; Uses norm value of
	 * perimeter edges.
	 * 
	 * @param successorEdge
	 * @return
	 */
	private double getThresholdMed(DirectedEdge successorEdge) {

		Node start = successorEdge.getFromNode();
		DirectedEdge dirEdge = successorEdge;
		ArrayList<Double> edgeLengths = new ArrayList<Double>();

		do {
			double dist = dirEdge.getFromNode().getCoordinate().distance(dirEdge.getToNode().getCoordinate());
			edgeLengths.add(dist);
			dirEdge = (DirectedEdge) dirEdge.getToNode().getData();
		} while (dirEdge.getFromNode() != start);
		Collections.sort(edgeLengths);
		int index = (int) Math.floor((edgeLengths.size() - 1) * _alpha);
		return edgeLengths.get(index);
	}

	private Perimeter findPerimeter(PlanarGraph graph) {
		int numExposed = 0;
		Node start, from, to;
		DirectedEdge successorEdge, inversePredecessorEdge, longestEdge;
		// Initialize
		start = findStart(graph); // use a "go west" algorithm
		longestEdge = null; // will be updated in loop 4 sure!
		double longest = Double.MIN_VALUE;
		from = start; // first node of perimeter

		// Special fake node/ege, to trick DirectedEdgeStar.getNextCWEdge
		Node fakeNode = new Node(new Coordinate(start.getCoordinate().x - 10, start.getCoordinate().y));
		DirectedEdge fakeEdge = new DirectedEdge(start, fakeNode, fakeNode.getCoordinate(), true);
		start.addOutEdge(fakeEdge); // edge pointing "west", used to trick the nextEdgeCW function
		inversePredecessorEdge = fakeEdge;
		successorEdge = nextEdgeCW(inversePredecessorEdge); // find next edge and perimeter node
		// no longer need fakeEdge, delete it again, phew, the hoops you have to jump
		// through sometimes..
		start.remove(fakeEdge);

		// do full clock wise roundtrip of perimeter
		// mark nodes and find longest edge
		do {

			// find new "to" node
			to = successorEdge.getToNode();

			// update stuff
			double len = from.getCoordinate().distance(to.getCoordinate());
			if (len > longest) {
				longest = len;
				longestEdge = successorEdge;
			}
			to.setMarked(true);
			numExposed++; // update number of exposed nodes
			from.setData(successorEdge); // set successor

			// move "from" to new position
			inversePredecessorEdge = successorEdge.getEdge().getDirEdge(to);
			from = to;
			successorEdge = nextEdgeCW(inversePredecessorEdge); // find next edge and perimeter node
		} while (from != start);

		return new Perimeter(numExposed, longestEdge);
	}

	/**
	 * Ignore edges marked as deleted when measuring degree of node
	 * 
	 * @param node
	 * @return
	 */
	private int getDegree(Node node) {
		int degree = 0;
		for (Object obj : node.getOutEdges().getEdges()) {
			DirectedEdge e = (DirectedEdge) obj;
			degree += e.isMarked() ? 0 : 1;
		}
		return degree;
	}

	private void markDeleted(DirectedEdge edge) {
		edge.setMarked(true);
		Edge parent = edge.getEdge();
		if (parent != null) {
			parent.setMarked(true);
			DirectedEdge other = parent.getDirEdge(edge.getToNode());
			if (other != null) {
				other.setMarked(true);
			}
		}
	}

	/**
	 * "Robot-arm" perimeter crawling step. Ignores edges in edge star marked as
	 * deleted
	 * 
	 * @param node
	 * @param inverseAngle
	 * @return
	 */
	private DirectedEdge nextEdgeCW(DirectedEdge outEdge) {
		// pick the out-edge of node, that is closest in clock-wise direction from
		// inverseAngle
		Node node = outEdge.getFromNode();
		DirectedEdgeStar star = node.getOutEdges();
		DirectedEdge successor = star.getNextCWEdge(outEdge);
		while (successor.isMarked()) {
			outEdge = successor;
			successor = star.getNextCWEdge(outEdge);
		}
		return successor;
	}

	/**
	 * Algorithm for finding an exposed node in a triangulation graph
	 * 
	 * @param graph
	 * @return
	 */
	private Node findStart(PlanarGraph graph) {
		// going west means going towards smaller values of x
		Node best = (Node) graph.nodeIterator().next();
		double mostWest = best.getCoordinate().x;
		DirectedEdgeStar star = best.getOutEdges();
		// move west for as long as possible
		while (true) {
			boolean improved = false;
			// examine outedges for node more west
			for (Object o : star.getEdges()) {
				DirectedEdge dirEdge = (DirectedEdge) o;
				Node to = dirEdge.getToNode();
				double toX = to.getCoordinate().x;
				if (toX < mostWest) {
					mostWest = toX;
					best = to;
					improved = true;
				}
			}
			if (!improved) {
				break;
			}
		}
		return best;
	}

	private Geometry toPolygon(Node node, int numExposed) {
		// the "data" of each perimeter node is a successor edge
		// number of exposed nodes recorded in _numExposed
		Coordinate[] shell = new Coordinate[numExposed + 1];
		// crawl perimeter
		for (int i = 0; i < numExposed; i++) {
			shell[i] = node.getCoordinate();
			DirectedEdge successorEdge = (DirectedEdge) node.getData();
			node = successorEdge.getToNode();
		}
		shell[numExposed] = shell[0];
		GeometryFactory fact = new GeometryFactory();
		// return just the concave hull
		return fact.createPolygon(fact.createLinearRing(shell), null);
	}

	public enum ThresholdHeuristic {
		// MST,

		/**
		 * Slower, but better
		 */
		AVG,
		/**
		 * Faster
		 */
		MED
	}

	private class Perimeter {
		int _numExposed;
		DirectedEdge _startEdge;

		public Perimeter(int numExposed, DirectedEdge startEdge) {
			super();
			this._numExposed = numExposed;
			this._startEdge = startEdge;
		}

		public int getNumExposed() {
			return _numExposed;
		}

		public DirectedEdge getStartEdge() {
			return _startEdge;
		}

	}
}