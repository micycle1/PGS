package org.locationtech.jts.operation.overlayng;

import java.util.Collection;
import java.util.List;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.noding.Noder;

/**
 * An OverlayNG optimised for situtations where the user wishes to compute
 * multiple different boolean operations on 2 shapes. Faster if a union AND
 * (intersection OR difference) is required. Possbily slower if just
 * intersection and differences are needed.
 * 
 * @author Michael Carleton
 *
 */
public class MultiOperationOverlayNG {

	protected final GeometryFactory geomFact;
	protected final Geometry geom0, geom1;
	protected InputGeometry inputGeom;
	protected PrecisionModel pm;
	protected Noder noder;

	protected OverlayGraph graph;
	protected OverlayLabeller labeller;

	protected boolean isOptimized = true;
	protected boolean isStrictMode = false;
	protected boolean isAreaResultOnly = false;

	/**
	 * 
	 * @param geom0 a geometry formed of edges (Lines or Polygons)
	 * @param geom1 a geometry formed of edges (Lines or Polygons)
	 */
	public MultiOperationOverlayNG(Geometry geom0, Geometry geom1) {
		this.geom0 = geom0;
		this.geom1 = geom1;
		inputGeom = new InputGeometry(geom0, geom1);
		geomFact = geom0.getFactory();
		pm = geom0.getFactory().getPrecisionModel();
	}

	protected void init() {
		graph = buildGraph(nodeEdges());
		labeller = new OverlayLabeller(graph, inputGeom);
		labeller.computeLabelling();
	}

	protected OverlayGraph buildGraph(Collection<Edge> edges) {
		OverlayGraph g = new OverlayGraph();
		for (Edge e : edges) {
			g.addEdge(e.getCoordinates(), e.createLabel());
		}
		return g;
	}

	public Geometry getResult(int opCode) {
		if (graph == null) {
			init();
		}

		labeller.markResultAreaEdges(opCode);
		labeller.unmarkDuplicateEdgesFromResultArea();
		return extractResult(graph.getResultAreaEdges(), opCode);
	}

	private List<Edge> nodeEdges() {
		/*
		 * Node the edges, using whatever noder is being used
		 */
		EdgeNodingBuilder nodingBuilder = new EdgeNodingBuilder(pm, noder);

		/*
		 * Optimize Intersection and Difference by clipping to the result extent, if
		 * enabled.
		 */
//		if (isOptimized) { // NOTE cannot be done in multi operation
//			Envelope clipEnv = OverlayUtil.clippingEnvelope(OverlayNG.DIFFERENCE, inputGeom, pm);
//			if (clipEnv != null) {
//				nodingBuilder.setClipEnvelope(clipEnv);
//			}
//		}

		List<Edge> mergedEdges = nodingBuilder.build(inputGeom.getGeometry(0), inputGeom.getGeometry(1));

		/*
		 * Record if an input geometry has collapsed. This is used to avoid trying to
		 * locate disconnected edges against a geometry which has collapsed completely.
		 */
		inputGeom.setCollapsed(0, !nodingBuilder.hasEdgesFor(0));
		inputGeom.setCollapsed(1, !nodingBuilder.hasEdgesFor(1));

		return mergedEdges;
	}

	private Geometry extractResult(List<OverlayEdge> resultAreaEdges, int opCode) {
		boolean isAllowMixedIntResult = !isStrictMode;

		// --- Build polygons
		PolygonBuilder polyBuilder = new PolygonBuilder(resultAreaEdges, geomFact);
		List<Polygon> resultPolyList = polyBuilder.getPolygons();
		boolean hasResultAreaComponents = !resultPolyList.isEmpty();

		List<LineString> resultLineList = null;
		List<Point> resultPointList = null;

		if (!isAreaResultOnly) {
			// handle lines (don't don't handle points)
			boolean allowResultLines = !hasResultAreaComponents || isAllowMixedIntResult || opCode == OverlayNG.SYMDIFFERENCE
					|| opCode == OverlayNG.UNION;
			if (allowResultLines) {
				LineBuilder lineBuilder = new LineBuilder(inputGeom, graph, hasResultAreaComponents, opCode, geomFact);
				lineBuilder.setStrictMode(isStrictMode);
				resultLineList = lineBuilder.getLines();
			}
		}

		Geometry resultGeom = OverlayUtil.createResultGeometry(resultPolyList, resultLineList, resultPointList, geomFact);

		// reset results for next call
		resultAreaEdges.forEach(e -> {
			e.setEdgeRing(null);
			e.setEdgeRingMax(null);
			e.setNextResult(null);
			e.setNextResultMax(null);
		});
		return resultGeom;
	}

}
