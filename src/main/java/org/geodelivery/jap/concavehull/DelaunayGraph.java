package org.geodelivery.jap.concavehull;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.operation.linemerge.LineMergeGraph;
import org.locationtech.jts.planargraph.PlanarGraph;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;

class DelaunayGraph {

	public DelaunayGraph() {
		super();
	}

	public PlanarGraph transform(Geometry geom) {

		DelaunayTriangulationBuilder triangulator = new DelaunayTriangulationBuilder();
		triangulator.setSites(geom);
		Geometry edges = triangulator.getEdges(new GeometryFactory());

		LineMergeGraph graph = new LineMergeGraph();
		for (int i = 0; i < edges.getNumGeometries(); i++) {
			LineString ls = (LineString) edges.getGeometryN(i);
			graph.addEdge(ls);
		}
		return graph;
	}

}