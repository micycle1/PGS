package org.geodelivery.jap.concavehull;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;

/**
 * Java implementation of ST_ConcaveHull from PostGis 2.0.
 * 
 * @author Pimin Konstantin Kefaloukos
 * 
 */
public class SnapHull {
	
	private SnapHull() {
	}

	// default values
	private static final double START_RANGE_RELATIVE = 0.03;
	private static final int MAX_HULL_POINTS = 2000;
	private static final int HULL_SEGMENT_FACTOR = 4;

	public static Geometry snapHull(Geometry src) {
		return snapHull(src, START_RANGE_RELATIVE, MAX_HULL_POINTS, HULL_SEGMENT_FACTOR);
	}

	public static Geometry snapHull(Geometry src, double hullSegmentFactor) {
		return snapHull(src, START_RANGE_RELATIVE, MAX_HULL_POINTS, hullSegmentFactor);
	}

	/**
	 * 
	 * @param src
	 * @param startRangeRelative 0.03
	 * @param maxHullPoints      any
	 * @param hullSegmentFactor  4
	 * @return
	 */
	public static Geometry snapHull(Geometry src, double startRangeRelative, int maxHullPoints,
			double hullSegmentFactor) {

		GeometryFactory gf = new GeometryFactory();
		Envelope srcEnv = src.getEnvelopeInternal();
		double srcDim = (srcEnv.getHeight() + srcEnv.getWidth()) / 2d;
		double startRange = startRangeRelative * srcDim;
		Coordinate[] coordinates = src.getCoordinates();
		ConvexHull ch = new ConvexHull(coordinates, gf);
		Geometry geometry = ch.getConvexHull();

		if (geometry instanceof Polygon) {
			// get the exterior ring
			LineString vexring = ((Polygon) geometry).getExteriorRing();
			float numVertices = (int) Math.min(vexring.getNumPoints() * hullSegmentFactor, maxHullPoints);
			double seglength = vexring.getLength() / numVertices;
			vexring = segmentize(vexring, seglength);
			Coordinate[] result = new Coordinate[vexring.getNumPoints()]; // upperbound on verts on boundary
			int bindex = 0;

			// build index of points
			STRtree index = new STRtree();
			for (Coordinate c : coordinates) {
				index.insert(new Envelope(c), c);
			}
			index.build();
			Coordinate previous = null;
			for (Coordinate c : vexring.getCoordinates()) {
				// This proceduce creates invalid polygons. Find better solution.
				Coordinate nearest = findNearest(c, startRange, srcDim, index);
				if (nearest != previous) {
					result[bindex++] = nearest;
					previous = nearest;
				}
			}
			Coordinate[] shell = new Coordinate[bindex];
			System.arraycopy(result, 0, shell, 0, bindex);
			Geometry p = gf.createPolygon(gf.createLinearRing(shell), null);

			if (!p.isValid()) {
				DouglasPeuckerSimplifier simp = new DouglasPeuckerSimplifier(p);
				p = simp.getResultGeometry();
			}
			return p;
		} else {
			return geometry; // linestring, point or empty, return convexhull
		}
	}

	private static Coordinate findNearest(Coordinate qc, double range, double maxRange, STRtree index) {
		while (range < maxRange) {
			Envelope searchEnv = new Envelope(qc.x - range, qc.x + range, qc.y - range, qc.y + range);
			@SuppressWarnings("rawtypes")
			List hits = index.query(searchEnv);
			if (hits.isEmpty()) {
				range *= 2;
			} else {
				Coordinate best = null;
				double bestDist = Double.MAX_VALUE;

				for (Object obj : hits) {
					Coordinate hit = (Coordinate) obj;
					double dist = qc.distance(hit);
					if (dist < bestDist) {
						bestDist = dist;
						best = hit;
					}
				}
				return best;
			}
		}
		return null;
	}

	private static LineString segmentize(LineString vexring, double seglength) {

		Coordinate[] vcoords = vexring.getCoordinates();
		GeometryFactory gf = new GeometryFactory();
		ArrayList<Coordinate> ext = new ArrayList<>();
		for (int i = 0; i < vcoords.length - 1; i++) {
			ext.add(vcoords[i]);
			// start debug
			double dist = vcoords[i].distance(vcoords[i + 1]); // OK
			double numSegments = Math.ceil(dist / seglength); // OK
			double actLength = dist / numSegments; // OK
			double factor = actLength / dist;
			double vectorX = factor * (vcoords[i + 1].x - vcoords[i].x);
			double vectorY = factor * (vcoords[i + 1].y - vcoords[i].y);
			Coordinate ins = vcoords[i];
			while (numSegments > 1) {
				ins = new Coordinate(ins.x + vectorX, ins.y + vectorY);
				ext.add(ins);
				numSegments--;
			}
			// end debug
			ext.add(vcoords[i + 1]);
		}
		return gf.createLinearRing(ext.toArray(new Coordinate[ext.size()]));
	}

}