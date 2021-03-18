package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;

import java.util.ArrayList;

import org.geodelivery.jap.concavehull.SnapHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.simplify.VWSimplifier;

import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.concavehull.ConcaveHull;
import uk.osgb.algorithm.concavehull.TriCheckerChi;
import uk.osgb.algorithm.minkowski_sum.Minkowski_Sum;

/**
 * Methods that affect the geometry or topology of shapes.
 * 
 * @author Michael Carleton
 *
 */
public class PTSMorphology {

	static {
		Minkowski_Sum.setGeometryFactory(PTS.GEOM_FACTORY);
	}

	/**
	 * 
	 * @param shape
	 * @param buffer extent/width of the buffer (may be positive or negative)
	 * @return
	 */
	public static PShape buffer(PShape shape, float buffer) {
		// TODO read
		// https://locationtech.github.io/jts/javadoc/org/locationtech/jts/operation/buffer/BufferOp.html
		return toPShape(fromPShape(shape).buffer(buffer, 4));
	}

	/**
	 * small edges are removed, while the general structure of the shape is
	 * preserved. This process is known as "opening" in computer vision.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape erosionDilation(PShape shape, float buffer) {
//		buffer = Math.abs(buffer);
		return toPShape(fromPShape(shape).buffer(-buffer).buffer(buffer));
	}

	/**
	 * Simplifies a PShape using the Douglas-Peucker algorithm.
	 * 
	 * @param shape
	 * @param distanceTolerance
	 * @return
	 * @see #topologySimplify(PShape, float)
	 */
	public static PShape simplify(PShape shape, float distanceTolerance) {
		return toPShape(DouglasPeuckerSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Simplifies a shape using the Visvalingam-Whyatt area-based algorithm.
	 * 
	 * @param shape
	 * @param distanceTolerance The simplification tolerance is specified as a
	 *                          distance.This is converted to an area tolerance by
	 *                          squaring it.
	 * @return
	 */
	public static PShape simplifyVW(PShape shape, float distanceTolerance) {
		return toPShape(VWSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Preserves topological structure (holes, etc.)
	 * 
	 * @param shape
	 * @param distanceTolerance
	 * @return
	 * @see #simplify(PShape, float)
	 */
	public static PShape topologySimplify(PShape shape, float distanceTolerance) {
		return toPShape(TopologyPreservingSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Computes the convex hull of multiple PShapes
	 * 
	 * @param shapes
	 * @return
	 */
	public static PShape convexHull(PShape... shapes) {
		Geometry g = fromPShape(shapes[0]);
		for (int i = 1; i < shapes.length; i++) {
			g = g.union(fromPShape(shapes[i]));
		}
		return toPShape(g.convexHull());
	}

	/**
	 * Uses Chi heurstic
	 * 
	 * @param points
	 * @param threshold euclidean distance threshold
	 * @return
	 * @see #concaveHull2(ArrayList, float)
	 */
	public static PShape concaveHull(ArrayList<PVector> points, float threshold) {
	
		// calls Ordnance Survey implementation
	
		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
			points.add(points.get(0)); // close geometry
		} else { // already closed
			coords = new Coordinate[points.size()];
		}
	
		for (int i = 0; i < coords.length; i++) {
			coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
		}
	
		Geometry g = PTS.GEOM_FACTORY.createPolygon(coords);
		ConcaveHull hull = new ConcaveHull(g);
		return toPShape(hull.getConcaveHullBFS(new TriCheckerChi(threshold), false, false).get(0));
	}

	/**
	 * Computes the concave hull of a shape using a different algorithm. Has a more
	 * "organic" structure compared to other concave method.
	 * 
	 * @param points
	 * @param threshold 0...1 (Normalized length parameter). Setting λP = 1 means
	 *                  that no edges will be removed from the Delaunay
	 *                  triangulation, so the resulting polygon will be the convex
	 *                  hull. Setting λP = 0 means that all edges that can be
	 *                  removed subject to the regularity constraint will be removed
	 *                  (however polygons that are eroded beyond the point where
	 *                  they provide a desirable characterization of the shape).
	 *                  Although the optimal parameter value varies for different
	 *                  shapes and point distributions, values of between 0.05–0.2
	 *                  typically produce optimal or near-optimal shape
	 *                  characterization across a wide range of point distributions.
	 * @return
	 * @see #concaveHull(ArrayList, float)
	 */
	public static PShape concaveHull2(ArrayList<PVector> points, float threshold) {
	
		/**
		 * (from https://doi.org/10.1016/j.patcog.2008.03.023) It is more convenient to
		 * normalize the threshold parameter with respect to a particular set of points
		 * P by using the maximum and minimum edge lengths of the Delaunay triangulation
		 * of P. Increasing l beyond the maximum edge length of the Delaunay
		 * triangulation cannot reduce the number of edges that will be removed (which
		 * will be zero anyway). Decreasing l beyond the minimum edge length of the
		 * Delaunay triangulation cannot increase the number of edges that will be
		 * removed.
		 */
	
		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
			points.add(points.get(0)); // close geometry
		} else { // already closed
			coords = new Coordinate[points.size()];
		}
	
		for (int i = 0; i < coords.length; i++) {
			coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
		}
	
		Geometry g = PTS.GEOM_FACTORY.createPolygon(coords);
	
		// TODO test AVG threshold heuristic
		org.geodelivery.jap.concavehull.ConcaveHull hull = new org.geodelivery.jap.concavehull.ConcaveHull(threshold);
	
		return toPShape(hull.transform(g));
	}

	/**
	 * Adjust segment factor to change between
	 * 
	 * @param shape
	 * @param segmentFactor
	 * @return
	 */
	public static PShape snapHull(PShape shape, float segmentFactor) {
		return toPShape(SnapHull.snapHull(fromPShape(shape), segmentFactor));
	}

	/**
		 * Minkowski addition a.k.a dilation
		 * 
		 * @return
		 */
		public static PShape minkSum(PShape source, PShape addition) {
	//		Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
			Geometry sum = Minkowski_Sum.compMinkSum(fromPShape(source), fromPShape(addition), false, false);
			return toPShape(sum);
		}

	/**
		 * TODO check, a.k.a erosion
		 * 
		 * @param source
		 * @param addition
		 * @return
		 */
		public static PShape minkDifference(PShape source, PShape addition) {
	//			Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
			Geometry sum = Minkowski_Sum.compMinkDiff(fromPShape(source), fromPShape(addition), false, false);
			return toPShape(sum);
		}

}
