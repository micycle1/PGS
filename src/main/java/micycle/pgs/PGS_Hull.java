package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.geodelivery.jap.concavehull.SnapHull;
import org.locationtech.jts.algorithm.hull.ConcaveHullOfPolygons;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.union.CascadedPolygonUnion;

import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.concavehull.ConcaveHull;
import uk.osgb.algorithm.concavehull.TriCheckerChi;

/**
 * Convex and concave hulls (amongst other variations) of polygons and point
 * sets.
 * <p>
 * or ... Find minimal coverings (hulls) of polygon arrangements and point/shape
 * sets.
 * 
 * @author Michael Carleton
 * @since 1.2.1
 */
public class PGS_Hull {

	private PGS_Hull() {
	}

	/**
	 * Computes the convex hull of multiple shapes.
	 * 
	 * @param shapes
	 * @return
	 * @see #convexHull(PShape...)
	 */
	public static PShape convexHull(List<PShape> shapes) {
		// TODO REMOVE
		Collection<Polygon> polygons = new ArrayList<>();
		shapes.forEach(s -> polygons.add((Polygon) fromPShape(s)));
		return toPShape(CascadedPolygonUnion.union(polygons).convexHull());
	}

	// https://github.com/locationtech/jts/pull/870
	/**
	 * Computes the convex hull of multiple shapes.
	 * 
	 * @param shapes varArgs
	 * @return
	 * @see #convexHull(List)
	 */
	public static PShape convexHull(PShape... shapes) { // TODO REMOVE
		return convexHull(Arrays.asList(shapes));
	}

	/**
	 * Convex hull of point set.
	 * 
	 * @param points
	 * @return
	 * @since 1.2.1
	 */
	public static PShape convexHull(Collection<PVector> points) {
		return null; // TODO
	}

	/**
	 * Computes concave hulls using a set of polygons as a constraint. The computed
	 * hull respects the polygon boundaries (i.e. the hull is guaranteed to contain
	 * the polygons). (In contrast the existing ConcaveHull class operates only on
	 * the vertices of the input, and the result hull may not contain all the area
	 * of polygons provided as input.)
	 * 
	 * @param shapes    a GROUP PShape, having multiple child PShapes, each of which
	 *                  is a polygon
	 * @param concavity a factor value between 0 and 1, specifying how concave the
	 *                  output is.
	 * @return concave hull of the input shapes
	 * @since 1.2.1
	 */
	public static PShape concaveHull(PShape shapes, double concavity) { // TODO
		ConcaveHullOfPolygons hull = new ConcaveHullOfPolygons(PGS_Conversion.fromPShape(shapes));
		hull.setHolesAllowed(true);
		hull.setMaximumEdgeLengthRatio(1 - concavity);
		return toPShape(hull.getHull());
	}

	/**
	 * Computes the concave hull of a point set using a breadth-first method.
	 * 
	 * @param points
	 * @param threshold euclidean distance threshold
	 * @return
	 * @since 1.1.0
	 * @see #concaveHullDFS(List, double)
	 * @see #concaveHullBFS2(List, double)
	 */
	public static PShape concaveHullBFS(List<PVector> points, double concavity) {
		concavity*=concavity; // square to make output change more linearly as concavity goes 0...1 
		// Geometry g = fromPShape(PGS_Conversion.toPointsPShape(points));
		Geometry g = prepareConcaveGeometry(points);
		org.locationtech.jts.algorithm.hull.ConcaveHull concaveHull = new org.locationtech.jts.algorithm.hull.ConcaveHull(g);
		concaveHull.setMaximumEdgeLengthRatio(concavity);
		return toPShape(concaveHull.getHull());
	}

	/**
	 * Computes the concave hull of a point set using a depth-first method. In
	 * contrast to the BFS method, the depth-first approach produces shapes that are
	 * more contiguous/less branching and spiral-like.
	 * 
	 * @param points
	 * @param threshold euclidean distance threshold
	 * @return
	 * @since 1.1.0
	 * @see #concaveHullBFS(List, double)
	 * @see #concaveHullBFS2(List, double)
	 */
	public static PShape concaveHullDFS(List<PVector> points, double threshold) {
		threshold*=threshold*threshold;
		List<PVector> closestList = PGS_Optimisation.farthestPointPair(points);
		threshold*=closestList.get(0).dist(closestList.get(1));
		ConcaveHull hull = new ConcaveHull(prepareConcaveGeometry(points));
		return toPShape(hull.getConcaveHullDFS(new TriCheckerChi(threshold)));
	}

	/**
	 * Computes the concave hull of a point set using a different algorithm. This
	 * approach has a more "organic" structure compared to other concaveBFS method.
	 * 
	 * @param points
	 * @param threshold 0...1 (Normalized length parameter).
	 *                  <ul>
	 *                  <li>Setting <code>threshold=1</code> means that no edges
	 *                  will be removed from the Delaunay triangulation, so the
	 *                  resulting polygon will be the convex hull.
	 *                  <li>Setting <code>threshold=0</code> means that all edges
	 *                  that can be removed subject to the regularity constraint
	 *                  will be removed (however polygons that are eroded beyond the
	 *                  point where they provide a desirable characterization of the
	 *                  shape).
	 *                  <li>Although the optimal parameter value varies for
	 *                  different shapes and point distributions, values between
	 *                  <code>0.05–0.2</code> typically produce optimal or
	 *                  near-optimal shape characterization across a wide range of
	 *                  point distributions.
	 * @return
	 * @see #concaveHullBFS(List, double)
	 */
	public static PShape concaveHullBFS2(List<PVector> points, double threshold) {
//		threshold*=threshold*threshold;
		/*-
		 * (from https://doi.org/10.1016/j.patcog.2008.03.023)
		 * It is more convenient to normalize the threshold parameter with respect to a
		 * particular set of points P by using the maximum and minimum edge lengths of
		 * the Delaunay triangulation of P. Increasing l beyond the maximum edge length
		 * of the Delaunay triangulation cannot reduce the number of edges that will be
		 * removed (which will be zero anyway). Decreasing l beyond the minimum edge
		 * length of the Delaunay triangulation cannot increase the number of edges that
		 * will be removed.
		 */
		org.geodelivery.jap.concavehull.ConcaveHull hull = new org.geodelivery.jap.concavehull.ConcaveHull(threshold);

		return toPShape(hull.transform(prepareConcaveGeometry(points)));
	}

	/**
	 * Prepares a multipoint geometry from a list of PVectors.
	 */
	private static Geometry prepareConcaveGeometry(List<PVector> points) {
		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
		} else { // already closed
			coords = new Coordinate[points.size()];
		}

		for (int i = 0; i < coords.length; i++) {
			if (i >= points.size()) {
				coords[i] = new Coordinate(points.get(0).x, points.get(0).y); // close geometry
			} else {
				coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
			}
		}

		return PGS.GEOM_FACTORY.createMultiPointFromCoords(coords);
	}

	/**
	 * Computes the "snap hull" for a shape, which is a convex hull that snaps to
	 * the shape. Adjust segment factor to change between
	 * 
	 * @param shape
	 * @param segmentFactor default = 4
	 * @return
	 */
	public static PShape snapHull(PShape shape, double segmentFactor) {
		// TODO JTS Concave Hull of Polygons behaves similarly, so remove this?
		return toPShape(SnapHull.snapHull(fromPShape(shape), segmentFactor));
	}

}
