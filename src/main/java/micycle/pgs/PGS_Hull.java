package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.Collection;
import java.util.List;

import org.locationtech.jts.algorithm.hull.ConcaveHullOfPolygons;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
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
	 * Computes the convex hull of a point set (the smallest convex polygon that
	 * contains all the points).
	 * 
	 * @param points a collection of points
	 * @return the minimum-area convex polygon containing the points
	 * @since 1.2.1
	 */
	public static PShape convexHull(Collection<PVector> points) {
		return toPShape(fromPShape(PGS_Conversion.toPointsPShape(points)).convexHull());
	}

	/**
	 * Computes the convex hull of the vertices from the input shape (the smallest
	 * <b>convex</b> polygon that contains all the shape's vertices).
	 * 
	 * @param shape a concave shape
	 * @return the minimum-area convex polygon containing the input's vertices
	 * @since 1.2.1
	 */
	public static PShape convexHull(PShape shape) {
		return toPShape(fromPShape(shape).convexHull());
	}

	/**
	 * Constructs a concave hull of a set of polygons, respecting the polygons as
	 * constraints (i.e. the hull is guaranteed to contain the polygons).
	 * 
	 * @param shapeSet  a GROUP PShape, having multiple child PShapes, each of which
	 *                  is a polygon
	 * @param concavity a factor value between 0 and 1, specifying how concave the
	 *                  output is.
	 * @param tight     sets whether the boundary of the hull polygon is kept tight
	 *                  to precisely the outer edges of the input polygons
	 * @return concave hull of the input shapes
	 * @since 1.2.1
	 */
	public static PShape concaveHull(PShape shapeSet, double concavity, boolean tight) {
		Geometry g = PGS_Conversion.fromPShape(shapeSet);
		if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)
				|| g.getGeometryType().equals(Geometry.TYPENAME_GEOMETRYCOLLECTION)) {
			g = g.union();
		}
		final ConcaveHullOfPolygons hull = new ConcaveHullOfPolygons(g);
		hull.setTight(tight);
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
		concavity *= concavity; // square to make output change more linearly as concavity goes 0...1
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
		threshold *= threshold * threshold;
		List<PVector> closestList = PGS_Optimisation.farthestPointPair(points);
		threshold *= closestList.get(0).dist(closestList.get(1));
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
	 * Computes a hull, having a variable level of convexity, of a shape.
	 * <p>
	 * When <code>convexity=0</code>, the original shape is reproduced (a pure
	 * concave hull); the hull tends towards a <i>convex hull</i> of the input as
	 * <code>convexity</code> goes to <code>1</code> (in other words, a hull with
	 * some level of "snapping" to the original shape).
	 * 
	 * @param shape
	 * @param convexity how convex the snap hull is, where 0 is the original shape
	 *                  (no convexity, fully concave) and 1 forms the convex hull of
	 *                  the shape's vertices
	 * @return
	 */
	public static PShape snapHull(PShape shape, double convexity) {
		convexity = Math.max(Math.min(convexity, 1), 0); // constrain 0...1
		ConcaveHullOfPolygons hull = new ConcaveHullOfPolygons(PGS_Conversion.fromPShape(shape)); // union in case multipolygon
		hull.setMaximumEdgeLengthRatio(convexity);
		return toPShape(hull.getHull());
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

}
