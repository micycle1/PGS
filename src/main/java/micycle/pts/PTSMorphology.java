package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.geodelivery.jap.concavehull.SnapHull;
import org.geotools.geometry.jts.JTS;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.simplify.VWSimplifier;

import micycle.pts.utility.PolygonDecomposition;
import processing.core.PApplet;
import processing.core.PConstants;
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
		// produces handled errors with geometries that have straight lines (like a
		// square)
		Geometry sum = Minkowski_Sum.compMinkSum(fromPShape(source), fromPShape(addition), true, true);
		return toPShape(sum);
	}

	/**
	 * Minkowski difference a.k.a erosion
	 * 
	 * @param source
	 * @param addition
	 * @return
	 */
	public static PShape minkDifference(PShape source, PShape addition) {
		Geometry sum = Minkowski_Sum.compMinkDiff(fromPShape(source), fromPShape(addition), true, true);
		return toPShape(sum);
	}

	/**
	 * Smoothes a geometry. The smoothing algorithm inserts new vertices which are
	 * positioned using Bezier splines.
	 * 
	 * @param shape
	 * @param fit   tightness of fit from 0 (loose) to 1 (tight)
	 * @return
	 */
	public static PShape smooth(PShape shape, float fit) {
		return toPShape(JTS.smooth(fromPShape(shape), fit));
	}

	/**
	 * Splits a shape into 4 equal quadrants
	 * 
	 * @param shape
	 * @return list containing the 4 split quadrants
	 */
	public static List<PShape> split(PShape shape) {
		// https://stackoverflow.com/questions/64252638/how-to-split-a-jts-polygon
		Geometry p = fromPShape(shape);
		ArrayList<PShape> ret = new ArrayList<>();

		final Envelope envelope = p.getEnvelopeInternal();
		double minX = envelope.getMinX();
		double maxX = envelope.getMaxX();
		double midX = minX + (maxX - minX) / 2.0;
		double minY = envelope.getMinY();
		double maxY = envelope.getMaxY();
		double midY = minY + (maxY - minY) / 2.0;

		Envelope llEnv = new Envelope(minX, midX, minY, midY);
		Envelope lrEnv = new Envelope(midX, maxX, minY, midY);
		Envelope ulEnv = new Envelope(minX, midX, midY, maxY);
		Envelope urEnv = new Envelope(midX, maxX, midY, maxY);
		Geometry UL = JTS.toGeometry(llEnv).intersection(p);
		Geometry UR = JTS.toGeometry(lrEnv).intersection(p);
		Geometry LL = JTS.toGeometry(ulEnv).intersection(p);
		Geometry LR = JTS.toGeometry(urEnv).intersection(p);
		ret.add(toPShape(UL));
		ret.add(toPShape(UR));
		ret.add(toPShape(LL));
		ret.add(toPShape(LR));

		return ret;
	}

	/**
	 * Splits a shape into 4^(1+recursions) rectangular partitions
	 * 
	 * @param shape
	 * @param splitDepth
	 */
	public static List<PShape> split(final PShape shape, int splitDepth) {
		splitDepth = Math.max(0, splitDepth);
		ArrayDeque<Geometry> stack = new ArrayDeque<>();
		stack.add(fromPShape(shape));

		ArrayList<PShape> ret = new ArrayList<>(); // add to when recursion depth reached
		ArrayList<Geometry> next = new ArrayList<>(); // add to when recursion depth reached

		int depth = 0;
		while (depth < splitDepth) {
			while (!stack.isEmpty()) {
				final Geometry slice = stack.pop();
				final Envelope envelope = slice.getEnvelopeInternal();
				final double minX = envelope.getMinX();
				final double maxX = envelope.getMaxX();
				final double midX = minX + (maxX - minX) / 2.0;
				final double minY = envelope.getMinY();
				final double maxY = envelope.getMaxY();
				final double midY = minY + (maxY - minY) / 2.0;

				Envelope llEnv = new Envelope(minX, midX, minY, midY);
				Envelope lrEnv = new Envelope(midX, maxX, minY, midY);
				Envelope ulEnv = new Envelope(minX, midX, midY, maxY);
				Envelope urEnv = new Envelope(midX, maxX, midY, maxY);
				Geometry UL = JTS.toGeometry(llEnv).intersection(slice);
				Geometry UR = JTS.toGeometry(lrEnv).intersection(slice);
				Geometry LL = JTS.toGeometry(ulEnv).intersection(slice);
				Geometry LR = JTS.toGeometry(urEnv).intersection(slice);
				next.add(UL);
				next.add(UR);
				next.add(LL);
				next.add(LR);
			}
			depth++;
			stack.addAll(next);
			next.clear();
		}

		stack.forEach(g -> ret.add(toPShape(g)));
		return ret;
	}

	/**
	 * Partitions a shape into simple polygons using Mark Bayazit's algorithm.
	 * 
	 * @param shape
	 * @return list of convex (simple) polygons comprising the original shape
	 */
	public static List<PShape> partition(PShape shape) {
		// https://mpen.ca/406/bayazit
		// retry GreedyPolygonSplitter()?

		Geometry g = fromPShape(shape);

		ArrayList<PShape> out = new ArrayList<>();

		for (int i = 0; i < g.getNumGeometries(); i++) {
			Geometry child = g.getGeometryN(i);
			if (child.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) { // skip any linestrings etc
				List<Polygon> decomposed = PolygonDecomposition.decompose((Polygon) child);
				for (Polygon polygon : decomposed) {
					out.add(toPShape(polygon));
				}
			}
		}

		return out;
	}

	/**
	 * Slice a shape using a line given by its start and endpoints.
	 * 
	 * @param shape
	 * @param p1    must be outside shape
	 * @param p2    must be outside shape
	 * @return a list containg two PShapes
	 */
	public static List<PShape> slice(PShape shape, PVector p1, PVector p2) {
		// adapted from https://gis.stackexchange.com/questions/189976/
		final Geometry poly = fromPShape(shape);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(poly);
		final LineSegment ls = new LineSegment(p1.x, p1.y, p2.x, p2.y);
		final LineString line = ls.toGeometry(PTS.GEOM_FACTORY);
		final Geometry nodedLinework = poly.getBoundary().union(line);
		final Geometry polys = polygonize(nodedLinework);

		final ArrayList<Polygon> leftSlices = new ArrayList<>();
		final ArrayList<Polygon> rightSlices = new ArrayList<>();

		for (int i = 0; i < polys.getNumGeometries(); i++) {
			final Polygon candpoly = (Polygon) polys.getGeometryN(i);
			if (cache.contains(candpoly.getInteriorPoint())) {
				if (ls.orientationIndex(candpoly.getCentroid().getCoordinate()) == Orientation.LEFT) {
					leftSlices.add(candpoly);
				} else {
					rightSlices.add(candpoly);
				}
			}
		}

		ArrayList<PShape> output = new ArrayList<>();
		output.add(toPShape(UnaryUnionOp.union(leftSlices)));
		output.add(toPShape(UnaryUnionOp.union(rightSlices)));
		return output;
	}

	/**
	 * Used by slice()
	 */
	@SuppressWarnings("unchecked")
	private static Geometry polygonize(Geometry geometry) {
		List<LineString> lines = LineStringExtracter.getLines(geometry);
		Polygonizer polygonizer = new Polygonizer();
		polygonizer.add(lines);
		Collection<Polygon> polys = polygonizer.getPolygons();
		Polygon[] polyArray = GeometryFactory.toPolygonArray(polys);
		return geometry.getFactory().createGeometryCollection(polyArray);
	}

}
