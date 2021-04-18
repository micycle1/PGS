package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import org.geotools.geometry.jts.JTS;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.simplify.VWSimplifier;

import micycle.pgs.utility.CornerRounding;
import micycle.pgs.utility.GaussianLineSmoothing;
import processing.core.PShape;
import uk.osgb.algorithm.minkowski_sum.Minkowski_Sum;

/**
 * Methods that affect the geometry or topology of shapes.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Morphology {

	static {
		Minkowski_Sum.setGeometryFactory(PGS.GEOM_FACTORY);
	}

	private PGS_Morphology() {
	}

	/**
	 * Computes a buffer area around the shape, having the given buffer width.
	 * 
	 * @param shape
	 * @param buffer extent/width of the buffer (may be positive or negative)
	 * @return
	 */
	public static PShape buffer(PShape shape, double buffer) {
		return toPShape(fromPShape(shape).buffer(buffer, 4));
	}

	/**
	 * small edges are removed, while the general structure of the shape is
	 * preserved. This process is known as "opening" in computer vision.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape erosionDilation(PShape shape, double buffer) {
		return toPShape(fromPShape(shape).buffer(-buffer).buffer(buffer));
	}

	/**
	 * Simplifies a shape using the Douglas-Peucker algorithm, reducing the
	 * complexity and number of vertices of the shape.
	 * <p>
	 * During the process shapes can be split, collapse to lines or disappear. Holes
	 * can be created or disappear.
	 * 
	 * @param shape
	 * @param distanceTolerance the tolerance to use
	 * @return simplifed copy of the shape
	 * @see #simplifyVW(PShape, double)
	 * @see #simplifyTopology(PShape, double)
	 */
	public static PShape simplify(PShape shape, double distanceTolerance) {
		return toPShape(DouglasPeuckerSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Simplifies a shape using the Visvalingam-Whyatt area-based algorithm,
	 * reducing the complexity and number of vertices of the shape.
	 * 
	 * @param shape
	 * @param distanceTolerance The simplification tolerance is specified as a
	 *                          distance.This is converted to an area tolerance by
	 *                          squaring it.
	 * @return simplifed copy of the shape
	 * @see #simplify(PShape, double)
	 * @see #simplifyTopology(PShape, double)
	 */
	public static PShape simplifyVW(PShape shape, double distanceTolerance) {
		return toPShape(VWSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Simplifies a shape, whilst preserving the topological structure of the shape
	 * (holes, etc.).
	 * 
	 * @param shape
	 * @param distanceTolerance the tolerance to use
	 * @return simplifed copy of the shape
	 * @see #simplify(PShape, double)
	 * @see #simplifyVW(PShape, double)
	 */
	public static PShape simplifyTopology(PShape shape, double distanceTolerance) {
		return toPShape(TopologyPreservingSimplifier.simplify(fromPShape(shape), distanceTolerance));
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
	 * @param subtract
	 * @return
	 */
	public static PShape minkDifference(PShape source, PShape subtract) {
		Geometry sum = Minkowski_Sum.compMinkDiff(fromPShape(source), fromPShape(subtract), true, true);
		return toPShape(sum);
	}

	/**
	 * Smoothes a shape. The smoothing algorithm inserts new vertices which are
	 * positioned using Bezier splines. The output shape tends to be a little larger
	 * than the input.
	 * 
	 * @param shape
	 * @param fit   tightness of fit from 0 (loose) to 1 (tight)
	 * @return smoothed copy of the shape
	 * @see #smoothGaussian(PShape, double)
	 */
	public static PShape smooth(PShape shape, double fit) {
		return toPShape(JTS.smooth(fromPShape(shape), fit));
	}

	/**
	 * Smoothes a shape by applying a gaussian filter to vertex coordinates. At
	 * larger values, this morphs the input shape much more visually than
	 * {@link #smooth(PShape, double)}.
	 * 
	 * @param shape
	 * @param sigma The standard deviation of the gaussian kernel. Larger values
	 *              provide more smoothing.
	 * @return smoothed copy of the shape
	 * @see #smooth(PShape, double)
	 */
	public static PShape smoothGaussian(PShape shape, double sigma) {
		Geometry g = fromPShape(shape);
		if (g.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			// TODO support holes
			Polygon p = (Polygon) g;
			return toPShape(GaussianLineSmoothing.get(p.getExteriorRing(), Math.max(sigma, 1)));
		}
		if (g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING)
				|| g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			LineString l = (LineString) g;
			return toPShape(GaussianLineSmoothing.get(l, Math.max(sigma, 1)));
		}
		System.err.println(g.getGeometryType() + " are not supported for this method (yet).");
		return new PShape(PShape.GROUP);
	}

	/**
	 * Rounds the corners of a shape by substituting a circular arc for the corners.
	 * Each corner is rounded in proportion to the smallest length of its 2
	 * constituent lines.
	 * 
	 * @param shape
	 * @param extent 0...1 (where 0 is no rounding, and 1 is the maximum rounding
	 *               whilst keeping shape valid). Values greater than 1 are allowed
	 *               by output undefined results.
	 * @return
	 */
	public static PShape round(PShape shape, double extent) {
		return CornerRounding.round(shape, extent);
	}

}
