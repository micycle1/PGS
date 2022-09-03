package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.Arrays;
import java.util.List;

import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.buffer.VariableBuffer;
import org.locationtech.jts.shape.CubicBezierCurve;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.simplify.VWSimplifier;

import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.PGS_Contour.OffsetStyle;
import micycle.pgs.color.RGB;
import micycle.pgs.commons.ChaikinCut;
import micycle.pgs.commons.CornerRounding;
import micycle.pgs.commons.DiscreteCurveEvolution;
import micycle.pgs.commons.GaussianLineSmoothing;
import micycle.pgs.commons.ShapeInterpolation;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.minkowski_sum.Minkowski_Sum;
import micycle.uniformnoise.UniformNoise;

/**
 * Methods that affect the geometry or topology of shapes.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Morphology {

	static {
		Minkowski_Sum.setGeometryFactory(PGS.GEOM_FACTORY);
	}

	private PGS_Morphology() {
	}

	/**
	 * Computes a rounded buffer area around the shape, having the given buffer
	 * width.
	 * 
	 * @param shape
	 * @param buffer extent/width of the buffer (which may be positive or negative)
	 * @return a polygonal shape representing the buffer region (which may be empty)
	 * @see #buffer(PShape, double, OffsetStyle)
	 */
	public static PShape buffer(PShape shape, double buffer) {
		return toPShape(fromPShape(shape).buffer(buffer, 8));
	}

	/**
	 * Computes a buffer area around the shape, having the given buffer width and
	 * buffer style (either round, miter, bevel).
	 * 
	 * @param shape
	 * @param buffer extent/width of the buffer (which may be positive or negative)
	 * @return a polygonal shape representing the buffer region (which may be empty)
	 * @see #buffer(PShape, double)
	 * @since 1.2.1
	 */
	public static PShape buffer(PShape shape, double buffer, OffsetStyle bufferStyle) {
		Geometry g = fromPShape(shape);
		BufferParameters bufParams = new BufferParameters(8, BufferParameters.CAP_FLAT, bufferStyle.style,
				BufferParameters.DEFAULT_MITRE_LIMIT);
		BufferOp b = new BufferOp(g, bufParams);
		return toPShape(b.getResultGeometry(buffer));
	}

	/**
	 * Buffers a shape with a varying buffer distance (interpolated between a start
	 * distance and an end distance) along the shape's perimeter.
	 * 
	 * @param shape         a single polygon or lineal shape
	 * @param startDistance the starting buffer amount
	 * @param endDistance   the terminating buffer amount
	 * @return a polygonal shape representing the variable buffer region (which may
	 *         beempty)
	 * @since 1.2.1
	 */
	public static PShape variableBuffer(PShape shape, double startDistance, double endDistance) {
		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing(); // variable buffer applies to linestrings only
		}
		return toPShape(VariableBuffer.buffer(g, startDistance, endDistance));
	}

	/**
	 * Applies a negative followed by a positive buffer (in a single operation). The
	 * effect of which is Small edges are removed, while the general structure of
	 * the shape is preserved. This process is known as "opening" in computer
	 * vision.
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
	 * Simplifies a shape via <i>Discrete Curve Evolution</i>.
	 * <p>
	 * Discrete curve evolution involves a stepwise elimination of kinks that are
	 * least relevant to the shape of the polygonal curve. The relevance of kinks is
	 * intended to reflect their contribution to the overall shape of the polygonal
	 * curve.
	 * 
	 * @param shape          a polygonal or lineal shape. GROUP shapes are not
	 *                       supported.
	 * @param removeFraction the fraction of least relevant kinks/vertices to
	 *                       remove. 0...1
	 * @return simplifed copy of the input shape
	 * @since 1.2.1
	 */
	public static PShape simplifyDCE(PShape shape, double removeFraction) {
		removeFraction = 1 - removeFraction; // since dce class is preserve-based, not remove-based
		Geometry g = fromPShape(shape);
		if (g instanceof Polygon) {
			LinearRing[] rings = new LinearRingIterator(g).getLinearRings();
			LinearRing[] dceRings = new LinearRing[rings.length];
			for (int i = 0; i < rings.length; i++) {
				LinearRing ring = rings[i];
				DiscreteCurveEvolution dce = new DiscreteCurveEvolution(
						Math.max(4, (int) Math.round(removeFraction * ring.getNumPoints())));
				dceRings[i] = PGS.GEOM_FACTORY.createLinearRing(dce.process(ring));
			}
			LinearRing[] holes = null;
			if (dceRings.length > 1) {
				holes = Arrays.copyOfRange(dceRings, 1, dceRings.length);
			}
			return toPShape(PGS.GEOM_FACTORY.createPolygon(dceRings[0], holes));
		} else if (g instanceof LineString) {
			LineString l = (LineString) g;
			DiscreteCurveEvolution dce = new DiscreteCurveEvolution(Math.max(4, (int) Math.round(removeFraction * l.getNumPoints())));
			return toPShape(PGS.GEOM_FACTORY.createLineString(dce.process(l)));
		} else {
			System.err.println(g.getGeometryType() + " are not supported for the simplifyDCE() method (yet).");
			return shape;
		}
	}

	/**
	 * Computes a <i>Minkowski sum</i> (a.k.a dilation) of the two source shapes.
	 * <p>
	 * To instill you with intuition of what a Minkowski sum looks like, here are a
	 * few examples:
	 * <ul>
	 * <li>The sum of any shape and a point is that shape translated by that point.
	 * <li>The sum of any shape and two points is two translated (possibly
	 * overlapping) copies of that shape.
	 * <li>The sum of two circles is a larger circle (sum the radii) with its centre
	 * at the sum of the centres of the smaller circles.
	 * <li>The sum of any shape and a line is that shape swept through that line.
	 * Think of placing your shape in sand, and dragging it along the line.
	 * <li>Similarly, the sum of a shape and any curve is what you’d get by sweeping
	 * the shape through the curve.
	 * <li>The sum of two parallel lines is a longer line.
	 * <li>For perpendicular lines, you get a square.</li>
	 * </ul>
	 * 
	 * @return shape representing the Minkowski sum of source+addition
	 * @see #minkDifference(PShape, PShape)
	 */
	public static PShape minkSum(PShape source, PShape addition) {
		// produces handled errors with geometries that have straight lines (like a
		// square)
		Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition), true, true);
		return toPShape(sum);
	}

	/**
	 * Computes a <i>Minkowski difference</i> (a.k.a erosion) of the two source
	 * shapes.
	 * 
	 * @return shape representing the Minkowski difference of source-subtract
	 * @see #minkSum(PShape, PShape)
	 */
	public static PShape minkDifference(PShape source, PShape subtract) {
		Geometry sum = Minkowski_Sum.minkDiff(fromPShape(source), fromPShape(subtract), true, true);
		return toPShape(sum);
	}

	/**
	 * Smoothes a shape. The smoothing algorithm inserts new vertices which are
	 * positioned using Bezier splines. The output shape tends to be a little larger
	 * than the input.
	 * 
	 * @param shape shape to smooth
	 * @param alpha curvedness parameter (0 is linear, 1 is round, >1 is
	 *              increasingly curved)
	 * @return smoothed copy of the shape
	 * @see #smoothGaussian(PShape, double)
	 */
	public static PShape smooth(PShape shape, double alpha) {
		Geometry curve = CubicBezierCurve.bezierCurve(fromPShape(shape), alpha);
		return toPShape(curve);
	}

	/**
	 * Smoothes a shape by applying a gaussian filter to vertex coordinates. At
	 * larger values, this morphs the input shape much more visually than
	 * {@link #smooth(PShape, double) smooth()}.
	 * 
	 * @param shape the shape to smooth
	 * @param sigma The standard deviation of the gaussian kernel. Larger values
	 *              provide more smoothing.
	 * @return smoothed copy of the shape
	 * @see #smooth(PShape, double)
	 */
	public static PShape smoothGaussian(PShape shape, double sigma) {
		Geometry g = fromPShape(shape);
		if (g.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			LinearRingIterator lri = new LinearRingIterator(g);
			LineString[] rings = lri.getLinearRings();
			LinearRing[] ringSmoothed = new LinearRing[rings.length];
			for (int i = 0; i < rings.length; i++) {
				Coordinate[] coords = GaussianLineSmoothing.get(rings[i], Math.max(sigma, 1), 1).getCoordinates();
				if (coords.length > 2) {
					ringSmoothed[i] = PGS.GEOM_FACTORY.createLinearRing(coords);
				} else {
					ringSmoothed[i] = PGS.GEOM_FACTORY.createLinearRing();
				}
			}

			LinearRing[] holes = null;
			if (ringSmoothed.length > 1) {
				holes = Arrays.copyOfRange(ringSmoothed, 1, ringSmoothed.length);
			}
			return toPShape(PGS.GEOM_FACTORY.createPolygon(ringSmoothed[0], holes));
		}
		if (g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) || g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			LineString l = (LineString) g;
			return toPShape(GaussianLineSmoothing.get(l, Math.max(sigma, 1), 1));
		}
		System.err.println(g.getGeometryType() + " are not supported for the smoothGaussian() method (yet).");
		return shape;
	}

	/**
	 * Rounds the corners of a shape by substituting a circular arc for each corner.
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

	/**
	 * Smoothes a shape via iterated corner cutting (chaikin cutting). More
	 * iterations results in more smoothing.
	 * 
	 * @param shape
	 * @param ratio      Between 0...1. Determines how far along each edge to
	 *                   perform the cuts. 0 is no cutting; 1 is maximal cutting
	 *                   (cut at the midpoint of each edge).
	 * @param iterations number of cutting iterations/recursions to perform. A value
	 *                   of 1 simply cuts the corners; higher values effectively
	 *                   smooth the cut. Values greater than ~10 generally have no
	 *                   additional effect.
	 * @return a cut copy of the input shape
	 * @since 1.1.0
	 */
	public static PShape chaikinCut(PShape shape, double ratio, int iterations) {
		ratio = Math.max(ratio, 0.0001);
		ratio = Math.min(ratio, 0.9999);
		ratio /= 2; // constrain to 0...0.5
		PShape cut = ChaikinCut.chaikin(shape, (float) ratio, iterations);
		PGS_Conversion.setAllFillColor(cut, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(cut, RGB.PINK, 3);
		return cut;
	}

	/**
	 * Warps/perturbs a shape by displacing vertices along a line between each
	 * vertex and the shape centroid.
	 * 
	 * <p>
	 * Inputs may be densified before warping.
	 * 
	 * @param shape      a polygonal shape
	 * @param magnitude  magnitude of the displacement. The value defines the
	 *                   maximum euclidean displacement of a vertex compared to the
	 *                   shape centroid
	 * @param warpOffset offset angle that determines at which angle to begin the
	 *                   displacement.
	 * @param densify    whether to densify the shape (using distance=1) before
	 *                   warping. When true, shapes with long edges will undergo
	 *                   warping along the whole edge (rather than only at the
	 *                   original vertices).
	 * @return
	 */
	public static PShape radialWarp(PShape shape, double magnitude, double warpOffset, boolean densify) {
		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			System.err.println("radialWarp() expects (single) polygon input. The geometry resolved to a " + g.getGeometryType());
			return shape;
		}

		final Point point = g.getCentroid();
		final PVector c = new PVector((float) point.getX(), (float) point.getY());

		final List<PVector> coords;

		if (densify) {
			final Densifier d = new Densifier(fromPShape(shape));
			d.setDistanceTolerance(1);
			d.setValidate(false);
			coords = PGS_Conversion.toPVector(toPShape(d.getResultGeometry()));
		} else {
			coords = PGS_Conversion.toPVector(shape);
		}

		final UniformNoise noise = new UniformNoise(1337);
		coords.forEach(coord -> {
			PVector heading = PVector.sub(coord, c); // vector from center to each vertex
			final double angle = heading.heading() + warpOffset;
			float perturbation = noise.uniformNoise(Math.cos(angle), Math.sin(angle));
			perturbation -= 0.5f; // [0...1] -> [-0.5...0.5]
			perturbation *= magnitude * 2;
			coord.add(heading.normalize().mult(perturbation)); // add perturbation to vertex
		});
		return PGS_Conversion.fromPVector(coords);
	}

	/**
	 * Warps/perturbs a shape by displacing vertices according to a 2D noise vector
	 * field.
	 * 
	 * <p>
	 * Inputs may be densified before warping.
	 * 
	 * @param shape      a polygonal shape
	 * @param magnitude  magnitude of the displacement (acting as noise value
	 *                   multiplier). The value defines the maximum displacement of
	 *                   a vertex in the both x and y axes.
	 * @param noiseScale the scale of the 2D noise vector field. This affects how of
	 *                   the coarseness of warping. Smaller values (~0.2) lead to
	 *                   more fine warping (at edges), whereas larger values (~2)
	 *                   affect the shape geometry at a larger scale.
	 * @param densify    whether to densify the shape (using distance=1) before
	 *                   warping. When true, shapes with long edges will undergo
	 *                   warping along the whole edge (rather than only at the
	 *                   original vertices).
	 * @return
	 * @see #fieldWarp(PShape, double, double, double, boolean, int)
	 */
	public static PShape fieldWarp(PShape shape, double magnitude, double noiseScale, boolean densify) {
		return fieldWarp(shape, magnitude, noiseScale, 0, densify, 1337);
	}

	/**
	 * Warps/perturbs a shape by displacing vertices according to a 2D noise vector
	 * field.
	 * 
	 * <p>
	 * Inputs may be densified before warping.
	 * 
	 * @param shape      a polygonal shape
	 * @param magnitude  magnitude of the displacement (acting as noise value
	 *                   multiplier). The value defines the maximum displacement of
	 *                   a vertex in the both x and y axes.
	 * @param noiseScale the scale of the 2D noise vector field. This affects how of
	 *                   the coarseness of warping. Smaller values (~0.2) lead to
	 *                   more fine warping (at edges), whereas larger values (~2)
	 *                   affect the shape geometry at a larger scale.
	 * @param time       used to offset the underlying noise field and hence animate
	 *                   the warping over time
	 * @param densify    whether to densify the shape (using distance=1) before
	 *                   warping. When true, shapes with long edges will undergo
	 *                   warping along the whole edge (rather than only at the
	 *                   original vertices).
	 * @param noiseSeed  a seed to pass to the underlying noise generator
	 * 
	 * @see #fieldWarp(PShape, double, double, boolean)
	 * @return
	 */
	public static PShape fieldWarp(PShape shape, double magnitude, double noiseScale, double time, boolean densify, int noiseSeed) {
		float scale = Math.max(1, (float) noiseScale * 500f);
		final boolean pointsShape = shape.getKind() == PConstants.POINTS;

		final PShape copy;
		if (densify && !pointsShape) {
			final Densifier d = new Densifier(fromPShape(shape));
			d.setDistanceTolerance(1);
			d.setValidate(false);
			copy = toPShape(d.getResultGeometry());
		} else {
			copy = toPShape(fromPShape(shape));
		}

		final UniformNoise noise = new UniformNoise(noiseSeed);

		if (copy.getChildCount() == 0) {
			// setVertex() will act on group shapes, so treat a single shape as group of 1
			copy.addChild(copy);
		}

		for (PShape child : copy.getChildren()) {
			for (int i = 0; i < child.getVertexCount(); i++) {
				final PVector coord = child.getVertex(i);
				float dx = noise.uniformNoise(coord.x / scale, coord.y / scale + time) - 0.5f;
				float dy = noise.uniformNoise(coord.x / scale + (101 + time), coord.y / scale + (101 + time)) - 0.5f;
				child.setVertex(i, coord.x + (dx * (float) magnitude * 2), coord.y + (dy * (float) magnitude * 2));
			}
		}

		if (pointsShape) {
			return copy;
		} else {
			if (copy.getChildCount() == 1) {
				return toPShape(GeometryFixer.fix(fromPShape(copy.getChild(0))));
			} else {
				// don't apply geometryFixer to GROUP shape, since fixing a multigeometry
				// appears to merge shapes. TODO apply .fix() to shapes individually
				return copy;
			}
		}
	}

	/**
	 * Generates an intermediate shape between two shapes by interpolating between
	 * them. This process has many names: shape morphing / blending / averaging /
	 * tweening / interpolation.
	 * <p>
	 * The underlying technique rotates one of the shapes to minimise the total
	 * distance between each shape's vertices, then performs linear interpolation
	 * between vertices. This performs well in practice but the outcome worsens as
	 * shapes become more concave; more sophisticated techniques would employ some
	 * level of rigidity preservation.
	 * 
	 * @param from                a single polygon; the shape we want to morph from
	 * @param to                  a single polygon; the shape we want to morph
	 *                            <code>from</code> into
	 * @param interpolationFactor between 0...1
	 * @return a polygonal PShape
	 * @since 1.2.0
	 * @see #interpolate(PShape, PShape, int)
	 */
	public static PShape interpolate(PShape from, PShape to, double interpolationFactor) {
		final Geometry toGeom = fromPShape(to);
		final Geometry fromGeom = fromPShape(from);
		if (toGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON) && fromGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			final ShapeInterpolation tween = new ShapeInterpolation(toGeom, fromGeom);
			return toPShape(PGS.GEOM_FACTORY.createPolygon(tween.tween(interpolationFactor)));
		} else {
			System.err.println("morph() accepts holeless single polygons only (for now).");
			return from;
		}
	}

	/**
	 * Generates intermediate shapes (frames) between two shapes by interpolating
	 * between them. This process has many names: shape morphing / blending /
	 * averaging / tweening / interpolation.
	 * <p>
	 * This method is faster than calling
	 * {@link #interpolate(PShape, PShape, double) interpolate()} repeatedly for
	 * different interpolation factors.
	 * 
	 * @param from   a single polygon; the shape we want to morph from
	 * @param to     a single polygon; the shape we want to morph <code>from</code>
	 *               into
	 * @param frames the number of frames (including first and last) to generate. >=
	 *               2
	 * @return a GROUP PShape, where each child shape is a frame from the
	 *         interpolation
	 * @since 1.2.1
	 * @see #interpolate(PShape, PShape, double)
	 */
	public static PShape interpolate(PShape from, PShape to, int frames) {
		final Geometry toGeom = fromPShape(to);
		final Geometry fromGeom = fromPShape(from);
		if (toGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON) && fromGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			final ShapeInterpolation tween = new ShapeInterpolation(toGeom, fromGeom);
			final float fraction = 1f / (frames - 1);
			PShape out = new PShape();
			for (int i = 0; i < frames; i++) {
				out.addChild(toPShape(PGS.GEOM_FACTORY.createPolygon(tween.tween(fraction * i))));
			}
			return out;
		} else {
			System.err.println("morph() accepts holeless single polygons only (for now).");
			return from;
		}
	}

}
