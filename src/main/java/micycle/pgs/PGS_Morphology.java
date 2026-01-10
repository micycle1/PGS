package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Lineal;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.buffer.VariableBuffer;
import org.locationtech.jts.precision.GeometryPrecisionReducer;
import org.locationtech.jts.shape.CubicBezierCurve;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.simplify.VWSimplifier;

import com.gihub.micycle1.malleo.Malleo;

import micycle.pgs.PGS_Contour.OffsetStyle;
import micycle.pgs.commons.ChaikinCut;
import micycle.pgs.commons.ContourRegularization;
import micycle.pgs.commons.ContourRegularization.Parameters;
import micycle.pgs.commons.CornerRounding;
import micycle.pgs.commons.CornerRounding.RoundingStyle;
import micycle.pgs.commons.DiscreteCurveEvolution;
import micycle.pgs.commons.DiscreteCurveEvolution.DCETerminationCallback;
import micycle.pgs.commons.EllipticFourierDesc;
import micycle.pgs.commons.FastAtan2;
import micycle.pgs.commons.GaussianLineSmoothing;
import micycle.pgs.commons.LaneRiesenfeldSmoothing;
import micycle.pgs.commons.NewtonThieleRingMorpher;
import micycle.pgs.commons.SchneiderBezierFitter;
import micycle.uniformnoise.UniformNoise;
import net.jafama.FastMath;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.minkowski_sum.MinkowskiSum;

/**
 * Morphological editing operations for {@link PShape} polygons.
 *
 * <p>
 * This class hosts algorithms that <em>reshape</em> geometry, typically by
 * offsetting, simplifying, smoothing, warping, or deforming outlines; often
 * changing vertex count and sometimes changing topology (splitting/merging
 * parts, creating/removing holes).
 *
 * @author Michael Carleton
 */
public final class PGS_Morphology {

	static {
		MinkowskiSum.setGeometryFactory(PGS.GEOM_FACTORY);
	}

	private PGS_Morphology() {
	}

	/**
	 * Returns a rounded buffer region of the given shape at the specified distance.
	 * <p>
	 * The distance is in the same coordinate units as the shape: positive values
	 * expand the shape, negative values contract it. The returned shape is a
	 * polygonal PShape and may be empty. The input shape is not modified.
	 *
	 * @param shape  the source shape to buffer
	 * @param buffer distance (extent/width) of the buffer; may be positive or
	 *               negative
	 * @return a polygonal PShape representing the buffer region (may be empty)
	 * @see #buffer(PShape, double, OffsetStyle)
	 */
	public static PShape buffer(PShape shape, double buffer) {
		return buffer(shape, buffer, OffsetStyle.ROUND);
	}

	/**
	 * Returns a buffer region of the given shape using the specified join style.
	 * <p>
	 * The distance is in the same coordinate units as the shape: positive values
	 * expand the shape, negative values contract it. The bufferStyle controls how
	 * corners are joined (e.g. ROUND, MITER, BEVEL). The returned shape is a
	 * polygonal PShape and may be empty. The input shape is not modified.
	 *
	 * @param shape       the source shape to buffer
	 * @param buffer      distance (extent/width) of the buffer; may be positive or
	 *                    negative
	 * @param bufferStyle how to join offset segments (ROUND, MITER, BEVEL)
	 * @return a polygonal PShape representing the buffer region (may be empty)
	 * @see #buffer(PShape, double)
	 * @since 1.3.0
	 */
	public static PShape buffer(PShape shape, double buffer, OffsetStyle bufferStyle) {
		return buffer(shape, buffer, bufferStyle, CapStyle.ROUND);
	}

	/**
	 * Returns a buffer region of the given shape using the specified join and cap
	 * styles.
	 * <p>
	 * The distance is in the same coordinate units as the shape: positive values
	 * expand the shape, negative values contract it. bufferStyle controls how
	 * corners are joined; capStyle controls the end-cap style used for open
	 * geometries. The input shape is not modified; the returned PShape preserves
	 * the user data from the original geometry. The result is a polygonal PShape
	 * and may be empty.
	 *
	 * @param shape       the source shape to buffer
	 * @param buffer      distance (extent/width) of the buffer; may be positive or
	 *                    negative
	 * @param bufferStyle how to join offset segments (ROUND, MITER, BEVEL)
	 * @param capStyle    how to draw end caps for open geometries (e.g. ROUND,
	 *                    FLAT)
	 * @return a polygonal PShape representing the buffer region (may be empty)
	 * @since 2.1
	 */
	public static PShape buffer(PShape shape, double buffer, OffsetStyle bufferStyle, CapStyle capStyle) {
		Geometry g = fromPShape(shape);
		BufferParameters bufParams = createBufferParams(buffer, 0.5, bufferStyle, capStyle);
		BufferOp b = new BufferOp(g, bufParams);
		var out = b.getResultGeometry(buffer);
		out.setUserData(g.getUserData());
		return toPShape(out);
	}

	/**
	 * Buffers a shape with a varying buffer distance (interpolated between a start
	 * distance and an end distance) along the shape's perimeter.
	 * 
	 * @param shape         a polygon, lineal shape, or GROUP containing such shapes
	 * @param startDistance the starting buffer amount
	 * @param endDistance   the terminating buffer amount
	 * @return a polygonal shape representing the variable buffer region (which may
	 *         be empty)
	 * @since 1.3.0
	 */
	public static PShape variableBuffer(PShape shape, double startDistance, double endDistance) {
		return PGS.applyToLinealGeometries(shape, line -> {
			var buffer = (Polygon) VariableBuffer.buffer(line, startDistance, endDistance);
			return buffer.getExteriorRing();
		});
	}

	/**
	 * Applies a variable buffer to a shape. The buffer width at each vertex is
	 * determined by a callback function that considers the vertex's properties and
	 * its relative position along the shape's boundary.
	 * <p>
	 * Example usage:
	 * 
	 * <pre>
	 * {@code
	 * PShape bufferedShape = bufferWithCallback(inputShape, (coordinate, fraction) -> {
	 * 	// Example logic: buffer width decreases linearly from 10 units at the
	 * 	// start to 1 unit at the end
	 * 	return 10 - (fraction * (10 - 1));
	 * });
	 * }
	 * </pre>
	 *
	 * @param shape          A single polygon, lineal shape, or GROUP containing
	 *                       such shapes
	 * @param bufferCallback A callback function that receives the vertex coordinate
	 *                       and a double representing fractional distance (0...1)
	 *                       of the vertex along the shape's boundary. The function
	 *                       may use properties of the vertex, or its position, to
	 *                       determine the buffer width at that point.
	 * @return A new shape representing the original shape buffered with variable
	 *         widths as specified by the callback function. The width at each
	 *         vertex is calculated independently.
	 * @since 2.0
	 */
	public static PShape variableBuffer(PShape shape, BiFunction<Coordinate, Double, Double> bufferCallback) {
		return PGS.applyToLinealGeometries(shape, line -> {
			final Coordinate[] coords = line.getCoordinates();
			if (coords.length == 0) {
				// return an "empty buffer" geometry consistent with VariableBuffer expectations
				return null;
			}

			final double totalLength = line.getLength();
			final double[] bufferDistances = new double[coords.length];

			// Guard against degenerate/zero-length lines (all points same).
			if (totalLength == 0) {
				final double d0 = bufferCallback.apply(coords[0], 0.0);
				for (int i = 0; i < bufferDistances.length; i++) {
					bufferDistances[i] = d0;
				}
			} else {
				bufferDistances[0] = bufferCallback.apply(coords[0], 0.0);

				double runningLength = 0;
				Coordinate prev = coords[0];

				for (int i = 1; i < coords.length; i++) {
					runningLength += prev.distance(coords[i]);
					final double fractionalDistance = runningLength / totalLength; // 0..1
					bufferDistances[i] = bufferCallback.apply(coords[i], fractionalDistance);
					prev = coords[i];
				}
			}

			final var vb = new VariableBuffer(line, bufferDistances);
			var buffer = (Polygon) vb.getResult();
			return buffer.getExteriorRing();
		});
	}

	/**
	 * Erodes (a negative buffer) a shape by a normalised amount (scaled to shape
	 * size).
	 * <p>
	 * {@code amount} is dimensionless: {@code amount == 1} corresponds to a full
	 * erosion (approximately to the maximum inscribed radius), often collapsing
	 * polygons to empty. {@code shape} may be a {@code GROUP}; each polygonal
	 * element is processed independently. The sign of {@code amount} is ignored
	 * (always erodes).
	 *
	 * @param shape  the source shape (polygonal or {@code GROUP})
	 * @param amount normalised erosion amount (dimensionless)
	 * @return a polygonal {@code PShape} of the eroded geometry (may be empty)
	 * @since 2.2
	 */
	public static PShape normalisedErosion(PShape shape, double amount) {
		double amt = -Math.abs(amount); // force erosion
		var polys = PGS.extractPolygons(fromPShape(shape));
		var buffered = polys.parallelStream().map(p -> {
			var mic = new MaximumInscribedCircle(p, 0.5);
			var r = mic.getRadiusLine().getLength() * (1 + 1e-3);
			var buffer = amt * r;
			var bufParams = createBufferParams(buffer, 0.5, OffsetStyle.ROUND, CapStyle.ROUND);
			BufferOp b = new BufferOp(p, bufParams);
			var out = b.getResultGeometry(buffer);
			return out;
		}).toList();

		return toPShape(buffered);
	}

	/**
	 * Applies a negative followed by a positive buffer (in a single operation), the
	 * effect of which is small edges/islands are removed, while the general
	 * structure of the shape is preserved.
	 * <p>
	 * This operation is known as "opening" in computer vision.
	 * 
	 * @param shape  polygonal shape
	 * @param buffer a positive number
	 * @return a polygonal {@code PShape} of the dilated geometry (may be empty)
	 * @see #dilationErosion(PShape, double)
	 */
	public static PShape erosionDilation(PShape shape, double buffer) {
		buffer = Math.abs(buffer);

		final int segments = (int) Math.ceil(BufferParameters.DEFAULT_QUADRANT_SEGMENTS + Math.sqrt(buffer));
		var in = fromPShape(shape);
		Geometry g = BufferOp.bufferOp(in, -buffer, segments);
		g = BufferOp.bufferOp(g, +buffer, segments);
		g.setUserData(in.getUserData());

		try {
			return toPShape(g);
		} catch (Exception e) {
			return toPShape(GeometryFixer.fix(g));
		}
	}

	/**
	 * Applies a positive followed by a negative buffer (in a single operation), the
	 * effect of which is small holes and gaps are filled in, while the general
	 * structure of the shape is preserved.
	 * <p>
	 * This operation is known as "closing" in computer vision.
	 * 
	 * @param shape  polygonal shape
	 * @param buffer a positive number
	 * @since 1.3.0
	 * @see #erosionDilation(PShape, double)
	 */
	public static PShape dilationErosion(PShape shape, double buffer) {
		buffer = Math.abs(buffer);

		final int segments = (int) Math.ceil(BufferParameters.DEFAULT_QUADRANT_SEGMENTS + Math.sqrt(buffer));
		var in = fromPShape(shape);
		Geometry g = BufferOp.bufferOp(in, buffer, segments);
		g = BufferOp.bufferOp(g, -buffer, segments);
		g.setUserData(in.getUserData());

		try {
			return toPShape(g);
		} catch (Exception e) {
			return toPShape(GeometryFixer.fix(g));
		}
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
	 * @see #simplifyVW(PShape, double) simplifyVW()
	 * @see #simplifyTopology(PShape, double) simplifyTopology()
	 * @see {@link PGS_Meshing#simplifyMesh(PShape, double, boolean) simplifyMesh()}
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
	 * @see #simplify(PShape, double) simplify()
	 * @see #simplifyTopology(PShape, double) simplifyTopology()
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
	 * @see #simplify(PShape, double) simplify()
	 * @see #simplifyVW(PShape, double) simplifyVW()
	 */
	public static PShape simplifyTopology(PShape shape, double distanceTolerance) {
		return toPShape(TopologyPreservingSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Simplifies a shape via <i>Discrete Curve Evolution</i>.
	 * <p>
	 * This algorithm simplifies a shape by iteratively removing kinks from the
	 * shape, <b>starting with those having the least shape-relevance</b>.
	 * <p>
	 * The simplification process terminates according to a user-specified
	 * {@link DCETerminationCallback#shouldTerminate(Coordinate, double, int)
	 * callback} that decides whether the DCE algorithm should terminate based on
	 * the current kink (having a candidate vertex), using its: coordinate,
	 * relevance score, and the number of vertices remaining in the simplified
	 * geometry. Implementations can use this method to provide custom termination
	 * logic which may depend on various factors, such as a threshold relevance
	 * score, a specific number of vertices to preserve, or other criteria.
	 * <p>
	 * Note: the termination callback is applied per ring (boundary, hole, line,
	 * etc) in the input.
	 * 
	 * @param shape               The input shape to be simplified, which can be a
	 *                            polygonal (inclusive of holes) or a lineal shape.
	 *                            GROUP shapes are supported.
	 * @param terminationCallback The callback used to determine when the
	 *                            simplification process should terminate.
	 *                            {@code true} if the DCE process should terminate;
	 *                            {@code false} otherwise.
	 * @return A new, simplified copy of the input shape, with the least significant
	 *         kinks or vertices removed according to the provided fraction.
	 * @since 2.0
	 */
	public static PShape simplifyDCE(PShape shape, DCETerminationCallback terminationCallback) {
		return PGS.applyToLinealGeometries(shape, ring -> {
			return DiscreteCurveEvolution.process(ring, terminationCallback);
		});
	}

	/**
	 * Simplify the shape using DCE, removing vertices with relevance < r.
	 * 
	 * @param shape              the input shape
	 * @param relevanceThreshold the relevance threshold; only vertices with
	 *                           relevance >= the threshold will be kept. 20 is a
	 *                           good starting value for generally imperceptible
	 *                           simplification.
	 * @return the simplified PShape
	 * @since 2.1
	 */
	public static PShape simplifyDCE(PShape shape, final double relevanceThreshold) {
		return simplifyDCE(shape, (currentVertex, relevance, verticesRemaining) -> relevance >= relevanceThreshold);
	}

	/**
	 * Creates a <a href="https://github.com/micycle1/Hobby-Curves"><i>Hobby
	 * Curve</i></a> from the vertices of the shape. This tends to simplify/round
	 * the <b>geometry</b> of shape, but may actually increase the number of
	 * vertices due to increased curvature.
	 * <p>
	 * You may want to consider simplifying a shape (reducing vertex count) with
	 * other methods before applying Hobby simplification.
	 * 
	 * @param shape   vertices to use as basis for the Hobby Curve
	 * @param tension a parameter that controls the tension of the curve (how
	 *                tightly it is "pulled" towards underlying vertices). Suitable
	 *                domain is [0.666...3].
	 * @return a Hobby Curve
	 * @since 1.4.0
	 */
	public static PShape simplifyHobby(PShape shape, double tension) {
		return PGS.applyToLinealGeometries(shape, ring -> {
			var points = PGS_Conversion.toPVector(toPShape(ring));
			if (ring.isClosed() && !points.get(0).equals(points.get(points.size() - 1))) {
				points.add(points.get(0).copy());
			}
			var g = fromPShape(PGS_Construction.createHobbyCurve(points, tension));
			if (g instanceof Polygon) {
				g = ((Polygon) g).getExteriorRing();
			}
			if (g instanceof Lineal) {
				return (LineString) g;
			}
			return null;
		});
	}

	/**
	 * Computes a <i>Minkowski sum</i> (a.k.a dilation) of the two source shapes.
	 * The <code>addition</code> shape should probably be centered on (0,0) for best
	 * results.
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
		Geometry sum = MinkowskiSum.minkSum(fromPShape(source), fromPShape(addition), true, true);
		return toPShape(sum);
	}

	/**
	 * Computes a <i>Minkowski difference</i> (a.k.a erosion) of the two source
	 * shapes. The <code>subtract</code> shape should probably be centered on (0,0)
	 * for best results.
	 * 
	 * @return shape representing the Minkowski difference of source-subtract
	 * @see #minkSum(PShape, PShape)
	 */
	public static PShape minkDifference(PShape source, PShape subtract) {
		Geometry sum = MinkowskiSum.minkDiff(fromPShape(source), fromPShape(subtract), true, true);
		return toPShape(sum);
	}

	/**
	 * Smoothes a shape. The smoothing algorithm inserts new vertices which are
	 * positioned using Bezier splines. The output shape tends to be a little larger
	 * than the input.
	 * <p>
	 * Note: this method effectively constructs a Bezier curve through the existing
	 * vertices. As a result, if the input geometry already has very dense / closely
	 * spaced vertices, the smoothing may have little or no perceptual effect. This
	 * differs from other smoothing approaches (e.g. Gaussian) that operate at a
	 * spatial scale and are therefore largely invariant to vertex density.
	 * </p>
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
	 * Smoothes a shape by <em>fitting</em> one or more cubic Bezier curve segments
	 * to each lineal component (polylines and polygon rings), then
	 * <em>resampling</em> the fitted Beziers to produce a new vertex sequence.
	 * <p>
	 * This method uses Philip J. Schneider’s curve fitting algorithm. Unlike
	 * {@link #smooth(PShape, double) smooth()}, which constructs a Bezier curve
	 * <em>through</em> the existing vertices, this method approximates the input
	 * within a user-specified tolerance and can substantially simplify noisy or
	 * densely-vertexed input while producing a visually smoother result.
	 * </p>
	 * <p>
	 * The {@code maxDeviation} parameter controls how closely the fitted Bezier(s)
	 * must follow the original polyline/ring: smaller values preserve the original
	 * shape more strictly (often producing more Bezier segments and/or more output
	 * vertices), while larger values allow a smoother, more generalised result.
	 * </p>
	 * <p>
	 * Implementation note: the fitted Bezier segments are sampled at a fixed
	 * spacing (currently 2 units in the coordinate system of the input geometry) to
	 * create the returned JTS geometry, which is then converted back to a
	 * {@link PShape}.
	 * </p>
	 *
	 * @param shape        shape whose lineal geometry (LineStrings and polygon
	 *                     rings) will be Bezier-fit and resampled
	 * @param maxDeviation maximum allowed deviation (error tolerance) between the
	 *                     input vertices and the fitted Bezier curve(s); must be
	 *                     {@code > 0}
	 * @return a smoothed copy of {@code shape} produced by piecewise cubic Bezier
	 *         fitting and resampling
	 *
	 * @since 2.2
	 * @see SchneiderBezierFitter
	 */
	public static PShape smoothBezierFit(PShape shape, double maxDeviation) {
		return PGS.applyToLinealGeometries(shape, ring -> {
			return SchneiderBezierFitter.fitAndSample(ring, maxDeviation, PGS_Conversion.BEZIER_SAMPLE_DISTANCE);
		});
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
	 * @see #smoothGaussianNormalised(PShape, double)
	 * @see #smooth(PShape, double)
	 */
	public static PShape smoothGaussian(PShape shape, double sigma) {
		return PGS.applyToLinealGeometries(shape, ring -> GaussianLineSmoothing.get(ring, sigma));
	}

	/**
	 * Applies Gaussian smoothing to each lineal geometry in a {@link PShape} using
	 * a normalised amount in [0..1], intended to be scale-invariant across children
	 * of different sizes. {@code amount=0} leaves geometry unchanged;
	 * {@code amount=1} collapses (per geometry) using the extreme-sigma fallback.
	 *
	 * @param shape  input shape
	 * @param amount normalised smoothing amount in [0..1]
	 * @return new shape with smoothed lineal components
	 * @see #smoothGaussian(PShape, double)
	 * @since 2.2
	 */
	public static PShape smoothGaussianNormalised(PShape shape, double amount) {
		return PGS.applyToLinealGeometries(shape, ring -> GaussianLineSmoothing.getNormalised(ring, amount));
	}

	/**
	 * Calculates the Elliptic Fourier Descriptors (EFD) of a specified shape,
	 * yielding a simplified/smoothed shape representation based on the specified
	 * number of descriptors.
	 * <p>
	 * The EFD technique is an approach for shape analysis and simplification that
	 * decomposes a shape into a sequence of elliptic harmonic components. These
	 * components encapsulate the contour details of the shape: lower-order
	 * harmonics capture the broad geometry of the shape, while higher-order
	 * harmonics register the detailed, high-frequency contour characteristics,
	 * analogous to Principal Component Analysis (PCA). This technique is
	 * particularly effective for generating condensed or smoother versions of
	 * complex shapes.
	 * 
	 * @param shape       A polygonal shape to be transformed using the EFD.
	 * @param descriptors The desired level of the EFD, denoting the quantity of
	 *                    harmonics to be retained in the output. The maximum value
	 *                    is half the total number of vertices in the shape, while
	 *                    the minimum allowable value is 2. As the number of
	 *                    harmonics is increased, the output tends towards the input
	 *                    shape.
	 * @return A new PShape, simplified through the application of the Elliptic
	 *         Fourier Descriptors up to the indicated order. This shape will always
	 *         have the same number of vertices as the original.
	 * @since 1.4.0
	 */
	public static PShape smoothEllipticFourier(PShape shape, int descriptors) {
		return PGS.applyToLinealGeometries(shape, ring -> {
			int descriptorz = Math.min(ring.getCoordinates().length / 2, descriptors); // max=#vertices/2
			descriptorz = Math.max(2, descriptorz); // min=2
			if (ring.isClosed()) {
				final EllipticFourierDesc efd = new EllipticFourierDesc((LinearRing) ring, descriptorz);
				Coordinate[] coords = efd.createPolygon();
				return PGS.GEOM_FACTORY.createLinearRing(coords);
			} else {
				return null; // open linestrings not supported
			}
		});
	}

	/**
	 * Smooths a shape using Lane-Riesenfeld curve subdivision with 4-point
	 * refinement to reduce contraction.
	 * 
	 * @param shape                 A shape having lineal geometries (polygons or
	 *                              linestrings). Can be a GROUP shape consisting of
	 *                              these.
	 * @param degree                The degree of the LR algorithm. Higher degrees
	 *                              influence the placement of vertices and the
	 *                              overall shape of the curve, but only slightly
	 *                              increase the number of vertices generated.
	 *                              Increasing the degree also increases the
	 *                              contraction of the curve toward its control
	 *                              points. The degree does not directly control the
	 *                              smoothness of the curve. A value of 3 or 4 is
	 *                              usually sufficient for most applications.
	 * @param subdivisions          The number of times the subdivision process is
	 *                              applied. More subdivisions result in finer
	 *                              refinement and visually smoother curves between
	 *                              vertices. A value of 3 or 4 is usually
	 *                              sufficient for most applications.
	 * @param antiContractionFactor The weight parameter for the 4-point refinement.
	 *                              Controls the interpolation strength. A value of
	 *                              0 effectively disables the contraction
	 *                              reduction. Generally suitable values are in
	 *                              [0...0.1]. Larger values may create
	 *                              self-intersecting geometry.
	 * @return A Shape having same structure as the input, whose geometries are now
	 *         smooth.
	 * @since 2.1
	 */
	public static PShape smoothLaneRiesenfeld(PShape shape, int degree, int subdivisions, double antiContractionFactor) {
		return PGS.applyToLinealGeometries(shape, lineal -> LaneRiesenfeldSmoothing.subdivide(lineal, degree, subdivisions, antiContractionFactor));
	}

	/**
	 * Rounds polygon corners by replacing each corner with a circular arc.
	 *
	 * <p>
	 * This processes only the linear content of the input <code>PShape</code> - the
	 * contour paths (closed contours for polygon exteriors and interior
	 * contours/holes, and open polylines). Non‑linear or unsupported children are
	 * ignored.
	 * </p>
	 * <p>
	 * The <code>radius</code> is nominal; it is clamped so it cannot exceed what
	 * adjacent edges can support (this prevents overlapping or invalid contours).
	 * </p>
	 *
	 * @param shape  a polygonal <code>PShape</code> or a <code>GROUP</code>
	 *               <code>PShape</code> containing polygonal children; holes
	 *               (interior contours) are supported
	 * @param radius nominal radius used to round corners; clamped by adjacent edge
	 *               lengths
	 * @return a non-null <code>PShape</code> with rounded corners; possibly an
	 *         empty <code>GROUP</code> when nothing remains
	 */
	public static PShape round(PShape shape, double radius) {
		return PGS.applyToLinealGeometries(shape, ring -> {
			var rounded = CornerRounding.roundCorners(toPShape(ring), radius, RoundingStyle.CIRCLE);
			var g = fromPShape(rounded);
			if (g instanceof Polygon) {
				g = ((Polygon) g).getExteriorRing();
			}
			if (g instanceof Lineal) {
				return (LineString) g;
			}
			return null; // pointal or other...
		});
	}

	/**
	 * Smoothes a shape by recursively cutting its corners, a technique introduced
	 * by George Chaikin in 1974.
	 * <p>
	 * This method can be used to generate smooth-looking curves from a limited
	 * number of points. More iterations result in more smoothing.
	 * 
	 * @param shape      The shape to be smoothed
	 * @param ratio      A ratio (between 0 and 1) determining how far along each
	 *                   edge to perform the two cuts. For example, a ratio of 0.5
	 *                   will cut the underlying edge twice, at 0.25x and 0.75x
	 *                   along its length. A value of 1 will cut each edge once,
	 *                   directly at its midpoint. It is recommended to use a value
	 *                   of 0.5 for this parameter.
	 * @param iterations The number of cutting iterations/recursions to perform. A
	 *                   value of 1 will simply cut the corners once, higher values
	 *                   will effectively smooth the cut. Values greater than ~10
	 *                   generally have no additional visual effect.
	 * @return A copy of the input shape with corners cut.
	 * @since 1.1.0
	 */
	public static PShape chaikinCut(PShape shape, double ratio, int iterations) {
		ratio = Math.max(ratio, 1e-6);
		ratio = Math.min(ratio, 1 - 1e-6);
		ratio /= 2; // constrain to 0...0.5
		float r = (float) ratio;

		return PGS.applyToLinealGeometries(shape, ring -> {
			var cut = ChaikinCut.chaikin(toPShape(ring), r, iterations);
			var g = fromPShape(cut);
			if (g instanceof Polygon) {
				g = ((Polygon) g).getExteriorRing();
			}
			if (g instanceof Lineal) {
				return (LineString) g;
			}
			return null; // pointal or other...
		});
	}

	/**
	 * Radially warps a polygon by moving each boundary vertex inward/outward along
	 * the ray from the polygon centroid to that vertex, creating a warping or
	 * perturbing effect.
	 * <p>
	 * Optionally, the input boundary can be densified before warping by inserting
	 * additional vertices at a spacing of ~1 unit. This causes long edges to warp
	 * smoothly along their full length rather than only at the original corner
	 * vertices.
	 * 
	 * @param shape      A polygonal {@link PShape} (or GROUP of polygons) to be
	 *                   distorted. The warp is applied to each polygon ring
	 *                   independently.
	 * @param magnitude  Controls the strength of the warp. Larger values produce
	 *                   larger radial displacements from the original boundary
	 *                   (i.e., larger inward/outward movement). A value of
	 *                   {@code 0} produces an unchanged shape.
	 * @param warpOffset An angular phase offset (in radians) added to each vertex's
	 *                   polar angle before sampling the noise field. Changing
	 *                   {@code warpOffset} does not change the warp magnitude; it
	 *                   rotates the noise pattern around the centroid (i.e., shifts
	 *                   where bulges/indentations occur along the boundary). This
	 *                   is useful for animation by incrementing {@code warpOffset}
	 *                   over time. The warp has a period of 2π. A typical/useful
	 *                   domain is {@code [0, 2*Math.PI)}.
	 * @param densify    A boolean parameter determining whether the shape should be
	 *                   densified (by inserting additional vertices at a distance
	 *                   of 1) before warping. If true, shapes with long edges will
	 *                   experience warping along their entire length, not just at
	 *                   the original vertices.
	 * @return A new PShape object that has been radially warped according to the
	 *         specified parameters.
	 */
	public static PShape radialWarp(PShape shape, double magnitude, double warpOffset, boolean densify) {
		final UniformNoise noise = new UniformNoise(1337);

		return PGS.applyToLinealGeometries(shape, line -> {

			// radialWarp is defined for polygon rings; if we get an open line, just return
			// it unchanged
			if (!line.isClosed()) {
				return line;
			}
			final Point centroid = line.getCentroid();
			final PVector c = new PVector((float) centroid.getX(), (float) centroid.getY());

			Geometry working = line;
			if (densify) {
				final Densifier d = new Densifier(line);
				d.setDistanceTolerance(1);
				d.setValidate(false);
				working = d.getResultGeometry();
			}

			final Coordinate[] coords = working.getCoordinates();
			if (coords.length == 0) {
				return line;
			}

			// Warp all unique vertices; then explicitly re-close
			final int n = coords.length;
			for (int i = 0; i < n - 1; i++) { // ignore last coordinate (closure); we re-close after warping
				final double x = coords[i].x;
				final double y = coords[i].y;

				double dx = x - c.x;
				double dy = y - c.y;

				final double len = Math.sqrt(dx * dx + dy * dy);
				if (len == 0) {
					continue; // vertex at centroid
				}

				final double angle = FastAtan2.atan2(dy, dx) + warpOffset;

				float perturbation = noise.uniformNoise(FastMath.cos(angle), FastMath.sin(angle));
				perturbation -= 0.5f; // [0..1] -> [-0.5..0.5]
				perturbation *= (float) (magnitude * 2.0);

				// normalize heading and displace
				dx /= len;
				dy /= len;

				coords[i].x = x + dx * perturbation;
				coords[i].y = y + dy * perturbation;
			}

			// ensure exact closure
			coords[n - 1].x = coords[0].x;
			coords[n - 1].y = coords[0].y;

			// preserve ring-ness if possible
			if (line instanceof LinearRing) {
				return PGS.GEOM_FACTORY.createLinearRing(coords);
			}
			return PGS.GEOM_FACTORY.createLineString(coords);
		});
	}

	/**
	 * Warps/perturbs a shape by displacing vertices (both positively and
	 * negatively) according to the magnitude of a sine wave which follows the shape
	 * perimeter at some frequency.
	 * 
	 * @param shape     single polygonal shape
	 * @param magnitude maximum perpendicular displacement along the shape perimeter
	 * @param frequency sine wave frequency. Values less than 1 will result in an
	 *                  offset that does not smoothly join up.
	 * @param phase     sine wave phase. corresponds to the fraction (0...1) around
	 *                  the shape perimeter where the wave starts (0 displacement).
	 * @return warped polygonal shape
	 * @since 1.3.0
	 */
	public static PShape sineWarp(PShape shape, double magnitude, double frequency, double phase) {
		Geometry g = fromPShape(shape);
		if (g instanceof Polygonal) {
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)) {
				g = g.getGeometryN(0);
			}
			g = ((Polygon) g).getExteriorRing();
		}
		final LengthIndexedLine l = new LengthIndexedLine(g);
		final double length = l.getEndIndex(); // perimeter length
		final CoordinateList coords = new CoordinateList();

		for (double distance = 0; distance < length; distance++) {
			final Coordinate coord = l.extractPoint(distance, Math.sin(Math.PI * 2 * frequency * distance / length + (Math.PI * 2 * phase)) * magnitude);
			coords.add(coord);
		}
		coords.closeRing();

		Geometry out = GeometryFixer.fix(PGS.GEOM_FACTORY.createPolygon(coords.toCoordinateArray()));
		return PGS_Conversion.toPShape(out);
	}

	/**
	 * Warps/perturbs a shape by displacing vertices according to a 2D noise vector
	 * field.
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
	 * @return warped polygonal shape
	 * @see #fieldWarp(PShape, double, double, double, boolean, int)
	 */
	public static PShape fieldWarp(PShape shape, double magnitude, double noiseScale, boolean densify) {
		return fieldWarp(shape, magnitude, noiseScale, 0, densify, 1337);
	}

	/**
	 * Warps/perturbs a shape by displacing vertices according to a 2D noise vector
	 * field.
	 * <p>
	 * Inputs may be densified before warping for more finely-grained warping.
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
	 * @return warped polygonal shape
	 */
	public static PShape fieldWarp(PShape shape, double magnitude, double noiseScale, double time, boolean densify, long noiseSeed) {
		float scale = Math.max(1, (float) noiseScale * 500f);
		final boolean pointsShape = shape.getKind() == PConstants.POINTS;

		final PShape copy;
		if (densify && !pointsShape) {
			final Densifier d = new Densifier(fromPShape(shape));
			d.setDistanceTolerance(PGS_Conversion.BEZIER_SAMPLE_DISTANCE);
			d.setValidate(false);
			copy = toPShape(d.getResultGeometry());
		} else {
			copy = PGS_Conversion.copy(shape);
		}

		final UniformNoise noise = new UniformNoise((int) (noiseSeed % Integer.MAX_VALUE));

		if (copy.getChildCount() == 0) {
			// setVertex() will act on group shapes, so treat a single shape as group of 1
			copy.addChild(copy);
		}

		for (PShape child : copy.getChildren()) {
			int vCount = child.getVertexCount();
			if (vCount == 0)
				continue;

			// Determine if the shape is closed.
			boolean isClosed = child.isClosed() || (vCount > 1 && child.getVertex(0).equals(child.getVertex(vCount - 1)));

			// If closed, we iterate up to N-1 and handle the last vertex separately to
			// ensure closure.
			int limit = isClosed ? vCount - 1 : vCount;

			for (int i = 0; i < limit; i++) {
				final PVector coord = child.getVertex(i);
				float dx = noise.uniformNoise(coord.x / scale, coord.y / scale + time) - 0.5f;
				float dy = noise.uniformNoise(coord.x / scale + (101 + time), coord.y / scale + (101 + time)) - 0.5f;
				child.setVertex(i, coord.x + (dx * (float) magnitude * 2), coord.y + (dy * (float) magnitude * 2));
			}

			// If the shape was closed, sync the last vertex with the newly warped first
			// vertex.
			if (isClosed && vCount > 1) {
				PVector firstV = child.getVertex(0);
				child.setVertex(vCount - 1, firstV.x, firstV.y);
			}
		}

		if (pointsShape) {
			return copy;
		} else {
			if (copy.getChildCount() == 1) {
				// Fix self-intersections or invalid geometries caused by warping
				return toPShape(GeometryFixer.fix(fromPShape(copy.getChild(0))));
			} else {
				// Return group as-is (fixing individual children would be safer but requires a
				// loop)
				return copy;
			}
		}
	}

	/**
	 * Applies a pinch warping effect to a shape, distorting vertices towards a
	 * specified point. The strength of the pinch effect decreases with distance
	 * from the pinch point, creating a localized warping effect.
	 *
	 * @param shape      The shape. The original shape is not modified; instead, a
	 *                   new warped shape is returned.
	 * @param pinchPoint The point in 2D space towards which the shape will be
	 *                   pinched. Vertices closer to this point will be affected
	 *                   more strongly.
	 * @param weight     The strength of the pinching effect. Higher values create
	 *                   stronger distortion.
	 * @return A new PShape object comprising the warped geometry.
	 * @since 2.0
	 */
	public static PShape pinchWarp(PShape shape, PVector pinchPoint, double weight) {
		return PGS.applyToLinealGeometries(shape, line -> {
			final var gf = line.getFactory();
			final var coords = line.getCoordinates();

			if (coords.length == 0) {
				return line;
			}

			final boolean closed = line.isClosed();

			for (int i = 0; i < coords.length; i++) {
				// if closed, we'll re-close explicitly after warping to avoid drift
				if (closed && i == coords.length - 1) {
					break;
				}

				final double x = coords[i].x;
				final double y = coords[i].y;

				final double dx = pinchPoint.x - x;
				final double dy = pinchPoint.y - y;

				final double distance = Math.sqrt(dx * dx + dy * dy);
				final double w = weight / (distance + 1.0);

				coords[i].x = x + dx * w;
				coords[i].y = y + dy * w;
			}

			if (closed) {
				coords[coords.length - 1].x = coords[0].x;
				coords[coords.length - 1].y = coords[0].y;
			}

			return gf.createLineString(coords);
		});
	}

	/**
	 * Generates an intermediate shape between two shapes by interpolating between
	 * their exterior rings. This process has many names: shape morphing / blending
	 * / averaging / tweening / interpolation.
	 * <p>
	 * Note the interpolated shape may self-intersect (this implementation is not
	 * "rigid").
	 * 
	 * @param from                a single polygon; the shape we want to morph from
	 * @param to                  a single polygon; the shape we want to morph
	 *                            <code>from</code> into
	 * @param interpolationFactor between 0...1
	 * @return a polygonal PShape
	 * @since 1.2.0
	 * @see #interpolate(PShape, PShape, int)
	 * @implNote Uses {@link NewtonThieleRingMorpher} for higher-quality
	 *           interpolation.
	 */
	public static PShape interpolate(PShape from, PShape to, double interpolationFactor) {
		return interpolate(List.of(from, to), interpolationFactor);
	}

	/**
	 * Generates an intermediate shape from a sequence of input shapes by
	 * interpolating (morphing) between their exterior rings.
	 * <p>
	 * This is a generalisation of {@link #interpolate(PShape, PShape, double)} to
	 * more than two shapes. The interpolation follows the order of {@code shapes}.
	 * <p>
	 * Note the interpolated shape may self-intersect (this implementation is not
	 * "rigid").
	 *
	 * @param shapes              a list of single-polygon {@link PShape}s; only the
	 *                            exterior ring is used.
	 * @param interpolationFactor interpolation parameter in the range
	 *                            {@code [0..1]}
	 * @return a polygonal {@link PShape} representing the interpolated shape
	 * @since 2.2.0
	 * @see #interpolate(PShape, PShape, double)
	 * @implNote Uses {@link NewtonThieleRingMorpher} for higher-quality
	 *           interpolation.
	 */
	public static PShape interpolate(List<PShape> shapes, double interpolationFactor) {
		var rings = shapes.stream().map(s -> ((Polygon) fromPShape(s)).getExteriorRing()).toArray(LinearRing[]::new);
		NewtonThieleRingMorpher m = new NewtonThieleRingMorpher(rings);
		var tween = m.interpolate(interpolationFactor);
		return toPShape(tween);
	}

	/**
	 * Generates intermediate shapes (frames) by interpolating (morphing) through a
	 * sequence of shapes. This process has many names: shape morphing / blending /
	 * averaging / tweening / interpolation.
	 * <p>
	 * The returned frames include both endpoints: the first frame corresponds to
	 * {@code t = 0} (the first shape in {@code shapes}) and the last frame
	 * corresponds to {@code t = 1} (the last shape in {@code shapes}). Intermediate
	 * frames are evenly spaced in {@code [0..1]} using {@code t = i/(frames-1)}.
	 * <p>
	 * This method is faster than calling {@link #interpolate(List, double)} (or
	 * {@link #interpolate(PShape, PShape, double)}) repeatedly for different
	 * interpolation factors.
	 *
	 * @param shapes a list of single-polygon {@link PShape}s, in the order they
	 *               should be morphed through; only the exterior ring is used.
	 * @param frames the number of frames (including first and last) to generate;
	 *               must be {@code >= 2}
	 * @return a GROUP {@link PShape} whose children are the generated frames
	 * @since 2.2.0
	 * @see #interpolate(List, double)
	 * @see #interpolate(PShape, PShape, double)
	 */
	public static PShape interpolate(List<PShape> shapes, int frames) {
		var rings = shapes.stream().map(s -> ((Polygon) fromPShape(s)).getExteriorRing()).toArray(LinearRing[]::new);
		NewtonThieleRingMorpher m = new NewtonThieleRingMorpher(rings);

		final double fraction = 1d / (frames - 1);
		PShape out = new PShape();
		for (int i = 0; i < frames; i++) {
			out.addChild(toPShape(m.interpolate(fraction * i)));
		}

		return out;
	}

	/**
	 * As-rigid-as-possible (ARAP) 2D deformation of a polygon {@link PShape} using
	 * point handles.
	 * <h2>Handle semantics</h2>
	 * <ul>
	 * <li>{@code handles} are points in the <em>rest</em> (original) shape's
	 * coordinate space.</li>
	 * <li>{@code handleTargets} are the desired positions for those same handles in
	 * the <em>deformed</em> shape.</li>
	 * <li>Both lists must have the same size and matching order (i.e., index
	 * {@code i} in {@code handles} maps to index {@code i} in
	 * {@code handleTargets}).</li>
	 * <li>ARAP typically requires at least 2 handles for a stable solve.</li>
	 * </ul>
	 *
	 * <h2>Performance notes</h2>
	 * <p>
	 * This method rebuilds and refines a triangulation on every call. For
	 * interactive dragging (re-solving every frame), prefer using {@link Malleo}
	 * directly: build the triangulation and call
	 * {@link Malleo#prepareHandles(List)} once, then repeatedly call
	 * {@link Malleo#solve(Malleo.CompiledHandles, List)} with updated targets.
	 *
	 * <h2>Output</h2>
	 * <p>
	 * Returns the deformed polygon boundary. The result may self-intersect
	 * depending on handle motion and mesh quality.
	 *
	 * @param shape         the rest shape to deform (expected to be a single
	 *                      polygon {@code PShape})
	 * @param handles       handle locations in rest-space
	 * @param handleTargets target locations for each handle, in the same order as
	 *                      {@code handles}
	 * @return a new {@code PShape} representing the deformed shape
	 * @since 2.2
	 */
	public static PShape arapDeform(PShape shape, List<PVector> handles, List<PVector> handleTargets) {
		var t = PGS_Triangulation.delaunayTriangulationMesh(shape);
		PGS_Triangulation.refine(t, 15, 50); // refine
		var g = PGS_Triangulation.toGeometry(t);

		Malleo m = new Malleo(g);
		var mHandles = Arrays.asList(PGS.toCoords(handles));
		var mTargets = Arrays.asList(PGS.toCoords(handleTargets));

		var compiledHandles = m.prepareHandles(mHandles);

		var deformed = m.solve(compiledHandles, mTargets);

		return toPShape(deformed);
	}

	/**
	 * Reduces the precision of a shape, whilst ensuring the output shape is valid.
	 * <p>
	 * This method effectively rounds vertices to the nearest value given by
	 * <code>precision</code>.
	 * 
	 * @param shape     shape to reduce
	 * @param precision the exact grid size with which to round shape vertices.
	 *                  should be non-zero and positive
	 * @return reduced copy of input
	 * @since 1.3.0
	 */
	public static PShape reducePrecision(PShape shape, double precision) {
		var pm = new PrecisionModel(-Math.max(Math.abs(precision), 1e-10));
		if (shape.getFamily() == PShape.GROUP) {
			// pointwise preserves polygon faces (doesn't merge)
			return toPShape(GeometryPrecisionReducer.reducePointwise(fromPShape(shape), pm));
		} else {
			return toPShape(GeometryPrecisionReducer.reduce(fromPShape(shape), pm));
		}
	}

	/**
	 * Regularises (straightens) the contour of a lineal {@link PShape} by snapping
	 * edges toward a small set of principal directions and simplifying the result.
	 * The prinicipal direction is derived from the shape's longest edge.
	 *
	 * @param shape     a lineal {@code PShape} to regularise (or a group containing
	 *                  lineal children)
	 * @param maxOffset maximum allowed offset. Used to constrain how far the
	 *                  regularised contour may deviate from the input; must be
	 *                  &gt;= 0
	 * @return a new {@code PShape} whose linework has been regularised
	 * @see #regularise(PShape, double, double)
	 * @since 2.2
	 */
	public static PShape regularise(PShape shape, double maxOffset) {
		var params = Parameters.builder().maximumOffset(maxOffset);
		return PGS.applyToLinealGeometries(shape, l -> {
			return ContourRegularization.regularize(l, params.build());
		});
	}

	/**
	 * Regularises (straightens) the contour of a lineal {@link PShape} by snapping
	 * edges toward principal directions and simplifying the result.
	 * <p>
	 * This overload lets you provide an explicit <em>principal axis
	 * orientation</em> (in degrees). Edges are snapped to be parallel to that axis
	 * or to its orthogonal (axis + 90°), subject to the {@code maxOffset}
	 * constraint.
	 *
	 * @param shape           a lineal {@code PShape} to regularize (or a group
	 *                        containing lineal children)
	 * @param maxOffset       maximum allowed offset used to constrain how far the
	 *                        regularised contour may deviate from the input; must
	 *                        be &gt;= 0
	 * @param axisOrientation principal axis direction, in degrees, expected in the
	 *                        range {@code [0,180)} (values outside this range are
	 *                        normalised)
	 * @return a new {@code PShape} whose linework has been regularised
	 * @see #regularise(PShape, double)
	 * @since 2.2
	 */
	public static PShape regularise(PShape shape, double maxOffset, double axisOrientation) {
		var d = new ContourRegularization.UserDefinedDirections(5, axisOrientation);
		var params = Parameters.builder().maximumOffset(maxOffset).directions(d);
		return PGS.applyToLinealGeometries(shape, l -> {
			return ContourRegularization.regularize(l, params.build());
		});
	}

	/**
	 * The end cap style to use. Cap style specifies the shape of the ends of
	 * buffered unclosed lines; it has no effect in polygons.
	 */
	public enum CapStyle {

		/**
		 * The usual round end caps.
		 */
		ROUND(BufferParameters.CAP_ROUND),
		/**
		 * End caps are truncated flat at the line ends.
		 */
		FLAT(BufferParameters.CAP_FLAT),
		/**
		 * End caps are squared off at the buffer distance beyond the line ends.
		 */
		SQUARE(BufferParameters.CAP_SQUARE);

		final int style;

		private CapStyle(int style) {
			this.style = style;
		}
	}

	private static BufferParameters createBufferParams(double r, double delta, OffsetStyle bufferStyle, CapStyle capStyle) {
		r = Math.abs(r);

		// compute the number of points for the full circle
		double ang = Math.acos(1.0 - delta / r);
		// if delta/r > 2 or so acos will fail – clamp it
		if (Double.isNaN(ang) || ang <= 0) {
			// in this degenerate case just fall back to a small number
			ang = Math.PI / 8.0;
		}

		// total points
		double nPtsDbl = Math.PI / ang;
		int nPts = (int) Math.ceil(nPtsDbl);

		// segments per quadrant
		int quadSeg = (int) Math.ceil(nPts / 4.0);

		// enforce a sensible minimum
		quadSeg = Math.max(quadSeg, BufferParameters.DEFAULT_QUADRANT_SEGMENTS);

		// cap style affects linestrings only
		return new BufferParameters(quadSeg, capStyle.style, bufferStyle.style, BufferParameters.DEFAULT_MITRE_LIMIT);
	}

}
