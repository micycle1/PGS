package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static processing.core.PConstants.GROUP;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;

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

import micycle.hobbycurves.HobbyCurve;
import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.PGS_Contour.OffsetStyle;
import micycle.pgs.color.Colors;
import micycle.pgs.commons.ChaikinCut;
import micycle.pgs.commons.CornerRounding;
import micycle.pgs.commons.DiscreteCurveEvolution;
import micycle.pgs.commons.EllipticFourierDesc;
import micycle.pgs.commons.GaussianLineSmoothing;
import micycle.pgs.commons.ShapeInterpolation;
import micycle.pgs.commons.DiscreteCurveEvolution.DCETerminationCallback;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.minkowski_sum.MinkowskiSum;
import micycle.uniformnoise.UniformNoise;

/**
 * Methods that affect the geometry or topology of shapes.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Morphology {

	static {
		MinkowskiSum.setGeometryFactory(PGS.GEOM_FACTORY);
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
		final int segments = (int) Math.ceil(BufferParameters.DEFAULT_QUADRANT_SEGMENTS + Math.sqrt(buffer));
		return toPShape(fromPShape(shape).buffer(buffer, segments));
	}

	/**
	 * Computes a buffer area around the shape, having the given buffer width and
	 * buffer style (either round, miter, bevel).
	 * 
	 * @param shape
	 * @param buffer extent/width of the buffer (which may be positive or negative)
	 * @return a polygonal shape representing the buffer region (which may be empty)
	 * @see #buffer(PShape, double)
	 * @since 1.3.0
	 */
	public static PShape buffer(PShape shape, double buffer, OffsetStyle bufferStyle) {
		Geometry g = fromPShape(shape);
		final int segments = (int) Math.ceil(BufferParameters.DEFAULT_QUADRANT_SEGMENTS + Math.sqrt(buffer));
		BufferParameters bufParams = new BufferParameters(segments, BufferParameters.CAP_FLAT, bufferStyle.style,
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
	 *         be empty)
	 * @since 1.3.0
	 */
	public static PShape variableBuffer(PShape shape, double startDistance, double endDistance) {
		Geometry g = fromPShape(shape);
		if (!g.getGeometryType().equals(Geometry.TYPENAME_LINEARRING) && !g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			g = ((Polygon) g).getExteriorRing(); // variable buffer applies to linestrings only
		}
		return toPShape(VariableBuffer.buffer(g, startDistance, endDistance));
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
	 * @param shape          A single polygon or lineal shape
	 * @param bufferCallback A callback function that receives the vertex coordinate
	 *                       and a double representing tractional distance (0...1)
	 *                       of the vertex along the shape's boundary. The function
	 *                       may use properties of the vertex, or its position, to
	 *                       determine the buffer width at that point.
	 * @return A new shape representing the original shape buffered with variable
	 *         widths as specified by the callback function. The width at each
	 *         vertex is calculated independently.
	 * @since 2.0
	 */
	public static PShape variableBuffer(PShape shape, BiFunction<Coordinate, Double, Double> bufferCallback) {
		final Geometry inputGeometry = fromPShape(shape);
		if (!(inputGeometry instanceof Lineal || inputGeometry instanceof Polygon)) {
			throw new IllegalArgumentException("The geometry must be linear or a non-multi polygonal shape.");
		}
		var coords = inputGeometry.getCoordinates();
		double[] bufferDistances = new double[coords.length];
		double totalLength = inputGeometry.getLength();
		double running_length = 0;
		Coordinate previousCoordinate = coords[0];

		for (int i = 1; i < coords.length; i++) {
			running_length += previousCoordinate.distance(coords[i]);
			double fractionalDistance = running_length / totalLength; // 0...1
			bufferDistances[i] = bufferCallback.apply(coords[i], fractionalDistance);
			previousCoordinate = coords[i];
		}

		bufferDistances[0] = bufferCallback.apply(coords[0], 0.0);

		VariableBuffer variableBuffer = new VariableBuffer(inputGeometry, bufferDistances);
		return toPShape(variableBuffer.getResult());
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
	 * @return
	 * @see #dilationErosion(PShape, double)
	 */
	public static PShape erosionDilation(PShape shape, double buffer) {
		buffer = Math.abs(buffer);

		final int segments = (int) Math.ceil(BufferParameters.DEFAULT_QUADRANT_SEGMENTS + Math.sqrt(buffer));
		Geometry g = BufferOp.bufferOp(fromPShape(shape), -buffer, segments);
		g = BufferOp.bufferOp(g, +buffer, segments);

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
	 * @return
	 * @since 1.3.0
	 * @see #erosionDilation(PShape, double)
	 */
	public static PShape dilationErosion(PShape shape, double buffer) {
		buffer = Math.abs(buffer);

		final int segments = (int) Math.ceil(BufferParameters.DEFAULT_QUADRANT_SEGMENTS + Math.sqrt(buffer));
		Geometry g = BufferOp.bufferOp(fromPShape(shape), buffer, segments);
		g = BufferOp.bufferOp(g, -buffer, segments);
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
	 * shape, starting with those having the least shape-relevance.
	 * <p>
	 * The simplification process terminates according to a user-specified
	 * {@link DCETerminationCallback#shouldTerminate(Coordinate, double, int)
	 * callback} that decides whether the DCE algorithm should terminate based on
	 * the current kink (having a candidate vertex), using its coordinates,
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
	 * @see PGS_Morphology#simplifyDCE(PShape, int)
	 */
	public static PShape simplifyDCE(PShape shape, DCETerminationCallback terminationCallback) {
		Geometry g = fromPShape(shape);
		switch (g.getGeometryType()) {
			case Geometry.TYPENAME_GEOMETRYCOLLECTION :
			case Geometry.TYPENAME_MULTIPOLYGON :
			case Geometry.TYPENAME_MULTILINESTRING :
				PShape group = new PShape(GROUP);
				for (int i = 0; i < g.getNumGeometries(); i++) {
					group.addChild(simplifyDCE(toPShape(g.getGeometryN(i)), terminationCallback));
				}
				return group;
			case Geometry.TYPENAME_LINEARRING :
			case Geometry.TYPENAME_POLYGON :
				// process each ring individually
				LinearRing[] rings = new LinearRingIterator(g).getLinearRings();
				LinearRing[] dceRings = new LinearRing[rings.length];
				for (int i = 0; i < rings.length; i++) {
					LinearRing ring = rings[i];
					Coordinate[] dce = DiscreteCurveEvolution.process(ring, terminationCallback);
					dceRings[i] = PGS.GEOM_FACTORY.createLinearRing(dce);
				}
				LinearRing[] holes = null;
				if (dceRings.length > 1) {
					holes = Arrays.copyOfRange(dceRings, 1, dceRings.length);
				}
				return toPShape(PGS.GEOM_FACTORY.createPolygon(dceRings[0], holes));
			case Geometry.TYPENAME_LINESTRING :
				LineString l = (LineString) g;
				Coordinate[] dce = DiscreteCurveEvolution.process(l, terminationCallback);
				return toPShape(PGS.GEOM_FACTORY.createLineString(dce));
			default :
				System.err.println(g.getGeometryType() + " are not supported for the simplifyDCE() method."); // pointal geoms
				return new PShape(); // return empty (so element is invisible if not processed)
		}
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
		tension = Math.max(tension, 0.668); // prevent degeneracy
		double[][] vertices = PGS_Conversion.toArray(shape, false);
		HobbyCurve curve = new HobbyCurve(vertices, tension, shape.isClosed(), 0.5, 0.5);
		List<PVector> points = new ArrayList<>();
		for (double[] b : curve.getBeziers()) {
			int i = 0;
			PVector p1 = new PVector((float) b[i++], (float) b[i++]);
			PVector cp1 = new PVector((float) b[i++], (float) b[i++]);
			PVector cp2 = new PVector((float) b[i++], (float) b[i++]);
			PVector p2 = new PVector((float) b[i++], (float) b[i]);
			PShape bezier = PGS_Conversion.fromCubicBezier(p1, cp1, cp2, p2);
			points.addAll(PGS_Conversion.toPVector(bezier));
		}

		return PGS_Conversion.fromPVector(points);
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

		switch (g.getGeometryType()) {
			case Geometry.TYPENAME_GEOMETRYCOLLECTION :
			case Geometry.TYPENAME_MULTIPOLYGON :
			case Geometry.TYPENAME_MULTILINESTRING :
				PShape group = new PShape(GROUP);
				for (int i = 0; i < g.getNumGeometries(); i++) {
					group.addChild(smoothGaussian(toPShape(g.getGeometryN(i)), sigma));
				}
				return group;
			case Geometry.TYPENAME_POLYGON :
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
			case Geometry.TYPENAME_LINEARRING :
			case Geometry.TYPENAME_LINESTRING :
				LineString l = (LineString) g;
				return toPShape(GaussianLineSmoothing.get(l, Math.max(sigma, 1), 1));
			default :
				System.err.println(g.getGeometryType() + " are not supported for the smoothGaussian() method."); // pointal geoms
				return new PShape(); // return empty (so element is invisible if not processed)
		}
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
		Geometry g = fromPShape(shape);
		switch (g.getGeometryType()) {
			case Geometry.TYPENAME_GEOMETRYCOLLECTION :
			case Geometry.TYPENAME_MULTIPOLYGON :
			case Geometry.TYPENAME_MULTILINESTRING :
				PShape group = new PShape(GROUP);
				for (int i = 0; i < g.getNumGeometries(); i++) {
					group.addChild(smoothEllipticFourier(toPShape(g.getGeometryN(i)), descriptors));
				}
				return group;
			case Geometry.TYPENAME_POLYGON :
				LinearRingIterator lri = new LinearRingIterator(g);
				LinearRing[] rings = lri.getLinearRings();
				LinearRing[] ringProcessed = new LinearRing[rings.length];
				for (int i = 0; i < rings.length; i++) {
					descriptors = Math.min(rings[i].getCoordinates().length / 2, descriptors); // max=#vertices/2
					descriptors = Math.max(2, descriptors); // min=2
					final EllipticFourierDesc efd = new EllipticFourierDesc(rings[i], descriptors);
					Coordinate[] coords = efd.createPolygon();
					ringProcessed[i] = PGS.GEOM_FACTORY.createLinearRing(coords);
				}

				LinearRing[] holes = null;
				if (ringProcessed.length > 1) {
					holes = Arrays.copyOfRange(ringProcessed, 1, ringProcessed.length);
				}
				return toPShape(PGS.GEOM_FACTORY.createPolygon(ringProcessed[0], holes));
			case Geometry.TYPENAME_LINEARRING :
				descriptors = Math.min(shape.getVertexCount() / 2, descriptors); // max=#vertices/2
				descriptors = Math.max(2, descriptors); // min=2
				LinearRing l = (LinearRing) g;
				final EllipticFourierDesc efd = new EllipticFourierDesc(l, descriptors);
				return toPShape(PGS.GEOM_FACTORY.createLinearRing(efd.createPolygon()));
			case Geometry.TYPENAME_LINESTRING :
			default :
				System.err.println(g.getGeometryType() + " are not supported for the smoothEllipticFourier() method."); // pointal/string
																														// geoms
				return new PShape(); // return empty (so element is invisible if not processed)
		}
	}

	/**
	 * Modifies the corners of a specified shape by replacing each angular corner
	 * with a smooth, circular arc. The radius of each arc is determined
	 * proportionally to the shorter of the two lines forming the corner.
	 * 
	 * @param shape  The original PShape object whose corners are to be rounded.
	 * @param extent Specifies the degree of corner rounding, with a range from 0 to
	 *               1. A value of 0 corresponds to no rounding, whereas a value of
	 *               1 yields maximum rounding while still maintaining the validity
	 *               of the shape. Values above 1 are accepted but may produce
	 *               unpredictable results.
	 * @return A new PShape object with corners rounded to the specified extent.
	 */
	public static PShape round(PShape shape, double extent) {
		return CornerRounding.round(shape, extent);
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
		PShape cut = ChaikinCut.chaikin(shape, (float) ratio, iterations);
		PGS_Conversion.setAllFillColor(cut, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(cut, Colors.PINK, 3);
		return cut;
	}

	/**
	 * Distorts a polygonal shape by radially displacing its vertices along the line
	 * connecting each vertex with the shape's centroid, creating a warping or
	 * perturbing effect.
	 * <p>
	 * The shape's input vertices can optionally be densified prior to the warping
	 * operation.
	 * 
	 * @param shape      A polygonal PShape object to be distorted.
	 * @param magnitude  The degree of the displacement, which determines the
	 *                   maximum Euclidean distance a vertex will be moved in
	 *                   relation to the shape's centroid.
	 * @param warpOffset An offset angle, which establishes the starting angle for
	 *                   the displacement process.
	 * @param densify    A boolean parameter determining whether the shape should be
	 *                   densified (by inserting additional vertices at a distance
	 *                   of 1) before warping. If true, shapes with long edges will
	 *                   experience warping along their entire length, not just at
	 *                   the original vertices.
	 * @return A new PShape object that has been radially warped according to the
	 *         specified parameters.
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
		if (!coords.get(0).equals(coords.get(coords.size() - 1))) {
			coords.add(coords.get(0));
		}
		return PGS_Conversion.fromPVector(coords);
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
	 * @return
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
			final Coordinate coord = l.extractPoint(distance,
					Math.sin(Math.PI * 2 * frequency * distance / length + (Math.PI * 2 * phase)) * magnitude);
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
	 * @return
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
	 * @return
	 */
	public static PShape fieldWarp(PShape shape, double magnitude, double noiseScale, double time, boolean densify, long noiseSeed) {
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

		final UniformNoise noise = new UniformNoise((int) (noiseSeed % Integer.MAX_VALUE));

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
		final Geometry fromGeom = fromPShape(from);
		final Geometry toGeom = fromPShape(to);
		if (toGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON) && fromGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			final ShapeInterpolation tween = new ShapeInterpolation(fromGeom, toGeom);
			return toPShape(PGS.GEOM_FACTORY.createPolygon(tween.tween(interpolationFactor)));
		} else {
			System.err.println("interpolate() accepts holeless single polygons only (for now).");
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
	 * @since 1.3.0
	 * @see #interpolate(PShape, PShape, double)
	 */
	public static PShape interpolate(PShape from, PShape to, int frames) {
		final Geometry fromGeom = fromPShape(from);
		final Geometry toGeom = fromPShape(to);
		if (toGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON) && fromGeom.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			final ShapeInterpolation tween = new ShapeInterpolation(fromGeom, toGeom);
			final float fraction = 1f / (frames - 1);
			PShape out = new PShape();
			for (int i = 0; i < frames; i++) {
				out.addChild(toPShape(PGS.GEOM_FACTORY.createPolygon(tween.tween(fraction * i))));
			}
			return out;
		} else {
			System.err.println("interpolate() accepts holeless single polygons only (for now).");
			return from;
		}
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
		return toPShape(GeometryPrecisionReducer.reduce(fromPShape(shape), new PrecisionModel(-Math.max(Math.abs(precision), 1e-10))));
	}

}
