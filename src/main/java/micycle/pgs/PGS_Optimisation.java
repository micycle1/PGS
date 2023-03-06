package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static processing.core.PConstants.GROUP;
import static micycle.pgs.PGS_Construction.createEllipse;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.algorithm.construct.LargestEmptyCircle;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.util.GeometricShapeFactory;

import almadina.rectpacking.RBPSolution;
import almadina.rectpacking.Rect;
import almadina.rectpacking.RectPacking.PackingHeuristic;
import micycle.pgs.color.RGB;
import micycle.pgs.commons.ClosestPointPair;
import micycle.pgs.commons.FarthestPointPair;
import micycle.pgs.commons.MaximumInscribedAARectangle;
import micycle.pgs.commons.MaximumInscribedRectangle;
import micycle.pgs.commons.MinimumBoundingEllipse;
import micycle.pgs.commons.MinimumBoundingTriangle;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Solve geometric optimisation problems, such as bounding volumes, inscribed
 * areas, optimal distances, etc.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Optimisation {

	private PGS_Optimisation() {
	}

	/**
	 * Computes the shape's envelope (bounding box). The vertices of the output
	 * PShape begin at the top-left corner of the envelope and are arranged
	 * counter-clockwise.
	 * 
	 * @param shape a rectangular shape that covers/bounds the input
	 * @return
	 */
	public static PShape envelope(PShape shape) {
		return toPShape(fromPShape(shape).getEnvelope());
	}

	/**
	 * The Maximum Inscribed Circle is determined by a point in the interior of the
	 * area which has the farthest distance from the area boundary, along with a
	 * boundary point at that distance.
	 * 
	 * @param shape
	 * @param tolerance the distance tolerance for computing the centre point
	 *                  (around 1)
	 */
	public static PShape maximumInscribedCircle(PShape shape, double tolerance) {
		MaximumInscribedCircle mic = new MaximumInscribedCircle(fromPShape(shape), tolerance);
		final double wh = mic.getRadiusLine().getLength() * 2;
		Polygon circle = createEllipse(PGS.coordFromPoint(mic.getCenter()), wh, wh);
		return toPShape(circle);
	}

	/**
	 * Return the maximum circle (at a given centerpoint inside/outside the circle)
	 * 
	 * @param shape
	 * @param centerPoint
	 * @return A circular PShape
	 */
	public static PShape maximumInscribedCircle(PShape shape, PVector centerPoint) {
		Geometry g = fromPShape(shape);
		Point p = PGS.pointFromPVector(centerPoint);
		Coordinate closestEdgePoint = DistanceOp.nearestPoints(g.getBoundary(), p)[0];
		double radius = PGS.coordFromPVector(centerPoint).distance(closestEdgePoint);
		Polygon circle = createEllipse(p.getCoordinate(), radius * 2, radius * 2);
		return toPShape(circle);
	}

	/**
	 * Finds an approximate largest area rectangle (of arbitrary orientation)
	 * contained within a polygon.
	 * 
	 * @param shape
	 * @return a rectangle shape
	 * @see #maximumInscribedAARectangle(PShape, boolean)
	 *      maximumInscribedAARectangle() - the largest axis-aligned rectangle
	 */
	public static PShape maximumInscribedRectangle(PShape shape) {
		Polygon polygon = (Polygon) fromPShape(shape);
		MaximumInscribedRectangle mir = new MaximumInscribedRectangle(polygon);
		return toPShape(mir.computeMIR());
	}

	/**
	 * Finds the rectangle with a maximum area whose sides are parallel to the
	 * x-axis and y-axis ("axis-aligned"), contained within a convex shape.
	 * <p>
	 * This method computes the MIR for convex shapes only; if a concave shape is
	 * passed in, the resulting rectangle will be computed based on its convex hull.
	 * <p>
	 * This method uses a brute force algorithm to perform an exhaustive search for
	 * a solution (therefore it is slow relative to other
	 * {@link micycle.pgs.PGS_Optimisation PGS_Optimisation} methods).
	 * 
	 * @param shape
	 * @param fast  whether to compute MIR based on a lower resolution input. When
	 *              true, processing is ~6 times faster but potentially a little
	 *              inaccurate
	 * @return a rectangle shape
	 * @since 1.3.0
	 * @see #maximumInscribedRectangle(PShape) maximumInscribedRectangle() -- the
	 *      largest rectangle of arbitrary orientation
	 */
	public static PShape maximumInscribedAARectangle(PShape shape, boolean fast) {
		double f = fast ? 5 : 2;
		final MaximumInscribedAARectangle mir = new MaximumInscribedAARectangle(fromPShape(shape), f);
		int[] r = mir.getInscribedRectangle();

		final GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setCentre(new Coordinate((r[0] + r[2] / 2d) * f, (r[1] + r[3] / 2d) * f));
		shapeFactory.setWidth(r[2] * f);
		shapeFactory.setHeight(r[3] * f);
		return toPShape(shapeFactory.createRectangle());
	}

	/**
	 * Finds the largest area <i>perimeter square</i> of the input. A <i>perimeter
	 * square</i> is a square whose 4 vertices each lie on the perimeter of the
	 * input shape (within the given tolerance).
	 * <p>
	 * If the input is convex, the output forms a fully inscribed square; if the
	 * input is concave the output is not necessarily inscribed.
	 * <p>
	 * Does not respect holes (for now...).
	 * 
	 * @param shape
	 * @param tolerance a value of 2-5 is usually suitable
	 * @return shape representing the maximum square
	 * @since 1.3.1
	 */
	public static PShape maximumPerimeterSquare(PShape shape, double tolerance) {
		shape = PGS_Morphology.simplify(shape, tolerance / 2);
		final Polygon p = (Polygon) PGS_Conversion.fromPShape(shape);
		Geometry buffer = p.getExteriorRing().buffer(tolerance / 2, 4);
		final Envelope e = buffer.getEnvelopeInternal();
		buffer = DouglasPeuckerSimplifier.simplify(buffer, tolerance / 2);
		final IndexedPointInAreaLocator pia = new IndexedPointInAreaLocator(buffer);
		shape = PGS_Processing.densify(shape, Math.max(0.5, tolerance)); // min of 0.5
		final List<PVector> points = PGS_Conversion.toPVector(shape);

		double maxDiagonal = 0;
		PVector[] maxAreaVertices = new PVector[0];
		for (int i = 0; i < points.size(); i++) {
			for (int j = 0; j < points.size(); j++) {
				final PVector a = points.get(i);
				final PVector b = points.get(j);
				double dist = PGS.distanceSq(a, b);

				if (dist < maxDiagonal) {
					continue;
				}

				final PVector m = PVector.add(a, b).div(2);
				final PVector n = new PVector(b.y - a.y, a.x - b.x).div(2);
				final PVector c = PVector.sub(m, n);

				final PVector d = PVector.add(m, n);
				// do envelope checks first -- slightly faster
				if (within(c, e) && within(d, e)) {
					if (pia.locate(new Coordinate(c.x, c.y)) != Location.EXTERIOR) {
						if (pia.locate(new Coordinate(d.x, d.y)) != Location.EXTERIOR) {
							maxDiagonal = dist;
							maxAreaVertices = new PVector[] { a, c, b, d, a }; // closed vertices
						}
					}
				}
			}
		}

		PShape out = PGS_Conversion.fromPVector(maxAreaVertices);
		out.setStroke(true);
		out.setStroke(micycle.pgs.color.RGB.PINK);
		out.setStrokeWeight(4);
		return out;
	}

	private static boolean within(PVector p, Envelope rect) {
		return p.x >= rect.getMinX() && p.x <= rect.getMaxX() && p.y <= rect.getMaxY() && p.y >= rect.getMinY();
	}

	/**
	 * Computes the Minimum Bounding Circle (MBC) for the points in a Geometry. The
	 * MBC is the smallest circle which covers all the vertices of the input shape
	 * (this is also known as the Smallest Enclosing Circle). This is equivalent to
	 * computing the Maximum Diameter of the input vertex set.
	 */
	public static PShape minimumBoundingCircle(PShape shape) {
		MinimumBoundingCircle mbc = new MinimumBoundingCircle(fromPShape(shape));
		final double wh = mbc.getRadius() * 2;
		Polygon circle = createEllipse(mbc.getCentre(), wh, wh);
		return toPShape(circle);
	}

	/**
	 * Computes the minimum bounding rectangle that encloses a shape. Unlike the
	 * envelope for a shape, the rectangle returned by this method can have any
	 * orientation (it's not axis-aligned).
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumBoundingRectangle(PShape shape) {
		Polygon md = (Polygon) MinimumDiameter.getMinimumRectangle(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Computes the minimum bounding ellipse that encloses a shape.
	 * 
	 * @param shape
	 * @param errorTolerance Mean-squared error tolerance (this value does not
	 *                       correspond to a pixel distance). 0.001 to 0.01
	 *                       recommended. Higher values are a looser (yet quicker)
	 *                       fit.
	 * @return
	 */
	public static PShape minimumBoundingEllipse(PShape shape, double errorTolerance) {
		final Geometry hull = fromPShape(shape).convexHull();
		final Coordinate[] coords = hull.getCoordinates();

		double[][] points = new double[coords.length][2];
		for (int i = 0; i < points.length; i++) {
			points[i][0] = coords[i].x;
			points[i][1] = coords[i].y;
		}

		final MinimumBoundingEllipse e = new MinimumBoundingEllipse(points, Math.max(errorTolerance, 0.001));
		double[][] eEoords = e.getBoundingCoordinates(100);

		final PShape ellipse = new PShape(PShape.PATH);
		ellipse.setFill(true);
		ellipse.setFill(RGB.WHITE);
		ellipse.beginShape();
		for (int i = 0; i < eEoords.length; i++) {
			ellipse.vertex((float) eEoords[i][0], (float) eEoords[i][1]);
		}
		ellipse.endShape();

		return ellipse;
	}

	/**
	 * Computes the minimum-area bounding triangle that encloses a shape.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumBoundingTriangle(PShape shape) {
		MinimumBoundingTriangle mbt = new MinimumBoundingTriangle(fromPShape(shape));
		return toPShape(mbt.getTriangle());
	}

	/**
	 * Computes the minimum diameter of a shape.
	 * <p>
	 * The minimum diameter is defined to be the width of the smallest band that
	 * contains the shape, where a band is a strip of the plane defined by two
	 * parallel lines. This can be thought of as the smallest hole that the geometry
	 * can be moved through, with a single rotation.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumDiameter(PShape shape) {
		LineString md = (LineString) MinimumDiameter.getMinimumDiameter(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Constructs the Largest Empty Circle for a set of obstacle geometries, up to a
	 * specified tolerance. Valid obstacles are point and line shapes (such as a
	 * POINTS PShape).
	 * <p>
	 * The Largest Empty Circle is the largest circle which has its center in the
	 * convex hull of the obstacles (the boundary), and whose interior does not
	 * intersect with any obstacle. The circle center is the point in the interior
	 * of the boundary which has the farthest distance from the obstacles (up to
	 * tolerance).
	 * 
	 * @param obstacles a shape representing the obstacles (points and lines)
	 * @param tolerance the distance tolerance for computing the circle center point
	 * @return
	 */
	public static PShape largestEmptyCircle(PShape obstacles, double tolerance) {
		LargestEmptyCircle lec = new LargestEmptyCircle(fromPShape(obstacles), Math.max(0.01, tolerance));
		double wh = lec.getRadiusLine().getLength() * 2;
		Polygon circle = createEllipse(PGS.coordFromPoint(lec.getCenter()), wh, wh);
		return toPShape(circle);
	}

	/**
	 * Covers a polygon with n circles such that no circle’s center lies outside the
	 * polygon. Circles will generally cover most of the shape and have some mutual
	 * overlap.
	 * 
	 * @param shape shape to cover
	 * @param n     number of circles to generate
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 * @see #circleCoverage(PShape, int, long)
	 * @since 1.3.1
	 */
	public static List<PVector> circleCoverage(PShape shape, int n) {
		return circleCoverage(shape, n, System.nanoTime());
	}

	/**
	 * Covers a polygon with n circles such that no circle’s center lies outside the
	 * polygon. Circles will generally cover most of the shape and have some mutual
	 * overlap.
	 * 
	 * @param shape shape to cover
	 * @param n     number of circles to generate
	 * @param seed  random seed
	 * @return A list of PVectors, each representing one circle: (.x, .y) represent
	 *         the center point and .z represents radius.
	 * @see #circleCoverage(PShape, int)
	 * @since 1.3.1
	 */
	public static List<PVector> circleCoverage(PShape shape, int n, long seed) {
		// same as 'Simple Methods to Represent Shapes with Sample Spheres'
		int nSeedPoints = (int) (PGS_ShapePredicates.area(shape) / 100); // ~one point every 10 units
		List<PVector> points = PGS_Processing.generateRandomGridPoints(shape, nSeedPoints, false, 0.5, seed);
		points.addAll(PGS_Conversion.toPVector(shape)); // incl. shape vertices

		List<PVector> circles = new ArrayList<>(n);
		PGS_PointSet.cluster(points, n, seed).forEach(group -> {
			if (group.size() < 2) { // unlikely
				return;
			}
			Geometry clusterPoints = PGS.GEOM_FACTORY.createMultiPointFromCoords(PGS.toCoords(group));
			MinimumBoundingCircle mbc = new MinimumBoundingCircle(clusterPoints);
			Coordinate mbcp = mbc.getCentre();
			circles.add(new PVector((float) mbcp.x, (float) mbcp.y, (float) mbc.getRadius()));
		});

		return circles;
	}

	public enum RectPackHeuristic {

		/**
		 * Packs rectangles such that the wasted/empty area within the bin is minimised.
		 * This heuristic tends to generate packings that are less dense but are better
		 * at covering the whole area (particularly if there is much spare area) than
		 * the other heuristics.
		 */
		BestAreaFit(PackingHeuristic.BestAreaFit),
		/**
		 * Packs rectangles such that the total touching perimeter length is maximised.
		 * In practice, rectangles are packed against the left-most and upper-most
		 * boundary of the bin first.
		 */
		TouchingPerimeter(PackingHeuristic.TouchingPerimeter),
		/**
		 * Packs rectangles such that the distance between the top-right corner of each
		 * rectangle and that of the bin is maximised.
		 */
		TopRightCornerDistance(PackingHeuristic.TopRightCornerDistance);

		private final PackingHeuristic h;

		private RectPackHeuristic(PackingHeuristic h) {
			this.h = h;
		}
	}

	/**
	 * Packs a collection of rectangles, according to the given packing heuristic,
	 * into rectangular 2D bin(s). Within each bin rectangles are packed flush with
	 * each other, having no overlap. Each rectangle is packed parallel to the edges
	 * of the plane.
	 * <p>
	 * When packed rectangles fill one bin, any remaining rectangles will be packed
	 * into additional bin(s).
	 * 
	 * @param rectangles a collection of rectangles (represented by PVectors),
	 *                   specifying their width (.x) and height (.y)
	 * @param binWidth   the width of each bin's area in which to pack the
	 *                   rectangles
	 * @param binHeight  the height of each bin's area in which to pack the
	 *                   rectangles
	 * @param heuristic  the packing heuristic to use. The heuristic determines
	 *                   rules for how every subsequent rectangle is placed
	 * @since 1.3.1
	 * @return a GROUP PShape, where each immediate child is a GROUP shape
	 *         corresponding to a bin; the child shapes of each bin are rectangles.
	 *         Bins are positioned at (0, 0).
	 */
	public static PShape rectPack(List<PVector> rectangles, int binWidth, int binHeight, RectPackHeuristic heuristic) {
		RBPSolution packer = new RBPSolution(binWidth, binHeight);
		List<Rect> rects = rectangles.stream().map(p -> Rect.of(Math.round(p.x), Math.round(p.y))).collect(Collectors.toList());

		packer.pack(rects, heuristic.h);

		PShape bins = new PShape(GROUP);
		packer.getBins().forEach(bin -> {
			PShape binGroup = new PShape(GROUP);
			bin.getPackedRects().forEach(r -> {
				binGroup.addChild(PGS.createRect(r.x, r.y, r.width, r.height));
			});
			bins.addChild(binGroup);
		});
		return bins;
	}

	/**
	 * Returns the nearest point of the shape to the given point. If the shape is
	 * has multiple children/geometries (a GROUP shape), the single closest point is
	 * returned.
	 * 
	 * @param shape
	 * @param point
	 * @return
	 * @see #closestPoints(PShape, PVector)
	 */
	public static PVector closestPoint(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		Coordinate coord = DistanceOp.nearestPoints(g, PGS.pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Returns the nearest point for each "island" / separate polygon in the GROUP
	 * input shape.
	 * 
	 * @param shape a GROUP shape
	 * @param point
	 * @return list of closest points for each child shape. Output is identical to
	 *         {@link #closestPoint(PShape, PVector)} if the input shape is a single
	 *         polygon
	 * @see #closestPoint(PShape, PVector)
	 */
	public static List<PVector> closestPoints(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		ArrayList<PVector> points = new ArrayList<>();
		for (int i = 0; i < g.getNumGeometries(); i++) {
			final Coordinate coord = DistanceOp.nearestPoints(g.getGeometryN(i), PGS.pointFromPVector(point))[0];
			points.add(PGS.toPVector(coord));
		}
		return points;
	}

	/**
	 * Computes the closest pair of points in a set of points. This method runs in
	 * O(n*log(n)), rather than the naive O(n*n) brute-force approach.
	 * 
	 * @param points a set of 2D points, represented by PVectors
	 * @return a List<PVector> containing exactly two elements which are the closest
	 *         pair of points among those in the set.
	 * @since 1.1.0
	 * @see #farthestPointPair(Collection)
	 */
	public static List<PVector> closestPointPair(Collection<PVector> points) {
		final ClosestPointPair closestPointPair = new ClosestPointPair(points);
		return closestPointPair.execute();
	}

	/**
	 * Computes the farthest pair of points (the "diametral pair") in a set of n
	 * points.
	 * <p>
	 * This method runs in O(n*log(n)), rather than the naive O(n*n) brute-force
	 * approach. However, it must first compute the convex hull of the point set, so
	 * there is more overhead; on small datasets, the brute-force approach is likely
	 * faster).
	 * 
	 * @param points a set of 2D points, represented by PVectors
	 * @return a List<PVector> containing exactly two elements which are the
	 *         farthest pair of points among those in the set.
	 * @since 1.1.0
	 * @see #closestPointPair(Collection)
	 * @see #closestPoints(PShape, PVector)
	 */
	public static List<PVector> farthestPointPair(Collection<PVector> points) {
		final FarthestPointPair fpp = new FarthestPointPair(points);
		final List<PVector> out = new ArrayList<>();
		out.add(fpp.either());
		out.add(fpp.other());
		return out;
	}

	/**
	 * Sorts the faces/child shapes of a GROUP shape according to hilbert curve
	 * index of each face's centroid coordinate. This ensures that nearby faces have
	 * a similar index in the list of children.
	 * 
	 * @param mesh group shape
	 * @return a copy of the input shape, having the same faces/child shapes in a
	 *         different order
	 * @since 1.3.0
	 */
	public static PShape hilbertSortFaces(PShape mesh) {
		Map<PVector, PShape> map = new HashMap<>(mesh.getChildCount());
		PGS_Conversion.getChildren(mesh).forEach(child -> {
			PVector centroid = PGS_ShapePredicates.centroid(child);
			map.put(centroid, child);
		});

		List<PVector> points = new ArrayList<>(map.keySet());
		return PGS_Conversion.flatten(PGS_PointSet.hilbertSort(points).stream().map(map::get).collect(Collectors.toList()));
	}

	/**
	 * Solves the Problem of Apollonius (finding a circle tangent to three other
	 * circles in the plane). Circles are represented by PVectors, where the z
	 * coordinate is interpreted as radius.
	 * 
	 * @param c1 One of the circles in the problem
	 * @param c2 One of the circles in the problem
	 * @param c3 One of the circles in the problem
	 * @param s1 An indication if the solution should be externally or internally
	 *           tangent (+1/-1) to c1
	 * @param s2 An indication if the solution should be externally or internally
	 *           tangent (+1/-1) to c2
	 * @param s3 An indication if the solution should be externally or internally
	 *           tangent (+1/-1) to c3
	 * @return The circle (as a PVector) that is tangent to c1, c2 and c3.
	 */
	public static PVector solveApollonius(PVector c1, PVector c2, PVector c3, int s1, int s2, int s3) {
		// https://github.com/DIKU-Steiner/ProGAL/blob/master/src/ProGAL/geom2d/ApolloniusSolver.java

		double x1 = c1.x;
		double y1 = c1.y;
		double r1 = c1.z;
		double x2 = c2.x;
		double y2 = c2.y;
		double r2 = c2.z;
		double x3 = c3.x;
		double y3 = c3.y;
		double r3 = c3.z;

		// Currently optimized for fewest multiplications. Should be optimized for
		// readability
		double v11 = 2 * x2 - 2 * x1;
		double v12 = 2 * y2 - 2 * y1;
		double v13 = x1 * x1 - x2 * x2 + y1 * y1 - y2 * y2 - r1 * r1 + r2 * r2;
		double v14 = 2 * s2 * r2 - 2 * s1 * r1;

		double v21 = 2 * x3 - 2 * x2;
		double v22 = 2 * y3 - 2 * y2;
		double v23 = x2 * x2 - x3 * x3 + y2 * y2 - y3 * y3 - r2 * r2 + r3 * r3;
		double v24 = 2 * s3 * r3 - 2 * s2 * r2;

		double w12 = v12 / v11;
		double w13 = v13 / v11;
		double w14 = v14 / v11;

		double w22 = v22 / v21 - w12;
		double w23 = v23 / v21 - w13;
		double w24 = v24 / v21 - w14;

		double P = -w23 / w22;
		double Q = w24 / w22;
		double M = -w12 * P - w13;
		double N = w14 - w12 * Q;

		double a = N * N + Q * Q - 1;
		double b = 2 * M * N - 2 * N * x1 + 2 * P * Q - 2 * Q * y1 + 2 * s1 * r1;
		double c = x1 * x1 + M * M - 2 * M * x1 + P * P + y1 * y1 - 2 * P * y1 - r1 * r1;

		// Find a root of a quadratic equation. This requires the circle centers not
		// to be e.g. colinear
		double D = b * b - 4 * a * c;
		double rs = (-b - Math.sqrt(D)) / (2 * a);
		double xs = M + N * rs;
		double ys = P + Q * rs;
		return new PVector((float) xs, (float) ys, (float) rs);
	}

}
