package micycle.pgs;

import static micycle.pgs.PGS_Construction.createCircle;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static processing.core.PConstants.GROUP;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.locationtech.jts.algorithm.MinimumAreaRectangle;
import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.algorithm.construct.LargestEmptyCircle;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
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
import micycle.pgs.color.Colors;
import micycle.pgs.commons.ClosestPointPair;
import micycle.pgs.commons.FarthestPointPair;
import micycle.pgs.commons.LargestEmptyCircles;
import micycle.pgs.commons.MaximumInscribedAARectangle;
import micycle.pgs.commons.MaximumInscribedRectangle;
import micycle.pgs.commons.MaximumInscribedTriangle;
import micycle.pgs.commons.MinimumBoundingEllipse;
import micycle.pgs.commons.MinimumBoundingTriangle;
import micycle.pgs.commons.Nullable;
import micycle.pgs.commons.SpiralIterator;
import micycle.pgs.commons.VisibilityPolygon;
import processing.core.PShape;
import processing.core.PVector;
import whitegreen.dalsoo.DalsooPack;

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
	 * @return polygonal shape having 4 coordinates
	 * @deprecated since 2.1; use {@link micycle.pgs.PGS_Hull#boundingBox(PShape)
	 *             boundingBox(PShape)} instead.
	 */
	@Deprecated
	public static PShape envelope(PShape shape) {
		return toPShape(fromPShape(shape).getEnvelope());
	}

	/**
	 * Computes the maximum inscribed circle within a given shape.
	 * <p>
	 * The Maximum Inscribed Circle (MIC) is defined as the largest possible circle
	 * that can be completely contained within the area of the input shape. It is
	 * determined by locating a point inside the shape that has the greatest
	 * distance from the shape's boundary (i.e., the center of the MIC), and
	 * returning a circle centered at this point with a radius equal to that
	 * distance.
	 * </p>
	 * <p>
	 * This method automatically selects a reasonable tolerance value for computing
	 * the center point of the MIC, balancing precision and computational
	 * efficiency.
	 * </p>
	 *
	 * @param shape the {@link PShape} representing the area within which to compute
	 *              the MIC
	 * @return a {@link PShape} instance representing the maximum inscribed circle
	 * @since 2.1
	 */
	public static PShape maximumInscribedCircle(PShape shape) {
		MaximumInscribedCircle mic = new MaximumInscribedCircle(fromPShape(shape));
		final double r = mic.getRadiusLine().getLength();
		Polygon circle = createCircle(PGS.coordFromPoint(mic.getCenter()), r);
		return toPShape(circle);
	}

	/**
	 * Computes the maximum inscribed circle within a given shape, using a specified
	 * tolerance.
	 * <p>
	 * The Maximum Inscribed Circle (MIC) is the largest possible circle that can be
	 * fully contained within the area of the input shape. The center of the MIC is
	 * the point in the interior that is farthest from the shape's boundary, and the
	 * radius is the distance from this center point to the closest boundary point.
	 * </p>
	 * 
	 * @param shape     the {@link PShape} representing the area within which to
	 *                  compute the MIC
	 * @param tolerance the distance tolerance for computing the center point; must
	 *                  be non-negative. Smaller values result in a more accurate
	 *                  circle but may require more computation (typical values are
	 *                  around 1).
	 * @return a {@link PShape} instance representing the maximum inscribed circle
	 */
	public static PShape maximumInscribedCircle(PShape shape, double tolerance) {
		MaximumInscribedCircle mic = new MaximumInscribedCircle(fromPShape(shape), tolerance);
		final double r = mic.getRadiusLine().getLength();
		Polygon circle = createCircle(PGS.coordFromPoint(mic.getCenter()), r);
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
		Polygon circle = createCircle(p.getCoordinate(), radius);
		return toPShape(circle);
	}

	/**
	 * Finds an approximate largest area rectangle (of arbitrary orientation)
	 * contained within a polygon.
	 * 
	 * @param shape a polygonal shape
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
	 * Finds an approximate largest area triangle (of arbitrary orientation)
	 * contained within a polygon.
	 * 
	 * @param shape a polygonal shape
	 * @return a triangular shape
	 * @since 2.1
	 */
	public static PShape maximumInscribedTriangle(PShape shape) {
		Polygon polygon = (Polygon) fromPShape(shape);
		var mit = new MaximumInscribedTriangle(polygon);
		return toPShape(mit.computeMIT());
	}

	/**
	 * Finds the rectangle with a maximum area whose sides are parallel to the
	 * x-axis and y-axis ("axis-aligned"), contained/insribed within a convex shape.
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
	 * The method does not respect holes (for now...).
	 * 
	 * @param shape
	 * @param tolerance a value of 2-5 is usually suitable
	 * @return shape representing the maximum square
	 * @since 1.4.0
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
		for (final PVector a : points) {
			for (final PVector b : points) {
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
		out.setStroke(micycle.pgs.color.Colors.PINK);
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
		final double r = mbc.getRadius();
		Polygon circle = createCircle(mbc.getCentre(), r);
		return toPShape(circle);
	}

	/**
	 * Computes the minimum-width bounding rectangle that encloses a shape. Unlike
	 * the envelope for a shape, the rectangle returned by this method can have any
	 * orientation (it's not axis-aligned).
	 * <p>
	 * The minimum-width enclosing rectangle does not necessarily have the minimum
	 * possible area. Use {@link #minimumAreaRectangle(PShape)
	 * minimumAreaRectangle()} to compute this.
	 * 
	 * @param shape The shape to compute the minimum bounding rectangle for.
	 * @return A PShape object representing the minimum bounding rectangle.
	 */
	public static PShape minimumWidthRectangle(PShape shape) {
		Geometry md = MinimumDiameter.getMinimumRectangle(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Computes the minimum-area rectangle that encloses a shape.
	 * <p>
	 * The minimum-area enclosing rectangle does not necessarily have the minimum
	 * possible width. Use {@link #minimumWidthRectangle(PShape)
	 * minimumBoundingRectangle()} to compute this.
	 * 
	 * @param shape The shape to compute the minimum-area rectangle for.
	 * @return A PShape object representing the minimum-area rectangle.
	 * @since 1.4.0
	 */
	public static PShape minimumAreaRectangle(PShape shape) {
		return toPShape(MinimumAreaRectangle.getMinimumRectangle(fromPShape(shape)));
	}

	/**
	 * Computes the minimum bounding ellipse that encloses a shape.
	 * 
	 * @param shape
	 * @param errorTolerance Mean-squared error tolerance (this value does not
	 *                       correspond to a pixel distance). 0.001 to 0.01
	 *                       recommended. Higher values are a looser (yet quicker)
	 *                       fit.
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
		ellipse.setFill(Colors.WHITE);
		ellipse.beginShape();
		for (double[] eEoord : eEoords) {
			ellipse.vertex((float) eEoord[0], (float) eEoord[1]);
		}
		ellipse.endShape();

		return ellipse;
	}

	/**
	 * Computes the minimum-area bounding triangle that encloses a shape.
	 * 
	 * @param shape
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
	 */
	public static PShape minimumDiameter(PShape shape) {
		LineString md = (LineString) MinimumDiameter.getMinimumDiameter(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Computes the minimum-width annulus (the donut-like region between two
	 * concentric circles with minimal width that encloses the given
	 * {@link PShape}).
	 * <p>
	 * The annulus is defined as the region between two concentric circles (with
	 * computed center and radii), such that all vertices of the input shape lie
	 * between the inner and outer circle, and the distance between these circles
	 * (the "width" of the annulus) is minimized.
	 * <p>
	 * The algorithm considers only the <b>vertices of the input shape</b> (not the
	 * filled area or edges) as points to be enclosed.
	 *
	 * @param shape the {@link PShape} whose <b>vertices</b> are to be enclosed by
	 *              the minimum-width annulus; must be non-null and contain at least
	 *              three points
	 * @return a {@link PShape} representing the minimum-width annulus (as a ring
	 *         shape)
	 * @since 2.1
	 */
	public static PShape minimumWidthAnnulus(PShape shape) {
		var points = PGS_Conversion.toPVector(shape);
		var d = PGS_ShapePredicates.diameter(shape);
		var convexHull = PGS_Conversion.toPVector(PGS_Hull.convexHull(points));
		var tree = PGS.makeKdtree(points);

		var bounds = new double[4];
		PGS_Hull.boundingBox(shape, bounds); // write to bounds
		bounds[0] -= d;
		bounds[1] -= d;
		bounds[2] += d;
		bounds[3] += d;

		var vd = PGS_Voronoi.innerVoronoi(points, bounds);
		var fpvd = PGS_Voronoi.farthestPointVoronoi(points);

		/*
		 * NOTE here we find inner edges/vertices in a generic way without geometric
		 * shortcuts given by VD/FPVD geometry (i.e. we know the circumcircle of each
		 * FPVD vertex -- its circumradius gives us outerR immediately).
		 */
		var a = PGS_Meshing.extractInnerEdgesAndVertices(vd);
		var b = PGS_Meshing.extractInnerEdgesAndVertices(fpvd);
		var vdVertices = a.getRight();
		var vdEdges = a.getLeft();
		var fpvdVertices = b.getRight();
		var fpvdEdges = b.getLeft();

		var overlayVerices = PGS_SegmentSet.intersections(vdEdges, fpvdEdges);

		/*
		 * Candidate centers for the smallest-width annulus have 3 sources. For each
		 * candidate we find the distance both to the closest and farthest vertex. The
		 * difference between these distances is the annulus width; we select the
		 * candidate having the smallest width.
		 */
		/*-
		 * The 3 centerpoint sources are:
		 * 		Vertices of the FPVD
		 * 		Vertices of the VD
		 * 		Vertices from the intersection between edges of VD and FPVD
		 */
		var candidates = new ArrayList<PVector>();
		candidates.addAll(vdVertices);
		candidates.addAll(fpvdVertices);
		candidates.addAll(overlayVerices);

		var result = candidates.parallelStream().map(v -> {
			// farthest point must lie on convex hull
			var far = PGS_Optimisation.farthestPoint(convexHull, v);
			var close = tree.query1nn(new double[] { v.x, v.y });
			double outerR = far.dist(v);
			double innerR = close.dist();
			return Triple.of(v, outerR, innerR);
		}).min(Comparator.comparingDouble(t -> t.getMiddle() - t.getRight())).get();

		return PGS_Construction.createRing(result.getLeft().x, result.getLeft().y, result.getMiddle(), result.getRight());
	}

	/**
	 * Computes the largest empty circle that does not intersect any obstacles (up
	 * to a specified tolerance).
	 * <p>
	 * Valid obstacles are point, line or polygonal shapes.
	 * <p>
	 * The circle center lies within the interior of the convex hull of the
	 * obstacles.
	 * 
	 * @param obstacles A PShape representing the obstacles.
	 * @param tolerance A double representing the tolerance for the circle
	 *                  computation.
	 * @return A PShape representing the largest empty circle that does not
	 *         intersect the obstacles and lies within the specified boundary.
	 */
	public static PShape largestEmptyCircle(PShape obstacles, double tolerance) {
		return largestEmptyCircle(obstacles, null, tolerance);
	}

	/**
	 * Computes the largest empty circle that does not intersect any obstacles and
	 * lies within the specified boundary (up to a specified tolerance).
	 * <p>
	 * Valid obstacles are point, line or polygonal shapes.
	 * <p>
	 * The circle center is the point in the interior of the boundary which has the
	 * farthest distance from the obstacles (up to the tolerance).
	 * 
	 * @param obstacles A PShape representing the obstacles.
	 * @param boundary  A PShape representing the polygonal boundary, or null if
	 *                  there is no boundary constraint.
	 * @param tolerance A double representing the tolerance for the circle
	 *                  computation.
	 * @return A PShape representing the largest empty circle that does not
	 *         intersect the obstacles and lies within the specified boundary.
	 */
	public static PShape largestEmptyCircle(PShape obstacles, @Nullable PShape boundary, double tolerance) {
		LargestEmptyCircle lec = new LargestEmptyCircle(fromPShape(obstacles), boundary == null ? null : fromPShape(boundary), Math.max(0.01, tolerance));
		double r = lec.getRadiusLine().getLength();
		Polygon circle = createCircle(PGS.coordFromPoint(lec.getCenter()), r);
		return toPShape(circle);
	}

	/**
	 * Computes the {@code n} largest empty circles that do not intersect any
	 * obstacles (nor each other) within an optional {@code boundary}.
	 * <p>
	 * The empty circles are found with a specified {@code tolerance} value, which
	 * limits the precision of the computation.
	 * <p>
	 * Valid obstacles are point, line or polygonal shapes.
	 *
	 * @param obstacles PShape containing the obstacles in the 2D space
	 * @param boundary  polygonal PShape defining the boundary of the space, or
	 *                  {@code null} if no boundary is defined (in which case the
	 *                  convex hull of obstacles is used as boundary).
	 * @param n         the number of largest empty circles to find
	 * @param tolerance the tolerance value for the computation
	 * @return a list of {@code PVector} objects representing the centers and radii
	 *         of the found largest empty circles as {@code PVector(x, y, r)}, where
	 *         {@code x} and {@code y} are the center coordinates, and {@code r} is
	 *         the radius
	 * @since 1.4.0
	 */
	public static List<PVector> largestEmptyCircles(PShape obstacles, @Nullable PShape boundary, int n, double tolerance) {
		tolerance = Math.max(0.01, tolerance);
		LargestEmptyCircles lecs = new LargestEmptyCircles(obstacles == null ? null : fromPShape(obstacles), boundary == null ? null : fromPShape(boundary),
				tolerance);

		final List<PVector> out = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			double[] c = lecs.findNextLEC();
			out.add(new PVector((float) c[0], (float) c[1], (float) c[2]));
		}

		return out;
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
	 * @since 1.4.0
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
	 * @since 1.4.0
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

	/**
	 * Various packing heuristics for
	 * {@link micycle.pgs.PGS_Optimisation#rectPack(List, int, int, RectPackHeuristic)
	 * rectpack()}.
	 */
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
	 * @since 1.4.0
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
	 * Packs a list of irregular polygonal shapes into (potentially multiple)
	 * rectangular containers (bins), while attempting to minimise the occupied
	 * space of the packing.
	 * <p>
	 * Every bin has the dimensions given by the width and height parameters; when
	 * packed shapes fill/overflow one bin, any remaining shapes will be packed into
	 * additional bin(s). Multiple bins are arranged in a grid, having the maximum
	 * number of columns specified by the <code>binColumns</code> parameter.
	 * <p>
	 * Bins are packed top-to-bottom vertically.
	 * 
	 * @param shapes     a list of PShapes to be packed within a bin(s)
	 * @param binWidth   the width of each bin/container to pack the shapes into
	 * @param binHeight  the height of each bin/container to pack the shapes into
	 * @param binColumns the number of columns to arrange the bins into (>= 1, only
	 *                   applies when there are multiple bins).
	 * @param spacing    the amount of spacing between each packed shape (>= 0).
	 * @return a new GROUP PShape object containing the packed shapes arranged in
	 *         columns
	 * @since 1.4.0
	 */
	public static PShape binPack(List<PShape> shapes, double binWidth, double binHeight, int binColumns, double spacing) {
		if (shapes.isEmpty()) {
			return new PShape();
		}
		binColumns = Math.max(1, binColumns); // enforce >= 1
		double[][][] polys = new double[shapes.size()][0][0];
		for (int i = 0; i < polys.length; i++) {
			polys[i] = PGS_Conversion.toArray(shapes.get(i), false);
		}

		PShape packing = new PShape(GROUP);
		DalsooPack pack = new DalsooPack(polys, spacing, null, 1, binWidth, binHeight, 0); // pack vertically
		pack.packAll(true, false); // use abey pack -- most efficient method
		pack.getPackedPolys(binColumns).forEach(p -> packing.addChild(PGS_Conversion.fromArray(p, true)));
		PGS_Conversion.disableAllStroke(packing);

		return packing;
	}

	/**
	 * Returns the closest vertex of a shape to a query point. For GROUP shapes, any
	 * child geometry's vertex may be returned.
	 *
	 * @param shape      the PShape to search for the closest vertex
	 * @param queryPoint the query PVector
	 * @return a new PVector at the position of the closest vertex (not a reference
	 *         to existing shape data)
	 * @since 2.1
	 */
	public static PVector closestVertex(PShape shape, PVector queryPoint) {
		List<PVector> vertices = PGS_Conversion.toPVector(shape);
		if (vertices.isEmpty()) {
			return null;
		}
		float minDistSq = Float.POSITIVE_INFINITY;
		PVector closest = null;
		for (PVector v : vertices) {
			float distSq = PVector.dist(v, queryPoint);
			if (distSq < minDistSq) {
				minDistSq = distSq;
				closest = v;
			}
		}
		return closest;
	}

	/**
	 * Returns the nearest point along the edges of the given shape to the specified
	 * query point.
	 * <p>
	 * This method computes the point on the perimeter (including all edges, not
	 * only the vertices) of the given shape that is closest to the given
	 * {@code point}. For composite shapes (such as GROUP shapes made of multiple
	 * child geometries), the single closest point across all children is returned.
	 * </p>
	 * <p>
	 * <strong>Note:</strong> The nearest location may be somewhere along an edge of
	 * the shape, not necessarily at one of the original shape's vertices.
	 * </p>
	 *
	 * @param shape the {@code PShape} to search for the closest boundary point.
	 * @param point the {@code PVector} point to which the nearest point is sought.
	 * @return a new {@code PVector} representing the exact coordinates of the
	 *         closest point on the shape's boundary or edge (not a reference to the
	 *         original coordinate).
	 * @see #closestPoints(PShape, PVector)
	 */
	public static PVector closestPoint(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		Coordinate coord = DistanceOp.nearestPoints(g, PGS.pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Finds the closest point in the collection to a specified point.
	 *
	 * @param points the collection of points to search within
	 * @param point  the point to find the closest neighbor for
	 * @return the closest point from the collection to the specified point
	 * @since 2.1
	 */
	public static PVector closestPoint(Collection<PVector> points, PVector point) {
		if (points == null || points.isEmpty()) {
			return null; // Handle empty or null collection
		}

		PVector closest = null;
		float minDistanceSq = Float.MAX_VALUE;

		for (PVector p : points) {
			float dx = p.x - point.x;
			float dy = p.y - point.y;
			float distanceSq = dx * dx + dy * dy;
			if (distanceSq < minDistanceSq) {
				minDistanceSq = distanceSq;
				closest = p;
			}
		}

		return closest;
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
	 * @return a List of PVectors containing exactly two elements which are the
	 *         closest pair of points among those in the set.
	 * @since 1.1.0
	 * @see #farthestPointPair(Collection)
	 */
	public static List<PVector> closestPointPair(Collection<PVector> points) {
		final ClosestPointPair closestPointPair = new ClosestPointPair(points);
		return closestPointPair.execute();
	}

	/**
	 * Returns the farthest vertex of a shape from a query point. For GROUP shapes,
	 * any child geometry's vertex may be returned.
	 *
	 * @param shape      the PShape to search for the farthest vertex
	 * @param queryPoint the query PVector
	 * @return a new PVector at the position of the farthest vertex (not a reference
	 *         to existing shape data)
	 * @since 2.1
	 */
	public static PVector farthestVertex(PShape shape, PVector queryPoint) {
		List<PVector> vertices = PGS_Conversion.toPVector(shape);
		if (vertices.isEmpty()) {
			return null;
		}
		float maxDistSq = Float.NEGATIVE_INFINITY;
		PVector farthest = null;
		for (PVector v : vertices) {
			float distSq = PVector.dist(v, queryPoint);
			if (distSq > maxDistSq) {
				maxDistSq = distSq;
				farthest = v;
			}
		}
		return farthest;
	}

	/**
	 * Finds the farthest point in the collection from a specified point.
	 *
	 * @param points the collection of points to search within
	 * @param point  the point from which the farthest neighbor is sought
	 * @return the farthest point from the collection to the specified point
	 * @since 2.1
	 */
	public static PVector farthestPoint(Collection<PVector> points, PVector point) {
		if (points == null || points.isEmpty()) {
			return null; // Handle empty or null collection
		}

		PVector farthest = null;
		float maxDistanceSq = Float.NEGATIVE_INFINITY;

		for (PVector p : points) {
			float dx = p.x - point.x;
			float dy = p.y - point.y;
			float distanceSq = dx * dx + dy * dy;
			if (distanceSq > maxDistanceSq) {
				maxDistanceSq = distanceSq;
				farthest = p;
			}
		}

		return farthest;
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
	 * @return a List of PVectors containing exactly two elements which are the
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
	 * Reorders the faces of a mesh into an anti-clockwise “spiral” (breadth-first
	 * rings) starting from a given face, then returns a new, flattened PShape
	 * containing exactly those faces in spiral order.
	 * 
	 * @param mesh      mesh-like GROUP PShape
	 * @param startFace One of the child‐faces of {@code mesh}. This face will
	 *                  appear first in the returned ordering; subsequent faces
	 *                  follow in concentric breadth‐first “rings” around it, sorted
	 *                  anti-clockwise.
	 * @return A new, flattened PShape whose set of faces equals the children of
	 *         {@code mesh}, but ordered in a spiral starting at
	 *         {@code startingFace}.
	 * @since 2.1
	 */
	public static PShape spiralSortFaces(PShape mesh, PShape startFace) {
		var faces = SpiralIterator.spiral(startFace, PGS_Conversion.getChildren(mesh));
		return PGS_Conversion.flatten(faces);
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

	/**
	 * Computes a visibility polygon / isovist, the area visible from a given point
	 * in a space, considering occlusions caused by obstacles. In this case,
	 * obstacles comprise the line segments of input shape.
	 * 
	 * @param obstacles shape representing obstacles, which may have any manner of
	 *                  polygon and line geometries.
	 * @param viewPoint view point from which to compute visibility. If the input if
	 *                  polygonal, the viewpoint may lie outside the polygon.
	 * @return a polygonal shape representing the visibility polygon.
	 * @since 1.4.0
	 * @see #visibilityPolygon(PShape, Collection)
	 */
	public static PShape visibilityPolygon(PShape obstacles, PVector viewPoint) {
		VisibilityPolygon vp = new VisibilityPolygon();
		vp.addGeometry(fromPShape(obstacles));
		return toPShape(vp.getIsovist(PGS.coordFromPVector(viewPoint), true));
	}

	/**
	 * Computes a visibility polygon / isovist, the area visible from a set of given
	 * points in space, considering occlusions caused by obstacles. In this case,
	 * obstacles comprise the line segments of input shape.
	 * 
	 * @param obstacles  shape representing obstacles, which may have any manner of
	 *                   polygon and line geometries.
	 * @param viewPoints viewpoints from which to compute visibility. If the input
	 *                   if polygonal, viewpoints may lie outside the polygon.
	 * @return a polygonal shape representing the visibility polygon (possibly a
	 *         GROUP shape of disjoint visibility polygons).
	 * @since 1.4.0
	 * @see #visibilityPolygon(PShape, PVector)
	 */
	public static PShape visibilityPolygon(PShape obstacles, Collection<PVector> viewPoints) {
		VisibilityPolygon vp = new VisibilityPolygon();
		vp.addGeometry(fromPShape(obstacles));
		return toPShape(vp.getIsovist(new CoordinateList(PGS.toCoords(viewPoints)), true));
	}

}
