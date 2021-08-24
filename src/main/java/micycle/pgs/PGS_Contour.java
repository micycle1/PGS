package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS.prepareLinesPShape;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.tinfour.common.Vertex;
import org.tinfour.contour.Contour;
import org.tinfour.contour.ContourBuilderForTin;
import org.tinfour.standard.IncrementalTin;
import org.twak.camp.Corner;
import org.twak.camp.Machine;
import org.twak.camp.Skeleton;
import org.twak.utils.collections.Loop;
import org.twak.utils.collections.LoopL;

import hageldave.jplotter.misc.Contours;
import hageldave.jplotter.renderables.Lines.SegmentDetails;
import micycle.medialAxis.MedialAxis;
import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.color.RGB;
import micycle.pgs.utility.SolubSkeleton;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Methods that produce a variety of methods for producing different kinds of
 * shape contours.
 *
 * <p>
 * A 2D contour is a closed sequence (a cycle) of 3 or more connected 2D
 * oriented straight line segments called contour edges. The endpoints of the
 * contour edges are called vertices. Each contour edge shares its endpoints
 * with at least two other contour edges.
 *
 * @author Michael Carleton
 *
 */
public final class PGS_Contour {

	/**
	 * TODO implement 'Base Point Split Algorithm for Generating Polygon Skeleton
	 * Lines'
	 */

	private PGS_Contour() {
	}

	/**
	 * Computes the medial axis of the given shape. The 3 parameters can be used to
	 * prune the medial axis according to different features (at the same time).
	 *
	 * @param shape
	 * @param axialThreshold    Prune edges based on their axial gradient. The axial
	 *                          gradient measures the change in the width of the
	 *                          shape per unit length of the axis (measured per edge
	 *                          segment). Between 0...1, where 0 is no pruning and 1
	 *                          is maximal pruning for this feature.
	 * @param distanceThreshold Prune edges based on the spatial distance between
	 *                          the medial axis root and edge's tail coordinate.
	 *                          Between 0...1, where 0 is no pruning and 1 is
	 *                          maximal pruning for this feature.
	 * @param areaThreshold     Prune edges based on the sum of each edge and its
	 *                          descendants underlying feature area. Between 0...1,
	 *                          where 0 is no pruning and 1 is maximal pruning for
	 *                          this feature.
	 * @return PShape of lines where lines represent medial axis edges
	 */
	public static PShape medialAxis(PShape shape, double axialThreshold, double distanceThreshold, double areaThreshold) {
		final Geometry g = fromPShape(shape);
		final MedialAxis m = new MedialAxis(g);

		final PShape lines = PGS.prepareLinesPShape(RGB.PINK, PShape.ROUND, 4);
		m.getPrunedEdges(axialThreshold, distanceThreshold, areaThreshold).forEach(e -> {
			lines.vertex((float) e.head.position.x, (float) e.head.position.y);
			lines.vertex((float) e.tail.position.x, (float) e.tail.position.y);
		});
		lines.endShape();
		return lines;
	}

	/**
	 * Roughly, it is the geometric graph whose edges are the traces of vertices of
	 * shrinking mitered offset curves of the polygon
	 *
	 * @param shape
	 * @param p
	 * @return shape with two children: one child contains bones; one contains
	 *         branches
	 * @see #straightSkeletonSolub(List)
	 */
	public static PShape straightSkeleton(PShape shape) {
		// https://github.com/Agent14zbz/ZTools/blob/main/src/main/java/geometry/ZSkeleton.java

		final Machine speed = new Machine(1); // every edge same speed

		Geometry g = fromPShape(shape);
		Polygon polygon;
		if (g.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			polygon = (Polygon) g;
			if (polygon.getCoordinates().length > 1000) {
				polygon = (Polygon) DouglasPeuckerSimplifier.simplify(polygon, 1);
			}
		} else {
			System.err.println("MultiPolygon not supported yet.");
			return new PShape();
		}

		HashSet<Coordinate> edgeCoordsSet = new HashSet<>();

		Skeleton skeleton;
		LoopL<org.twak.camp.Edge> loopL = new LoopL<>(); // list of loops
		ArrayList<Corner> corners = new ArrayList<>();
		Loop<org.twak.camp.Edge> loop = new Loop<>();

		LinearRing exterior = polygon.getExteriorRing();
		Coordinate[] coords = exterior.getCoordinates();
		if (!Orientation.isCCW(coords)) {
			reverse(coords); // exterior should be CCW
		}

		for (int j = 0; j < coords.length - 1; j++) {
			double a = coords[j].x;
			double b = coords[j].y;
			corners.add(new Corner(a, b));
			edgeCoordsSet.add(coords[j]);
		}
		corners.add(new Corner(coords[0].x, coords[0].y)); // close loop
		edgeCoordsSet.add(coords[0]); // close loop

		for (int j = 0; j < corners.size() - 1; j++) {
			org.twak.camp.Edge edge = new org.twak.camp.Edge(corners.get(j), corners.get((j + 1) % (corners.size() - 1)));
			edge.machine = speed;
			loop.append(edge);
		}
		loopL.add(loop);

		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			corners = new ArrayList<>();
			LinearRing hole = polygon.getInteriorRingN(i);
			coords = hole.getCoordinates();
			if (Orientation.isCCW(coords)) {
				reverse(coords);
			}

			for (int j = 0; j < coords.length - 1; j++) {
				double a = coords[j].x;
				double b = coords[j].y;
				corners.add(new Corner(a, b));
				edgeCoordsSet.add(coords[j]);
			}
			corners.add(new Corner(coords[0].x, coords[0].y)); // close loop
			edgeCoordsSet.add(coords[0]); // close loop

			loop = new Loop<>();
			for (int j = 0; j < corners.size() - 1; j++) {
				org.twak.camp.Edge edge = new org.twak.camp.Edge(corners.get(j), corners.get((j + 1) % (corners.size() - 1)));
				edge.machine = speed;
				loop.append(edge);
			}
			loopL.add(loop);
		}

		PShape lines = new PShape();
		lines.setFamily(PConstants.GROUP);
		PShape bones = prepareLinesPShape(null, null, 4);
		PShape branches = prepareLinesPShape(RGB.composeColor(40, 235, 180, 128), null, null);
		try {
			skeleton = new Skeleton(loopL, true);
			skeleton.skeleton();

			skeleton.output.edges.map.values().forEach(e -> {
				boolean a = edgeCoordsSet.contains(new Coordinate(e.start.x, e.start.y));
				boolean b = edgeCoordsSet.contains(new Coordinate(e.end.x, e.end.y));
				if (a ^ b) { // branch (xor)
					branches.vertex((float) e.start.x, (float) e.start.y);
					branches.vertex((float) e.end.x, (float) e.end.y);
				} else {
					if (!a) { // bone
						bones.vertex((float) e.start.x, (float) e.start.y);
						bones.vertex((float) e.end.x, (float) e.end.y);
					}
				}
			});
		} catch (Exception ignore) {
			// hide init or collision errors from console
		}

		bones.endShape();
		branches.endShape();
		lines.addChild(branches);
		lines.addChild(bones);
		return lines;
	}

	/**
	 * An alternative straight skeleton implementation. Not robust, but fast. Does
	 * not support holes.
	 *
	 * @param shape a hull
	 * @return shape with two children: one child contains bones; one contains
	 *         branches
	 * @see #straightSkeleton(PShape) straightSkeleton(PShape). Fully robust, but
	 *      slower.
	 */
	public static PShape straightSkeletonSolub(PShape shape) {
		ArrayList<PVector> points = new ArrayList<>();
		Polygon polygon = (Polygon) fromPShape(shape);

		Coordinate[] coords = polygon.getExteriorRing().getCoordinates();
		if (Orientation.isCCW(coords)) {
			reverse(coords); // exterior should be CW
		}
		for (Coordinate coordinate : coords) {
			points.add(new PVector((float) coordinate.x, (float) coordinate.y));
		}
		points.remove(0); // remove closing point

		SolubSkeleton skeleton = new SolubSkeleton(points, 0);
		skeleton.run();

		PShape lines = new PShape();
		lines.setFamily(PConstants.GROUP);
		PShape branches = prepareLinesPShape(RGB.composeColor(40, 235, 180, 128), null, null);
		PShape bones = prepareLinesPShape(null, null, 4);

		skeleton.branches.forEach(branch -> {
			branches.vertex(branch.sp1.x, branch.sp1.y);
			branches.vertex(branch.sp2.x, branch.sp2.y);
		});
		skeleton.bones.forEach(bone -> {
			bones.vertex(bone.sp1.x, bone.sp1.y);
			bones.vertex(bone.sp2.x, bone.sp2.y);
		});

		bones.endShape();
		branches.endShape();
		lines.addChild(branches);
		lines.addChild(bones);
		return lines;
	}

	/**
	 * Generates a topographic-like isoline contour map from the shape's vertices.
	 * The "elevation" (or z value) of points is the euclidean distance between a
	 * point in the shape and the given "high" point.
	 * <p>
	 * Assigns each point feature a number equal to the distance between geometry's
	 * centroid and the point.
	 *
	 * @param shape
	 * @param highPoint       position of "high" point within the shape
	 * @param intervalSpacing distance between successive isolines
	 * @return
	 */
	public static PShape isolines(PShape shape, PVector highPoint, double intervalSpacing) {

		/*
		 * Also See:
		 * https://github.com/hageldave/JPlotter/blob/master/jplotter/src/main/java/
		 * hageldave/jplotter/misc/Contours.java
		 * https://blog.bruce-hill.com/meandering-triangles
		 * http://indiemaps.com/blog/2008/06/isolining-package-for-actionscript-3/
		 */

		Geometry g = fromPShape(shape);
		if (g.getCoordinates().length > 2000) {
			g = DouglasPeuckerSimplifier.simplify(g, 1);
		}
		final int buffer = (int) Math.max(10, Math.round(intervalSpacing) + 1);
		PreparedGeometry cache = PreparedGeometryFactory.prepare(g.buffer(10));

		final List<Vertex> tinVertices = new ArrayList<>(200);
		double maxDist = 0;

		/**
		 * Poisson a little faster, but isolines are more rough
		 */
//		ArrayList<PVector> randomPoints = pd.generate(e[0].x - buffer, e[0].y - buffer, e[3].x + buffer,
//				e[1].y + buffer, intervalSpacing, 6);
//		PoissonDistribution pd = new PoissonDistribution(0);
		Coordinate[] e = g.getEnvelope().getCoordinates(); // envelope/bounding box of shape
		ArrayList<PVector> randomPoints = generateGrid(e[0].x - buffer, e[0].y - buffer, e[3].x + buffer, e[1].y + buffer, intervalSpacing,
				intervalSpacing);

		for (PVector v : randomPoints) {
			/**
			 * Major bottleneck of method is isoline computation so reduce points to only
			 * those needed.
			 */
			if (cache.covers(PGS.pointFromPVector(v))) {
				double d = highPoint.dist(v);
				maxDist = Math.max(d, maxDist);
				tinVertices.add(new Vertex(v.x, v.y, d));
			}
//			if (g.isWithinDistance(PTS.pointFromPVector(v), 10)) {
//				double d = highPoint.dist(v);
//				maxDist = Math.max(d, maxDist);
//				tinVertices.add(new Vertex(v.x, v.y, d, 0));
//			}
		}

		final IncrementalTin tin = new IncrementalTin(intervalSpacing);
		tin.add(tinVertices, null); // insert point set; points are triangulated upon insertion

		double[] intervals = generateDoubleSequence(0, maxDist, intervalSpacing);

		/**
		 * A null valuator tells the builder to just use the z values from the vertices
		 * rather than applying any adjustments to their values.
		 */
		final ContourBuilderForTin builder = new ContourBuilderForTin(tin, null, intervals, true);

		List<org.tinfour.contour.Contour> contours = builder.getContours();

		PShape parent = new PShape(PConstants.GROUP);
		parent.setKind(PConstants.GROUP);

		LineDissolver ld = new LineDissolver();
		for (org.tinfour.contour.Contour contour : contours) {
			Coordinate[] coords = new Coordinate[contour.getCoordinates().length / 2];
			for (int i = 0; i < contour.getCoordinates().length; i += 2) {
				float vx = (float) contour.getCoordinates()[i];
				float vy = (float) contour.getCoordinates()[i + 1];
				coords[i / 2] = new Coordinate(vx, vy);
			}
			ld.add(GEOM_FACTORY.createLineString(coords));
		}

		PShape out = new PShape();
		try {
			// NOTE need to use intersection() rather than checkling whether vertices are
			// contained within the shape (faster) because vertices of longer (straight)
			// line segments may lie within the shape when the segment extends outside the
			// shape
			out = toPShape(DouglasPeuckerSimplifier.simplify(ld.getResult(), 1).intersection(g));
			PGS_Conversion.disableAllFill(out);
		} catch (Exception e2) {
			// catch non-noded intersection
		}
		return out;
	}

	/**
	 * Generates a topographic-like isoline contour map from the given points. This
	 * method uses the Z value of each PVector point as the "elevation" of that
	 * location in the map.
	 *
	 * @param points               List of PVectors: the z coordinate for each
	 *                             PVector defines the contour height at that
	 *                             location
	 * @param intervalValueSpacing contour height distance represented by successive
	 *                             isolines (e.g. a value of 1 will generate
	 *                             isolines at each 1 unit of height)
	 * @param isolineMin           minimum value represented by isolines
	 * @param isolineMax           maximum value represented by isolines
	 * @return a map of {isoline -> height of the isoline}
	 */
	public static Map<PShape, Float> isolines(Collection<PVector> points, double intervalValueSpacing, double isolineMin,
			double isolineMax) {
		// lines = max-min/spacing
		final IncrementalTin tin = new IncrementalTin(10);
		points.forEach(point -> tin.add(new Vertex(point.x, point.y, point.z)));

		double[] intervals = generateDoubleSequence(isolineMin, isolineMax, intervalValueSpacing);

		final ContourBuilderForTin builder = new ContourBuilderForTin(tin, null, intervals, false);
		List<Contour> contours = builder.getContours();

		Map<PShape, Float> isolines = new HashMap<>(contours.size());

		for (Contour contourLine : contours) {
			final PShape isoline = new PShape();
			isoline.setFamily(PShape.PATH);
			isoline.setStroke(true);
			isoline.setStrokeWeight(2);
			isoline.setStroke(RGB.PINK);
			isoline.beginShape();

			final double[] coords = contourLine.getCoordinates(); // [x1, y1, x2, y2, ...]
			for (int i = 0; i < coords.length; i += 2) {
				float vx = (float) coords[i];
				float vy = (float) coords[i + 1];
				isoline.vertex(vx, vy);
			}

			isoline.endShape();
			isolines.put(isoline, (float) contourLine.getZ());
		}

		return isolines;
	}

	/**
	 * Generates isolines from a 2D regular grid of elevation values. This method
	 * generates isolines having the exact "height" of the provided iso value,
	 * unlike the other methods that generate isolines at each iso value step.
	 *
	 * @param values   z-coordinates for the grid points
	 * @param isoValue the iso value for which the contour (iso) lines should be
	 *                 computed
	 * @param scaleX   factor to scale contour lines by in X direction. 1 returns
	 *                 the lines up to the same dimensions of the underlying matrix.
	 * @param scaleY   factor to scale contour lines by in Y direction. 1 returns
	 *                 the lines up to the same dimensions of the underlying matrix.
	 * @return
	 */
	public static PShape isolinesFromGrid(double[][] values, double isoValue, double scaleX, double scaleY) {
		PShape lines = prepareLinesPShape(null, null, null);

		List<SegmentDetails> segments = Contours.computeContourLines(values, isoValue, 0);
		segments.forEach(s -> {
			lines.vertex((float) (s.p0.getX() * scaleX), (float) (s.p0.getY() * scaleY));
			lines.vertex((float) (s.p1.getX() * scaleX), (float) (s.p1.getY() * scaleY));
		});
		lines.endShape();
		return lines;
	}

	/**
	 * Specifies the join style for offset curves.
	 */
	public enum OffsetStyle {

		MITER(BufferParameters.JOIN_MITRE), BEVEL(BufferParameters.JOIN_BEVEL), ROUND(BufferParameters.JOIN_ROUND);

		private final int style;

		private OffsetStyle(int style) {
			this.style = style;
		}
	}

	/**
	 * Produces inwards offset curves from the shape. Curves will be generated until
	 * they collapse.
	 *
	 * @param shape   a single polygon or multipolygon (GROUP PShape)
	 * @param spacing Spacing between successive offset curves. Should be >=1.
	 * @return A GROUP PShape, where each child shape is one curve
	 * @see #offsetCurvesOutward(PShape, double, int)
	 */
	public static PShape offsetCurvesInward(PShape shape, OffsetStyle style, double spacing) {
		return offsetCurves(shape, style, spacing, 0, false);
	}

	/**
	 * Produces offset curves that emanate outwards from the shape.
	 *
	 * @param shape   a single polygon or multipolygon (GROUP PShape)
	 * @param spacing Spacing between successive offset curves. Should be >=1.
	 * @param curves  number of offset curves (including the original shape outline)
	 * @return A GROUP PShape, where each child shape is one curve
	 * @see #offsetCurvesInward(PShape, double)
	 */
	public static PShape offsetCurvesOutward(PShape shape, OffsetStyle style, double spacing, final int curves) {
		return offsetCurves(shape, style, spacing, curves, true);
	}

	/**
	 * Generic method for offset curves.
	 */
	private static PShape offsetCurves(PShape shape, OffsetStyle style, double spacing, final int curves, boolean outwards) {
		Geometry g = fromPShape(shape);

		if (g.getCoordinates().length > 2000) {
			g = DouglasPeuckerSimplifier.simplify(g, 1);
		}

		final BufferParameters bufParams = new BufferParameters(8, BufferParameters.CAP_FLAT, style.style,
				BufferParameters.DEFAULT_MITRE_LIMIT);
//		bufParams.setSimplifyFactor(5); // can produce "poor" yet interesting results

		spacing = Math.max(1, Math.abs(spacing)); // ensure positive and >=1
		spacing = outwards ? spacing : -spacing;

		final PShape parent = new PShape(PConstants.GROUP);
		int currentCurves = 0;
		while ((outwards && currentCurves < curves) || (!outwards && !g.isEmpty())) {
			// current geometry's ring(s) to PATH PShape(s)
			for (LinearRing ring : new LinearRingIterator(g)) { // iterate over rings individually
				final Coordinate[] coords = ring.getCoordinates();
				PShape lines = new PShape(PShape.PATH);
				lines.setStroke(true);
				lines.setStrokeWeight(2);
				lines.setStroke(RGB.PINK);
				lines.beginShape();

				for (int i = 0; i < coords.length; i++) {
					Coordinate coord = coords[i];
					lines.vertex((float) coord.x, (float) coord.y);
				}
				lines.endShape();
				parent.addChild(lines);
			}

			BufferOp b = new BufferOp(g, bufParams);
			g = b.getResultGeometry(spacing);
			if (style != OffsetStyle.MITER) {
				// simplify because rounded buffers produce LOADS of dense vertices
				g = DouglasPeuckerSimplifier.simplify(g, outwards ? 0.1 : 0.5);
			}
			currentCurves++;
		}

		return parent;
	}

	/**
	 * Start inclusive; end exclusive
	 */
	private static double[] generateDoubleSequence(double start, double end, double step) {
		double[] sequence = new double[(int) Math.ceil((end - start) / step)];
		for (int i = 0; i < sequence.length; i++) {
			sequence[i] = start + i * step;
		}
		return sequence;
	}

	/**
	 * Generates a grid of points
	 *
	 * @param minX
	 * @param minY
	 * @param maxX
	 * @param maxY
	 * @param spacingX
	 * @param spacingY
	 * @return
	 */
	private static ArrayList<PVector> generateGrid(double minX, double minY, double maxX, double maxY, double spacingX, double spacingY) {
		ArrayList<PVector> grid = new ArrayList<>();
		double[] y = generateDoubleSequence(minY, maxY, spacingY);
		double[] x = generateDoubleSequence(minX, maxX, spacingX);

		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < y.length; j++) {
				grid.add(new PVector((float) x[i], (float) y[j]));
			}
		}
		return grid;
	}

	private static <T> void reverse(T[] a) {
		// used in straightSkeleton()
		int l = a.length;
		for (int j = 0; j < l / 2; j++) {
			T temp = a[j];
			a[j] = a[l - j - 1];
			a[l - j - 1] = temp;
		}
	}

}
