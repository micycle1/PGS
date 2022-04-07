package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS.prepareLinesPShape;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.joml.Vector2d;
import org.joml.Vector2dc;
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
import org.tinfour.utils.SmoothingFilter;
import org.twak.camp.Corner;
import org.twak.camp.Edge;
import org.twak.camp.Machine;
import org.twak.camp.Skeleton;
import org.twak.utils.collections.Loop;
import org.twak.utils.collections.LoopL;

import hageldave.jplotter.misc.Contours;
import hageldave.jplotter.renderables.Lines.SegmentDetails;
import kendzi.math.geometry.skeleton.SkeletonConfiguration;
import kendzi.math.geometry.skeleton.SkeletonOutput;
import micycle.medialAxis.MedialAxis;
import micycle.pgs.PGS.GeometryIterator;
import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.color.RGB;
import micycle.pgs.commons.PEdge;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Methods for producing different kinds of shape contours.
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
	 * Computes the straight skeleton for a shape.
	 * <p>
	 * A straight skeleton is a skeletal structure similar to the medial axis,
	 * consisting of straight-line segments only. Roughly, it is the geometric graph
	 * whose edges are the traces of vertices of shrinking mitered offset curves of
	 * the polygon.
	 *
	 * @param shape a single polygon (that can contain holes), or a multi polygon
	 *              (whose polygons can contain holes)
	 * @return when the input is a single polygon, returns a GROUP PShape containing
	 *         3 children: child 1 = GROUP PShape of skeleton faces; child 2 = LINES
	 *         PShape of branches (lines that connect skeleton to edge); child 3 =
	 *         LINES PShape of bones (the pure straight skeleton). For
	 *         multi-polygons, a master GROUP shape of skeleton GROUP shapes
	 *         (described above) is returned.
	 */
	public static PShape straightSkeleton(PShape shape) {
		final Geometry g = fromPShape(shape);
		if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOLYGON)) {
			PShape group = new PShape(PConstants.GROUP);
			GeometryIterator gi = new GeometryIterator(g);
			gi.forEach(p -> group.addChild(straightSkeleton((Polygon) p)));
			return group;
		} else if (g.getGeometryType().equals(Geometry.TYPENAME_POLYGON)) {
			return straightSkeleton((Polygon) g);
		}
		return shape;
	}

	/**
	 * 
	 * @param polygon a single polygon that can contain holes
	 * @return
	 */
	private static PShape straightSkeleton(Polygon polygon) {
		/*
		 * Kenzi implementation (since PGS 1.2.1) is much faster (~50x!) but can fail on
		 * more complicated inputs. Therefore try Kenzi implementation first, but fall
		 * back to Twak implementation if it fails.
		 */
		try {
			return straightSkeletonKendzi(polygon);
		} catch (Exception e) {
			return straightSkeletonTwak(polygon);
		}
	}

	private static PShape straightSkeletonTwak(Polygon polygon) {
		if (polygon.getCoordinates().length > 1000) {
			polygon = (Polygon) DouglasPeuckerSimplifier.simplify(polygon, 2);
		}

		final Set<Coordinate> edgeCoordsSet = new HashSet<>();
		final Skeleton skeleton;
		final LoopL<Edge> loops = new LoopL<>(); // list of loops
		final Machine speed = new Machine(1); // every edge same speed

		final LinearRing[] rings = new LinearRingIterator(polygon).getLinearRings();
		for (int i = 0; i < rings.length; i++) {
			loops.add(ringToLoop(rings[i], i > 0, edgeCoordsSet, speed));
		}

		final PShape lines = new PShape(PConstants.GROUP);
		final PShape faces = new PShape(PConstants.GROUP);
		/*
		 * Create PEdges first to prevent lines being duplicated in output shapes since
		 * faces share branches and bones.
		 */
		final Set<PEdge> branchEdges = new HashSet<>();
		final Set<PEdge> boneEdges = new HashSet<>();
		try {
			skeleton = new Skeleton(loops, true);
			skeleton.skeleton(); // compute skeleton

			skeleton.output.faces.values().forEach(f -> {
				final List<Point3d> vertices = f.getLoopL().iterator().next().asList();
				List<PVector> faceVertices = new ArrayList<>();

				for (int i = 0; i < vertices.size(); i++) {
					final Point3d p1 = vertices.get(i);
					final Point3d p2 = vertices.get((i + 1) % vertices.size());
					faceVertices.add(new PVector((float) p1.x, (float) p1.y));
					final boolean a = edgeCoordsSet.contains(new Coordinate(p1.x, p1.y)); // NOTE Coordinate()
					final boolean b = edgeCoordsSet.contains(new Coordinate(p2.x, p2.y));
					if (a ^ b) { // branch (xor)
						branchEdges.add(new PEdge(p1.x, p1.y, p2.x, p2.y));
					} else {
						if (!a) { // bone
							boneEdges.add(new PEdge(p1.x, p1.y, p2.x, p2.y));
						}
					}
				}

				PShape face = PGS_Conversion.fromPVector(faceVertices);
				face.setStroke(true);
				face.setStrokeWeight(2);
				face.setStroke(RGB.composeColor(147, 112, 219));
				faces.addChild(face);
			});
		} catch (Exception ignore) {
			// hide init or collision errors from console
		}

		final PShape bones = prepareLinesPShape(null, null, 4);
		boneEdges.forEach(e -> {
			bones.vertex(e.a.x, e.a.y);
			bones.vertex(e.b.x, e.b.y);
		});
		bones.endShape();

		final PShape branches = prepareLinesPShape(RGB.composeColor(40, 235, 180), null, null);
		branchEdges.forEach(e -> {
			branches.vertex(e.a.x, e.a.y);
			branches.vertex(e.b.x, e.b.y);
		});
		branches.endShape();

		lines.addChild(faces);
		lines.addChild(branches);
		lines.addChild(bones);

		return lines;
	}

	private static PShape straightSkeletonKendzi(Polygon polygon) {
		final LinearRing[] rings = new LinearRingIterator(polygon).getLinearRings();
		Set<Vector2dc> edgeCoordsSet = new HashSet<>();
		final List<Vector2dc> points = ringToVec(rings[0], edgeCoordsSet);
		final List<List<Vector2dc>> holes = new ArrayList<>();
		for (int i = 1; i < rings.length; i++) {
			holes.add(ringToVec(rings[i], edgeCoordsSet));
		}

		final SkeletonOutput so = kendzi.math.geometry.skeleton.Skeleton.skeleton(points, holes, new SkeletonConfiguration());
		final PShape skeleton = new PShape(PConstants.GROUP);
		final PShape faces = new PShape(PConstants.GROUP);
		/*
		 * Create PEdges first to prevent lines being duplicated in output shapes since
		 * faces share branches and bones.
		 */
		final Set<PEdge> branchEdges = new HashSet<>();
		final Set<PEdge> boneEdges = new HashSet<>();
		so.getFaces().forEach(f -> {
			/*
			 * q stores the index of second vertex of the face that is a shape vertex. This
			 * is used to rotate f.getPoints() so that the vertices of every face PShape
			 * begin at the shape edge.
			 */
			int q = 0;
			for (int i = 0; i < f.getPoints().size(); i++) {
				final Vector2dc p1 = f.getPoints().get(i);
				final Vector2dc p2 = f.getPoints().get((i + 1) % f.getPoints().size());
				final boolean a = edgeCoordsSet.contains(p1);
				final boolean b = edgeCoordsSet.contains(p2);
				if (a ^ b) { // branch (xor)
					branchEdges.add(new PEdge(p1.x(), p1.y(), p2.x(), p2.y()));
					q = i;
				} else {
					if (!a) { // bone
						boneEdges.add(new PEdge(p1.x(), p1.y(), p2.x(), p2.y()));
					} else {
						q = i;
					}
				}
			}

			List<PVector> faceVertices = new ArrayList<>(f.getPoints().size());
			Collections.rotate(f.getPoints(), -q + 1);
			f.getPoints().forEach(p -> faceVertices.add(new PVector((float) p.x(), (float) p.y())));

			PShape face = PGS_Conversion.fromPVector(faceVertices);
			face.setStroke(true);
			face.setStrokeWeight(2);
			face.setStroke(RGB.composeColor(147, 112, 219));
			faces.addChild(face);
		});

		final PShape bones = prepareLinesPShape(null, null, 4);
		boneEdges.forEach(e -> {
			bones.vertex(e.a.x, e.a.y);
			bones.vertex(e.b.x, e.b.y);
		});
		bones.endShape();

		final PShape branches = prepareLinesPShape(RGB.composeColor(40, 235, 180), null, null);
		branchEdges.forEach(e -> {
			branches.vertex(e.a.x, e.a.y);
			branches.vertex(e.b.x, e.b.y);
		});
		branches.endShape();

		skeleton.addChild(faces);
		skeleton.addChild(branches);
		skeleton.addChild(bones);

		return skeleton;
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
			/*
			 * Need to use intersection() rather than checkling whether vertices are
			 * contained within the shape (faster) because vertices of longer (straight)
			 * line segments may lie within the shape when the segment extends outside the
			 * shape
			 */
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
		return isolines(points, intervalValueSpacing, isolineMin, isolineMax, 0);
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
	 * @param smoothing            Number of contour smoothing passes to perform.
	 *                             The best choice for this value depends on the
	 *                             requirements of the application. Values in the
	 *                             range 5 to 40 are good candidates for
	 *                             investigation.
	 * @return a map of {isoline -> height of the isoline}
	 */
	public static Map<PShape, Float> isolines(Collection<PVector> points, double intervalValueSpacing, double isolineMin, double isolineMax,
			int smoothing) {
		final IncrementalTin tin = new IncrementalTin(10);
		points.forEach(point -> tin.add(new Vertex(point.x, point.y, point.z)));

		double[] intervals = generateDoubleSequence(isolineMin, isolineMax, intervalValueSpacing);

		SmoothingFilter filter = null;
		if (smoothing > 0) {
			filter = new SmoothingFilter(tin, smoothing);
		}
		final ContourBuilderForTin builder = new ContourBuilderForTin(tin, filter, intervals, false);
		List<Contour> contours = builder.getContours();

		Map<PShape, Float> isolines = new HashMap<>(contours.size());

		for (Contour contourLine : contours) {
			final double[] coords = contourLine.getXY(); // [x1, y1, x2, y2, ...]
			final PShape isoline = new PShape();
			isoline.setFamily(PShape.PATH);
			isoline.setStroke(true);
			isoline.setStrokeWeight(2);
			isoline.setStroke(RGB.PINK);

			isoline.beginShape();
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
	 * @see #offsetCurvesOutward(PShape, OffsetStyle, double, int)
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
	 * @see #offsetCurvesInward(PShape, OffsetStyle, double)
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

	private static Loop<Edge> ringToLoop(LinearRing ring, boolean hole, Set<Coordinate> edgeCoordsSet, Machine speed) {
		Coordinate[] coords = ring.getCoordinates();
		if (!hole && !Orientation.isCCW(coords)) {
			reverse(coords); // exterior should be CCW
		}
		if (hole && Orientation.isCCW(coords)) {
			reverse(coords); // holes should be CW
		}

		List<Corner> corners = new ArrayList<>();
		Loop<Edge> loop = new Loop<>();

		for (int j = 0; j < coords.length; j++) {
			corners.add(new Corner(coords[j].x, coords[j].y));
			edgeCoordsSet.add(coords[j]);
		}

		for (int j = 0; j < corners.size() - 1; j++) {
			Edge edge = new Edge(corners.get(j), corners.get((j + 1) % (corners.size() - 1)));
			edge.machine = speed;
			loop.append(edge);
		}

		return loop;
	}

	private static List<Vector2dc> ringToVec(LinearRing ring, Set<Vector2dc> edgeCoordsSet) {
		final List<Vector2dc> points = new ArrayList<>();
		Coordinate[] coords = ring.getCoordinates();
		/*
		 * Kendzi polygons are unclosed (cannot start and end with the same point),
		 * unlike a LinearRing.
		 */
		for (int i = 0; i < coords.length - 1; i++) { // note - 1
			final Vector2dc p = new Vector2d(coords[i].x, coords[i].y);
			points.add(p);
			edgeCoordsSet.add(p);
		}
		return points;
	}

}
