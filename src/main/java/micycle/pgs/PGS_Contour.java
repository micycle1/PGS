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
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.vecmath.Point3d;

import org.jgrapht.alg.interfaces.ShortestPathAlgorithm;
import org.jgrapht.alg.shortestpath.BFSShortestPath;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.joml.Vector2d;
import org.joml.Vector2dc;
import org.locationtech.jts.algorithm.Angle;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.buffer.OffsetCurve;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.SimpleTriangle;
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

import com.google.common.collect.Lists;

import kendzi.math.geometry.skeleton.SkeletonConfiguration;
import kendzi.math.geometry.skeleton.SkeletonOutput;
import micycle.medialAxis.MedialAxis;
import micycle.medialAxis.MedialAxis.MedialDisk;
import micycle.pgs.PGS.GeometryIterator;
import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.color.ColorUtils;
import micycle.pgs.color.Colors;
import micycle.pgs.commons.PEdge;
import net.jafama.FastMath;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Methods for producing different kinds of shape contours.
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

	/*
	 * TODO implement 'Base Point Split Algorithm for Generating Polygon Skeleton
	 * Lines'
	 */

	private PGS_Contour() {
	}

	/**
	 * Computes the medial axis of the given shape, providing a characterization of
	 * the skeleton of a shape.
	 * <p>
	 * The method employs three parameters to prune (simplify) the resulting medial
	 * axis according to distinct features.
	 *
	 * @param shape             The PShape to calculate the medial axis for.
	 * @param axialThreshold    A value in the range of 0 to 1 that determines the
	 *                          level of pruning based on the axial gradient. The
	 *                          axial gradient evaluates the variation in the
	 *                          shape's width per unit length along the axis,
	 *                          measured per edge segment. A value of 0 results in
	 *                          no pruning, while a value of 1 leads to the maximum
	 *                          possible pruning.
	 * @param distanceThreshold A value between 0 and 1 that determines the level of
	 *                          pruning based on the spatial distance from the root
	 *                          of the medial axis to the tail coordinate of each
	 *                          edge. A value of 0 results in no pruning, while a
	 *                          value of 1 results in maximum possible pruning.
	 * @param areaThreshold     A value between 0 and 1 that determines the level of
	 *                          pruning based on the aggregate feature area of each
	 *                          edge and its descendants in the medial axis. A value
	 *                          of 0 results in no pruning, while a value of 1
	 *                          results in the maximum possible pruning.
	 * @return A GROUP PShape containing maximal-length lines representing the
	 *         pruned edges of the medial axis.
	 */
	public static PShape medialAxis(PShape shape, double axialThreshold, double distanceThreshold, double areaThreshold) {
		final Geometry g = fromPShape(shape);
		final MedialAxis m = new MedialAxis(g);
		return PGS_SegmentSet.dissolve(m.getPrunedEdges(axialThreshold, distanceThreshold, areaThreshold).stream()
				.map(e -> new PEdge(e.head.position.x, e.head.position.y, e.tail.position.x, e.tail.position.y))
				.collect(Collectors.toList()));
	}

	/**
	 * Computes the chordal axis of a shape, which provides a characterization of
	 * the skeleton of a shape.
	 * <p>
	 * In its primitive form, the chordal axis is constructed by joining the
	 * midpoints of the chords and the centroids of junction and terminal triangles
	 * of the delaunay trianglution of a shape.
	 * <p>
	 * It can be considered a more useful alternative to the medial axis for
	 * obtaining skeletons of discrete shapes.
	 * 
	 * @param shape polygonal shape
	 * @return a GROUP PShape, where each group is a single maximum-length line
	 *         segment (possibly >2 vertices)
	 * @since 1.3.0
	 */
	@SuppressWarnings("unchecked")
	public static PShape chordalAxis(PShape shape) {
		/*-
		 * See 'Rectification of the Chordal Axis Transform and a New Criterion for
		 * Shape Decomposition' for CAT extension.
		 * See 'Morphological Analysis of Shapes' for CAT pruning techniques.
		 * See 'Shape Matching By Part Alignment Using Extended Chordal Axis Transform'.
		 */
		final IIncrementalTin triangulation = PGS_Triangulation.delaunayTriangulationMesh(shape);
		final SimpleGraph<SimpleTriangle, DefaultEdge> graph = PGS_Triangulation.toDualGraph(triangulation);

		final List<PEdge> edges = new ArrayList<>(graph.vertexSet().size());

		for (SimpleTriangle t : graph.vertexSet()) {
			/*
			 * For each triangle, determine how many edges it has in common with the shape
			 * boundary, and use this number [1,2,3] to classify it. Below, triangles are
			 * classified based on the number of neighbor triangles.
			 */
			switch (graph.outDegreeOf(t)) {
				case 1 : // Terminal triangle (2 edges in perimeter)
					final IQuadEdge interiorEdge; // one edge is interior
					if (t.getEdgeA().isConstrainedRegionBorder()) {
						if (t.getEdgeB().isConstrainedRegionBorder()) {
							interiorEdge = t.getEdgeC();
						} else {
							interiorEdge = t.getEdgeB();
						}
					} else {
						interiorEdge = t.getEdgeA();
					}
					PVector centroid = centroid(t);
					edges.add(new PEdge(centroid.x, centroid.y, midpoint(interiorEdge).x, midpoint(interiorEdge).y));
					break;
				case 2 : // Sleeve triangle (one edge in perimeter)
					final IQuadEdge interiorEdgeA; // 2 edges are interior
					final IQuadEdge interiorEdgeB;
					if (t.getEdgeA().isConstrainedRegionBorder()) {
						interiorEdgeA = t.getEdgeB();
						interiorEdgeB = t.getEdgeC();
					} else if (t.getEdgeB().isConstrainedRegionBorder()) {
						interiorEdgeA = t.getEdgeA();
						interiorEdgeB = t.getEdgeC();
					} else {
						interiorEdgeA = t.getEdgeA();
						interiorEdgeB = t.getEdgeB();
					}
					PVector midpoint1 = midpoint(interiorEdgeA);
					PVector midpoint2 = midpoint(interiorEdgeB);
					edges.add(new PEdge(midpoint1.x, midpoint1.y, midpoint2.x, midpoint2.y));
					break;
				case 3 : // Junction triangle (no edge in perimeter)
					/*
					 * Connect the midpoints of the two shortest sides to the midpoint of the
					 * longest side.
					 */
					double maxLength = t.getEdgeA().getLength();
					IQuadEdge longest = t.getEdgeA();
					IQuadEdge shortA = t.getEdgeB(), shortB = t.getEdgeC();
					if (t.getEdgeB().getLength() > maxLength) {
						maxLength = t.getEdgeB().getLength();
						shortA = t.getEdgeA();
						shortB = t.getEdgeC();
						longest = t.getEdgeB();
					}
					if (t.getEdgeC().getLength() > maxLength) {
						shortA = t.getEdgeA();
						shortB = t.getEdgeB();
						longest = t.getEdgeC();
					}
					final PVector midpointL = midpoint(longest);
					final PVector midpointA = midpoint(shortA);
					final PVector midpointB = midpoint(shortB);
					edges.add(new PEdge(midpointA.x, midpointA.y, midpointL.x, midpointL.y)); // A<->L
					edges.add(new PEdge(midpointB.x, midpointB.y, midpointL.x, midpointL.y)); // B<->L
					break;
				default :
					break;
			}
		}

		return PGS_SegmentSet.dissolve(edges);
	}

	/**
	 * Computes the straight skeleton for a shape.
	 * <p>
	 * A straight skeleton is a skeletal structure similar to the medial axis,
	 * consisting of straight-line segments only. Roughly, it is the geometric graph
	 * whose edges are the traces of vertices of shrinking mitered offset curves of
	 * the polygon.
	 * <p>
	 * For a single polygon, this method returns a GROUP PShape containing three
	 * children:
	 * <ul>
	 * <li>Child 0: GROUP PShape consisting of skeleton faces.</li>
	 * <li>Child 1: LINES PShape representing branches, which are lines connecting
	 * the skeleton to the polygon's edge.</li>
	 * <li>Child 2: LINES PShape composed of bones, depicting the pure straight
	 * skeleton of the polygon.</li>
	 * </ul>
	 * <p>
	 * For multi-polygons, the method returns a master GROUP PShape. This master
	 * shape includes multiple skeleton GROUP shapes, each corresponding to a single
	 * polygon and structured as described above.
	 * 
	 * @param shape a single polygon (that can contain holes), or a multi polygon
	 *              (whose polygons can contain holes)
	 * 
	 * @return PShape based on the input polygon structure, either as a single or
	 *         multi-polygon skeleton representation.
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

	private static PShape straightSkeleton(Polygon polygon) {
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
				List<Point3d> vertices = f.getLoopL().iterator().next().stream().toList();
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
				face.setStroke(ColorUtils.composeColor(147, 112, 219));
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

		final PShape branches = prepareLinesPShape(ColorUtils.composeColor(40, 235, 180), null, null);
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

	/**
	 * Generates a topographic-like isoline contour map from the shape's vertices
	 * and a given "high point". Isolines represent the "elevation", or euclidean
	 * distance, between a location in the shape and the "high point".
	 * <p>
	 * Assigns each point feature a number equal to the distance between geometry's
	 * centroid and the point.
	 *
	 * @param shape           the bounds in which to draw isolines
	 * @param highPoint       position of "high" point within the shape
	 * @param intervalSpacing distance between successive isolines
	 * @return PShape containing isolines linework
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

		/*
		 * "A null valuator tells the builder to just use the z values from the vertices
		 * rather than applying any adjustments to their values."
		 */
		final ContourBuilderForTin builder = new ContourBuilderForTin(tin, null, intervals, true);

		List<Contour> contours = builder.getContours();

		PShape parent = new PShape(PConstants.GROUP);
		parent.setKind(PConstants.GROUP);

		LineDissolver ld = new LineDissolver();
		for (Contour contour : contours) {
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
			out.setStrokeCap(PConstants.SQUARE);
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
			isoline.setStroke(Colors.PINK);

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
	 * Generates a contour map based on a distance field of a shape.
	 * <p>
	 * A distance field maps each point within the shape to the shortest distance
	 * between that point and the shape boundary.
	 * 
	 * @param shape   polygonal shape
	 * @param spacing distance represented by successive contour lines
	 * @return GROUP shape, where each child is a closed contour line or contour
	 *         line partition
	 * @since 1.3.0
	 */
	public static PShape distanceField(PShape shape, double spacing) {
		Geometry g = fromPShape(shape);
		MedialAxis m = new MedialAxis(g);

		List<PVector> disks = new ArrayList<>();
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		for (MedialDisk d : m.getDisks()) {
			disks.add(new PVector((float) d.position.x, (float) d.position.y, (float) d.distance));
			min = Math.min(d.distance, min);
			max = Math.max(d.distance, max);
		}

		PShape out = PGS_Conversion.flatten(PGS_Contour.isolines(disks, spacing, min, max).keySet());
		PShape i = PGS_ShapeBoolean.intersect(shape, out);
		PGS_Conversion.disableAllFill(i); // since some shapes may be polygons
		return i;
	}

	/**
	 * Generates a tree structure representing the shortest paths from a given start
	 * point to all other vertices in the provided mesh. The paths are computed
	 * using the existing connectivity of the mesh edges, ensuring that the
	 * shortest-path tree respects the original mesh structure. The tree is
	 * constructed using a Breadth-First Search (BFS) algorithm.
	 * <p>
	 * The shortest-path tree represents the minimal set of mesh edges required to
	 * connect the start point to all other vertices in the mesh, following the
	 * mesh's inherent connectivity. This ensures that the paths are constrained by
	 * the mesh's topology rather than creating arbitrary connections between
	 * vertices.
	 * <p>
	 * If the provided start point does not exactly match a vertex in the mesh, the
	 * closest vertex in the mesh to the start point is used as the actual starting
	 * point for the shortest-path computation.
	 *
	 * @param mesh    A GROUP shape representing a mesh from which the graph is
	 *                constructed. The mesh defines the connectivity between
	 *                vertices via its edges.
	 * @param source  The starting point from which the shortest paths are
	 *                calculated. If this point does not exactly match a vertex in
	 *                the mesh, the closest vertex in the mesh will be used as the
	 *                starting point.
	 * @param flatten Determines the format of the output shortest-path tree.
	 *                <p>
	 *                If {@code true}, the method returns a flattened representation
	 *                of the shortest-path tree as a single set of edges. This
	 *                removes duplicate edges and combines all paths into a single
	 *                structure.
	 *                <p>
	 *                If {@code false}, the method returns a GROUP shape of
	 *                individual paths, where each path is a separate line from the
	 *                start point to each vertex in the mesh. This representation
	 *                retains the structure of the shortest-path tree as a
	 *                collection of distinct paths.
	 *                <p>
	 * @return A PShape object representing the tree of shortest paths from the
	 *         start point to all other vertices in the mesh. The paths are
	 *         constrained by the mesh's edge connectivity.
	 * @since 2.1
	 */
	public static PShape distanceTree(PShape mesh, PVector source, boolean flatten) {
		var g = PGS_Conversion.toGraph(mesh);
		ShortestPathAlgorithm<PVector, PEdge> spa = new BFSShortestPath<>(g);

		final PVector sourceActual = PGS_Optimisation.closestPoint(g.vertexSet(), source);
		var paths = spa.getPaths(sourceActual);

		PShape out;
		if (flatten) {
			var edges = g.vertexSet().stream().filter(v -> !v.equals(sourceActual)) // Exclude the source vertex
					.flatMap(v -> paths.getPath(v).getEdgeList().stream()) // Flatten the edge lists into a single stream
					.collect(Collectors.toSet()); // Collect the edges into a Set to remove duplicates
			out = PGS_SegmentSet.toPShape(edges);
		}
		else {
			var pathLines = g.vertexSet().stream() //
					.filter(v -> v != sourceActual) // Exclude the source vertex
					.map(v -> PGS_Conversion.fromPVector(paths.getPath(v).getVertexList())).toList();
			out = PGS_Conversion.flatten(pathLines);
			PGS_Conversion.setAllStrokeColor(out, ColorUtils.setAlpha(Colors.PINK, 50), 4);			
		}

		return out;
	}

	/**
	 * Calculates the longest center line passing through a given shape (using
	 * default straightness weighting and smoothing parameters).
	 * <p>
	 * The center line is determined based on the medial axis of the shape; it
	 * endpoints are always leaf vertices of the medial axis, and the line will pass
	 * through the center coordinate of the shape's largest inscribed circle.
	 * 
	 * @param shape The non-GROUP PShape representing the input shape for which the
	 *              longest center line is to be calculated.
	 * @return A new PShape representing the smoothed longest center line of the
	 *         input shape.
	 * @see #centerLine(PShape, double, double)
	 * @since 1.4.0
	 */
	public static PShape centerLine(PShape shape) {
		return centerLine(shape, 0.7, 50);
	}

	/**
	 * Calculates the longest center line passing through a given shape. Note that
	 * shapes with sparse vertices will need densified beforehand.
	 * <p>
	 * The center line is determined based on the medial axis of the shape; it
	 * endpoints are always leaf vertices of the medial axis, and the line will pass
	 * through the center coordinate of the shape's largest inscribed circle.
	 * <p>
	 * The method can promote paths that are straighter, even if they are somewhat
	 * shorter than the longest possible line, by weighting them according to
	 * <code>straightnessWeighting</code>.
	 * <p>
	 * This method can be used to find the placement position for overlayed curved
	 * text label.
	 * 
	 * @param shape                 The non-GROUP PShape representing the input
	 *                              shape for which the longest center line is to be
	 *                              calculated.
	 * @param straightnessWeighting A value in [0...1] to determine how straighter
	 *                              paths should be weighted (preferred over the
	 *                              longest possible center line). 0 is no
	 *                              additional weighting - the longest path is
	 *                              chosen despite how concave it is - and 1 is
	 *                              maximum weighting. A good starting value to use
	 *                              is ~0.7.
	 * @param smoothing             Gaussian smoothing parameter. A good starting
	 *                              value to use is ~50.
	 * @return A new PShape representing the smoothed longest center line of the
	 *         input shape.
	 * @since 1.4.0
	 */
	public static PShape centerLine(PShape shape, double straightnessWeighting, double smoothing) {
		/*
		 * For general polygons, the medial root always trifurcates into 3 distinct
		 * paths/sub-trees. I assume the longest center line passes through the root, so
		 * we know that one vertex (origin) of the longest line belongs to the one of
		 * the paths and the other vertex (terminus) belongs to one of the two other
		 * paths. So we look at the euclidean distance between pairs of leaf vertices,
		 * where each leaf in the pair belongs to a distinct path, and select the pair
		 * that maximises this distance to be the longest center line. One additional
		 * factor is we multiply the euclidean distance by the angle the line makes with
		 * the root node in order to promote paths that are straighter, even if they are
		 * somewhat shorter.
		 */
		MedialAxis m = new MedialAxis(fromPShape(shape));

		List<micycle.medialAxis.MedialAxis.Edge> longestPath = new ArrayList<>();
		List<MedialDisk> subTree1 = m.getDescendants(m.getRoot().children.get(0)).stream().filter(d -> d.degree == 0)
				.collect(Collectors.toList());
		List<MedialDisk> subTree2 = m.getDescendants(m.getRoot().children.get(1)).stream().filter(d -> d.degree == 0)
				.collect(Collectors.toList());
		if (m.getRoot().children.size() == 2) {
			// special case of elliptical (etc.) shapes
			longestPath = new ArrayList<>(m.getEdges());
		} else {
			List<MedialDisk> subTree3 = m.getDescendants(m.getRoot().children.get(2)).stream().filter(d -> d.degree == 0)
					.collect(Collectors.toList());

			MedialDisk longestPathD1 = null;
			MedialDisk longestPathD2 = null;

			final List<List<MedialDisk>> diskPairs = new ArrayList<>();
			diskPairs.addAll(Lists.cartesianProduct(subTree1, subTree2));
			diskPairs.addAll(Lists.cartesianProduct(subTree1, subTree3));
			diskPairs.addAll(Lists.cartesianProduct(subTree2, subTree3));

			double maxWeight = Double.NEGATIVE_INFINITY;
			for (List<MedialDisk> diskPair : diskPairs) {
				final MedialDisk d1 = diskPair.get(0);
				final MedialDisk d2 = diskPair.get(1);
				double angle = Angle.angleBetween(d1.position, m.getRoot().position, d2.position);
				// could also use euclid distanc between d1 and d2
				angle = FastMath.pow(1 + angle, straightnessWeighting);
				double pathweight = (d1.distance + d2.distance) * Math.max((angle - 1), 1);
				if (pathweight > maxWeight) {
					maxWeight = pathweight;
					longestPathD1 = d1;
					longestPathD2 = d2;
				}
			}

			longestPath.addAll(m.getEdgesToRoot(longestPathD1));
			longestPath.addAll(m.getEdgesToRoot(longestPathD2));
		}

		List<LineString> strings = longestPath.stream().map(e -> e.lineString).collect(Collectors.toList());
		MultiLineString stringsGeometry = GEOM_FACTORY.createMultiLineString(GeometryFactory.toLineStringArray(strings));
		PShape longestPathShape = toPShape(LineDissolver.dissolve(stringsGeometry));
		longestPathShape = PGS_Morphology.simplify(longestPathShape, 1);
		longestPathShape = PGS_Morphology.smoothGaussian(longestPathShape, smoothing);

		return longestPathShape;
	}

	/**
	 * Specifies the join style for offset curves.
	 */
	public enum OffsetStyle {

		MITER(BufferParameters.JOIN_MITRE), BEVEL(BufferParameters.JOIN_BEVEL), ROUND(BufferParameters.JOIN_ROUND);

		final int style;

		private OffsetStyle(int style) {
			this.style = style;
		}
	}

	/**
	 * Generates inward-facing offset curves from a shape. Curves are generated
	 * until they collapse.
	 * 
	 * @param shape   a path, polygon or multipolygon (GROUP) shape
	 * @param style   the type of curve join (BEVEL, MITER, or ROUND)
	 * @param spacing the distance between each curve, must be >=1
	 * @return a GROUP PShape where each child is a single curve or a group of
	 *         curves created at the same step
	 * @see {@link #offsetCurvesOutward(PShape, OffsetStyle, double, int)
	 *      offsetCurvesOutward()} for outward-facing curves.
	 */
	public static PShape offsetCurvesInward(PShape shape, OffsetStyle style, double spacing) {
		return offsetCurves(shape, style, spacing, 0, false);
	}

	/**
	 * Generates N inward-facing offset curves from a shape.
	 * 
	 * @param shape   a path, polygon or multipolygon (GROUP) shape
	 * @param style   the type of curve join (BEVEL, MITER, or ROUND)
	 * @param spacing the distance between each curve, must be >=1
	 * @param curves  the number of curves to generate (including the original shape
	 *                outline)
	 * @return a GROUP PShape where each child is a single curve or a group of
	 *         curves created at the same step
	 * @see {@link #offsetCurvesOutward(PShape, OffsetStyle, double, int)
	 *      offsetCurvesOutward()} for outward-facing curves.
	 */
	public static PShape offsetCurvesInward(PShape shape, OffsetStyle style, double spacing, int curves) {
		return offsetCurves(shape, style, spacing, curves, false);
	}

	/**
	 * Generates N outward-facing offset curves from a shape.
	 *
	 * @param shape   a path, polygon or multipolygon (GROUP) shape
	 * @param style   the type of curve join (BEVEL, MITER, or ROUND)
	 * @param spacing the distance between each curve, must be >=1
	 * @param curves  the number of curves to generate (including the original shape
	 *                outline)
	 * @return a GROUP PShape where each child is a single curve or a group of
	 *         curves created at the same step
	 * @see #offsetCurvesInward(PShape, OffsetStyle, double)
	 */
	public static PShape offsetCurvesOutward(PShape shape, OffsetStyle style, double spacing, final int curves) {
		return offsetCurves(shape, style, spacing, curves, true);
	}

	/**
	 * Generic method for producing offset curves from shapes.
	 */
	private static PShape offsetCurves(PShape shape, OffsetStyle style, double spacing, final int curves, boolean outwards) {
		Geometry g = fromPShape(shape);

		// handle non-polygonal / path shapes
		if (g.getGeometryType().equals(Geometry.TYPENAME_LINESTRING)) {
			List<Geometry> strings = new ArrayList<>(curves);
			for (int i = 0; i < curves; i++) {
				strings.add(
						OffsetCurve.getCurve(g, spacing * (outwards ? 1 : -1) * i, 8, style.style, BufferParameters.DEFAULT_MITRE_LIMIT));
			}
			return toPShape(strings);
		}

		if (g.getCoordinates().length > 2000) {
			g = DouglasPeuckerSimplifier.simplify(g, 0.25);
		}

		final BufferParameters bufParams = new BufferParameters(8, BufferParameters.CAP_FLAT, style.style,
				BufferParameters.DEFAULT_MITRE_LIMIT);
//		bufParams.setSimplifyFactor(5); // can produce "poor" yet interesting results

		spacing = Math.max(1, Math.abs(spacing)); // ensure positive and >=1
		spacing = outwards ? spacing : -spacing;

		final PShape parent = new PShape(PConstants.GROUP);
		int currentCurves = 0;
		while ((outwards && currentCurves < curves) || (!outwards && !g.isEmpty() && curves == 0)
				|| (!outwards && currentCurves < curves)) {
			LinearRing[] rings = new LinearRingIterator(g).getLinearRings();
			if (rings.length == 1) {
				PShape curve = toPShape(rings[0]);
				curve.setFill(false);
				curve.setStrokeWeight(2);
				parent.addChild(curve);
			} else if (rings.length > 1) { // curves created at the same step are grouped together
				PShape curveParent = new PShape(PConstants.GROUP);
				for (LinearRing ring : rings) {
					PShape curve = toPShape(ring);
					curve.setFill(false);
					curve.setStrokeWeight(2);
					curveParent.addChild(curve);
				}
				parent.addChild(curveParent);
			}

			BufferOp b = new BufferOp(g, bufParams);
			g = b.getResultGeometry(spacing);
			if (style != OffsetStyle.MITER) {
				// simplify because rounded buffers produce LOADS of dense vertices
				g = DouglasPeuckerSimplifier.simplify(g, 0.05);
			}
			currentCurves++;
		}

		return parent;
	}

	private static PVector midpoint(IQuadEdge e) {
		return new PVector((float) (e.getA().x + e.getB().x) / 2, (float) (e.getA().y + e.getB().y) / 2);
	}

	private static PVector centroid(final SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new PVector((float) x, (float) y);
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
	 * Generates a grid of points.
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

		for (double element : x) {
			for (double element2 : y) {
				grid.add(new PVector((float) element, (float) element2));
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

		for (Coordinate coord : coords) {
			corners.add(new Corner(coord.x, coord.y));
			edgeCoordsSet.add(coord);
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
