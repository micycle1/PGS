package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static processing.core.PConstants.TRIANGLES;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.function.Consumer;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;

import earcut4j.Earcut;
import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.color.RGB;
import micycle.pgs.utility.PoissonDistribution;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Delaunay and earcut triangulation of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Triangulation {

	private PGS_Triangulation() {
	}

	/**
	 * Generates a Delaunay Triangulation from the given shape. The triangulation
	 * can be both constrained (meaning the triangulation is masked by the original
	 * shape) and refined (meaning additional points are inserted, usually leading
	 * to more uniform triangle shapes and sizes).
	 * 
	 * @param shape         the shape whose vertices to generate a triangulation of
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to maintain the Delaunay nature when
	 *                      constraining the triangulation, and whether to check
	 *                      that centroid locations lie within the shape during
	 *                      refinement. When pretty=true, triangles in the
	 *                      triangulation may be slightly more regular in
	 *                      shape/size. There is a small performance overhead which
	 *                      becomes more considerable at higher refinement levels.
	 *                      When constrain=false and refinements=0, this argument
	 *                      has no effect.
	 * @return a TRIANGLES PShape
	 * @see #delaunayTriangulationPoints(PShape, List, boolean, int, boolean)
	 * @see #delaunayTriangulationMesh(PShape, List, boolean, int, boolean)
	 */
	public static PShape delaunayTriangulation(PShape shape, List<PVector> steinerPoints, boolean constrain, int refinements,
			boolean pretty) {
		IncrementalTin tin = delaunayTriangulationMesh(shape, steinerPoints, constrain, refinements, pretty);

		PShape triangulation = new PShape(PShape.GEOMETRY);

		PGS_Conversion.setAllFillColor(triangulation, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(triangulation, RGB.PINK, 2);
		triangulation.beginShape(TRIANGLES);

		Consumer<Vertex[]> triangleVertexConsumer = t -> {
			triangulation.vertex((float) t[0].x, (float) t[0].y);
			triangulation.vertex((float) t[1].x, (float) t[1].y);
			triangulation.vertex((float) t[2].x, (float) t[2].y);
		};
		if (constrain) {
			TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		} else {
			TriangleCollector.visitTriangles(tin, triangleVertexConsumer);
		}

		triangulation.endShape();
		return triangulation;
	}

	/**
	 * Generates a Delaunay Triangulation from the given shape. The triangulation
	 * can be both constrained (meaning the triangulation is masked by the original
	 * shape) and refined (meaning additional points are inserted, usually leading
	 * to more uniform triangle shapes and sizes).
	 * <p>
	 * This method returns the triangulation as a list of points, rather than a
	 * PShape.
	 * 
	 * @param shape         the shape whose vertices to generate a triangulation of
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to maintain the Delaunay nature when
	 *                      constraining the triangulation, and whether to check
	 *                      that centroid locations lie within the shape during
	 *                      refinement. When pretty=true, triangles in the
	 *                      triangulation may be slightly more regular in
	 *                      shape/size. There is a small performance overhead which
	 *                      becomes more considerable at higher refinement levels.
	 *                      When constrain=false and refinements=0, this argument
	 *                      has no effect.
	 * @return List of PVector coordinates, where each consecutive triplet of
	 *         coordinates are the 3 vertices belonging to one triangle
	 * @see #delaunayTriangulationPoints(PShape, List, boolean, int, boolean)
	 * @see #delaunayTriangulationMesh(PShape, List, boolean, int, boolean)
	 */
	public static List<PVector> delaunayTriangulationPoints(PShape shape, List<PVector> steinerPoints, boolean constrain, int refinements,
			boolean pretty) {
		IncrementalTin tin = delaunayTriangulationMesh(shape, steinerPoints, constrain, refinements, pretty);

		ArrayList<PVector> triangles = new ArrayList<>();
		Consumer<Vertex[]> triangleVertexConsumer = t -> {
			triangles.add(toPVector(t[0]));
			triangles.add(toPVector(t[1]));
			triangles.add(toPVector(t[2]));
		};
		if (constrain) {
			TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		} else {
			TriangleCollector.visitTriangles(tin, triangleVertexConsumer);
		}
		return triangles;
	}

	/**
	 * Generates a Delaunay Triangulation from the given shape. The triangulation
	 * can be both constrained (meaning the triangulation is masked by the original
	 * shape) and refined (meaning additional points are inserted, usually leading
	 * to more uniform triangle shapes and sizes).
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 * 
	 * @param shape         the shape whose vertices to generate a triangulation of
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement passes to
	 *                      perform. Each pass inserts the centroids of every
	 *                      existing triangle into the triangulation. Should be 0 or
	 *                      greater (probably no more than 5).
	 * @param pretty        Whether to maintain the Delaunay nature when
	 *                      constraining the triangulation, and whether to check
	 *                      that centroid locations lie within the shape during
	 *                      refinement. When pretty=true, triangles in the
	 *                      triangulation may be slightly more regular in
	 *                      shape/size. There is a small performance overhead which
	 *                      becomes more considerable at higher refinement levels.
	 *                      When constrain=false and refinements=0, this argument
	 *                      has no effect.
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulation(PShape, List, boolean, int, boolean)
	 * @see #delaunayTriangulationPoints(PShape, List, boolean, int, boolean)
	 */
	public static IncrementalTin delaunayTriangulationMesh(PShape shape, Collection<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		final Geometry g = fromPShape(shape);
		final IncrementalTin tin = new IncrementalTin(10);

		final ArrayList<Vertex> vertices = new ArrayList<>();
		Coordinate[] coords = g.getCoordinates();
		for (int i = 0; i < coords.length; i++) {
			vertices.add(new Vertex(coords[i].x, coords[i].y, 0));
		}

		tin.add(vertices, null); // initial triangulation

		if (steinerPoints != null) {
			steinerPoints.forEach(v -> tin.add(new Vertex(v.x, v.y, 0))); // add steiner points
		}

		if (refinements > 0) {

			final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
			final ArrayList<Vertex> refinementVertices = new ArrayList<>();

			/**
			 * A possible optimisation is to recursely split within each triangle upto the
			 * refinement depth (in one pass), so perform many less location checks. Another
			 * is to rasterise the PShape and check pixel[] array.
			 */
			for (int i = 0; i < refinements; i++) {
				refinementVertices.clear();
				TriangleCollector.visitSimpleTriangles(tin, t -> {
					if (t.getArea() > 50) { // don't refine small triangles
						final Coordinate center = centroid(t); // use centroid rather than circumcircle center
						if (pretty || pointLocator.locate(center) != Location.EXTERIOR) {
							refinementVertices.add(new Vertex(center.x, center.y, 0));
						}
					}
				});
				tin.add(refinementVertices, null); // add refinement (steiner) points
			}
		}

		if (constrain) {
			List<IConstraint> constraints = new ArrayList<>();
			for (int n = 0; n < g.getNumGeometries(); n++) {
				boolean hole = false;

				LinearRingIterator lri = new LinearRingIterator(g.getGeometryN(n));
				for (LinearRing ring : lri) {
					ArrayList<Vertex> points = new ArrayList<>();
					Coordinate[] c = ring.getCoordinates();
					if (Orientation.isCCW(c) && !hole) {
						for (int i = 0; i < c.length; i++) {
							points.add(new Vertex(c[i].x, c[i].y, 0));
						}
					} else {
						// holes are CW; region to keep lies to left the of constraints
						for (int i = c.length - 1; i >= 0; i--) {
							points.add(new Vertex(c[i].x, c[i].y, 0));
						}
					}
					constraints.add(new PolygonConstraint(points));
					hole = true; // all rings except the first are holes
				}
			}
			tin.addConstraints(constraints, pretty); // true/false is negligible?
		}

		return tin;
	}

	/**
	 * Creates a Delaunay triangulation of the shape where additional steiner
	 * points, populated by poisson sampling, are included.
	 * 
	 * @param shape
	 * @param spacing (Minimum) spacing between poisson points
	 * @return list of PVectors, where each successive triplet of PVectors
	 *         correspond to the 3 vertices of one triangle
	 */
	public static List<PVector> poissonTriangulation(PShape shape, double spacing) {

		final Geometry g = fromPShape(shape).buffer(-spacing);
		Envelope e = g.getEnvelopeInternal();

		PoissonDistribution pd = new PoissonDistribution(0);
		List<PVector> poissonPoints = pd.generate(e.getMinX(), e.getMinY(), e.getMinX() + e.getWidth(), e.getMinY() + e.getHeight(),
				spacing, 8);
//		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
//		List<PVector> pp = poissonPoints.parallelStream()
//				.filter(p -> pointLocator.locate(PGS.coordFromPVector(p)) != Location.EXTERIOR)
//				.collect(Collectors.toList());

		IncrementalTin tin = delaunayTriangulationMesh(shape, poissonPoints, true, 0, true);

		ArrayList<PVector> triangles = new ArrayList<>();
		Consumer<Vertex[]> triangleVertexConsumer = t -> {
			triangles.add(toPVector(t[0]));
			triangles.add(toPVector(t[1]));
			triangles.add(toPVector(t[2]));
		};
		TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		return triangles;
	}

	/**
	 * Computes a triangulation of the shape according to the ear clipping
	 * ("earcut") method. The triangulation is constrained to the shape by default.
	 * Does not support holes (for now...).
	 * 
	 * @param shape
	 * @return a TRIANGLES PShape
	 * @since 1.1.0
	 */
	public static PShape earCutTriangulation(PShape shape) {
		return earCutTriangulation(PGS_Conversion.toPVector(shape));
	}

	/**
	 * Computes a triangulation of the given points according to the ear clipping
	 * ("earcut") method.
	 * 
	 * @param points
	 * @return a TRIANGLES PShape
	 */
	public static PShape earCutTriangulation(List<PVector> points) {
		double[] arrCoords = new double[points.size() * 2];

		for (int i = 0; i < points.size(); i++) {
			arrCoords[2 * i] = points.get(i).x;
			arrCoords[2 * i + 1] = points.get(i).y;
		}

		List<Integer> triangles = Earcut.earcut(arrCoords, null, 2);

		PShape triangulation = new PShape();
		triangulation.setFamily(PShape.GEOMETRY);
		triangulation.setStroke(true);
		triangulation.setStrokeWeight(2);
		triangulation.setStroke(RGB.PINK);
		triangulation.setFill(true);
		triangulation.setFill(micycle.pgs.color.RGB.composeColor(255, 255, 255, 255));

		triangulation.beginShape(TRIANGLES);
		for (int i = 0; i < triangles.size(); i += 3) {
			final int v1 = 2 * triangles.get(i);
			final int v2 = 2 * triangles.get(i + 1);
			final int v3 = 2 * triangles.get(i + 2);
			triangulation.vertex((float) arrCoords[v1], (float) arrCoords[v1 + 1]);
			triangulation.vertex((float) arrCoords[v2], (float) arrCoords[v2 + 1]);
			triangulation.vertex((float) arrCoords[v3], (float) arrCoords[v3 + 1]);
		}
		triangulation.endShape();

		return triangulation;
	}

	/**
	 * Generates a shape consisting of polygonal faces of an Urquhart graph. An
	 * Urquhart graph is obtained by removing the longest edge from each triangle in
	 * a triangulation.
	 * <p>
	 * In practice this is a way to tessellate a shape into polygons (with the
	 * resulting tessellation being in between a
	 * {@link #delaunayTriangulation(PShape, List, boolean, int, boolean)
	 * triangulation} and a {@link micycle.pgs.PGS_Processing#partition(PShape)
	 * partition}).
	 * <p>
	 * Note that this method processes a Delaunay triangulation. Process a shape
	 * using
	 * {@link #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * delaunayTriangulationTin()} first and then feed it to this method.
	 * 
	 * @param triangulation        a triangulation mesh
	 * @param discardOpenTriangles whether to discard "open" triangles from the
	 *                             output or include these individual triangles as
	 *                             faces. A triangle is open if it lies on the
	 *                             border of the triangulation and is not part of
	 *                             any larger face (which occurs if its longest edge
	 *                             has no neighbouring triangle)
	 * @return a GROUP PShape where each child shape is a single face
	 * @Since 1.1.0
	 */
	public static PShape urquhartFaces(final IncrementalTin triangulation, final boolean discardOpenTriangles) {
		// TODO functionality to merge small groups (area < x) into a neighbouring group

		/*-
		 * Algorithm:
		 * 0) Iterate over every triangle, finding its longest edge; add these into a set.
		 * 1) Attempt to collapse each longest edge.
		 * 2) Collapsing:
		 * 		if the neighbour triangle shares this edge and neither are in a group, add both triangles to a new face group
		 * 		if one triangle is part of a group, add the other triangle to the same group
		 * 		if both triangles are already in different groups, merge the groups together
		 * 3) Merge the triangles in each group into a polygon.
		 */

		/*
		 * Build a map of edges-> triangles. In combination with use getDual(), this is
		 * used to find a triangle edge's neighbouring triangle.
		 */
		final HashMap<IQuadEdge, SimpleTriangle> map = new HashMap<>();
		final HashSet<IQuadEdge> uniqueLongestEdges = new HashSet<IQuadEdge>();

		TriangleCollector.visitSimpleTriangles(triangulation, t -> {
			final IConstraint constraint = t.getContainingRegion();
			if (constraint != null && constraint.definesConstrainedRegion()) {
				map.put(t.getEdgeA(), t);
				map.put(t.getEdgeB(), t);
				map.put(t.getEdgeC(), t);
				uniqueLongestEdges.add(findLongestEdge(t));
			}
		});

		final HashMap<SimpleTriangle, HashSet<SimpleTriangle>> triangleGroups = new HashMap<>(); // map of triangles->face group ID

		uniqueLongestEdges.forEach(e -> {
			final SimpleTriangle t = map.get(e);
			final SimpleTriangle neighbour = map.get(e.getDual());

			if (neighbour == null) {
				// there is no neighbouring triangle which means this face is "open"
				// TODO maybe: Merge the open triangle into neighbouring group
				if (!discardOpenTriangles && !triangleGroups.containsKey(t)) {
					final HashSet<SimpleTriangle> group = new HashSet<>(); // create a new group
					group.add(t); // add this triangle to face group
					triangleGroups.put(t, group); // create mapping for this triangle
				}
				return;
			}

			if (triangleGroups.containsKey(neighbour)) {
				if (triangleGroups.containsKey(t)) {
					// this edge will join two existing groups, so merge the groups (if different)
					if (triangleGroups.get(t) != triangleGroups.get(neighbour)) {
						final HashSet<SimpleTriangle> mergeGroup = triangleGroups.get(t); // merge this triangle into neighbour group
						triangleGroups.get(neighbour).addAll(mergeGroup);
						mergeGroup.forEach(member -> {
							// repoint old group members to new group
							triangleGroups.put(member, triangleGroups.get(neighbour)); // overwrite mapping for each old group member
						});
					}

				} else {
					triangleGroups.get(neighbour).add(t); // add this triangle to an existing face
					triangleGroups.put(t, triangleGroups.get(neighbour));// create mapping for this triangle
				}
			} else {
				if (triangleGroups.containsKey(t)) {
					triangleGroups.get(t).add(neighbour);
					triangleGroups.put(neighbour, triangleGroups.get(t)); // create mapping for neighbour
				} else {
					// neither this triangle or neighbour are part of a group, so create a new one
					// containing both
					final HashSet<SimpleTriangle> group = new HashSet<>(); // create a new group
					group.add(t); // add this triangle to face group
					group.add(neighbour); // add neighbouring triangle to face group too
					triangleGroups.put(t, group); // create mapping for this triangle
					triangleGroups.put(neighbour, group); // create mapping for neighbour
				}
			}
		});

		final PShape out = new PShape(PShape.GROUP);

		new HashSet<>(triangleGroups.values()).forEach(group -> { // intermediate hashset to remove duplicate values

			final ArrayList<Polygon> triangles = new ArrayList<>(group.size());

			for (SimpleTriangle triangle : group) {
				final Coordinate[] coords = new Coordinate[] { toCoord(triangle.getVertexA()), toCoord(triangle.getVertexB()),
						toCoord(triangle.getVertexC()), toCoord(triangle.getVertexA()) };
				triangles.add(PGS.GEOM_FACTORY.createPolygon(coords));
			}

			/*
			 * Use .buffer(0) because it's faster than CascadedPolygonUnion. convexHull() is
			 * much faster still but not it's deterministic (polygons wiggle around).
			 */
			final PShape face = toPShape(PGS.GEOM_FACTORY.buildGeometry(triangles).buffer(0));
			face.setStrokeWeight(3);
			out.addChild(face);
		});

		return out;
	}

	private static PVector toPVector(Vertex v) {
		return new PVector((float) v.getX(), (float) v.getY());
	}

	private static Coordinate centroid(SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new Coordinate(x, y);
	}

	private static Coordinate toCoord(Vertex v) {
		return new Coordinate(v.x, v.y);
	}

	/**
	 * Calculate the longest edge of a given triangle.
	 */
	private static IQuadEdge findLongestEdge(SimpleTriangle t) {
		if (t.getEdgeA().getLength() > t.getEdgeB().getLength()) {
			if (t.getEdgeC().getLength() > t.getEdgeA().getLength()) {
				return t.getEdgeC();
			} else {
				return t.getEdgeA();
			}
		} else {
			if (t.getEdgeC().getLength() > t.getEdgeB().getLength()) {
				return t.getEdgeC();
			} else {
				return t.getEdgeB();
			}
		}
	}
}
