package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.tinfour.common.IConstraint;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;

import earcut4j.Earcut;
import micycle.pgs.PGS.LinearRingIterator;
import micycle.pgs.color.RGB;
import micycle.pgs.utility.Nullable;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Delaunay and earcut triangulation of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Triangulation {

	private PGS_Triangulation() {
	}

	/**
	 * Generates a constrained Delaunay Triangulation from the given shape.
	 * 
	 * @param shape the shape whose vertices to generate a triangulation from
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 */
	public static PShape delaunayTriangulation(PShape shape) {
		return delaunayTriangulation(shape, null, true, 0, true);
	}

	/**
	 * Generates a Delaunay Triangulation from the given shape. The triangulation
	 * can be both constrained (meaning the triangulation is masked by the original
	 * shape) and refined (meaning additional points are inserted, usually leading
	 * to more uniform triangle shapes and sizes).
	 * 
	 * @param shape         the shape whose vertices to generate a triangulation
	 *                      from
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
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #delaunayTriangulationPoints(PShape, List, boolean, int, boolean)
	 * @see #delaunayTriangulationMesh(PShape, List, boolean, int, boolean)
	 */
	public static PShape delaunayTriangulation(PShape shape, @Nullable Collection<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		final IncrementalTin tin = delaunayTriangulationMesh(shape, steinerPoints, constrain, refinements, pretty);

		final PShape triangulation = new PShape(PConstants.GROUP);

		final Consumer<Vertex[]> triangleVertexConsumer = t -> {
			final PShape triangle = new PShape(PShape.PATH);
			triangle.beginShape();
			triangle.vertex((float) t[0].x, (float) t[0].y);
			triangle.vertex((float) t[1].x, (float) t[1].y);
			triangle.vertex((float) t[2].x, (float) t[2].y);
			triangle.endShape(PConstants.CLOSE);
			triangulation.addChild(triangle);
		};

		if (constrain) {
			TriangleCollector.visitTrianglesConstrained(tin, triangleVertexConsumer);
		} else {
			TriangleCollector.visitTriangles(tin, triangleVertexConsumer);
		}

		PGS_Conversion.setAllFillColor(triangulation, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(triangulation, RGB.PINK, 2);

		return triangulation;
	}

	/**
	 * Generates a Delaunay Triangulation from a collection of points.
	 * 
	 * @param points the point collection to triangulate
	 * @return a TRIANGLES PShape
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 * @since 1.1.0
	 */
	public static PShape delaunayTriangulation(Collection<PVector> points) {
		return delaunayTriangulation(null, points, false, 0, false);
	}

	/**
	 * Generates a constrained Delaunay Triangulation from a collection of points.
	 * <p>
	 * This method returns the triangulation as a list of points, rather than a
	 * PShape.
	 * 
	 * @param shape the shape whose vertices to generate a triangulation from
	 * @return List of PVector coordinates, where each consecutive triplet of
	 *         coordinates are the 3 vertices belonging to one triangle
	 */
	public static List<PVector> delaunayTriangulationPoints(PShape shape) {
		return delaunayTriangulationPoints(shape, null, true, 0, true);
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
	public static List<PVector> delaunayTriangulationPoints(PShape shape, @Nullable Collection<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		final IncrementalTin tin = delaunayTriangulationMesh(shape, steinerPoints, constrain, refinements, pretty);

		final ArrayList<PVector> triangles = new ArrayList<>();
		final Consumer<Vertex[]> triangleVertexConsumer = t -> {
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
	 * Generates a Delaunay Triangulation from a collection of points.
	 * <p>
	 * This method returns the triangulation as a list of points, rather than a
	 * PShape.
	 * 
	 * @param points the point collection to triangulate
	 * @return List of PVector coordinates, where each consecutive triplet of
	 *         coordinates are the 3 vertices belonging to one triangle
	 * @see #delaunayTriangulationPoints(PShape, Collection, boolean, int, boolean)
	 * @since 1.1.0
	 */
	public static List<PVector> delaunayTriangulationPoints(Collection<PVector> points) {
		return delaunayTriangulationPoints(null, points, false, 0, false);
	}

	/**
	 * Generates a constrained Delaunay Triangulation from the given shape.
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 * 
	 * @param shape the shape whose vertices to generate a triangulation from
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 */
	public static IncrementalTin delaunayTriangulationMesh(PShape shape) {
		return delaunayTriangulationMesh(shape, null, true, 0, true);
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
	 * @param shape         the shape whose vertices to generate a triangulation
	 *                      from. <b>Can be null</b>.
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement/subdivision
	 *                      passes to perform. Each pass inserts the centroids of
	 *                      every existing triangle into the triangulation. Should
	 *                      be 0 or greater (probably no more than 5).
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
	public static IncrementalTin delaunayTriangulationMesh(@Nullable PShape shape, @Nullable Collection<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		Geometry g = shape == null ? PGS.GEOM_FACTORY.createEmpty(2) : fromPShape(shape);
		final IncrementalTin tin = new IncrementalTin(10);

		final ArrayList<Vertex> vertices = new ArrayList<>();
		final Coordinate[] coords = g.getCoordinates();
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

			/*
			 * A possible optimisation is to recursely split within each triangle upto the
			 * refinement depth (in one pass), so perform many less location checks. Another
			 * is to rasterise the PShape and check pixel[] array. TODO See 'sqrt(3)
			 * Subdivision' by Leif Kobbelt
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
			// If geom is a point set, constrain tin using its concave hull.
			if (g.getGeometryType().equals(Geometry.TYPENAME_MULTIPOINT)) {
				g = fromPShape(PGS_Processing.concaveHull2(PGS_Conversion.toPVector(shape), 0.3));
			}
			List<IConstraint> constraints = new ArrayList<>();
			for (int n = 0; n < g.getNumGeometries(); n++) {
				boolean hole = false;

				LinearRingIterator lri = new LinearRingIterator(g.getGeometryN(n));
				for (LinearRing ring : lri) {
					ArrayList<Vertex> points = new ArrayList<>();
					Coordinate[] c = ring.getCoordinates();
					if (c.length == 0) {
						continue;
					}
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
			if (!constraints.isEmpty()) {
				tin.addConstraints(constraints, pretty); // true/false is negligible?
			}
		}

		return tin;
	}

	/**
	 * Generates a Delaunay Triangulation from a collection of points.
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 * 
	 * @param points the point collection to triangulate
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 * @since 1.1.0
	 */
	public static IncrementalTin delaunayTriangulationMesh(Collection<PVector> points) {
		return delaunayTriangulationMesh(null, points, false, 0, false);
	}

	/**
	 * Creates a Delaunay triangulation of the shape where additional steiner
	 * points, populated by poisson sampling, are included.
	 * 
	 * @param shape
	 * @param spacing (Minimum) spacing between poisson points
	 * @return a GROUP PShape, where each child shape is one triangle
	 * @see #poissonTriangulationPoints(PShape, double)
	 */
	public static PShape poissonTriangulation(PShape shape, double spacing) {
		final Envelope e = fromPShape(shape).getEnvelopeInternal();

		final List<PVector> poissonPoints = PGS_PointSet.poisson(e.getMinX(), e.getMinY(), e.getMinX() + e.getWidth(),
				e.getMinY() + e.getHeight(), spacing, 0);

		final IncrementalTin tin = delaunayTriangulationMesh(shape, poissonPoints, true, 0, false);

		final PShape triangulation = new PShape(PConstants.GROUP);

		TriangleCollector.visitTrianglesConstrained(tin, t -> {
			final PShape triangle = new PShape(PShape.PATH);
			triangle.beginShape();
			triangle.vertex((float) t[0].x, (float) t[0].y);
			triangle.vertex((float) t[1].x, (float) t[1].y);
			triangle.vertex((float) t[2].x, (float) t[2].y);
			triangle.endShape(PConstants.CLOSE);
			triangulation.addChild(triangle);
		});

		PGS_Conversion.setAllFillColor(triangulation, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(triangulation, RGB.PINK, 2);

		return triangulation;
	}

	/**
	 * Creates a Delaunay triangulation of the shape where additional steiner
	 * points, populated by poisson sampling, are included.
	 * 
	 * @param shape
	 * @param spacing (Minimum) spacing between poisson points
	 * @return list of PVectors, where each successive triplet of PVectors
	 *         correspond to the 3 vertices of one triangle
	 * @see #poissonTriangulation(PShape, double)
	 */
	public static List<PVector> poissonTriangulationPoints(PShape shape, double spacing) {
		final Envelope e = fromPShape(shape).getEnvelopeInternal();

		final List<PVector> poissonPoints = PGS_PointSet.poisson(e.getMinX(), e.getMinY(), e.getMinX() + e.getWidth(),
				e.getMinY() + e.getHeight(), spacing, 0);

		final IncrementalTin tin = delaunayTriangulationMesh(shape, poissonPoints, true, 0, false);

		final ArrayList<PVector> triangles = new ArrayList<>();
		TriangleCollector.visitTrianglesConstrained(tin, t -> {
			triangles.add(toPVector(t[0]));
			triangles.add(toPVector(t[1]));
			triangles.add(toPVector(t[2]));
		});
		return triangles;
	}

	/**
	 * Computes a triangulation of the shape according to the ear clipping
	 * ("earcut") method. The triangulation is constrained to the shape by default.
	 * Does not support holes (for now...).
	 * 
	 * @param shape shape whose vertices to triangulate
	 * @return a GROUP PShape, where each child shape is one triangle
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
	 * @return a GROUP PShape, where each child shape is one triangle
	 */
	public static PShape earCutTriangulation(List<PVector> points) {
		final double[] arrCoords = new double[points.size() * 2];

		for (int i = 0; i < points.size(); i++) {
			arrCoords[2 * i] = points.get(i).x;
			arrCoords[2 * i + 1] = points.get(i).y;
		}

		final List<Integer> triangles = Earcut.earcut(arrCoords, null, 2);

		final PShape triangulation = new PShape(PConstants.GROUP);

		for (int i = 0; i < triangles.size(); i += 3) {
			final int v1 = 2 * triangles.get(i);
			final int v2 = 2 * triangles.get(i + 1);
			final int v3 = 2 * triangles.get(i + 2);

			final PShape triangle = new PShape(PShape.PATH);
			triangle.beginShape();
			triangle.vertex((float) arrCoords[v1], (float) arrCoords[v1 + 1]);
			triangle.vertex((float) arrCoords[v2], (float) arrCoords[v2 + 1]);
			triangle.vertex((float) arrCoords[v3], (float) arrCoords[v3 + 1]);
			triangle.endShape(PConstants.CLOSE);
			triangulation.addChild(triangle);
		}

		PGS_Conversion.setAllFillColor(triangulation, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(triangulation, RGB.PINK, 2);

		return triangulation;
	}

	static PVector toPVector(final Vertex v) {
		return new PVector((float) v.getX(), (float) v.getY());
	}

	/**
	 * Computes the centroid/barycentre of a triangle.
	 */
	private static Coordinate centroid(final SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new Coordinate(x, y);
	}
}
