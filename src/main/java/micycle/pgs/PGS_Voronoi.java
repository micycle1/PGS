package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.HilbertSort;
import org.tinfour.voronoi.BoundedVoronoiBuildOptions;
import org.tinfour.voronoi.BoundedVoronoiDiagram;
import org.tinfour.voronoi.ThiessenPolygon;

import micycle.pgs.color.RGB;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Voronoi Diagrams of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Voronoi {

	private PGS_Voronoi() {
	}

	/**
	 * Generates a Voronoi diagram for a single shape, where shape vertices are
	 * voronoi point sites. In this method each voronoi cell designates the area
	 * closest to some vertex.
	 * <p>
	 * Note: If the input shape is polygonal, the output is sensitive to how densely
	 * populated lines are in the input. Consider processing a shape with
	 * {@link micycle.pgs.PGS_Processing#densify(PShape, double)
	 * densify(density=~10)} method first before using this method on a polygon.
	 * 
	 * @param shape     a shape whose vertices to use as Voronoi sites
	 * @param constrain whether to constrain the diagram lines to the shape (if
	 *                  polygonal). When true, output voronoi cells are
	 *                  cropped/constrained to the shape outline.
	 * 
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(Collection)
	 */
	public static PShape innerVoronoi(final PShape shape, final boolean constrain) {
		return innerVoronoi(shape, constrain, null);
	}

	/**
	 * Generates a Voronoi diagram for a single shape, where shape vertices are
	 * voronoi point sites. In this method each voronoi cell designates the area
	 * closest to some vertex.
	 * <p>
	 * Note: If the input shape is polygonal, the output is sensitive to how densely
	 * populated lines are in the input. Consider processing a shape with
	 * {@link micycle.pgs.PGS_Processing#densify(PShape, double)
	 * densify(density=~10)} method first before using this method on a polygon.
	 * 
	 * @param shape     a shape whose vertices to use as Voronoi sites
	 * @param constrain whether to constrain the diagram lines to the shape (if
	 *                  polygonal). When true, output voronoi cells are
	 *                  cropped/constrained to the shape outline.
	 * @param bounds    an array of the form [minX, minY, maxX, maxY] defining the
	 *                  boundary of the voronoi diagram. the boundary must fully
	 *                  contain the shape.
	 * 
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(Collection)
	 */
	public static PShape innerVoronoi(final PShape shape, final boolean constrain, double[] bounds) {
		final IncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape, null, false, 0, false);

		final Geometry g = fromPShape(shape);
		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		final double x, y, w, h;
		if (bounds == null) {
			final Envelope envelope = g.getEnvelopeInternal();
			x = envelope.getMinX();
			y = envelope.getMinY();
			w = envelope.getMaxX() - envelope.getMinX();
			h = envelope.getMaxY() - envelope.getMinY();
		} else {
			x = bounds[0];
			y = bounds[1];
			w = bounds[2] - bounds[0];
			h = bounds[3] - bounds[1];
		}
		options.setBounds(new Rectangle2D.Double(x, y, w, h));

		final BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(tin.getVertices(), options);

		List<Geometry> faces = v.getPolygons().stream().map(PGS_Voronoi::toPolygon).collect(Collectors.toList());
		if (constrain && g instanceof Polygonal) {
			faces = faces.parallelStream().map(f -> OverlayNG.overlay(f, g, OverlayNG.INTERSECTION)).collect(Collectors.toList());
		}

		return PGS_Conversion.toPShape(faces);
	}

	/**
	 * Generates a Voronoi diagram for a set of points. In this method each voronoi
	 * cell designates the area closest to some point.
	 * 
	 * @param points the set of points to use as Voronoi sites
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(PShape)
	 */
	public static PShape innerVoronoi(Collection<PVector> points) {
		return innerVoronoi(PGS_Conversion.toPointsPShape(points), false);
	}

	/**
	 * Generates a Voronoi diagram for a set of points. In this method each voronoi
	 * cell designates the area closest to some point.
	 * 
	 * @param points the set of points to use as Voronoi sites
	 * @param bounds an array of the form [minX, minY, maxX, maxY] defining the
	 *               boundary of the voronoi diagram. the boundary must fully
	 *               contain the points.
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(PShape)
	 */
	public static PShape innerVoronoi(Collection<PVector> points, double[] bounds) {
		return innerVoronoi(PGS_Conversion.toPointsPShape(points), false, bounds);
	}

	/**
	 * Generates a Voronoi diagram for a set of disjoint shapes. In this method each
	 * voronoi cell designates the area closest to some individual shape.
	 * <p>
	 * Note: Each geometry primitive in a <code>POINTS</code> or <code>LINES</code>
	 * shape is treated as a distinct voronoi site (rather than a singular site
	 * representing the full mass of points or lines).
	 * 
	 * @param shape a GROUP PShape consisting of any number of non-intersecting
	 *              polygonal, lineal, or points child shapes
	 * @return a GROUP PShape, where each child shape is a Voronoi cell, bounded by
	 *         the envelope all shapes
	 * @since 1.3.0
	 * @return GROUP shape consisting of voronoi cells; each cell corresponds to an
	 *         area around a line segment for which the closest line segment to any
	 *         point in that area is the line segment
	 */
	public static PShape compoundVoronoi(PShape shape) {
		return compoundVoronoi(shape, null);
	}

	/**
	 * Generates a Voronoi diagram for a set of disjoint shapes. In this method each
	 * voronoi cell designates the area closest to some individual shape.
	 * <p>
	 * Note: Each geometry primitive in a <code>POINTS</code> or <code>LINES</code>
	 * shape is treated as a distinct voronoi site (rather than a singular site
	 * representing the full mass of points or lines).
	 * 
	 * @param shape  a GROUP PShape consisting of any number of non-intersecting
	 *               polygonal, lineal, or points child shapes
	 * @param bounds an array of the form [minX, minY, maxX, maxY] defining the
	 *               boundary of the voronoi diagram. the boundary must fully
	 *               contain the shape.
	 * @return a GROUP PShape, where each child shape is a Voronoi cell, bounded by
	 *         the envelope all shapes
	 * @since 1.3.0
	 * @return GROUP shape consisting of voronoi cells; each cell corresponds to an
	 *         area around a line segment for which the closest line segment to any
	 *         point in that area is the line segment
	 */
	public static PShape compoundVoronoi(PShape shape, double[] bounds) {
		Geometry g = fromPShape(shape);
		Geometry densified = Densifier.densify(g, 2);

		List<Vertex> vertices = new ArrayList<>();
		final List<List<Vertex>> segmentVertexGroups = new ArrayList<>();

		for (int i = 0; i < densified.getNumGeometries(); i++) {
			Geometry geometry = densified.getGeometryN(i);
			List<Vertex> featureVertices;
			switch (geometry.getGeometryType()) {
				case Geometry.TYPENAME_LINEARRING :
				case Geometry.TYPENAME_POLYGON :
				case Geometry.TYPENAME_LINESTRING :
				case Geometry.TYPENAME_POINT :
					featureVertices = toVertex(geometry.getCoordinates());
					if (!featureVertices.isEmpty()) {
						segmentVertexGroups.add(featureVertices);
						vertices.addAll(featureVertices);
					}
					break;
				case Geometry.TYPENAME_MULTILINESTRING :
				case Geometry.TYPENAME_MULTIPOINT :
				case Geometry.TYPENAME_MULTIPOLYGON : // nested multi polygon
					for (int j = 0; j < geometry.getNumGeometries(); j++) {
						featureVertices = toVertex(geometry.getGeometryN(j).getCoordinates());
						if (!featureVertices.isEmpty()) {
							segmentVertexGroups.add(featureVertices);
							vertices.addAll(featureVertices);
						}
					}
					break;
				default :
					break;
			}
		}

		if (vertices.size() > 2500) {
			HilbertSort hs = new HilbertSort();
			hs.sort(vertices);
		}
		final IncrementalTin tin = new IncrementalTin(2);
		tin.add(vertices, null); // initial triangulation
		if (!tin.isBootstrapped()) {
			return new PShape(); // shape probably empty
		}

		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		final double x, y, w, h;
		if (bounds == null) {
			final Envelope envelope = g.getEnvelopeInternal();
			x = envelope.getMinX();
			y = envelope.getMinY();
			w = envelope.getMaxX() - envelope.getMinX();
			h = envelope.getMaxY() - envelope.getMinY();
		} else {
			x = bounds[0];
			y = bounds[1];
			w = bounds[2] - bounds[0];
			h = bounds[3] - bounds[1];
		}
		options.setBounds(new Rectangle2D.Double(x, y, w, h));

		final BoundedVoronoiDiagram voronoi = new BoundedVoronoiDiagram(tin);

		// Map densified vertices to the voronoi cell they define.
		final HashMap<Vertex, ThiessenPolygon> vertexCellMap = new HashMap<>();
		voronoi.getPolygons().forEach(p -> vertexCellMap.put(p.getVertex(), p));

		PShape voronoiCells = new PShape();

		/*
		 * There is a voronoi cell for each densified vertex. We first group densified
		 * vertices by their source geometry and then union/dissolve the cells belonging
		 * to each vertex group.
		 */
		segmentVertexGroups.forEach(vertexGroup -> {
			PShape cellSegments = new PShape(PConstants.GROUP);
			vertexGroup.forEach(segmentVertex -> {
				ThiessenPolygon thiessenCell = vertexCellMap.get(segmentVertex);
				if (thiessenCell != null) { // null if degenerate input
					PShape cellSegment = new PShape(PShape.PATH);
					cellSegment.beginShape();
					for (IQuadEdge e : thiessenCell.getEdges()) {
						cellSegment.vertex((float) e.getA().x, (float) e.getA().y);
					}
					cellSegment.endShape(PConstants.CLOSE);
					cellSegments.addChild(cellSegment);
				}
			});
			voronoiCells.addChild(PGS_ShapeBoolean.unionMesh(cellSegments));
		});

		PGS_Conversion.setAllFillColor(voronoiCells, RGB.WHITE);
		PGS_Conversion.setAllStrokeColor(voronoiCells, RGB.PINK, 2);

		return voronoiCells;
	}

	private static Polygon toPolygon(ThiessenPolygon polygon) {
		Coordinate[] coords = new Coordinate[polygon.getEdges().size() + 1];
		int i = 0;
		for (IQuadEdge e : polygon.getEdges()) {
			coords[i++] = new Coordinate(e.getA().x, e.getA().y);
		}
		coords[i] = new Coordinate(polygon.getEdges().get(0).getA().x, polygon.getEdges().get(0).getA().y); // close polygon
		return PGS.GEOM_FACTORY.createPolygon(coords);
	}

	private static List<Vertex> toVertex(Coordinate[] coords) {
		final boolean closed = coords[0].equals2D(coords[coords.length - 1]) && coords.length > 1;
		List<Vertex> vertexes = new ArrayList<>(coords.length - (closed ? 1 : 0));
		for (int i = 0; i < coords.length - (closed ? 1 : 0); i++) {
			Coordinate coord = coords[i];
			vertexes.add(new Vertex(coord.x, coord.y, 0));
		}
		return vertexes;
	}
}
