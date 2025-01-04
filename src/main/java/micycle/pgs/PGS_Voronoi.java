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

import micycle.pgs.color.Colors;
import micycle.pgs.commons.Nullable;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Voronoi Diagrams of shapes and point sets. Supports polygonal constraining
 * and relaxation to generate centroidal Voronoi.
 * 
 * @author Michael Carleton
 *
 */
@SuppressWarnings("squid:S3776")
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
	 * @param shape     A shape whose vertices to use as Voronoi sites
	 * @param constrain A flag indicating whether or not to constrain the resulting
	 *                  diagram to the original shape (if it is polygonal).
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(Collection)
	 * @see #innerVoronoi(PShape, boolean, double[], Collection, int)
	 */
	public static PShape innerVoronoi(final PShape shape, final boolean constrain) {
		return innerVoronoi(shape, constrain, null, null, 0);
	}

	/**
	 * Generates an inner Voronoi diagram of a given shape with a specified number
	 * of relaxations.
	 * 
	 * @param shape       The shape to generate the inner Voronoi diagram for.
	 * @param relaxations The number of times to relax the diagram.
	 * @return The generated inner Voronoi diagram as a GROUP PShape, where each
	 *         child shape is a Voronoi cell
	 */
	public static PShape innerVoronoi(final PShape shape, final int relaxations) {
		return innerVoronoi(shape, true, null, null, relaxations);
	}

	/**
	 * Generates an inner Voronoi diagram of a given shape with additional sites.
	 * 
	 * @param shape           The shape to generate the inner Voronoi diagram for.
	 * @param additionalSites A collection of PVector points representing additional
	 *                        sites to be used in the diagram.
	 * @return The generated inner Voronoi diagram as a GROUP PShape, where each
	 *         child shape is a Voronoi cell.
	 * @see #innerVoronoi(PShape, boolean, double[], Collection, int)
	 */
	public static PShape innerVoronoi(final PShape shape, Collection<PVector> additionalSites) {
		return innerVoronoi(shape, true, null, additionalSites, 0);
	}

	/**
	 * Generates an inner Voronoi diagram of a given shape with additional sites and
	 * relaxation.
	 * 
	 * @param shape           The shape to generate the inner Voronoi diagram for.
	 * @param additionalSites A collection of PVector points representing additional
	 *                        sites to be used in the diagram.
	 * @param relaxations     The number of times to relax the diagram.
	 * @return The generated inner Voronoi diagram as a GROUP PShape, where each
	 *         child shape is a Voronoi cell.
	 * @see #innerVoronoi(PShape, boolean, double[], Collection, int)
	 */
	public static PShape innerVoronoi(final PShape shape, Collection<PVector> additionalSites, int relaxations) {
		return innerVoronoi(shape, true, null, additionalSites, relaxations);
	}

	/**
	 * Generates an inner Voronoi diagram for the specified shape, optionally
	 * constrained to the shape's bounds, with additional options for including
	 * Steiner points and specifying the number of Lloyd's relaxations to apply.
	 * This method allows for the creation of a more customized Voronoi diagram by
	 * providing additional control over the generation process.
	 * 
	 * @param shape         The shape to generate the inner Voronoi diagram for.
	 * @param constrain     Indicates whether the resulting Voronoi diagram should
	 *                      be constrained to the bounds of the shape. If true,
	 *                      Voronoi cells will be clipped to the shape's boundary.
	 * @param bounds        An optional array of doubles specifying the bounding box
	 *                      within which the Voronoi diagram is generated. If null,
	 *                      no bounding box is applied beyond the shape's bounds
	 *                      when constraining. The array should be of the form
	 *                      [minX, minY, maxX, maxY].
	 * @param steinerPoints An optional collection of PVector points representing
	 *                      additional sites (Steiner points) to be used in the
	 *                      diagram. These points are in addition to the vertices of
	 *                      the input shape.
	 * @param relaxations   The number of Lloyd's relaxations to apply to the
	 *                      Voronoi diagram. This process helps to create more
	 *                      evenly sized cells by adjusting the position of sites
	 *                      based on their cell centroids.
	 * @return A GROUP PShape where each child shape represents a Voronoi cell.
	 * @since 2.0
	 * @see #innerVoronoi(PShape, boolean)
	 * @see #innerVoronoi(PShape, Collection)
	 * @see #innerVoronoi(PShape, Collection, int)
	 */
	public static PShape innerVoronoi(final PShape shape, final boolean constrain, @Nullable final double[] bounds,
			@Nullable final Collection<PVector> steinerPoints, final int relaxations) {
		final Geometry g = fromPShape(shape);
		BoundedVoronoiDiagram v = innerVoronoiRaw(shape, constrain, bounds, steinerPoints, relaxations);
		List<Geometry> faces = new ArrayList<>();
		if (v != null && v.getPolygons() != null) {
			faces = v.getPolygons().stream().filter(p -> p.getEdges().size() > 1).map(PGS_Voronoi::toPolygon).collect(Collectors.toList());
		}
		if (constrain && g instanceof Polygonal) {
			faces = faces.parallelStream().map(f -> OverlayNG.overlay(f, g, OverlayNG.INTERSECTION)).collect(Collectors.toList());
			faces.removeIf(f -> f.getNumPoints() == 0); // (odd, artifacts of intersection?)
		}

		PShape facesShape = PGS_Conversion.toPShape(faces);
		for (int i = 0; i < faces.size(); i++) {
			if (faces.get(i).getUserData() != null) {
				facesShape.getChild(i).setName(Integer.toString((int) faces.get(i).getUserData()));
			}
		}
		return facesShape;
	}

	/**
	 * Generates a Voronoi diagram of a given shape, where shape vertices are
	 * voronoi point sites. In this method each voronoi cell designates the area
	 * closest to some vertex.
	 * <p>
	 * Note: If the input shape is polygonal, the output is sensitive to how densely
	 * populated lines are in the input. It may be desirable to first process a
	 * polygonal shape with
	 * {@link micycle.pgs.PGS_Processing#densify(PShape, double)
	 * densify(density=~10)} before using this method.
	 * <p>
	 * The diagram may be "relaxed" into a <i>Centroidal Voronoi Diagram</i>. The
	 * relaxation process is a technique used to improve the quality of the Voronoi
	 * diagram. It involves moving the vertices of the diagram slightly to reduce
	 * the maximum distance between a vertex and the centroid of its associated
	 * cell. The process is repeated for a specified number of
	 * <code>relaxations</code> iterations. This process aims to reduce the number
	 * of irregular shaped polygons in the Voronoi diagram and produce a smoother
	 * and more evenly distributed diagram.
	 * 
	 * @param shape         The shape to generate the inner Voronoi diagram for
	 *                      (using its vertices for Voronoi sites).
	 * @param constrain     A flag indicating whether or not to constrain the
	 *                      resulting diagram to the original shape (if it is
	 *                      polygonal).
	 * @param bounds        an optional array of the form [minX, minY, maxX, maxY]
	 *                      representing the bounds of the diagram. The boundary
	 *                      must fully contain the shape (but needn't contain all
	 *                      steiner points).
	 * @param steinerPoints an optional collection of PVector points representing
	 *                      Steiner points to be used as additional sites in the
	 *                      diagram.
	 * @param relaxations   the number of times to relax the diagram. 0 or greater.
	 * 
	 * @return a GROUP PShape, where each child shape is a Voronoi cell. The
	 *         <code>.name</code> value of each cell is set to the integer index of
	 *         its vertex site.
	 * @see #innerVoronoi(Collection)
	 */
	public static BoundedVoronoiDiagram innerVoronoiRaw(final PShape shape, final boolean constrain, @Nullable final double[] bounds,
			@Nullable final Collection<PVector> steinerPoints, final int relaxations) {
		final Geometry g = fromPShape(shape);
		final List<Vertex> vertices = new ArrayList<>();
		final Coordinate[] coords = g.getCoordinates();
		if (coords.length < 3) { // at least 3 vertices are required
			return null;
		}

		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		final Rectangle2D boundsRect;
		if (bounds == null) {
			final Envelope e = g.getEnvelopeInternal();
			boundsRect = new Rectangle2D.Double(e.getMinX(), e.getMinY(), e.getWidth(), e.getHeight());
		} else {
			boundsRect = new Rectangle2D.Double(bounds[0], bounds[1], bounds[2] - bounds[0], bounds[3] - bounds[1]);
		}
		options.setBounds(boundsRect);
		options.enableAutomaticColorAssignment(false);

		for (int i = 0; i < coords.length; i++) {
			Coordinate p = coords[i];
			if (boundsRect.contains(p.x, p.y)) {
				vertices.add(new Vertex(p.x, p.y, Double.NaN, i));
			}
		}
		if (steinerPoints != null) {
			steinerPoints.forEach(p -> {
				if (boundsRect.contains(p.x, p.y)) {
					vertices.add(new Vertex(p.x, p.y, vertices.size()));
				}
			});
		}

		HilbertSort hs = new HilbertSort();
		hs.sort(vertices); // prevent degenerate insertion

		BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(vertices, options);

		for (int i = 0; i < relaxations; i++) {
			double maxDistDelta = 0;
			List<Vertex> newSites = new ArrayList<>(vertices.size());
			for (ThiessenPolygon p : v.getPolygons()) {
				final Vertex newSite;
				final PVector centroid = computeCentroid(p);
				if (p.getVertex().getIndex() == 0 || steinerPoints == null) {
					PVector site = new PVector((float) p.getVertex().x, (float) p.getVertex().y);
					site.add(PVector.sub(centroid, site).mult(1.5f)); // over-relax (1.5x)
					newSite = new Vertex(site.x, site.y, p.getIndex());
					maxDistDelta = Math.max(maxDistDelta, p.getVertex().getDistance(centroid.x, centroid.y));
				} else {
					newSite = p.getVertex();
				}
				if (boundsRect.contains(newSite.x, newSite.y)) {
					newSites.add(newSite);
				}
			}
			if (maxDistDelta < 0.05) {
				break; // sufficiently converged, exit relaxation early
			}
			v = new BoundedVoronoiDiagram(newSites, options);
		}

		return v;
	}

	/**
	 * Generates a Voronoi diagram for a set of points. In this method each voronoi
	 * cell designates the area closest to some point.
	 * 
	 * @param points the set of points to use as Voronoi sites
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(PShape, boolean, double[], Collection, int)
	 */
	public static PShape innerVoronoi(Collection<PVector> points) {
		return innerVoronoi(PGS_Conversion.toPointsPShape(points), false);
	}

	/**
	 * Generates a Voronoi diagram for a set of points, with relaxation. In this
	 * method each voronoi cell designates the area closest to some point.
	 * 
	 * @param points      the set of points to use as Voronoi sites
	 * @param relaxations the number of times to relax the diagram. 0 or greater.
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(PShape, boolean, double[], Collection, int)
	 */
	public static PShape innerVoronoi(Collection<PVector> points, int relaxations) {
		return innerVoronoi(PGS_Conversion.toPointsPShape(points), false, null, null, relaxations);
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
	 */
	public static PShape innerVoronoi(Collection<PVector> points, double[] bounds) {
		return innerVoronoi(PGS_Conversion.toPointsPShape(points), false, bounds, null, 0);
	}

	/**
	 * Generates a boundary-constrained Voronoi diagram for a set of points, with
	 * relaxation. In this method each voronoi cell designates the area closest to
	 * some point.
	 * 
	 * @param points      the set of points to use as Voronoi sites
	 * @param bounds      an array of the form [minX, minY, maxX, maxY] representing
	 *                    the bounds of the diagram. The boundary must fully contain
	 *                    the shape.
	 * @param relaxations the number of times to relax the diagram. 0 or greater.
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #innerVoronoi(PShape, boolean, double[], Collection, int)
	 */
	public static PShape innerVoronoi(Collection<PVector> points, double[] bounds, int relaxations) {
		return innerVoronoi(PGS_Conversion.toPointsPShape(points), false, bounds, null, relaxations);
	}

	/**
	 * Converts a collection of points into a Raw Voronoi diagram object within
	 * specified bounds, optionally applying Lloyd's relaxation to improve the
	 * diagram's properties.
	 *
	 * @param points      a Collection of {@link PVector} objects representing the
	 *                    points to be used as the sites of the Voronoi diagram.
	 * @param bounds      an array of four doubles specifying the rectangular
	 *                    boundary of the Voronoi diagram in the form [minX, minY,
	 *                    maxX, maxY]. The bounds must be large enough to completely
	 *                    enclose all the points.
	 * @param relaxations an integer specifying the number of Lloyd's relaxation
	 *                    iterations to perform on the Voronoi diagram. This process
	 *                    can help make the cells more uniform in size and shape. A
	 *                    value of 0 indicates no relaxation.
	 * @return a {@link PShape} object representing the Voronoi diagram, where each
	 *         child shape corresponds to a Voronoi cell associated with each input
	 *         point. The PShape is of type GROUP, allowing for collective
	 *         manipulation of the cells.
	 * @since 2.0
	 * @see #innerVoronoi(PShape, boolean, double[], Collection, int) for a more
	 *      detailed method that allows for additional customization of the Voronoi
	 *      diagram generation process.
	 */

	public static BoundedVoronoiDiagram innerVoronoiRaw(Collection<PVector> points, double[] bounds, int relaxations) {
		return innerVoronoiRaw(PGS_Conversion.toPointsPShape(points), false, bounds, null, relaxations);
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
	 * @return GROUP shape consisting of voronoi cells; each cell corresponds to an
	 *         area around a line segment for which the closest line segment to any
	 *         point in that area is the line segment
	 * @since 1.3.0
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
	 * @return GROUP shape consisting of voronoi cells; each cell corresponds to an
	 *         area around a line segment for which the closest line segment to any
	 *         point in that area is the line segment
	 * @since 1.3.0
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

		/*
		 * There is a voronoi cell for each densified vertex. We first group densified
		 * vertices by their source geometry and then union/dissolve the cells belonging
		 * to each vertex group.
		 */
		final List<PShape> faces = segmentVertexGroups.parallelStream().map(vertexGroup -> {
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
			return PGS_ShapeBoolean.unionMesh(cellSegments);
		}).collect(Collectors.toList());

		PShape voronoiCells = PGS_Conversion.flatten(faces);
		PGS_Conversion.setAllFillColor(voronoiCells, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(voronoiCells, Colors.PINK, 2);

		return voronoiCells;
	}

	static Polygon toPolygon(ThiessenPolygon polygon) {
		Coordinate[] coords = new Coordinate[polygon.getEdges().size() + 1];
		int i = 0;
		for (IQuadEdge e : polygon.getEdges()) {
			coords[i++] = new Coordinate(e.getA().x, e.getA().y);
		}
		coords[i] = new Coordinate(polygon.getEdges().get(0).getA().x, polygon.getEdges().get(0).getA().y); // close polygon

		Polygon p = PGS.GEOM_FACTORY.createPolygon(coords);
		p.setUserData(polygon.getIndex()); // preserve polygon index
		return p;
	}

	private static PVector computeCentroid(ThiessenPolygon polygon) {
		double xSum = 0;
		double ySum = 0;
		int n = 0;
		for (IQuadEdge e : polygon.getEdges()) {
			xSum += e.getA().x;
			ySum += e.getA().y;
			n++;
		}
		return new PVector((float) xSum / n, (float) ySum / n);
	}

	private static List<Vertex> toVertex(Coordinate[] coords) {
		final boolean closed = coords[0].equals2D(coords[coords.length - 1]) && coords.length > 1;
		List<Vertex> vertices = new ArrayList<>(coords.length - (closed ? 1 : 0));
		for (int i = 0; i < coords.length - (closed ? 1 : 0); i++) {
			Coordinate coord = coords[i];
			vertices.add(new Vertex(coord.x, coord.y, 0));
		}
		return vertices;
	}
}
