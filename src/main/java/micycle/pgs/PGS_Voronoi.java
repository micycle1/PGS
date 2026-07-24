package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static micycle.pgs.PGS.GEOM_FACTORY;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import org.locationtech.jts.coverage.CoverageSimplifier;
import org.locationtech.jts.coverage.CoverageUnion;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.TopologyException;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.operation.overlay.snap.GeometrySnapper;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.relateng.RelateNG;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.utils.HilbertSort;
import org.tinfour.voronoi.BoundedVoronoiBuildOptions;
import org.tinfour.voronoi.BoundedVoronoiDiagram;
import org.tinfour.voronoi.ThiessenPolygon;

import com.github.micycle1.geoblitz.HilbertParallelPolygonUnion;
import com.github.quickhull3d.PowerDiagram2D;
import com.github.quickhull3d.PowerDiagram2D.Rect;

import micycle.pgs.color.Colors;
import micycle.pgs.commons.AdditivelyWeightedVoronoi;
import micycle.pgs.commons.DiscreteCurveEvolution;
import micycle.pgs.commons.DiscreteCurveEvolution.DCETerminationCallback;
import micycle.pgs.commons.FarthestPointVoronoi;
import micycle.pgs.commons.ManhattanVoronoi;
import micycle.pgs.commons.MultiplicativelyWeightedVoronoi;
import micycle.pgs.commons.Nullable;
import micycle.pgs.commons.PEdge;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Voronoi diagram utilities for 2D point sets and {@link PShape} polygons.
 *
 * <p>
 * This class generates several variants of Voronoi diagrams, including:
 * standard (unweighted) diagrams, additively/multiplicatively weighted
 * diagrams, farthest-point Voronoi, and polygon-constrained (“inner”) Voronoi.
 *
 * <h2>Centroidal Voronoi (relaxation)</h2>
 * <p>
 * Several {@code innerVoronoi(...)} overloads support Lloyd-style relaxation by
 * repeatedly rebuilding the diagram and moving sites toward cell centroids,
 * producing centroidal Voronoi tessellations (CVTs) inside a boundary polygon.
 *
 * @author Michael Carleton
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
		BoundedVoronoiDiagram v = innerVoronoiRaw(shape, bounds, steinerPoints, relaxations);
		List<Geometry> faces = new ArrayList<>();
		if (v != null && v.getPolygons() != null) {
			faces = v.getPolygons().stream().filter(p -> p.getEdges().size() > 1).map(PGS_Voronoi::toPolygon).collect(Collectors.toList());
		}
		if (constrain) {
			final Geometry g = fromPShape(shape);
			if (g instanceof Polygonal) {
				final var index = RelateNG.prepare(g);
				faces = faces.parallelStream().map(f -> {
					final var relation = index.evaluate(f);
					if (relation.isContains()) {
						return f;
					} else if (relation.isDisjoint()) {
						return g.getFactory().createEmpty(2);
					}
					return OverlayNG.overlay(f, g, OverlayNG.INTERSECTION);
				}).collect(Collectors.toList());
				faces.removeIf(f -> f.isEmpty() || f.getNumPoints() == 0);
			}
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
	 * @param bounds        an optional array of the form [minX, minY, maxX, maxY]
	 *                      representing the bounds of the diagram. The boundary
	 *                      must fully contain the shape (but needn't contain all
	 *                      steiner points).
	 * @param steinerPoints an optional collection of PVector points representing
	 *                      Steiner points to be used as additional sites in the
	 *                      diagram.
	 * @param relaxations   the number of times to relax the diagram. 0 or greater.
	 * @return a GROUP PShape, where each child shape is a Voronoi cell. The
	 *         <code>.name</code> value of each cell is set to the integer index of
	 *         its vertex site.
	 * @see #innerVoronoi(Collection)
	 */
	public static BoundedVoronoiDiagram innerVoronoiRaw(final PShape shape, @Nullable final double[] bounds, @Nullable final Collection<PVector> steinerPoints,
			final int relaxations) {
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
			if (maxDistDelta < 1e-6) {
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

	public static BoundedVoronoiDiagram innerVoronoiRaw(Collection<PVector> points, @Nullable double[] bounds, int relaxations) {
		return innerVoronoiRaw(PGS_Conversion.toPointsPShape(points), bounds, null, relaxations);
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

		List<Vertex> vertices = new ArrayList<>(Math.max(16, densified.getNumPoints()));
		List<List<Vertex>> segmentVertexGroups = new ArrayList<>(Math.max(16, densified.getNumGeometries()));

		collectVertexGroups(densified, segmentVertexGroups, vertices);

		if (vertices.size() > 2500) {
			HilbertSort hs = new HilbertSort();
			hs.sort(vertices);
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

		final BoundedVoronoiDiagram voronoi = new BoundedVoronoiDiagram(vertices, options);

		// Map densified vertices to the voronoi cell they define.
		final HashMap<Vertex, ThiessenPolygon> vertexCellMap = new HashMap<>();
		voronoi.getPolygons().forEach(p -> vertexCellMap.put(p.getVertex(), p));

		/*
		 * There is a voronoi cell for each densified vertex. We first group densified
		 * vertices by their source geometry and then union/dissolve the cells belonging
		 * to each vertex group.
		 */
		final List<Geometry> faces = segmentVertexGroups.parallelStream().map(vertexGroup -> {
			var cells = new ArrayList<Geometry>(vertexGroup.size());
			vertexGroup.forEach(segmentVertex -> {
				ThiessenPolygon thiessenCell = vertexCellMap.get(segmentVertex);
				if (thiessenCell != null) { // null if degenerate input
					cells.add(toPolygon(thiessenCell));
				}
			});

			try {
				return CoverageUnion.union(cells.toArray(new Geometry[0]));
			} catch (TopologyException e) {
				var gf = cells.get(0).getFactory();
				var valid = GeometryFixer.fix(gf.createGeometryCollection(cells.toArray(new Geometry[0])));
				return HilbertParallelPolygonUnion.union(valid);
			}

		}).toList();

		PShape voronoiCells = toPShape(faces);
		PGS_Conversion.setAllFillColor(voronoiCells, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(voronoiCells, Colors.PINK, 2);
		return voronoiCells;
	}

	/**
	 * Generates an <b>additively weighted Voronoi diagram</b> (AWVD) for a set of
	 * weighted point sites, clipped to the provided bounding rectangle.
	 * <p>
	 * AWVDs are a generalisation of standard Voronoi diagrams where each site has
	 * an additive weight. Distances are compared using an adjusted metric of the
	 * form:
	 * 
	 * <pre>
	 *   d(p, s) = ||p - s|| - w
	 * </pre>
	 * 
	 * where {@code s} is the site location and {@code w} is its weight. Increasing
	 * a site's weight tends to expand its cell; decreasing it tends to shrink the
	 * cell. Unlike standard Voronoi diagrams, AWVD cell boundaries are generally
	 * <i>curved</i> (hyperbolic arcs), and some sites may end up with empty cells
	 * depending on weights and configuration.
	 * <p>
	 * Each input {@link PVector} encodes one site where:
	 * <ul>
	 * <li>{@code (.x, .y)} is the site coordinate</li>
	 * <li>{@code .z} is the site's weight (in the same units as {@code x/y})</li>
	 * </ul>
	 * <p>
	 * Post-processing:
	 * <ul>
	 * <li>If {@code forceConforming} is {@code true}, additional meshing/coverage
	 * operations are applied to remove tiny gaps between adjacent cells and
	 * simplify the interior boundaries.</li>
	 * </ul>
	 *
	 * @param weightedSites a collection of weighted sites encoded as PVectors:
	 *                      {@code (.x, .y)} position and {@code .z} weight
	 * @param bounds        an array of the form {@code [minX, minY, maxX, maxY]}
	 *                      defining the clipping bounds of the diagram; must fully
	 *                      contain all sites
	 * @return a GROUP {@link PShape} where each child shape is a (possibly curved)
	 *         AWVD cell polygon clipped to {@code bounds}
	 * @since 2.2
	 */
	public static PShape additivelyWeightedVoronoi(Collection<PVector> weightedSites, double[] bounds) {
		var sites = weightedSites.stream().map(s -> PGS.coordFromPVector(s)).toList();
		var e = new Envelope(bounds[0], bounds[2], bounds[1], bounds[3]); // x,x,y,y

		AdditivelyWeightedVoronoi vd = new AdditivelyWeightedVoronoi(GEOM_FACTORY, 0.25);
		List<? extends Geometry> cells = vd.computeCells(sites, e);

		// Produces rather dense output, so simplify using conservative DCE relevance
		final DCETerminationCallback dceCallback = (currentVertex, relevance, verticesRemaining) -> relevance >= 10;
		cells = cells.stream().map(cell -> {
			var ring = DiscreteCurveEvolution.process((LineString) cell.getBoundary(), dceCallback);
			return ring;
		}).toList();

		var awvd = toPShape(cells);
		PGS_Conversion.setAllFillColor(awvd, Colors.WHITE);
		PGS_Conversion.setAllStrokeColor(awvd, Colors.PINK, 2);

		return awvd;
	}

	/**
	 * Generates a Multiplicatively Weighted Voronoi Diagrams diagram for a set of
	 * weighted sites.
	 * <p>
	 * MWVDs are a generalisation of Voronoi diagrams where each site has a weight
	 * associated with it. These weights influence the boundaries between cells in
	 * the diagram. Instead of being equidistant from generator points, the
	 * boundaries are defined by the <b>ratio</b> of distances to the weighted
	 * generator points. This results in characteristically curved cell boundaries,
	 * unlike the straight line boundaries seen in standard Voronoi diagrams.
	 * 
	 * @param weightedSites A list of PVectors, each representing one site:
	 *                      <code>(.x, .y)</code> represent the coordinate and
	 *                      <b><code>.z</code> represents weight</b>.
	 * @param bounds        an array of the form [minX, minY, maxX, maxY]
	 *                      representing the bounds of the diagram. The boundary
	 *                      must cover all points.
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @since 2.0
	 */
	public static PShape multiplicativelyWeightedVoronoi(Collection<PVector> weightedSites, double[] bounds) {
		return multiplicativelyWeightedVoronoi(weightedSites, bounds, false);
	}

	/**
	 * Generates a Multiplicatively Weighted Voronoi Diagrams diagram for a set of
	 * weighted sites.
	 * <p>
	 * MWVDs are a generalisation of Voronoi diagrams where each site has a weight
	 * associated with it. These weights influence the boundaries between cells in
	 * the diagram. Instead of being equidistant from generator points, the
	 * boundaries are defined by the <b>ratio</b> of distances to the weighted
	 * generator points. This results in characteristically curved cell boundaries,
	 * unlike the straight line boundaries seen in standard Voronoi diagrams.
	 * 
	 * @param weightedSites   A list of PVectors, each representing one site:
	 *                        <code>(.x, .y)</code> represent the coordinate and
	 *                        <b><code>.z</code> represents weight</b>.
	 * @param bounds          an array of the form [minX, minY, maxX, maxY]
	 *                        representing the bounds of the diagram. The boundary
	 *                        must cover all points.
	 * @param forceConforming Whether to apply additional processing to ensure the
	 *                        MWVD creates a conforming mesh. A conforming mesh has
	 *                        no tiny gaps between adjacent cells.
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @since 2.0
	 */
	public static PShape multiplicativelyWeightedVoronoi(Collection<PVector> weightedSites, double[] bounds, boolean forceConforming) {
		var faces = MultiplicativelyWeightedVoronoi.getMWVFromPVectors(weightedSites.stream().toList(), bounds);
		Geometry geoms = GEOM_FACTORY.createGeometryCollection(faces.toArray(new Geometry[] {}));
		if (forceConforming) {
			geoms = GeometrySnapper.snapToSelf(geoms, 1e-5, true); // slow
		}
		var s = PGS_Conversion.toPShape(geoms);
//		if (forceConforming) {
//			s = PGS_Meshing.fixBreaks(s, 1e-4); // faster than GeometrySnapper, less robust
//		}
		return s;
	}

	/**
	 * Generates the <b>farthest-point Voronoi diagram</b> (FPVD) for a set of
	 * sites.
	 * <p>
	 * The farthest-point Voronoi diagram partitions the plane into regions such
	 * that each region consists of all points for which a particular site is the
	 * <b>farthest</b> among all provided sites (not the nearest). Only sites that
	 * are convex hull vertices have regions in the FPVD.
	 * <p>
	 * The resulting diagram is not clipped to a bounding box and may extend well
	 * beyond the convex hull of the input sites, but it is still represented by a
	 * finite set of edges.
	 *
	 * @param sites a collection of {@link PVector} representing the sites; only
	 *              convex hull vertices have regions
	 * @return a {@link PShape} representing the farthest-point Voronoi diagram as a
	 *         set of edges
	 * @see #farthestPointVoronoi(Collection, double[])
	 * @since 2.1
	 */
	public static PShape farthestPointVoronoi(Collection<PVector> sites) {
		FarthestPointVoronoi fpvd = new FarthestPointVoronoi();
		fpvd.setSites(sites.stream().map(s -> PGS.coordFromPVector(s)).toList());

		var edges = fpvd.getDCEL().getEdges().stream().map(e -> {
			var a = PGS.toPVector(e.origVertex);
			var b = PGS.toPVector(e.destVertex);
			return new PEdge(a, b);
		}).toList();

		return PGS_SegmentSet.toPShape(edges);
	}

	/**
	 * Generates a <b>farthest-point Voronoi diagram</b> (FPVD) for a given set of
	 * sites and a bounding box.
	 * <p>
	 * The <i>farthest-point Voronoi diagram</i> is a variant of the Voronoi diagram
	 * in which the region for each site <code>p</code> consists of all points in
	 * the plane for which <code>p</code> is the <b>farthest</b> site among all
	 * sites. In contrast to a regular (nearest-point) Voronoi diagram, where each
	 * region surrounds and contains its generating site, in the FPVD, regions do
	 * <b>not</b> hug their site; instead, the generator of a region is {typically
	 * distant and not even contained within} its region.
	 * <p>
	 * <b>Properties:</b>
	 * <ul>
	 * <li>Only sites that are vertices of the convex hull have non-empty regions in
	 * the FPVD, since only those can be farthest from some location in the
	 * plane.</li>
	 * <li>A useful interpretation: all points in a given FPVD region share the same
	 * farthest generator site. However, the generator site is not visually apparent
	 * from the region itself, as it is not located within or even near the region.
	 * </ul>
	 * 
	 * @param sites  A collection of {@link PVector}s representing sites; only the
	 *               convex hull vertices will have corresponding regions in the
	 *               output diagram.
	 * @param bounds A double array of form <code>[minX, minY, maxX, maxY]</code>
	 *               representing the axis-aligned bounding box for clipping the
	 *               diagram.
	 * @return A {@link PShape} representing the farthest-point Voronoi diagram.
	 *         Each cell corresponds to the region for one convex hull vertex site.
	 * @since 2.1
	 */
	public static PShape farthestPointVoronoi(Collection<PVector> sites, double[] bounds) {
		FarthestPointVoronoi fpvd = new FarthestPointVoronoi();
		Envelope e = new Envelope(bounds[0], bounds[2], bounds[1], bounds[3]); // x,x,y,y
		fpvd.setClipEnvelope(e);
		fpvd.setSites(sites.stream().map(s -> PGS.coordFromPVector(s)).toList());

		return toPShape(fpvd.getDiagram());
	}

	/**
	 * Computes a <b>power diagram</b> (a.k.a. <i>Laguerre–Voronoi</i> diagram) for
	 * a set of <b>weighted</b> sites, with no clipping bounds.
	 * <p>
	 * Each site is given as a {@link PVector} where {@code (.x, .y)} is the site
	 * location and {@code .z} is its weight.
	 * <h3>Intuition</h3> A power diagram is the weighted analogue of a standard
	 * Voronoi diagram, but it still produces <b>straight-edged (polygonal)
	 * cells</b>. Increasing a site's weight can allow it to “win” territory even
	 * when it is farther away in ordinary Euclidean distance.
	 * <p>
	 * Unlike an <b>additively-weighted Voronoi diagram</b> (Apollonius diagram),
	 * which typically yields <b>curved</b> boundaries, power diagrams use <i>power
	 * distance</i> (squared distance with a weight offset), which keeps boundaries
	 * <b>linear</b>.
	 *
	 * @param weightedSites collection of sites encoded as PVectors:
	 *                      {@code (.x, .y)} = position, {@code .z} = weight
	 * @return a GROUP {@link PShape} whose children are the (closed) polygonal
	 *         cells of the power diagram; empty/degenerate cells are omitted
	 * @see #powerDiagram(Collection, double[])
	 * @since 2.2
	 */
	public static PShape powerDiagram(Collection<PVector> weightedSites) {
		return powerDiagram(weightedSites, null);
	}

	/**
	 * Computes a <b>power diagram</b> (a.k.a. <i>Laguerre–Voronoi</i> diagram) for
	 * a set of <b>weighted</b> sites.
	 * <p>
	 * Each site is given as a {@link PVector} where {@code (.x, .y)} is the site
	 * location and {@code .z} is its weight.
	 * <h3>Intuition</h3> A power diagram is the weighted analogue of a standard
	 * Voronoi diagram, but it still produces <b>straight-edged (polygonal)
	 * cells</b>. Conceptually, each site has an associated “strength” (its weight)
	 * that offsets distance: a site with a larger weight can “win” territory even
	 * if it is farther away in ordinary Euclidean terms. Power cells may be empty
	 * (i.e. fewer cells than sites) and may not contain the site.
	 * <p>
	 * Unlike an <b>additively-weighted Voronoi diagram</b> (a.k.a. Apollonius
	 * diagram), where distance is modified by <i>subtracting</i> a radius/weight
	 * and cell boundaries are typically <b>curved</b> (circular arcs), the power
	 * diagram uses <i>power distance</i> (squared distance with a weight offset),
	 * which keeps boundaries <b>linear</b> and cells convex.
	 * <p>
	 * Note: in practice, weights often need to differ substantially in magnitude
	 * (roughly on the order of ~100×) before the effect is visually obvious.
	 *
	 * @param weightedSites collection of sites encoded as PVectors:
	 *                      {@code (.x, .y)} = position, {@code .z} = weight
	 * @param bounds        optional clipping bounds as
	 *                      {@code [minX, minY, maxX, maxY]}. If {@code null}, the
	 *                      diagram is left unclipped.
	 * @return a GROUP {@link PShape} whose children are the (closed) polygonal
	 *         cells of the power diagram; empty/degenerate cells are omitted
	 * @since 2.2
	 * @see #powerDiagram(Collection)
	 */
	public static PShape powerDiagram(Collection<PVector> weightedSites, @Nullable double[] bounds) {
		// NOTE r^2
		var sites = weightedSites.stream().map(z -> new PowerDiagram2D.Site(z.x, z.y, z.z * z.z)).toList();
		final Rect r = bounds == null ? null : new Rect(bounds[0], bounds[1], bounds[2], bounds[3]);
		var cells = PowerDiagram2D.computeCells(sites, r);
		var faces = cells.stream().map(cell -> {
			var points = cell.polygon().stream().map(q -> new PVector((float) q.x(), (float) q.y())).collect(Collectors.toList());
			if (!points.get(0).equals(points.get(points.size() - 1))) {
				points.add(points.get(0)); // unclosed by default, so close
			}
			return PGS_Conversion.fromPVector(points);
		}).filter(Objects::nonNull).toList();

		return PGS_Conversion.flatten(faces);
	}

	/**
	 * Computes a <b>Manhattan (L1) Voronoi diagram</b> for a set of sites,
	 * optionally clipped to an axis-aligned bounding box.
	 * <p>
	 * In a Manhattan Voronoi diagram, distance is measured using the <i>L1</i>
	 * (a.k.a. “city-block” or “taxicab”) metric. Each output cell contains the
	 * points for which a given site is the <b>nearest</b> site under this metric
	 * (ties may occur along cell boundaries).
	 * <p>
	 * If {@code bounds} is {@code null}, clipping bounds are computed automatically
	 * from the input sites using their axis-aligned envelope (i.e. the min/max of
	 * {@code x} and {@code y}). Note that this envelope is often a tight fit; if
	 * you want visible “infinite” outer cells, pass an expanded bounding box.
	 * <p>
	 * Compared to a standard (Euclidean/L2) Voronoi diagram, Manhattan Voronoi
	 * cells tend to align with the coordinate axes and produce characteristic
	 * 45°/axis- aligned edges.
	 *
	 * @param sites  collection of {@link PVector} sites (only {@code x} and
	 *               {@code y} are used)
	 * @param bounds optional clipping bounds as {@code [minX, minY, maxX, maxY]}
	 *               defining the axis-aligned rectangle to which the diagram is
	 *               restricted. If {@code null}, bounds are derived from the sites'
	 *               envelope.
	 * @return a {@link PShape} representing the (optionally clipped) Manhattan
	 *         Voronoi cells (a GROUP shape whose children are polygonal regions)
	 * @since 2.2
	 */
	public static PShape manhattanVoronoi(Collection<PVector> sites, @Nullable double[] bounds) {
		var coords = sites.stream().map(PGS::coordFromPVector).toList();
		Envelope e;
		if (bounds == null) {
			var mp = GEOM_FACTORY.createMultiPointFromCoords(coords.toArray(Coordinate[]::new));
			e = mp.getEnvelopeInternal();
		} else {
			e = new Envelope(bounds[0], bounds[2], bounds[1], bounds[3]);
		}
		var vSites = ManhattanVoronoi.generate(coords, e, false);

		var cells = vSites.stream().map(s -> s.toPolygon(GEOM_FACTORY)).toList();
		return toPShape(cells);
	}

	static Polygon toPolygon(ThiessenPolygon polygon) {
		Coordinate[] coords = new Coordinate[polygon.getEdges().size() + 1];
		int i = 0;
		for (IQuadEdge e : polygon.getEdges()) {
			coords[i++] = new Coordinate(e.getA().x, e.getA().y);
		}
		coords[i] = new Coordinate(polygon.getEdges().get(0).getA().x, polygon.getEdges().get(0).getA().y); // close polygon

		Polygon p = GEOM_FACTORY.createPolygon(coords);
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

	/**
	 * Collects coordinate sets into groups, handling nested GeometryCollections
	 * uniformly.
	 */
	private static void collectVertexGroups(Geometry geom, List<List<Vertex>> groups, List<Vertex> allVertices) {
		if (geom == null || geom.isEmpty()) {
			return;
		}

		// GeometryCollection covers MultiPoint/MultiLineString/MultiPolygon and more.
		if (geom instanceof GeometryCollection gc) {
			for (int i = 0; i < gc.getNumGeometries(); i++) {
				collectVertexGroups(gc.getGeometryN(i), groups, allVertices);
			}
			return;
		}

		// For Polygon/LineString/LinearRing/Point etc.
		List<Vertex> featureVertices = toVertex(geom.getCoordinates());
		if (!featureVertices.isEmpty()) {
			groups.add(featureVertices);
			allVertices.addAll(featureVertices);
		}
	}
}
