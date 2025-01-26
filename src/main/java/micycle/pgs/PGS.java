package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static processing.core.PConstants.GROUP;
import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.function.UnaryOperator;

import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.GeometryFilter;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.Noder;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.snapround.SnapRoundingNoder;
import org.locationtech.jts.operation.linemerge.LineMerger;
import org.locationtech.jts.operation.polygonize.Polygonizer;

import micycle.pgs.color.Colors;
import micycle.pgs.commons.Nullable;
import micycle.pgs.commons.PEdge;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * This class houses functions used by the library internally.
 * 
 * @author Michael Carleton
 */
final class PGS {

	/** Defines number of vertices that comprise constructed geometries. */
	static final int SHAPE_SAMPLES = 80;

	/**
	 * PGS global geometry factory (uses 32 bit float precision).
	 */
	public static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	private PGS() {
	}

	/**
	 * Create a LINES PShape, ready for vertices (shape.vertex(x, y) calls).
	 * 
	 * @param strokeColor  nullable (default = {@link Colors#PINK})
	 * @param strokeCap    nullable (default = <code>ROUND</code>)
	 * @param strokeWeight nullable (default = <code>2</code>)
	 * @return LINES PShape ready for vertex calls
	 */
	static final PShape prepareLinesPShape(@Nullable Integer strokeColor, @Nullable Integer strokeCap, @Nullable Integer strokeWeight) {
		if (strokeColor == null) {
			strokeColor = Colors.PINK;
		}
		if (strokeCap == null) {
			strokeCap = ROUND;
		}
		if (strokeWeight == null) {
			strokeWeight = 2;
		}
		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(strokeCap);
		lines.setStroke(true);
		lines.setStrokeWeight(strokeWeight);
		lines.setStroke(strokeColor);
		lines.beginShape(LINES);
		return lines;
	}

	/**
	 * Euclidean distance between two points
	 */
	static final double distance(Point a, Point b) {
		double deltaX = a.getX() - b.getX();
		double deltaY = a.getY() - b.getY();
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	static final double distanceSq(PVector a, PVector b) {
		float dx = a.x - b.x;
		float dy = a.y - b.y;
		return (dx * dx + dy * dy);
	}

	static final LineString createLineString(PVector a, PVector b) {
		return GEOM_FACTORY.createLineString(new Coordinate[] { coordFromPVector(a), coordFromPVector(b) });
	}

	static final SegmentString createSegmentString(PVector a, PVector b) {
		return new NodedSegmentString(new Coordinate[] { PGS.coordFromPVector(a), PGS.coordFromPVector(b) }, null);
	}

	static final Point createPoint(double x, double y) {
		return GEOM_FACTORY.createPoint(new Coordinate(x, y));
	}

	/**
	 * Creates a stroked rectangle.
	 */
	static final PShape createRect(double x, double y, double w, double h) {
		final PShape rect = new PShape(PShape.PATH);
		rect.setFill(true);
		rect.setFill(255);
		rect.setStroke(true);
		rect.setStrokeWeight(4);
		rect.setStroke(Colors.PINK);
		rect.beginShape();
		rect.vertex((float) x, (float) y);
		rect.vertex((float) (x + w), (float) y);
		rect.vertex((float) (x + w), (float) (y + h));
		rect.vertex((float) x, (float) (y + h));
		rect.endShape(PConstants.CLOSE);
		return rect;
	}

	static final Point pointFromPVector(final PVector p) {
		return GEOM_FACTORY.createPoint(new Coordinate(p.x, p.y));
	}

	static final Coordinate coordFromPoint(final Point p) {
		return new Coordinate(p.getX(), p.getY());
	}

	static final Coordinate coordFromPVector(final PVector p) {
		return new Coordinate(p.x, p.y);
	}

	static final Coordinate[] toCoords(final Collection<PVector> points) {
		CoordinateList coords = new CoordinateList();
		points.forEach(p -> coords.add(coordFromPVector(p)));
		return coords.toCoordinateArray();
	}

	static final PVector toPVector(Coordinate c) {
		return new PVector((float) c.x, (float) c.y);
	}

	/**
	 * Reflection-based workaround to get the fill color of a PShape (this field is
	 * usually inaccessible).
	 */
	static final int getPShapeFillColor(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("fillColor");
			f.setAccessible(true);
			return f.getInt(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * Reflection-based workaround to get the stroke color of a PShape (this field
	 * is usually inaccessible).
	 */
	static final int getPShapeStrokeColor(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("strokeColor");
			f.setAccessible(true);
			return f.getInt(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * Reflection-based workaround to get the stroke strokeWeight of a PShape (this
	 * field is usually inaccessible).
	 */
	static final float getPShapeStrokeWeight(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("strokeWeight");
			f.setAccessible(true);
			return f.getFloat(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * Requires a closed hole
	 * 
	 * @param points
	 * @return
	 */
	static final boolean isClockwise(List<PVector> points) {
		boolean closed = true;
		if (points.get(0).equals(points.get(points.size() - 1))) {
			closed = false;
			points.add(points.get(0)); // mutate list
		}
		double area = 0;

		for (int i = 0; i < (points.size()); i++) {
			int j = (i + 1) % points.size();
			area += points.get(i).x * points.get(j).y;
			area -= points.get(j).x * points.get(i).y;
		}

		if (!closed) {
			points.remove(points.size() - 1); // revert mutation
		}

		return (area < 0);
	}

	/**
	 * Nodes (optional) then polygonizes a set of line segments.
	 * 
	 * @param segments list of segments (noded or non-noded)
	 * @param node     whether to node the segments before polygonization. If the
	 *                 segments constitute a conforming mesh, then set this as
	 *                 false; otherwise true.
	 * @return a GROUP PShape, where each child shape represents a polygon face
	 *         formed by the given edges
	 */
	static final PShape polygonizeSegments(Collection<SegmentString> segments, boolean node) {
		if (node) {
			segments = nodeSegmentStrings(segments);
		}
		final List<PEdge> meshEdges = new ArrayList<>(segments.size());
		segments.forEach(ss -> { // ss is not necessarily a single edge (can be many connected edges)
			for (int i = 0; i < ss.size() - 1; i++) {
				meshEdges.add(new PEdge(toPVector(ss.getCoordinate(i)), toPVector(ss.getCoordinate(i + 1))));
			}
		});
		Collections.shuffle(meshEdges);
		return polygonizeEdges(meshEdges);
	}

	/**
	 * Polygonizes a set of edges.
	 * 
	 * @param edges a collection of NODED (i.e. non intersecting / must only meet at
	 *              their endpoints) edges. The collection can contain duplicates.
	 * @return a GROUP PShape, where each child shape represents a polygon face
	 *         formed by the given edges
	 */
	static final PShape polygonizeEdges(Collection<PEdge> edges) {
		return polygonizeEdgesRobust(edges);
	}

	/**
	 * Polygonizes a set of edges using JTS Polygonizer (occasionally
	 * FastPolygonizer is not robust enough).
	 * 
	 * @param edges a collection of NODED (i.e. non intersecting / must onlymeet at
	 *              their endpoints) edges. The collection can containduplicates.
	 * @return a GROUP PShape, where each child shape represents a polygon face
	 *         formed by the given edges
	 */
	@SuppressWarnings("unchecked")
	private static final PShape polygonizeEdgesRobust(Collection<PEdge> edges) {
		final Set<PEdge> edgeSet = new HashSet<>(edges);
		final Polygonizer polygonizer = new Polygonizer();
		polygonizer.setCheckRingsValid(false);
		edgeSet.forEach(ss -> {
			/*
			 * NOTE: If the same LineString is added more than once to the polygonizer, the
			 * string is "collapsed" and not counted as an edge. Therefore a set is used to
			 * ensure strings are added once only to the polygonizer. A PEdge is used to
			 * determine this (since LineString hashcode doesn't work).
			 */
			final LineString l = createLineString(ss.a, ss.b);
			polygonizer.add(l);
		});
		return PGS_Conversion.toPShape(polygonizer.getPolygons());
	}

	/**
	 * Computes a robust noding for a collection of SegmentStrings.
	 * 
	 * @param segments
	 * @return
	 */
	@SuppressWarnings("unchecked")
	static final Collection<SegmentString> nodeSegmentStrings(Collection<SegmentString> segments) {
		/*
		 * Other noder implementations do not node correctly (fail to detect
		 * intersections) on many inputs; furthermore, using a very small tolerance
		 * (i.e. ~1e-10) on SnappingNoder noder on a small tolerance misses
		 * intersections too (hence 0.01 chosen as suitable). "Noding robustness issues
		 * are generally caused by nearly coincident line segments, or by very short
		 * line segments. Snapping mitigates both of these situations.".
		 */
		Noder noder = new SnapRoundingNoder(new PrecisionModel(-5e-3));
		noder.computeNodes(segments);
		return noder.getNodedSubstrings();
	}

	static final <T> HashSet<T> makeHashSet(int expectedSize) {
		// required capacity = actual_capacity / fill_factor + 1 (to avoid rehashing)
		return new HashSet<>((int) ((expectedSize) / 0.75 + 1));
	}

	/**
	 * Computes an <b>ordered</b> list of <b>vertices</b> that make up the boundary
	 * of a polygon from an <b>unordered</b> collection of <b>edges</b>. The
	 * underlying approach is around ~10x faster than JTS .buffer(0) and ~3x faster
	 * than {@link LineMerger}.
	 * <p>
	 * For now, this method does not properly support multi-shapes, nor unclosed
	 * edge collections (that form unclosed linestrings).
	 * <p>
	 * Notably, unlike {@link LineMerger} this approach does not merge successive
	 * boundary segments that together form a straight line into a single longer
	 * segment.
	 * 
	 * @param edges unordered/random collection of edges (containing no duplicates),
	 *              that together constitute the boundary of a single polygon / a
	 *              closed ring
	 * @return list of sequential vertices belonging to the polygon that follow some
	 *         constant winding (may wind clockwise or anti-clockwise). Note: this
	 *         vertex list is not closed (having same start and end vertex) by
	 *         default!
	 */
	static List<PVector> fromEdges(Collection<PEdge> edges) {
		// NOTE same as org.locationtech.jts.operation.linemerge.LineSequencer ?
		// map of vertex to the 2 edges that share it
		final HashMap<PVector, HashSet<PEdge>> vertexEdges = new HashMap<>((int) ((edges.size()) / 0.75 + 1));

		/*
		 * Build up map of vertex->edge to later find edges sharing a given vertex in
		 * O(1). When the input is valid (edges form a closed loop) every vertex is
		 * shared by 2 edges.
		 */
		for (PEdge e : edges) {
			if (vertexEdges.containsKey(e.a)) {
				vertexEdges.get(e.a).add(e);
			} else {
				HashSet<PEdge> h = new HashSet<>();
				h.add(e);
				vertexEdges.put(e.a, h);
			}
			if (vertexEdges.containsKey(e.b)) {
				vertexEdges.get(e.b).add(e);
			} else {
				HashSet<PEdge> h = new HashSet<>();
				h.add(e);
				vertexEdges.put(e.b, h);
			}
		}

		List<PVector> vertices = new ArrayList<>(edges.size() + 1); // boundary vertices

		// begin by choosing a random edge
		final PEdge startingEdge = edges.iterator().next();
		vertices.add(startingEdge.a);
		vertices.add(startingEdge.b);
		vertexEdges.get(startingEdge.a).remove(startingEdge);
		vertexEdges.get(startingEdge.b).remove(startingEdge);

		while (vertices.size() < edges.size()) {
			final PVector lastVertex = vertices.get(vertices.size() - 1);
			Set<PEdge> connectedEdges = vertexEdges.get(lastVertex);

			if (connectedEdges.isEmpty()) {
				/*
				 * This will be hit if the input is malformed (contains multiple disjoint shapes
				 * for example), and break when the first loop is closed. On valid inputs the
				 * while loop will break before this statement can be hit.
				 */
				break;
			}

			final PEdge nextEdge = connectedEdges.iterator().next();
			if (nextEdge.a.equals(lastVertex)) {
				vertices.add(nextEdge.b);
				vertexEdges.get(nextEdge.b).remove(nextEdge);
			} else {
				vertices.add(nextEdge.a);
				vertexEdges.get(nextEdge.a).remove(nextEdge);
			}
			connectedEdges.remove(nextEdge); // remove this edge from vertex mapping
			if (connectedEdges.isEmpty()) {
				vertexEdges.remove(lastVertex); // have used both edges connected to this vertex -- now remove!
			}
		}

		return vertices;
	}

	static SimpleWeightedGraph<PVector, PEdge> makeCompleteGraph(List<PVector> points) {
		SimpleWeightedGraph<PVector, PEdge> graph = new SimpleWeightedGraph<>(PEdge.class);

		// Add all vertices before starting the edge creation process
		for (PVector vertex : points) {
			graph.addVertex(vertex);
		}

		// Create edges between all pairs of vertices
		for (int i = 0; i < points.size(); i++) {
			PVector a = points.get(i);
			for (int j = i + 1; j < points.size(); j++) {
				PVector b = points.get(j);
				PEdge e = new PEdge(a, b);
				graph.addEdge(a, b, e);
				graph.setEdgeWeight(e, e.length());
			}
		}
		return graph;
	}

	/**
	 * Extracts all the polygons from a given geometry. If the geometry instance is
	 * a MultiPolygon, each individual polygon is extracted and added to the result
	 * list. Other geometry types contained within the input geometry are ignored.
	 */
	static List<Polygon> extractPolygons(Geometry g) {
		List<Polygon> polygons = new ArrayList<>(g.getNumGeometries());
		g.apply((GeometryFilter) geom -> {
			if (geom instanceof Polygon) {
				polygons.add((Polygon) geom);
			} else if (geom instanceof MultiPolygon) {
				for (int i = 0; i < geom.getNumGeometries(); i++) {
					polygons.add((Polygon) geom.getGeometryN(i));
				}
			}
		});
		return polygons;
	}

	/**
	 * Extracts all the LinearRings from a given polygon. This includes both the
	 * exterior ring and all interior rings (if any).
	 */
	static List<LinearRing> extractLinearRings(Polygon polygon) {
		List<LinearRing> rings = new ArrayList<>(1 + polygon.getNumInteriorRing());
		rings.add(polygon.getExteriorRing());
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			rings.add(polygon.getInteriorRingN(i));
		}

		return rings;
	}

	/**
	 * Provides convenient iteration of the child geometries of a JTS MultiGeometry.
	 * This iterator does not recurse all geometries (as does
	 * {@link org.locationtech.jts.geom.GeometryCollectionIterator
	 * GeometryCollectionIterator}), but returns the first level geometries only.
	 * 
	 * @author Michael Carleton
	 */
	static final class GeometryIterator implements Iterable<Geometry> {

		private final Geometry g;

		public GeometryIterator(Geometry g) {
			this.g = g;
		}

		@Override
		public Iterator<Geometry> iterator() {
			return new Iterator<>() {
				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < g.getNumGeometries();
				}

				@Override
				public Geometry next() {
					if (!hasNext()) {
						throw new NoSuchElementException();
					}
					return g.getGeometryN(currentIndex++);
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
		}
	}

	/**
	 * Provides convenient iteration of exterior and linear rings (if any) of a
	 * polygonal JTS geometry. Supports MultiGeometries.
	 * 
	 * @author Michael Carleton
	 */
	static final class LinearRingIterator implements Iterable<LinearRing> {

		private LinearRing[] array;

		/**
		 * Constructs the iterator for the given geometry. The first ring returned by
		 * the iterator is the exterior ring; all other rings (if any) are interior
		 * rings.
		 * 
		 * @param g input geometry
		 */
		public LinearRingIterator(Geometry g) {
			ArrayList<LinearRing> rings = new ArrayList<>(g.getNumGeometries());
			for (int i = 0; i < g.getNumGeometries(); i++) {
				Polygon poly = (Polygon) g.getGeometryN(i);
//				if (poly.getNumPoints() == 0) {
//					continue;
//				}
				rings.add(poly.getExteriorRing());
				for (int j = 0; j < poly.getNumInteriorRing(); j++) {
					rings.add(poly.getInteriorRingN(j));
				}
			}
			array = rings.toArray(new LinearRing[rings.size()]);
		}

		public LinearRing[] getLinearRings() {
			return array;
		}

		@Override
		public Iterator<LinearRing> iterator() {
			return new Iterator<>() {

				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < array.length;
				}

				@Override
				public LinearRing next() {
					if (!hasNext()) {
						throw new NoSuchElementException();
					}
					return array[currentIndex++];
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
		}
	}

	/**
	 * Applies a given function to all lineal elements (LineString and LinearRing)
	 * in the input PShape. The function processes each LineString or LinearRing and
	 * returns a modified LineString. The method preserves the structure of the
	 * input geometry, including exterior-hole relationships in polygons and
	 * multi-geometries.
	 *
	 * @param shape    The input PShape containing lineal elements to process.
	 * @param function A UnaryOperator that takes a LineString or LinearRing and
	 *                 returns a modified LineString.
	 * @return A new PShape with the processed lineal elements. Returns an empty
	 *         PShape for unsupported geometry types.
	 * @since 2.1
	 */
	static PShape applyToLinealGeometries(PShape shape, UnaryOperator<LineString> function) {
		Geometry g = fromPShape(shape);
		switch (g.getGeometryType()) {
			case Geometry.TYPENAME_GEOMETRYCOLLECTION :
			case Geometry.TYPENAME_MULTIPOLYGON :
			case Geometry.TYPENAME_MULTILINESTRING :
				PShape group = new PShape(GROUP);
				for (int i = 0; i < g.getNumGeometries(); i++) {
					group.addChild(applyToLinealGeometries(toPShape(g.getGeometryN(i)), function));
				}
				return group;
			case Geometry.TYPENAME_LINEARRING :
			case Geometry.TYPENAME_POLYGON :
				// Preserve exterior-hole relations
				LinearRing[] rings = new LinearRingIterator(g).getLinearRings();
				LinearRing[] processedRings = new LinearRing[rings.length];
				for (int i = 0; i < rings.length; i++) {
					LinearRing ring = rings[i];
					LineString out = function.apply(ring);
					if (out.isClosed()) {
						processedRings[i] = GEOM_FACTORY.createLinearRing(out.getCoordinates());
					} else {
						System.err.println("Output LineString is not closed. Closing it automatically.");
						Coordinate[] closedCoords = Arrays.copyOf(out.getCoordinates(), out.getNumPoints() + 1);
						closedCoords[closedCoords.length - 1] = closedCoords[0]; // Close the ring
						processedRings[i] = GEOM_FACTORY.createLinearRing(closedCoords);
					}
				}
				if (processedRings.length == 0) {
					return new PShape();
				}
				LinearRing[] holes = null;
				if (processedRings.length > 1) {
					holes = Arrays.copyOfRange(processedRings, 1, processedRings.length);
				}
				return toPShape(GEOM_FACTORY.createPolygon(processedRings[0], holes));
			case Geometry.TYPENAME_LINESTRING :
				LineString l = (LineString) g;
				LineString out = function.apply(l);
				return toPShape(out);
			default :
				return new PShape(); // Return empty (so element is invisible if not processed)
		}
	}

}
