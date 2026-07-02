package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.function.UnaryOperator;

import org.apache.commons.lang3.ArrayUtils;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.GeometryFilter;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.util.GeometryTransformer;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.Noder;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.snapround.SnapRoundingNoder;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.tinspin.index.IndexConfig;
import org.tinspin.index.kdtree.KDTree;

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
	 * Precision model that's suitable to guarantee conformity under float
	 * coordinates.
	 */
	public static final PrecisionModel PM = new PrecisionModel(1024); // grid = 1/1024 == Math.ulp(1e4f)

	/**
	 * PGS global geometry factory.
	 */
	public static final GeometryFactory GEOM_FACTORY = new GeometryFactory(PM);

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
		return new Coordinate(p.x, p.y, p.z);
	}

	static final Coordinate[] toCoords(final Collection<PVector> points) {
		CoordinateList coords = new CoordinateList();
		points.forEach(p -> coords.add(coordFromPVector(p)));
		return coords.toCoordinateArray();
	}

	static final PVector toPVector(Coordinate c) {
		return new PVector((float) c.x, (float) c.y, (float) c.z);
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
	 * For Processing's y-axis down system, a negative area from the shoelace
	 * formula in the function's logic will actually correspond to what is visually
	 * counter-clockwise in Processing.
	 */
	static boolean isClockwise(final List<PVector> points) {
		if (points == null || points.size() < 3) {
			throw new IllegalArgumentException("Polygon must have at least 3 points.");
		}

		double area = 0;
		int n = points.size();

		for (int i = 0; i < n; i++) {
			PVector p1 = points.get(i);
			PVector p2 = points.get((i + 1) % n);
			area += (p1.x * p2.y - p2.x * p1.y);
		}

		return area < 0; // negative area means clockwise (in standard y-up system)
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
	static final PShape polygonizeSegments(Collection<? extends SegmentString> segments, boolean node) {
		if (node) {
			segments = nodeSegmentStrings(segments);
		}
		final List<PEdge> meshEdges = new ArrayList<>(segments.size());
		segments.forEach(ss -> { // ss is not necessarily a single edge (can be many connected edges)
			for (int i = 0; i < ss.size() - 1; i++) {
				meshEdges.add(new PEdge(toPVector(ss.getCoordinate(i)), toPVector(ss.getCoordinate(i + 1))));
			}
		});
//		Collections.shuffle(meshEdges);
		return polygonizeNodedEdges(meshEdges);
	}

	/**
	 * Given a (possibly non‐noded) set of PEdge’s, optionally nodes all
	 * intersections, and polygonizes. Does NOT convert back into PEdge; it goes
	 * straight from noded SegmentStrings → JTS LineStrings → Polygonizer → PShape.
	 *
	 * @param edges the input edges
	 * @param node  if true, splits at all interior intersections
	 * @return a GROUP PShape whose children are each polygon face
	 */
	static final PShape polygonizeEdges(Collection<PEdge> edges) {
		List<NodedSegmentString> segStrs = new HashSet<>(edges).stream().map(e -> {
			Coordinate c0 = coordFromPVector(e.a);
			Coordinate c1 = coordFromPVector(e.b);
			return new NodedSegmentString(new Coordinate[] { c0, c1 }, null);
		}).toList();

		Collection<SegmentString> noded = nodeSegmentStrings(segStrs);

		Polygonizer polygonizer = new Polygonizer();
		for (SegmentString ss : noded) {
			var coords = ss.getCoordinates();
			for (int i = 0; i < coords.length - 1; i++) {
				Coordinate p0 = coords[i];
				Coordinate p1 = coords[i + 1];
				LineString ls = GEOM_FACTORY.createLineString(new Coordinate[] { p0, p1 });
				polygonizer.add(ls);
			}
		}

		@SuppressWarnings("unchecked")
		Collection<Polygon> polys = polygonizer.getPolygons();
		return PGS_Conversion.toPShape(polys);
	}

	/**
	 * Polygonizes a set of pre-noded edges.
	 * 
	 * @param edges a collection of NODED (i.e. non intersecting / must only meet at
	 *              their endpoints) edges. The collection can contain duplicates.
	 * @return a GROUP PShape, where each child shape represents a polygon face
	 *         formed by the given edges
	 */
	static final PShape polygonizeNodedEdges(Collection<PEdge> edges) {
		return polygonizeEdgesRobust(edges);
	}

	/**
	 * Polygonizes a set of edges using JTS Polygonizer (occasionally
	 * FastPolygonizer is not robust enough).
	 * 
	 * @param edges a collection of NODED (i.e. non intersecting / must only meet at
	 *              their endpoints) edges. The collection can contain duplicates.
	 * @return a GROUP PShape, where each child shape represents a polygon face
	 *         formed by the given edges
	 */
	@SuppressWarnings("unchecked")
	private static final PShape polygonizeEdgesRobust(Collection<PEdge> edges) {
		final Set<PEdge> edgeSet = new HashSet<>(edges);
		final Polygonizer polygonizer = new Polygonizer();
		// polygonizer.setCheckRingsValid(false);
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

		List<Polygon> polys = (List<Polygon>) polygonizer.getPolygons();
		// polygonizer preserves poly order but not vertex order
		polys.forEach(poly -> poly.normalize());
		return PGS_Conversion.toPShape(polys);
	}

	/**
	 * Post-processes polygons produced by JTS
	 * {@link org.locationtech.jts.operation.polygonize.Polygonizer Polygonizer} so
	 * that nested rings are interpreted as holes of their enclosing polygon.
	 * <p>
	 * This is necessary because {@code Polygonizer} returns <em>all</em> bounded
	 * faces implied by the input linework. For example, when a ring lies inside
	 * another ring, {@code Polygonizer} will typically produce both the enclosing
	 * polygon-with-hole <em>and</em> the inner “hole face” as a standalone polygon.
	 * This method removes those hole faces by classifying faces by nesting depth
	 * (odd depth = hole, even depth = filled).
	 * </p>
	 *
	 * @param polygonizerFaces polygons returned by {@code Polygonizer}.
	 * @param dissolve         if {@code true}, unions the kept faces into a
	 *                         dissolved geometry; if {@code false}, returns a
	 *                         {@link GeometryCollection} of the kept faces.
	 * @return a geometry containing only the “filled” faces, with holes inferred
	 *         from nesting.
	 */
	static Geometry dropHolePolygons(List<Polygon> polygonizerFaces, boolean dissolve) {
		// method could be optimised, but shouldn't need used much
		if (polygonizerFaces == null || polygonizerFaces.isEmpty()) {
			return new GeometryFactory().createGeometryCollection();
		}

		record Face(Polygon face, Polygon shellOnly, double area) {
		}

		// Build faces with "shell-only" geometry (exterior ring only)
		List<Face> faces = new ArrayList<>(polygonizerFaces.size());
		for (Polygon p : polygonizerFaces) {
			if (p == null || p.isEmpty()) {
				continue;
			}
			LinearRing shell = p.getExteriorRing();
			Polygon shellOnly = GEOM_FACTORY.createPolygon(shell, null);
			faces.add(new Face(p, shellOnly, shellOnly.getArea()));
		}

		// Sort by increasing area to find smallest containing parent efficiently (still
		// O(n^2))
		faces.sort(Comparator.comparingDouble(Face::area));

		int n = faces.size();
		int[] parent = new int[n];
		Arrays.fill(parent, -1);

		// Parent of i = smallest-area shell that covers a point inside i's shell
		for (int i = 0; i < n; i++) {
			Point testPt = faces.get(i).shellOnly().getInteriorPoint();
			for (int j = i + 1; j < n; j++) {
				if (faces.get(j).shellOnly().covers(testPt)) {
					parent[i] = j;
					break;
				}
			}
		}

		// Keep even-depth faces (filled), drop odd-depth faces (holes)
		List<Polygon> kept = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			int depth = 0;
			for (int p = parent[i]; p != -1; p = parent[p]) {
				depth++;
			}
			if ((depth & 1) == 0) {
				kept.add(faces.get(i).face());
			}
		}

		if (!dissolve) {
			return GEOM_FACTORY.createGeometryCollection(kept.toArray(Geometry[]::new));
		}
		return kept.isEmpty() ? GEOM_FACTORY.createGeometryCollection() : UnaryUnionOp.union(kept);
	}

	/**
	 * Computes a robust noding for a collection of SegmentStrings.
	 * 
	 * @param segments
	 * @return
	 */
	@SuppressWarnings("unchecked")
	static final Collection<SegmentString> nodeSegmentStrings(Collection<? extends SegmentString> segments) {
		/*
		 * Other noder implementations do not node correctly (fail to detect
		 * intersections) on many inputs; furthermore, using a very small tolerance
		 * (i.e. ~1e-10) on SnappingNoder noder on a small tolerance misses
		 * intersections too (hence 1/1024 chosen as suitable). "Noding robustness
		 * issues are generally caused by nearly coincident line segments, or by very
		 * short line segments. Snapping mitigates both of these situations.".
		 */
		Noder noder = new SnapRoundingNoder(PGS.PM);
		noder.computeNodes(segments);
		return noder.getNodedSubstrings();
	}

	static final <T> HashSet<T> makeHashSet(int expectedSize) {
		// required capacity = actual_capacity / fill_factor + 1 (to avoid rehashing)
		return new HashSet<>((int) ((expectedSize) / 0.75 + 1));
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
		List<Polygon> polygons = new ArrayList<>();
		g.apply((GeometryFilter) geom -> {
			if (geom instanceof Polygon) {
				polygons.add((Polygon) geom);
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
	 * Creates a 2D KDTree populated with <code>points</code>.
	 */
	static KDTree<PVector> makeKdtree(Collection<PVector> points) {
		KDTree<PVector> tree = KDTree.create(IndexConfig.create(2).setDefensiveKeyCopy(false));
		points.forEach(p -> {
			tree.insert(new double[] { p.x, p.y }, p);
		});
		return tree;
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
				// if (poly.getNumPoints() == 0) {
				// continue;
				// }
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
	 * Apply a transformation to every lineal element in a {@code PShape},
	 * preserving geometry structure and polygon/hole relationships, and return a
	 * non-null result.
	 *
	 * <p>
	 * The geometry encoded by {@code shape} (via {@code fromPShape}) is traversed
	 * and {@code function} is applied to each lineal component: {@code LineString}
	 * and polygon rings ({@code LinearRing}, passed to the function as a
	 * {@code LineString}). The function may return a replacement
	 * {@code LineString}, or {@code null} to drop that element.
	 *
	 * <h3>Structure preservation</h3>
	 * <ul>
	 * <li><b>GeometryCollection / MultiPolygon / MultiLineString</b>
	 * <ul>
	 * <li>Children are processed in index order; the relative order of surviving
	 * children is preserved.</li>
	 * <li>Children for which the function yields {@code null} (or that become
	 * empty) are omitted from the result.</li>
	 * <li>If the <i>input</i> encodes a multi/collection geometry, the returned
	 * {@code PShape} is always of kind {@code GROUP} (it may be empty if nothing
	 * survives), even if only a single child remains after filtering.</li>
	 * </ul>
	 * </li>
	 *
	 * <li><b>Polygon</b>
	 * <ul>
	 * <li>Rings are visited shell-first (exterior, then holes in interior-ring
	 * index order), preserving exterior–hole relationships.</li>
	 * <li>If the exterior ring is dropped or becomes invalid, the entire polygon is
	 * dropped.</li>
	 * <li>Holes that are dropped or become invalid are omitted; remaining holes
	 * retain their original order.</li>
	 * <li>Ring orientation is enforced: exterior is clockwise (CW); holes are
	 * counter-clockwise (CCW).</li>
	 * </ul>
	 * </li>
	 *
	 * <li><b>LinearRing</b>
	 * <ul>
	 * <li>If a {@code LinearRing} is encountered outside a polygon, it is treated
	 * as an exterior ring for closure/orientation rules.</li>
	 * </ul>
	 * </li>
	 * </ul>
	 *
	 * <h3>Additional behavior</h3>
	 * <ul>
	 * <li>Non-closed ring outputs are closed when possible (if at least two points
	 * exist).</li>
	 * <li>Rings must have at least 4 coordinates (including repeated first/last)
	 * after closing; otherwise they are dropped.</li>
	 * <li>{@code LineString} elements return the transformed line, or are dropped
	 * if {@code function} returns {@code null}.</li>
	 * <li>Unsupported geometry types are ignored (dropped). If the root geometry is
	 * unsupported, an empty {@code PShape} is returned.</li>
	 * <li>No full topology validation is performed; run JTS validators if
	 * needed.</li>
	 * </ul>
	 *
	 * <h3>Return contract</h3>
	 * <ul>
	 * <li>This method never returns {@code null}. If no geometry survives, an empty
	 * {@code PShape} is returned (for multi/collection inputs, an empty
	 * {@code GROUP} {@code PShape}).</li>
	 * </ul>
	 *
	 * @param shape    input {@code PShape} encoding geometries to transform (must
	 *                 be convertible via {@code fromPShape})
	 * @param function operator applied to each {@code LineString}; polygon rings
	 *                 are passed as {@code LineString}. Returning {@code null}
	 *                 drops that element.
	 * @return a non-null {@code PShape} representing the transformed geometry
	 * @since 2.1
	 */
	static PShape applyToLinealGeometries(PShape shape, UnaryOperator<LineString> fn) {
		final Geometry in = fromPShape(shape);

		if (in instanceof Point || in instanceof MultiPoint) {
			return new PShape();
		}

		final boolean rootIsMultiPolygon = in instanceof MultiPolygon;
		final boolean rootIsMultiLineString = in instanceof MultiLineString;
		final boolean rootIsGeomCollection = (in instanceof GeometryCollection) && !rootIsMultiPolygon && !rootIsMultiLineString && !(in instanceof MultiPoint);

		final Object rootUserData = in.getUserData();

		Geometry out = new PGS_Transformer(fn).transform(in);

		// Never return null; match empty policies
		if (out == null || out.isEmpty()) {
			return new PShape(PConstants.GROUP);
		}

		// Preserve "GROUP-ness" for multi/collection roots even if only one child
		// survives
		if (rootIsMultiPolygon && out instanceof Polygon p) {
			out = GEOM_FACTORY.createMultiPolygon(new Polygon[] { p });
		} else if (rootIsMultiLineString && out instanceof LineString ls && !(out instanceof MultiLineString)) {
			out = GEOM_FACTORY.createMultiLineString(new LineString[] { ls });
		} else if (rootIsGeomCollection && !(out instanceof GeometryCollection)) {
			out = GEOM_FACTORY.createGeometryCollection(new Geometry[] { out });
		}

		out.setUserData(rootUserData);
		return toPShape(out);
	}

	static boolean isEmptyShape(PShape s) {
		if (s == null) {
			return true;
		}
		if (s.getChildCount() > 0) {
			return false;
		}
		if (s.getVertexCount() > 0) {
			return false;
		}
		return true;
	}

	private static class PGS_Transformer extends GeometryTransformer {

		private final UnaryOperator<LineString> fn;

		PGS_Transformer(UnaryOperator<LineString> fn) {
			this.fn = fn;
		}

		@Override
		protected Geometry transformPolygon(Polygon p, Geometry parent) {
			// Own the polygon traversal order: shell first, then holes by index.
			LinearRing shell = processRing(p.getExteriorRing(), false);
			if (shell == null) {
				return null; // drop whole polygon
			}

			List<LinearRing> holes = new ArrayList<>(p.getNumInteriorRing());
			for (int i = 0; i < p.getNumInteriorRing(); i++) {
				LinearRing h = processRing(p.getInteriorRingN(i), true);
				if (h != null) {
					holes.add(h);
				}
			}

			Polygon out = GEOM_FACTORY.createPolygon(shell, holes.toArray(LinearRing[]::new));
			out.setUserData(p.getUserData());
			return out;
		}

		@Override
		protected Geometry transformLinearRing(LinearRing ring, Geometry parent) {
			// Standalone rings: treat as exterior policy (CW)
			LinearRing out = processRing(ring, false);
			if (out != null) {
				out.setUserData(ring.getUserData());
			}
			return out;
		}

		@Override
		protected Geometry transformLineString(LineString ls, Geometry parent) {
			// Note: GeometryTransformer may route rings here too; ensure we handle them as
			// rings.
			if (ls instanceof LinearRing r) {
				return transformLinearRing(r, parent);
			}

			LineString res = fn.apply(ls);
			if (res == null || res.isEmpty()) {
				return null;
			}

			LineString out = GEOM_FACTORY.createLineString(res.getCoordinateSequence());
			out.setUserData(ls.getUserData());
			return out;
		}

		@Override
		protected Geometry transformGeometryCollection(GeometryCollection gc, Geometry parent) {
			// Preserve order; filter null/empty; preserve container type
			List<Geometry> kept = new ArrayList<>(gc.getNumGeometries());
			for (int i = 0; i < gc.getNumGeometries(); i++) {
				Geometry t = transform(gc.getGeometryN(i));
				if (t != null && !t.isEmpty()) {
					kept.add(t);
				}
			}

			if (gc instanceof MultiPolygon) {
				List<Polygon> polys = new ArrayList<>();
				for (Geometry g : kept) {
					if (g instanceof Polygon p) {
						polys.add(p);
					}
				}
				return GEOM_FACTORY.createMultiPolygon(polys.toArray(Polygon[]::new));
			}

			if (gc instanceof MultiLineString) {
				List<LineString> lines = new ArrayList<>();
				for (Geometry g : kept) {
					if (g instanceof LineString ls) {
						lines.add(ls);
					}
				}
				return GEOM_FACTORY.createMultiLineString(lines.toArray(LineString[]::new));
			}

			return GEOM_FACTORY.createGeometryCollection(kept.toArray(Geometry[]::new));
		}

		private LinearRing processRing(LinearRing ring, boolean isHole) {
			// Apply fn to ring (passed as LineString)
			LineString res = fn.apply(ring);
			if (res == null || res.isEmpty()) {
				return null;
			}

			Coordinate[] coords = res.getCoordinates();

			// Ensure closed when possible
			if (coords.length >= 2 && !coords[0].equals2D(coords[coords.length - 1])) {
				coords = Arrays.copyOf(coords, coords.length + 1);
				coords[coords.length - 1] = coords[0];
			}

			// Need at least 4 coordinates for a valid ring
			if (coords.length < 4) {
				return null;
			}

			// Enforce orientation: exterior CW, holes CCW
			boolean ccw = Orientation.isCCWArea(coords);
			if (isHole && !ccw) {
				ArrayUtils.reverse(coords);
			}
			if (!isHole && ccw) {
				ArrayUtils.reverse(coords);
			}

			return GEOM_FACTORY.createLinearRing(coords);
		}
	}

}
