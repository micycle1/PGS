package micycle.pgs;

import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS.prepareLinesPShape;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;
import static micycle.pgs.color.RGB.composeColor;

import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateArrays;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.util.LinearComponentExtracter;
import org.locationtech.jts.geomgraph.Edge;
import org.locationtech.jts.geomgraph.EdgeIntersection;
import org.locationtech.jts.geomgraph.GeometryGraph;
import org.locationtech.jts.geomgraph.index.EdgeSetIntersector;
import org.locationtech.jts.geomgraph.index.MonotoneChainIndexer;
import org.locationtech.jts.geomgraph.index.SegmentIntersector;
import org.locationtech.jts.geomgraph.index.SimpleMCSweepLineIntersector;
import org.locationtech.jts.geomgraph.index.SimpleSweepLineIntersector;
import org.locationtech.jts.index.chain.MonotoneChainBuilder;
import org.locationtech.jts.index.chain.MonotoneChainOverlapAction;
import org.locationtech.jts.index.hprtree.HPRtree;
import org.locationtech.jts.noding.SegmentIntersectionDetector;
import org.locationtech.jts.operation.linemerge.LineMerger;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;
import org.tinfour.common.IIncrementalTinNavigator;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.voronoi.BoundedVoronoiBuildOptions;
import org.tinfour.voronoi.BoundedVoronoiDiagram;
import org.tinfour.voronoi.ThiessenPolygon;
import org.tinspin.index.rtree.Entry;
import org.tinspin.index.rtree.RTree;
import org.tinspin.index.rtree.RTreeIterator;
import org.tinspin.index.rtree.RTreeQueryKnn;

import de.alsclo.voronoi.Voronoi;
import micycle.pgs.color.RGB;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import valenpe7.bentley_ottmann.BentleyOttmann;
import valenpe7.bentley_ottmann.Segment;

/**
 * 
 * Voronoi Diagrams of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_Voronoi {

	private PGS_Voronoi() {

	}

	/**
	 * 
	 * @param shape     the shape whose vertices to use as vornoi sites
	 * @param constrain constrain the diagram lines to the shape? when true, lines
	 *                  are cropped
	 * @param p
	 * @return
	 */
	public static PShape test(PShape shape, boolean constrain) {
		// straight-walk through triangulation to find intersected edge?
		// TODO args for drawing interal and external and crossing
		// TODO circles diagram from PVectors where Z is radius
		final IncrementalTin tin = PGS_Triangulation.delaunayTriangulationTin(shape, null, constrain, 0, false);
		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		options.setBounds(new Rectangle2D.Double(-1000, -1000, 5000, 4000)); // should be appropriate

		final BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(tin.getVertices(), options);
		Intersection I = new Intersection();
		for (int i = 0; i < shape.getVertexCount() - 1; i++) {
			PVector a = shape.getVertex(i);
			PVector b = shape.getVertex(i + 1);
			I.insert(a.x, a.y, b.x, b.y);
		}

		final IIncrementalTinNavigator navigator = tin.getNavigator();

		PShape lines = PGS.prepareLinesPShape(RGB.PINK, PConstants.SQUARE, 2);
		if (constrain) {
			v.getPolygons().forEach(poly -> {
				poly.getEdges().forEach(e -> {

					// is point A within the shape?
					final boolean inA = tin
							.getRegionConstraint(navigator.getNeighborEdge(e.getA().x, e.getA().y)) != null;
//					final boolean inB = tin
//							.getRegionConstraint(navigator.getNeighborEdge(e.getB().x, e.getB().y)) != null;
					final List<PVector> intersections = I.testMultiple(e.getA().x, e.getA().y, e.getB().x, e.getB().y);

					/**
					 * This can skip some lines on very concave dense shapes (imagine a maze-like
					 * shapes) because the line may cross the shape twice and due to robustness
					 * errors between JTS and the intersection result, I can't reliably check which
					 * segments overlap the shape. If the line crosses a convex part, then ideally
					 * one segment needs output -- the segment within the part; if the line crosses
					 * a concave part then two segments need output, where each exist inside the
					 * shape, skipping the concave gap. For now however, lines which intersect twice
					 * are ignored.
					 */
//					if (intersections.size() == 2) {
//						if (cache.contains(PGS.createLineString(intersections.get(0), intersections.get(1)))) {
//							lines.vertex(intersections.get(0).x, intersections.get(0).y);
//							lines.vertex(intersections.get(1).x, intersections.get(1).y);
//						} else {...}}

					if (intersections.size() == 1) {
						lines.vertex(intersections.get(0).x, intersections.get(0).y); // vertex of intersection
						if (inA) {
							lines.vertex((float) e.getA().x, (float) e.getA().y);
						} else {
							lines.vertex((float) e.getB().x, (float) e.getB().y);
						}
					} else {
						if (intersections.size() == 0) {
							if (inA) {
								lines.vertex((float) e.getA().x, (float) e.getA().y);
								lines.vertex((float) e.getB().x, (float) e.getB().y);
							}
						}
					}
				});
			});
		} else { // display all voronoi polygons/lines
			v.getPolygons().forEach(poly -> {
				poly.getEdges().forEach(e -> {
					final boolean inA = tin
							.getRegionConstraint(navigator.getNeighborEdge(e.getA().x, e.getA().y)) != null;

					if (!inA) {
						final boolean inB = tin
								.getRegionConstraint(navigator.getNeighborEdge(e.getB().x, e.getB().y)) != null;
						if (!inB) {
							lines.vertex((float) e.getA().x, (float) e.getA().y);
							lines.vertex((float) e.getB().x, (float) e.getB().y);
						}
					}

				});
			});
		}

		lines.endShape();
		return lines;
	}

	/**
	 * 
	 * @param circles
	 * @param samples 50 is suitable
	 * @param p
	 * @return
	 */
	public static PShape voronoicirclesites(Iterable<PVector> circles, int samples, PApplet p) {
		final IncrementalTin tin = new IncrementalTin(5);

		final double angleInc = Math.PI * 2 / samples;

		p.fill(RGB.WHITE);
		p.noStroke();

		/**
		 * R-Tree insert circle bound, then use euclid distance - radius to find nearest
		 * circle
		 */
		RTree<PVector> tree = RTree.createRStar(2);

		HashMap<Double, PVector> map = new HashMap<>();
		circles.forEach(c -> {
//			p.ellipse(c.x, c.y, c.z * 2, c.z * 2); // TODO remove
			map.put(PGS.cantorPairing((int) c.x, (int) c.y), c);
			tree.insert(new double[] { c.x - c.z, c.y - c.z }, new double[] { c.x + c.z, c.y + c.z }, c);
			double angle = 0;
			while (angle < Math.PI * 2) {
				if (c.z > 0) {
					tin.add(new Vertex(c.z * Math.cos(angle) + c.x, c.z * Math.sin(angle) + c.y, 0));
				}
				angle += angleInc;
			}
		});

		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		options.setBounds(new Rectangle2D.Double(-100, -100, 3000, 3000)); // should be appropriate
		final BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(tin.getVertices(), options);

		final PShape lines = PGS.prepareLinesPShape(RGB.PINK, PConstants.SQUARE, 3);

		p.stroke(RGB.composeColor(50, 255, 100));
		p.strokeWeight(2);
		final HashSet<Integer> seen = new HashSet<>();
		for (ThiessenPolygon poly : v.getPolygons()) {
			poly: for (IQuadEdge e : poly.getEdges()) {

//				RTreeQueryKnn<PVector> closestA = tree.queryKNN(new double[] { e.getA().x, e.getB().x }, 2);

				final PVector a = new PVector((float) e.getA().x, (float) e.getA().y);
				final PVector b = new PVector((float) e.getB().x, (float) e.getB().y);
				final int hash = a.hashCode() + b.hashCode();

				// reduces edges to check by ~2/3rds
				if (seen.contains(hash)) {
					continue; // edge already seen
				} else {
					seen.add(hash);
				}

				for (PVector c : circles) {
					if (dist2D(c, a, c.z)) { // A is within circle
						if (!dist2D(c, b, c.z)) { // B is within circle too -- overlap
							lines.vertex(b.x, b.y);
							PVector intersect = PVector.sub(b, a).normalize().mult(c.z).add(c);
							lines.vertex(intersect.x, intersect.y);
						}
						continue poly;
					}
					if (dist2D(c, b, c.z)) { // B is within circle
						if (!dist2D(c, a, c.z)) {
							lines.vertex(a.x, a.y);
							PVector intersect = PVector.sub(a, b).normalize().mult(c.z).add(c);
							lines.vertex(intersect.x, intersect.y);
							continue poly;
						}

					}
				}
				lines.vertex((float) e.getA().x, (float) e.getA().y);
				lines.vertex((float) e.getB().x, (float) e.getB().y);
			}
		}

//		System.out.println("n " + n);
//		System.out.println(ld.getResult().getCoordinates().length);
//		System.out.println(lm.getMergedLineStrings().size());
//		v.getPolygons().forEach(poly -> {
//			poly.getEdges().forEach(e -> {
//
//
//			});
//		});

		lines.endShape();
		return lines;
	}

	/**
	 * TODO set clip envelope?
	 * 
	 * @param shape     the shape whose vertices to use as vornoi sites
	 * @param tolerance snapping used in underlying triangulation algorithm
	 * @return
	 */
	public static PShape voronoiDiagram(PShape shape, float tolerance) {
		Geometry g = fromPShape(shape);
		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setTolerance(tolerance);
		v.setSites(g);
		v.setClipEnvelope(new Envelope(0, 1000, 0, 1000)); // speeds up when lots of edges
//		v.setSites(new ArrayList<Coordinate>(Arrays.asList(g.getCoordinates())));
		Geometry out = v.getDiagram(GEOM_FACTORY);
		return toPShape(out.intersection(g)); //
	}

	/**
	 * Voronoi diagram of circle sites (rather than point sites) approximation.
	 * 
	 * https://sci-hub.do/https://www.sciencedirect.com/science/article/abs/pii/S1049965283710394
	 * 
	 * @return
	 */
	public static PShape voronoiCirclesDiagram(PShape shape, float tolerance) {

		/**
		 * Use kdtree instead? Insert kd sites at circle centers and check radius to
		 * determine if line inside
		 * 
		 * Represent circles/voronoi polygons as MonotoneChain?
		 */
		final Geometry g = fromPShape(shape);
//		final PreparedGeometry cache = PreparedGeometryFactory.prepare(g);
		final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);

		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setTolerance(tolerance);
		v.setSites(g);
//		v.setClipEnvelope(new Envelope(0, 1500, 0, 1000)); // TODO

		final Geometry out = v.getDiagram(GEOM_FACTORY);

		final LineDissolver ld = new LineDissolver();
		ld.add(out);
		final Geometry d = ld.getResult();
//		System.out.println("diss " + d.getNumGeometries());

		PShape lines = prepareLinesPShape(composeColor(100, 150, 200), null, null);

		for (int i = 0; i < d.getNumGeometries(); i++) {
			final LineString l = (LineString) d.getGeometryN(i);
//			if (child.getGeometryType() == "LineString") {
//				final LineString l = (LineString) child;

			if (pointLocator.locate(l.getStartPoint().getCoordinate()) == Location.EXTERIOR
					&& pointLocator.locate(l.getEndPoint().getCoordinate()) == Location.EXTERIOR) {
				lines.vertex((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY());
				lines.vertex((float) l.getEndPoint().getX(), (float) l.getEndPoint().getY());
			}

//			if (!cache.contains(l.getStartPoint()) && !cache.contains(l.getEndPoint())) {
//				lines.vertex((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY());
//				lines.vertex((float) l.getEndPoint().getX(), (float) l.getEndPoint().getY());
////				}
//			}
		}

		lines.endShape();
		return lines;
	}

	/**
	 * Doesn't use JTS so that voronoi can be applied to (sub-divided) circles
	 */
	public static PShape voronoiDiagram2(ArrayList<PVector> points) {
		ArrayList<de.alsclo.voronoi.graph.Point> graphIn = new ArrayList<>();
		points.forEach(point -> graphIn.add(new de.alsclo.voronoi.graph.Point(point.x, point.y)));
		Voronoi voronoi = new Voronoi(graphIn);

		PShape lines = prepareLinesPShape(null, null, null);
		voronoi.getGraph().edgeStream().forEach(edge -> {
			if (edge.getA() != null && edge.getB() != null) {
				lines.vertex((float) edge.getA().getLocation().x, (float) edge.getA().getLocation().y);
				lines.vertex((float) edge.getB().getLocation().x, (float) edge.getB().getLocation().y);
			}

		});
		lines.endShape();
		return lines;

	}

	// compare euclidean dist squared
	private static boolean dist2D(PVector a, PVector b, float d) {
		double deltaX = a.y - b.y;
		double deltaY = a.x - b.x;
		return deltaX * deltaX + deltaY * deltaY < d * d;

	}

	private static class IntersectionJTS {

		final ArrayList<Edge> polygonEdges;
		final ArrayList<Edge> edges;

		final SimpleMCSweepLineIntersector i;
		private final SegmentIntersector si;

		public IntersectionJTS(Polygon p) {
			i = new SimpleMCSweepLineIntersector();
			si = new SegmentIntersector(new RobustLineIntersector(), true, false);

			polygonEdges = new ArrayList<>();
			polygonEdges.add(new Edge(p.getCoordinates()));

			edges = new ArrayList<>();
		}

		public void putEdge(Edge e) {
			edges.add(e);
		}

		/**
		 * Computes intersections between the polygon and given (disjoint) edges
		 * 
		 * @return A map of edges->the coordinate on the edge that intersects. Only
		 *         intersecting edges are contained in the map.
		 */
		HashMap<Edge, Coordinate> compute() {
			HashMap<Edge, Coordinate> intersections = new HashMap<>();
			i.computeIntersections(polygonEdges, edges, si);
			edges.forEach(e -> {
				final Iterator<EdgeIntersection> iter2 = e.getEdgeIntersectionList().iterator();
				if (iter2.hasNext()) { // ONLY RETURN A SINGLE INTERSECTION
					intersections.put(e, iter2.next().coord);
				}
			});
			return intersections;
		}

	}

	// fast line intersection, backed by RTree to find candidate lines to test
	// intersection
	// build up based on one set, then test other set
	// see https://github.com/cardawid/balabanalgorithm
	// designed for voronoi polygon / polygon intersection
	// built to provide a much faster alternative to PTS' geometry.intersect() and
	// easily
	// works with any triangulation/voronoi implementation
	public static class Intersection {

		// TODO decide on PVectors / double[] etc

		private final RTree<E> rtree;

		public Intersection() {
			rtree = RTree.createRStar(2);
			// TODO look into JTS HPRtree
		}

		// insert a segment from the set to test against
		public void insert(double x1, double y1, double x2, double y2) {
			E e = new E(x1, y1, x2, y2);
			rtree.insert(new double[] { Math.min(x1, x2), Math.min(y1, y2) },
					new double[] { Math.max(x1, x2), Math.max(y1, y2) }, e);
		}

		// doesn't check whether for co-linear points -- assume not
		// breaks when the first intersection point is found
		public PVector test(double x1, double y1, double x2, double y2) {
			RTreeIterator<E> iterator = rtree.queryIntersect(new double[] { Math.min(x1, x2), Math.min(y1, y2) },
					new double[] { Math.max(x1, x2), Math.max(y1, y2) });
			while (iterator.hasNext()) {
				Entry<E> r = iterator.next();
				final PVector p1 = r.value().p1;
				final PVector p2 = r.value().p2;
				final PVector p3 = new PVector((float) x1, (float) y1);
				final PVector p4 = new PVector((float) x2, (float) y2);

				if (segmentsIntersect(p1, p2, p3, p4)) {
					double m1 = (p2.y - p1.y) / (p2.x - p1.x);
					double m2 = (p4.y - p3.y) / (p4.x - p3.x);
					double b1 = p1.y - m1 * p1.x;
					double b2 = p3.y - m2 * p3.x;
					double x = (b2 - b1) / (m1 - m2);
					double y = (m1 * b2 - m2 * b1) / (m1 - m2);
					return new PVector((float) x, (float) y);
				}
			}
			return null;
		}

		// test for up to 2 intersection points
		public List<PVector> testMultiple(double x1, double y1, double x2, double y2) {
			RTreeIterator<E> iterator = rtree.queryIntersect(new double[] { Math.min(x1, x2), Math.min(y1, y2) },
					new double[] { Math.max(x1, x2), Math.max(y1, y2) }); // usually no more than a few elements
			List<PVector> intersections = new ArrayList<>();
			while (iterator.hasNext()) {
				Entry<E> r = iterator.next();
				final PVector p1 = r.value().p1;
				final PVector p2 = r.value().p2;
				final PVector p3 = new PVector((float) x1, (float) y1);
				final PVector p4 = new PVector((float) x2, (float) y2);

				if (segmentsIntersect(p1, p2, p3, p4)) {
					double m1 = (p2.y - p1.y) / (p2.x - p1.x);
					double m2 = (p4.y - p3.y) / (p4.x - p3.x);
					double b1 = p1.y - m1 * p1.x;
					double b2 = p3.y - m2 * p3.x;
					double x = (b2 - b1) / (m1 - m2);
					double y = (m1 * b2 - m2 * b1) / (m1 - m2);
					intersections.add(new PVector((float) x, (float) y));
					if (intersections.size() == 2) {
						// can't be more than 2 intersections in a voronoi lines diagram?
						return intersections;
					}
				}
			}
			return intersections;
		}

		private static boolean segmentsIntersect(PVector p1, PVector p2, PVector p3, PVector p4) {

			// Get the orientation of points p3 and p4 in relation
			// to the line segment (p1, p2)
			int o1 = orientation(p1, p2, p3);
			int o2 = orientation(p1, p2, p4);
			int o3 = orientation(p3, p4, p1);
			int o4 = orientation(p3, p4, p2);

			// If the points p1, p2 are on opposite sides of the infinite
			// line formed by (p3, p4) and conversly p3, p4 are on opposite
			// sides of the infinite line formed by (p1, p2) then there is
			// an intersection.
			return (o1 != o2 && o3 != o4);
		}

		// Finds the orientation of point 'c' relative to the line segment (a, b)
		// Returns 0 if all three points are collinear.
		// Returns -1 if 'c' is clockwise to segment (a, b), i.e right of line formed by
		// the segment.
		// Returns +1 if 'c' is counter clockwise to segment (a, b), i.e left of line
		// formed by the segment.
		private static int orientation(PVector a, PVector b, PVector c) {
			double value = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
			// if (Math.abs(value) < 1e-7)
			// return 0;
			return (value > 0) ? -1 : +1;
		}

		private static class E {
			final PVector p1, p2;

			public E(double x1, double y1, double x2, double y2) {
				p1 = new PVector((float) x1, (float) y1);
				p2 = new PVector((float) x2, (float) y2);
			}
		}
	}
}
