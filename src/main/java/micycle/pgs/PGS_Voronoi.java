package micycle.pgs;

import static micycle.pgs.PGS_Conversion.fromPShape;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geomgraph.Edge;
import org.locationtech.jts.geomgraph.EdgeIntersection;
import org.locationtech.jts.geomgraph.index.EdgeSetIntersector;
import org.locationtech.jts.geomgraph.index.SegmentIntersector;
import org.locationtech.jts.geomgraph.index.SimpleMCSweepLineIntersector;
import org.tinfour.common.IIncrementalTinNavigator;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.voronoi.BoundedVoronoiBuildOptions;
import org.tinfour.voronoi.BoundedVoronoiDiagram;
import org.tinfour.voronoi.ThiessenPolygon;
import org.tinspin.index.PointDistanceFunction;
import org.tinspin.index.PointEntryDist;
import org.tinspin.index.kdtree.KDTree;
import org.tinspin.index.rtree.Entry;
import org.tinspin.index.rtree.RTree;
import org.tinspin.index.rtree.RTreeIterator;

import micycle.pgs.PGS.PEdge;
import micycle.pgs.color.RGB;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * 
 * Voronoi Diagrams of shapes and point sets.
 * 
 * @author Michael Carleton
 *
 */
public final class PGS_Voronoi {

	private PGS_Voronoi() {
	}

	/**
	 * Generates a Voronoi diagram from a shape, where shape vertices are Voronoi
	 * point sites. This method outputs the Voronoi diagram as lines.
	 * 
	 * @param shape     the shape whose vertices to use as Voronoi sites
	 * @param constrain whether to constrain the diagram lines to the shape. When
	 *                  true, the output includes only voronoi line segments within
	 *                  the shape.
	 * @return a PShape consisting of voronoi lines
	 */
	public static PShape voronoiDiagram(PShape shape, boolean constrain) {
		final IncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape, null, constrain, 0, false);

		final Geometry g = fromPShape(shape);
		
		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		options.setBounds(tin.getBounds());
		final BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(tin.getVertices(), options);

		final IIncrementalTinNavigator navigator = tin.getNavigator();

		final PShape voronoi = new PShape(PConstants.GROUP);
		final PShape axis = PGS.prepareLinesPShape(RGB.PINK, PConstants.ROUND, 4);
		final PShape lines = PGS.prepareLinesPShape(RGB.PINK, PConstants.SQUARE, 2);

		if (constrain) { // constrain: include only inner voronoi line segments
			final Coordinate[] coords = g.getCoordinates();

			final SweepLineSegmentIntersection intersection = new SweepLineSegmentIntersection(coords);
			final HashSet<Integer> seen = new HashSet<>();
			final ArrayList<Edge> edges = new ArrayList<>();

			for (ThiessenPolygon poly : v.getPolygons()) {
				for (IQuadEdge e : poly.getEdges()) {
					final Coordinate c1 = new Coordinate(e.getA().x, e.getA().y);
					final Coordinate c2 = new Coordinate(e.getB().x, e.getB().y);
					final int hash = c1.hashCode() ^ c2.hashCode(); // order-invariant hash
					if (seen.add(hash)) { // only process unique edges
						/*
						 * It's a little faster to filter edges to intersection-check here by first
						 * checking if one if its vertices is inside and the other outside the shape.
						 * The tradeoff is that we miss segments that cross outside the shape yet both
						 * vertices are within the shape (which can happen on very concave shapes).
						 */
						final boolean inA = tin.getRegionConstraint(navigator.getNeighborEdge(c1.x, c1.y)) != null;
						final boolean inB = tin.getRegionConstraint(navigator.getNeighborEdge(c2.x, c2.y)) != null;
						if (inA ^ inB) {
							edges.add(new Edge(new Coordinate[] { c1, c2 }));
						} else if (inA) { // edges are completely inside the shape (main medial axis)
							axis.vertex((float) c1.x, (float) c1.y);
							axis.vertex((float) c2.x, (float) c2.y);
						}
					}
				}
			}

			final HashMap<Edge, Coordinate> intersections = intersection.compute(edges);
			intersections.forEach((e, c) -> {
				/*
				 * When segments intersects >1 time, if it crosses a convex part, then ideally
				 * one segment needs output -- the segment within the part; if the segment
				 * crosses a concave part then two segments need output, where each exist inside
				 * the shape, skipping the concave gap. For now however, lines which intersect
				 * >1 time are ignored.
				 */
				if (e.isIsolated()) { // ignore lines with >1 intersection
					final Coordinate c1 = e.getCoordinates()[0];
					lines.vertex((float) c.x, (float) c.y);

					final boolean inA = tin.getRegionConstraint(navigator.getNeighborEdge(c1.x, c1.y)) != null;
					if (inA) {
						lines.vertex((float) c1.x, (float) c1.y);
					} else {
						final Coordinate c2 = e.getCoordinates()[1];
						lines.vertex((float) c2.x, (float) c2.y);
					}
				}
			});
			axis.endShape();
		} else { // no constraining: display all voronoi polygons/lines
			v.getPolygons().forEach(poly -> poly.getEdges().forEach(e -> {
				final boolean inA = tin.getRegionConstraint(navigator.getNeighborEdge(e.getA().x, e.getA().y)) != null;
				if (!inA) {
					final boolean inB = tin.getRegionConstraint(navigator.getNeighborEdge(e.getB().x, e.getB().y)) != null;
					if (!inB) {
						lines.vertex((float) e.getA().x, (float) e.getA().y);
						lines.vertex((float) e.getB().x, (float) e.getB().y);
					}
				}

			}));
		}
		lines.endShape();

		voronoi.addChild(lines);
		voronoi.addChild(axis);
		return voronoi;
	}

	/**
	 * Generates a Voronoi diagram from a set of points. This method outputs the
	 * Voronoi diagram as lines.
	 * 
	 * @param points    the set of points to use as Voronoi sites
	 * @param constrain whether to constrain the diagram's lines to the concave hull
	 *                  of the point set
	 * @return
	 */
	public static PShape voronoiDiagram(Collection<PVector> points, boolean constrain) {
		return voronoiDiagram(PGS_Conversion.toPointsPShape(points), constrain);
	}

	/**
	 * Generates a Voronoi diagram from a shape, where shape vertices are Voronoi
	 * point sites. This method outputs the Voronoi diagram as polygonal cells.
	 * 
	 * @param shape the shape whose vertices to use as Voronoi sites
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #voronoiCells(List)
	 */
	public static PShape voronoiCells(PShape shape) {
		final IncrementalTin tin = PGS_Triangulation.delaunayTriangulationMesh(shape, null, false, 0, false);

		final Envelope envelope = fromPShape(shape).getEnvelopeInternal();
		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		options.setBounds(new Rectangle2D.Double(envelope.getMinX(), envelope.getMinY(), envelope.getMaxX() - envelope.getMinX(),
				envelope.getMaxY() - envelope.getMinY()));

		final BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(tin.getVertices(), options);

		final PShape voronoi = new PShape(PConstants.GROUP);

		for (ThiessenPolygon poly : v.getPolygons()) {
			final PShape cell = new PShape(PShape.GEOMETRY);
			cell.setFill(true);
			cell.setFill(RGB.WHITE);
			cell.setStroke(true);
			cell.setStroke(RGB.PINK);
			cell.setStrokeWeight(3);
			cell.beginShape();
			for (IQuadEdge e : poly.getEdges()) {
				cell.vertex((float) e.getA().x, (float) e.getA().y);
			}
			cell.endShape(PShape.CLOSE);
			voronoi.addChild(cell);
		}

		return voronoi;
	}

	/**
	 * Generates a Voronoi diagram from a set of points. This method outputs the
	 * Voronoi diagram as polygonal cells.
	 * 
	 * @param points the set of points to use as Voronoi sites
	 * @return a GROUP PShape, where each child shape is a Voronoi cell
	 * @see #voronoiCells(PShape)
	 */
	public static PShape voronoiCells(Collection<PVector> points) {
		return voronoiCells(PGS_Conversion.toPointsPShape(points));
	}

	/**
	 * Generates a Voronoi diagram from circle sites (rather than point sites).
	 * <p>
	 * Circle sites are modelled by PVectors, where x, y correspond to the center of
	 * the site and z corresponds to the radius of the site.
	 * 
	 * @param circles       list of PVectors to use as circle sites
	 * @param circleSamples defines how many samples from each circle's
	 *                      circumference should be used to compute the voronoi
	 *                      diagram. 25-50 is a suitable range
	 * @param drawBranches  whether to the draw/output branches from the coming from
	 *                      the center of each circle
	 * @return
	 */
	public static PShape voronoiCirclesDiagram(Collection<PVector> circles, int circleSamples, boolean drawBranches) {
		final IncrementalTin tin = new IncrementalTin(5);

		final PointDistanceFunction pdf2D = (p1, p2) -> {
			final double deltaX = p1[0] - p2[0];
			final double deltaY = p1[1] - p2[1];
			return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
		};

		final KDTree<PVector> sites = KDTree.create(2, pdf2D); // TODO use vpTree?
		final List<PVector> sitesList = new ArrayList<>();
		final double angleInc = Math.PI * 2 / circleSamples;
		circles.forEach(c -> {
			if (c.z >= 0) {
				sitesList.add(c);
				double angle = 0;
				while (angle < Math.PI * 2) {
					tin.add(new Vertex(c.z * Math.cos(angle) + c.x, c.z * Math.sin(angle) + c.y, 0));
					angle += angleInc;
				}
			}
		});
		Collections.shuffle(sitesList); // shuffle vertices for more balanced KDTree
		sitesList.forEach(c -> sites.insert(new double[] { c.x, c.y }, c));

		final BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
		options.setBounds(new Rectangle2D.Double(-1000, -1000, 3000, 3000)); // should be enough
		final BoundedVoronoiDiagram v = new BoundedVoronoiDiagram(tin.getVertices(), options);

		final PShape lines = PGS.prepareLinesPShape(RGB.PINK, PConstants.SQUARE, 3);
		final HashSet<PEdge> seen = new HashSet<>();
		for (ThiessenPolygon poly : v.getPolygons()) {
			for (IQuadEdge e : poly.getEdges()) {
				final PVector a = new PVector((float) e.getA().x, (float) e.getA().y);
				final PVector b = new PVector((float) e.getB().x, (float) e.getB().y);
				final PEdge edge = new PEdge(a, b);

				if (!seen.add(edge)) { // reduces edges to check by ~2/3rds
					continue;
				}

				PointEntryDist<PVector> nearestSite = sites.query1NN(new double[] { e.getA().x, e.getA().y });

				if (nearestSite.dist() < nearestSite.value().z) {
					if (drawBranches && distGreater(a, b, nearestSite.value().z)) {
						PVector intersect = PVector.sub(b, a).normalize().mult(nearestSite.value().z).add(nearestSite.value());
						lines.vertex(b.x, b.y);
						lines.vertex(intersect.x, intersect.y);
					}
					continue;
				}

				nearestSite = sites.query1NN(new double[] { e.getB().x, e.getB().y });
				if (nearestSite.dist() < nearestSite.value().z) {
					if (drawBranches && distGreater(a, b, nearestSite.value().z)) {
						PVector intersect = PVector.sub(a, b).normalize().mult(nearestSite.value().z).add(nearestSite.value());
						lines.vertex(a.x, a.y);
						lines.vertex(intersect.x, intersect.y);
					}
					continue;
				}

				lines.vertex((float) e.getA().x, (float) e.getA().y);
				lines.vertex((float) e.getB().x, (float) e.getB().y);
			}
		}

		lines.endShape();
		return lines;
	}

	private static boolean distGreater(PVector a, PVector b, float d) {
		final double deltaX = a.x - b.x;
		final double deltaY = a.y - b.y;
		return deltaX * deltaX + deltaY * deltaY > (d * d);

	}

	/**
	 * Wrapper/helper class to run JTS SweepLineIntersector compute intersections
	 * between two sets of line segments.
	 * 
	 * @author Michael Carleton
	 *
	 */
	private static class SweepLineSegmentIntersection {

		final List<Edge> polygonEdges;

		private final EdgeSetIntersector i;
		private final SegmentIntersector si;

		/**
		 * 
		 * @param polygonCoords coords from polygon to test against
		 */
		SweepLineSegmentIntersection(Coordinate[] polygonCoords) {
			i = new SimpleMCSweepLineIntersector();
			si = new SegmentIntersector(new RobustLineIntersector(), true, false);

			polygonEdges = new ArrayList<>();
			polygonEdges.add(new Edge(polygonCoords)); // list of a single edge where the edge has many segments
		}

		/**
		 * Computes intersections between the polygon and the given (disjoint) edges
		 * 
		 * @return A map of given edges->the coordinate on the edge that intersects.
		 *         Only intersecting edges are contained in the map.
		 */
		@SuppressWarnings("unchecked")
		HashMap<Edge, Coordinate> compute(List<Edge> edges) {
			HashMap<Edge, Coordinate> intersections = new HashMap<>();
			i.computeIntersections(polygonEdges, edges, si);
			edges.forEach(e -> {
				final Iterator<EdgeIntersection> iter2 = e.getEdgeIntersectionList().iterator();
				if (iter2.hasNext()) {
					intersections.put(e, iter2.next().coord); // ONLY RETURN A SINGLE INTERSECTION
					if (iter2.hasNext()) {
						e.setIsolated(false); // mark lines with 2 intersections to ignore later
					}
				}
			});
			return intersections;
		}

	}

	/**
	 * Computes points of intersection between two sets of line segments. The
	 * segment set to test against is backed by an RTree to find quickly find
	 * possible intersecting candidates. The idea is you should insert all the
	 * segments from one set first, then test segments from the other set as you
	 * please (one-by-one).
	 * 
	 * <p>
	 * Built to provide a much faster alternative to PTS' geometry.intersect() for
	 * voronoi diagram cropping.
	 * 
	 * @author Michael Carleton
	 * @deprecated it's a little slower than IntersectionJTS
	 */
	@SuppressWarnings("unused")
	@Deprecated
	private static class SegmentIntersection {

		// TODO decide on PVectors / double[] etc
		// TODO look into JTS HPRtree
		private final RTree<E> rtree;

		SegmentIntersection() {
			rtree = RTree.createRStar(2);
		}

		/**
		 * Inserts a segment from the set to test against
		 */
		public void insert(double x1, double y1, double x2, double y2) {
			E e = new E(x1, y1, x2, y2);
			rtree.insert(new double[] { Math.min(x1, x2), Math.min(y1, y2) }, new double[] { Math.max(x1, x2), Math.max(y1, y2) }, e);
		}

		/**
		 * Computes a single intersection between the provided line segment and the test
		 * set (populated via insert()). If a point of intersection if found, this
		 * method returns early (even if there are multiple intersections in total).
		 * Doesn't check for co-linear points.
		 * 
		 * @return returns the first point of intersection; null if no intersection
		 */
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

		/**
		 * * Computes upto two intersections between the provided line segment and the
		 * test set (populated via insert()).
		 * 
		 * @return list of points that intersecting on the line segment. If there are no
		 *         intersections, the list is empty.
		 */
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
