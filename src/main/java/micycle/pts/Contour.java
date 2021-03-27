package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import static micycle.pts.PTS.GEOM_FACTORY;
import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.linearref.LinearIterator;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.linemerge.LineMergeEdge;
import org.locationtech.jts.operation.linemerge.LineMergeGraph;
import org.locationtech.jts.operation.linemerge.LineMerger;
import org.locationtech.jts.planargraph.Node;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.contour.ContourBuilderForTin;
import org.tinfour.standard.IncrementalTin;
import org.twak.camp.Corner;
import org.twak.camp.Machine;
import org.twak.camp.Skeleton;
import org.twak.utils.collections.Loop;
import org.twak.utils.collections.LoopL;

import hageldave.jplotter.misc.Contours;
import hageldave.jplotter.renderables.Lines.SegmentDetails;
import micycle.medialAxis.MedialAxis;
import micycle.medialAxis.MedialAxis.Branch;
import micycle.pts.PTS.LinearRingIterator;
import micycle.pts.color.RGB;
import micycle.pts.utility.PoissonDistribution;
import micycle.pts.utility.SolubSkeleton;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Houses a variety of methods for producing different kinds of shape contours.
 * 
 * <p>
 * A 2D contour is a closed sequence (a cycle) of 3 or more connected 2D
 * oriented straight line segments called contour edges. The endpoints of the
 * contour edges are called vertices. Each contour edge shares its endpoints
 * with at least two other contour edges.
 * 
 * @author Michael Carleton
 *
 */
public class Contour {

	/**
	 * TODO implement 'Base Point Split Algorithm for Generating Polygon Skeleton
	 * Lines'
	 */

	/**
	 * Set of points in space equidistant to 2 or more points on the surface. As
	 * density of boundary points goes to infinity, a voronoi diagram converges to a
	 * medial axis.
	 * 
	 * @param shape
	 * @param density          distance tolerance for boundary densification
	 *                         (smaller values more converge towards a more accurate
	 *                         axis, but slower), 5-20 is appropriate
	 * @param maximumCloseness the SQUARE of the
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public static PShape medialAxis(PShape shape, float density, float maximumCloseness, PApplet p) {
		// TODO remove the duplicate Voronoi edges ?
		final Geometry g = fromPShape(shape);
		final Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(density);
		d.setValidate(false); // don't perform validation processing (a little faster)
		final Geometry dense = d.getResultGeometry();

		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setSites(dense);
		Geometry voronoi = v.getDiagram(GEOM_FACTORY);

		final Geometry small = dense.buffer(-maximumCloseness);
		PreparedGeometry cache = PreparedGeometryFactory.prepare(small); // provides MUCH faster contains() check

		// inline lines creation
		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(ROUND);
		lines.setStroke(true);
		lines.setStrokeWeight(2);
		lines.setStroke(-1232222);
		lines.beginShape(LINES);

		ArrayList<LineString> axis = new ArrayList<LineString>();

		// TODO getCoordinates() call slow on cell too?

		// TODO compare both points at once?
		int z = 0;
		for (int i = 0; i < voronoi.getNumGeometries(); i++) {
			Polygon cell = (Polygon) voronoi.getGeometryN(i); // TODO .get(0) prevent occasional crash
			for (int j = 0; j < cell.getCoordinates().length - 1; j++) {
				Coordinate a = cell.getCoordinates()[j];
//				CoordinateSequence seq = geometryFactory.getCoordinateSequenceFactory().create(new Coordinate[] { a });
				if (cache.covers(GEOM_FACTORY.createPoint(a))) {
					Coordinate b = cell.getCoordinates()[j + 1];
//					seq = geometryFactory.getCoordinateSequenceFactory().create(new Coordinate[] { b });
					if (cache.covers(GEOM_FACTORY.createPoint(b))) {
//						lines.vertex((float) a.x, (float) a.y);
//						lines.vertex((float) b.x, (float) b.y);
//						p.stroke((float) a.x % 255, (float) b.y % 255, (z * 5) % 255);
//						p.line((float) a.x, (float) a.y, (float) b.x, (float) b.y);
						axis.add(GEOM_FACTORY.createLineString(new Coordinate[] { a, b }));
						z++;
//						lm.add(geometryFactory.createLineString(new Coordinate[] { a, b }));
					}
				}
			}
		}

//		System.out.println("axis lines" + axis.size());
		LineDissolver ld = new LineDissolver();
		ld.add(axis);
		Geometry medialAxisLines = ld.getResult();

		HashSet<Double> seen = new HashSet<Double>();

		LineMergeGraph graph = new LineMergeGraph();

		// now use LineMergeGraph to sew edges

		z = 0;
		for (int i = 0; i < medialAxisLines.getNumGeometries(); i++) {
			LineString l = (LineString) medialAxisLines.getGeometryN(i);
			graph.addEdge(l);
			p.strokeWeight(3);
			p.stroke((float) l.getStartPoint().getX() % 255, (float) l.getEndPoint().getY() % 255, (z * 5) % 255);
//			p.line((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY(), (float) l.getEndPoint().getX(),
//					(float) l.getEndPoint().getY());
			double o = cantorPairing(l.getStartPoint().getX(), l.getStartPoint().getY());
			if (!seen.add(o)) {
				p.strokeWeight(10);
//				p.point((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY());
			}
			z++;
		}
//		System.out.println("axis lines" + z);

		/**
		 * Use kdtree to get nearest 4 points of a site; then use tin naviagator to get
		 * triangles containing points?
		 */

//		PrimMinimumSpanningTree<Vertex, QuadEdge> pst = )

		Iterator<Node> ni = (Iterator<Node>) graph.nodeIterator();

		((HashSet<LineMergeEdge>) graph.getEdges()).forEach(e -> {
			if ((e.getDirEdge(0).getFromNode().getDegree() > 1) && (e.getDirEdge(0).getToNode().getDegree() > 1)) {
				final LineString l = e.getLine();
				p.strokeWeight(3);
				p.stroke((float) l.getStartPoint().getX() % 255, (float) l.getEndPoint().getY() % 255, 125);
				p.line((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY(),
						(float) l.getEndPoint().getX(), (float) l.getEndPoint().getY());
			}
		});

		p.stroke(123, 76, 81);
		p.strokeWeight(10);
		ni.forEachRemaining(n -> {
//			n.get
			if (n.getDegree() > 1) { // ignore linestring endpoints
				p.point((float) n.getCoordinate().x, (float) n.getCoordinate().y);
			}
		});

//		for (LineString l : ((List<LineString>) lm.getMergedLineStrings())) {
//			z++;
//			lines.vertex((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY());
//			lines.vertex((float) l.getEndPoint().getX(), (float) l.getEndPoint().getY());
//		}
//		System.out.println(z);

		lines.endShape();
//		return toPShape());
		return lines;
	}

	/**
	 * Returns a simplified medial axis graph leaving only maximal-length lines in
	 * which every unique segment appears once only. The output lines run between
	 * node vertices of the input, which are vertices which have either degree 1, or
	 * degree 3 or greater.
	 * 
	 * @return
	 */
	public static PShape dissolvedMedialAxis() {
		// use JTS GeometryGraph?
		// TODO break out dissolver+linemergegraph into here
		return null;
	}

	public static void medialAxis4(PShape shape, PApplet p) {
		Geometry g = fromPShape(shape);
		MedialAxis m = new MedialAxis(g);
//		m.getDissolvedGeometry();
//		p.shape(toPShape(
//				DouglasPeuckerSimplifier.simplify(m.getDissolvedGeometry(), p.map(p.mouseX, 0, p.width, 0, 100))));
		PShape s = toPShape(m.getDissolvedGeometry());

		s = toPShape(DouglasPeuckerSimplifier.simplify(m.getDissolvedGeometry(), 15));

		s.setStrokeWeight(15);
		for (int i = 0; i < s.getChildCount(); i++) {
			s.getChild(i).setStrokeWeight(6);
			s.getChild(i).setStrokeCap(ROUND);
			s.getChild(i).setStroke(p.color((i * 23.3f) % 255, 255 - (i * 33.3f) % 255, (i * 16.3f) % 255));
		}
		p.shape(s);

		p.beginShape(LINES);
		m.getEdges().forEach(e -> {
//			p.line((float) e.head.position.x, (float) e.head.position.y, (float) e.tail.position.x,
//					(float) e.tail.position.y);
//			p.vertex((float) e.head.position.x, (float) e.head.position.y);
//			p.vertex((float) e.tail.position.x, (float) e.tail.position.y);
		});
		p.endShape(LINES);
//		m.getLineMergeGraph()

		p.stroke(14, 150, 14);
//		m.getDisks().forEach(d -> {
//			p.ellipse((float) d.position.x, (float) d.position.y, (float) d.radius * 2, (float) d.radius * 2);
//		});
	}

	/**
	 * TODO remove...temporary method to test MedialAxis lib.
	 * 
	 * @param shape
	 * @param p
	 */
	public static void medialAxis3(PShape shape, PApplet p) {
		Geometry g = fromPShape(shape);
		MedialAxis m = new MedialAxis(g);
//		m.drawVDM(p);
//		m.drawVDMPrune(p, p.map(p.mouseX, 0, p.width, 0, 1));
//		m.drawVDMPrune(p, 500);
//		p.noFill();
		p.noStroke();
		p.fill(0);
//		m.getLeaves().forEach(l -> {
//			p.ellipse((float) l.position.x, (float) l.position.y, 5, 5);
//		});

//		m.getAncestors(m.getLeaves().get(0)).forEach(v -> {
//
//		});
//		List<micycle.medialAxis.MedialAxis.VD> ancestors = m
//				.getAncestors(m.getLeaves().get((int) (p.frameCount * 0.1) % m.getLeaves().size()));
//
//		micycle.medialAxis.MedialAxis.VD base = ancestors.get(0);
//		for (int i = 1; i < ancestors.size(); i++) {
//			p.line((float) base.position.x, (float) base.position.y, (float) ancestors.get(i).position.x,
//					(float) ancestors.get(i).position.y);
//			base = ancestors.get(i);
//		}

//		for

		p.fill(0);
		p.stroke(0, 255, 0);
		p.strokeWeight(8);
//		m.voronoiDisks.forEach(v -> {
//			final Coordinate pa = v.position;
//			p.point((float) pa.getX(), (float) pa.getY());
//		});

		p.beginShape(PConstants.TRIANGLES);
		p.noStroke();
		m.getDisks().forEach(d -> {
			final SimpleTriangle t = d.t;
//			double r 
			if (d.radius > 0) {
//				p.beginShape();
				p.noFill();
				if (d.distance < p.map(p.mouseX, 0, p.width, 0, (float) m.furthestNode.distance)) {
					p.fill(255, 0, 0);
				} else {
					p.fill(0, 255, 0);
				}
//				p.fill(150, p.map((float) d.distance, 0, (float) m.furthestNode.distance, 0, 255), 200);
//				p.strokeWeight(2);
//				p.stroke(255, 0, 0);
				Vertex v = t.getVertexA();
//				p.noStroke();
				p.vertex((float) v.getX(), (float) v.getY());
//				p.point((float) v.getX(), (float) v.getY());
				v = t.getVertexB();
				p.vertex((float) v.getX(), (float) v.getY());
//				p.point((float) v.getX(), (float) v.getY());
				v = t.getVertexC();
				p.vertex((float) v.getX(), (float) v.getY());
//				p.point((float) v.getX(), (float) v.getY());
//				p.endShape(p.CLOSE);

//				p.strokeWeight(2);
//				p.stroke(0, 150, 0);
//				System.out.println(t.area());
				p.point((float) d.position.x, (float) d.position.y);
//				p.ellipse((float) d.position.x, (float) d.position.y, (float) d.radius * 2, (float) d.radius * 2);
				p.strokeWeight(7);
//				TriangulationPoint va = t.points[0];
//				p.point((float) va.getX(), (float) va.getY());
//				va = t.points[1];
//				p.point((float) va.getX(), (float) va.getY());
//				va = t.points[2];
//				p.point((float) va.getX(), (float) va.getY());
			}
		});
		p.endShape();

//		m.validTris.forEach(t -> {
//
//			double[] c = MedialAxis.circumcircle(t);
//			if (c[2] > 100) {
//				p.beginShape();
//				p.noFill();
//				p.strokeWeight(2);
//				p.stroke(255, 0, 0);
//				TriangulationPoint v = t.points[0];
////				p.noStroke();
//				p.vertex((float) v.getX(), (float) v.getY());
////				p.point((float) v.getX(), (float) v.getY());
//				v = t.points[1];
//				p.vertex((float) v.getX(), (float) v.getY());
////				p.point((float) v.getX(), (float) v.getY());
//				v = t.points[2];
//				p.vertex((float) v.getX(), (float) v.getY());
////				p.point((float) v.getX(), (float) v.getY());
//				p.endShape(p.CLOSE);
//
//				p.strokeWeight(2);
//				p.stroke(0, 150, 0);
//				System.out.println(t.area());
//				p.point((float) c[0], (float) c[1]);
//				p.ellipse((float) c[0], (float) c[1], (float) c[2] * 2, (float) c[2] * 2);
//				p.strokeWeight(7);
//				TriangulationPoint va = t.points[0];
//				p.point((float) va.getX(), (float) va.getY());
//				va = t.points[1];
//				p.point((float) va.getX(), (float) va.getY());
//				va = t.points[2];
//				p.point((float) va.getX(), (float) va.getY());
//			}
//		});

//		p.stroke(0);
//		m.branches.forEach(b -> {
//			System.out.println(b.size());
//			b.forEach(vd -> {
//
//				p.point((float) vd.circumcircle[0], (float) vd.circumcircle[1]);
//			});
//		});

//		System.out.println(m.branches.size());

		int q = 0;
		p.strokeWeight(3);
//		for (List<MedialAxis.VD> disks : m.branches) {
////			System.out.println("len" + disks.size());
//			p.stroke((int) disks.get(0).circumcircle[0] % 255, q % 255, disks.size() * 71 % 255);
//
//			micycle.medialAxis.MedialAxis.VD last = disks.get(0); // use parent
//			for (MedialAxis.VD vd : disks) {
//				p.line((float) last.circumcircle[0], (float) last.circumcircle[1], (float) vd.circumcircle[0],
//						(float) vd.circumcircle[1]);
//				last = vd;
//			}
////			p.line((float) last.circumcircle[0], (float) last.circumcircle[1],
////					(float) last.children.get(0).circumcircle[0], (float) last.children.get(0).circumcircle[1]);
//			q += 37;
//		}

//		for (ArrayList<micycle.medialAxis.MedialAxis.VD> branch : m.getBranches()) {
//			micycle.medialAxis.MedialAxis.VD last = branch.get(0); // use parent
//			p.stroke((int) branch.get(0).circumcircle[0] % 255, q % 255, branch.size() * 71 % 255);
//			for (MedialAxis.VD vd : branch) {
//				p.line((float) last.circumcircle[0], (float) last.circumcircle[1], (float) vd.circumcircle[0],
//						(float) vd.circumcircle[1]);
//				last = vd;
//			}
//			q += 37;
//		}

		micycle.medialAxis.MedialAxis.MedialDisk d;
		d = m.nearestDisk(p.mouseX, p.mouseY);
//		System.err.println(d.featureArea);

		d = m.rootNode;
		p.noFill();
		p.stroke(50, 150, 250);
		float dist = p.map(p.mouseX, 0, p.width, 0, (float) m.furthestNode.distance) * 2;
//		dist = (float) d.radius * 2;
		p.ellipse((float) d.position.x, (float) d.position.y, dist, dist);
//		p.point((float) d.position.x, (float) d.position.y);
		p.fill(255);
//		p.text(d.depthBF, 10, 80);
//		p.textAlign(PConstants.LEFT);
//		p.text((float) d.featureArea, 0, 100); // feature area

		d = m.rootNode;
		p.noFill();
		p.stroke(50, 150, 250);
//		p.ellipse((float) d.position.x, (float) d.position.y, (float) d.radius * 2, (float) d.radius * 2);
//		p.point((float) d.position.x, (float) d.position.y);

		p.stroke(255, 0, 0);
//		TriangulationPoint v = d.t.points[0];
//		p.point((float) v.getX(), (float) v.getY());
//		v = d.t.points[1];
//		p.point((float) v.getX(), (float) v.getY());
//		v = d.t.points[2];
//		p.point((float) v.getX(), (float) v.getY());
//		p.point((float) d.position.x, (float) d.position.y);

		p.strokeWeight(4);
		p.colorMode(PConstants.HSB, 1, 1, 1, 1);
		float x = 0;
		for (Branch s : m.getBranches()) {
			p.stroke((x += 1.618034) % 1, 1, 1);
			s.edges.forEach(e -> {
				p.line((float) e.head.position.x, (float) e.head.position.y, (float) e.tail.position.x,
						(float) e.tail.position.y);
			});
		}
		p.colorMode(PConstants.RGB, 255, 255, 255, 255);

//		for (micycle.pts.MedialAxis.VD vd : m.getBifurcations()) {
//			p.noFill();
//			p.stroke(50, 150, 250);
////				p.ellipse((float) vd.circumcircle[0], (float) vd.circumcircle[1], (float) vd.circumcircle[2] * 2,
////						(float) vd.circumcircle[2] * 2);
//			p.stroke(255, 0, 0);
//			p.strokeWeight(15);
//			p.point((float) vd.circumcircle[0], (float) vd.circumcircle[1]);
//		}

//			for (int i = 0; i < 300; i++) {
//				SimpleTriangle t = m.voronoiDisks.get(i).t;
//
//				if (m.voronoiDisks.get(i).area > 1000) {
//			v = vd.t.getVertexA();
//			p.beginShape();
//			p.noStroke();
//			p.fill(100);
//			p.vertex((float) v.x, (float) v.y);
//			v = vd.t.getVertexB();
//			p.vertex((float) v.x, (float) v.y);
//			v = vd.t.getVertexC();
//			p.vertex((float) v.x, (float) v.y);
//			p.endShape();
//				}
//
//			}

	}

//	public static PShape gishurMedialAxis(PShape shape) {
//		gishur.x.XPolygon poly = Conversion.toXPolygon(shape);
//		gishur.x.voronoi.Skeleton s = new gishur.x.voronoi.Skeleton(poly);
//		s.execute();
////		s.checkEdges();
//		ListView<XPoint> points = s.getEdges();
//
////		for (int i = 0; i < s.chainCount(); ++i) {
////			SimpleList PL = s.getChain(i);
////		}
//
//		SimpleList<XPoint> L = s.getPoints();
//		SimpleListItem i = (SimpleListItem) L.first();
//
//		PShape lines = new PShape();
//		lines.setFamily(PShape.GEOMETRY);
//		lines.setStrokeCap(ROUND);
//		lines.setStroke(true);
//		lines.setStrokeWeight(3);
//		lines.setStroke(-1232222);
//		lines.beginShape();
//
//		while (i.next() != null) {
//			XPoint x = (XPoint) i.value();
//			lines.vertex((float) x.x, (float) x.y);
//			i = (SimpleListItem) i.next();
//		}
//		lines.endShape();
//
//
////		XPoint x = (XPoint) L.next(L.first()).value();
////		x.x;
////		L.get
//
////		SimpleList l = s.getLines(true);
////		l.nex
//		return lines;
//	}

	/**
	 * 
	 * @param shape
	 * @return
	 */
	public static SolubSkeleton solubSkeleton(List<PVector> points) {

		final Coordinate[] coords;
		if (!points.get(0).equals(points.get(points.size() - 1))) {
			coords = new Coordinate[points.size() + 1];
			points.add(points.get(0)); // close geometry
		} else { // already closed
			coords = new Coordinate[points.size()];
		}

		for (int i = 0; i < coords.length; i++) {
			coords[i] = new Coordinate(points.get(i).x, points.get(i).y);
		}

		Polygon p = GEOM_FACTORY.createPolygon(coords); // reverse
		points.clear();

		for (Coordinate coordinate : p.getExteriorRing().getCoordinates()) {
			points.add(new PVector((float) coordinate.x, (float) coordinate.y));
		}
		points.remove(points.size() - 1); // remove closing point

		SolubSkeleton skeleton = new SolubSkeleton(points, 20);
		return skeleton;
	}

	/**
	 * Straight skeleton. Not robust, but fast. Does not support holes
	 * 
	 * @param shape     a hull
	 * @param tolerance minimum closeness that skeleton "bone" is to nearest vertex
	 * @return SS object
	 * @see #straightSkeleton(PShape) straightSkeleton(PShape). Fully robust, but
	 *      slower.
	 */
	public static SolubSkeleton solubSkeleton(PShape shape, float tolerance) {

		ArrayList<PVector> points = new ArrayList<>();

		Polygon p = (Polygon) fromPShape(shape);

		// exterior ring is clockwise, so reverse() to get anti-clockwise
		for (Coordinate coordinate : p.getExteriorRing().reverse().getCoordinates()) {
			points.add(new PVector((float) coordinate.x, (float) coordinate.y));
		}
		points.remove(0); // remove closing point
		points.remove(0); // remove closing point

		SolubSkeleton skeleton = new SolubSkeleton(points, tolerance);
		return skeleton;
	}

	/**
	 * Roughly, it is the geometric graph whose edges are the traces of vertices of
	 * shrinking mitered offset curves of the polygon
	 * 
	 * @param shape
	 * @param p
	 * @return
	 */
	public static PShape straightSkeleton(PShape shape, PApplet p) {
		// https://github.com/Agent14zbz/ZTools/blob/main/src/main/java/geometry/ZSkeleton.java

		Machine speed = new Machine(1); // every edge same speed

		Geometry g = fromPShape(shape);
		Polygon polygon;
		if (g.getGeometryType() == Geometry.TYPENAME_POLYGON) {
			polygon = (Polygon) g;
			if (polygon.getCoordinates().length > 1000) {
				polygon = (Polygon) DouglasPeuckerSimplifier.simplify(polygon, 1);
			}
		} else {
			System.out.println("MultiPolygon not supported yet.");
			return new PShape();
		}

		HashSet<Double> edgeCoordsSet = new HashSet<>();

		Skeleton skeleton;

		LoopL<org.twak.camp.Edge> loopL = new LoopL<>(); // list of loops

		ArrayList<Corner> corners = new ArrayList<>();
		Loop<org.twak.camp.Edge> loop = new Loop<>();

		LinearRing exterior = polygon.getExteriorRing();
		if (!Orientation.isCCW(exterior.getCoordinates())) {
			exterior = exterior.reverse();
		}
		for (int j = 0; j < exterior.getCoordinates().length - 1; j++) {
			double a = exterior.getCoordinates()[j].x;
			double b = exterior.getCoordinates()[j].y;
			corners.add(new Corner(a, b));
			edgeCoordsSet.add(cantorPairing(a, b));
		}
		for (int j = 0; j < corners.size() - 1; j++) {
			org.twak.camp.Edge edge = new org.twak.camp.Edge(corners.get(j),
					corners.get((j + 1) % (corners.size() - 1)));
			edge.machine = speed;
			loop.append(edge);
		}
		loopL.add(loop);

		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			corners = new ArrayList<>();
			// holes should be clockwise
			LinearRing hole = polygon.getInteriorRingN(i).reverse();
//				System.out.println("hole:" + Orientation.isCCW(hole.getCoordinates()));
			for (int j = 0; j < hole.getNumPoints() - 1; j++) {
				corners.add(new Corner(hole.getCoordinates()[j].x, hole.getCoordinates()[j].y));
			}
			loop = new Loop<>();
			for (int j = 0; j < corners.size() - 1; j++) {
				org.twak.camp.Edge edge = new org.twak.camp.Edge(corners.get(j),
						corners.get((j + 1) % (corners.size() - 1)));
				edge.machine = speed;
				loop.append(edge);
			}
			loopL.add(loop);
		}

//		}

		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(ROUND);
		lines.setStroke(true);
		lines.setStrokeWeight(4);
		lines.setStroke(-1232222);
		lines.beginShape(LINES);
		try {
			skeleton = new Skeleton(loopL, true);
			skeleton.skeleton();
			skeleton.output.edges.map.values().forEach(e -> {
				boolean a = edgeCoordsSet.contains(cantorPairing(e.start.x, e.start.y));
				boolean b = edgeCoordsSet.contains(cantorPairing(e.end.x, e.end.y));

				if (a ^ b) { // branch
					lines.vertex((float) e.start.x, (float) e.start.y);
					lines.vertex((float) e.end.x, (float) e.end.y);
				} else {
					if (a) { // edge
//						p.stroke(0, 50, 180);
//						p.line((float) e.start.x, (float) e.start.y, (float) e.end.x, (float) e.end.y);

					} else { // bone
						lines.vertex((float) e.start.x, (float) e.start.y);
						lines.vertex((float) e.end.x, (float) e.end.y);
//						p.stroke(0);
//						p.strokeWeight(5);
//						p.point((float) e.start.x, (float) e.start.y);
					}
				}
			});

//			skeleton.output.faces.values().forEach(face -> {
//				face.edges.forEach(edge -> {
//					int n = 0;
//					for (SharedEdge l : edge) {
//						p.stroke(12, 56, 200);
//						if (n == edge.count() - 2) { // last edge belongs to the shape exterior
//							p.stroke(120, 156, 100);
//						}
//						p.line((float) l.start.x, (float) l.start.y, (float) l.end.x, (float) l.end.y);
//						n++;
//					}
//				});
//			});

//			skeleton.output.faces.values().forEach(f -> {
//				f.getLoopL().forEach(l -> {
//					l.forEach(v -> {
//						lines.vertex((float) v.x, (float) v.y);
//					});
//				});
//				final org.twak.camp.Edge e = f.edge;
//				lines.vertex((float) e.start.x, (float) e.start.y);
//				lines.vertex((float) e.end.x, (float) e.end.y);
//				f.topSE.forEach(e2 -> {
//					lines.vertex((float) e2.start.x, (float) e2.start.y);
//					lines.vertex((float) e2.end.x, (float) e2.end.y);
//				});
//			});
		} catch (Exception ignore) {
			// hide init or collision errors from console
		}

		lines.endShape();
		return lines;
	}

	/**
	 * TODO replace contained in straight skeleton within straightSkeleton(PShape,
	 * PApplet)
	 * 
	 * @param g
	 * @return
	 */
	private Skeleton skeletonise(Geometry g) {

		Polygon polygon;
		if (g.getGeometryType() == Geometry.TYPENAME_POLYGON) {
			polygon = (Polygon) g;
			if (polygon.getCoordinates().length > 1000) {
				polygon = (Polygon) DouglasPeuckerSimplifier.simplify(g, 1);
			}
		} else {
			System.out.println("MultiPolygon not supported yet.");
			return new Skeleton();
		}

		Machine speed = new Machine(1); // every edge same speed

		HashSet<Double> edgeCoordsSet = new HashSet<>();

		Skeleton skeleton;

		LoopL<org.twak.camp.Edge> loopL = new LoopL<>(); // list of loops

		ArrayList<Corner> corners = new ArrayList<>();
		Loop<org.twak.camp.Edge> loop = new Loop<>();

		LinearRing exterior = polygon.getExteriorRing();
		if (!Orientation.isCCW(exterior.getCoordinates())) {
			exterior = exterior.reverse();
		}
		Coordinate[] coords = exterior.getCoordinates();
		for (int j = 0; j < coords.length - 1; j++) {
			double a = coords[j].x;
			double b = coords[j].y;
			corners.add(new Corner(a, b));
			edgeCoordsSet.add(cantorPairing(a, b));
		}
		for (int j = 0; j < corners.size() - 1; j++) {
			org.twak.camp.Edge edge = new org.twak.camp.Edge(corners.get(j),
					corners.get((j + 1) % (corners.size() - 1)));
			edge.machine = speed;
			loop.append(edge);
		}
		loopL.add(loop);

		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			corners = new ArrayList<>();
			// holes should be clockwise
			LinearRing hole = polygon.getInteriorRingN(i).reverse();
//				System.out.println("hole:" + Orientation.isCCW(hole.getCoordinates()));
			for (int j = 0; j < hole.getNumPoints() - 1; j++) {
				corners.add(new Corner(hole.getCoordinates()[j].x, hole.getCoordinates()[j].y));
			}
			loop = new Loop<>();
			for (int j = 0; j < corners.size() - 1; j++) {
				org.twak.camp.Edge edge = new org.twak.camp.Edge(corners.get(j),
						corners.get((j + 1) % (corners.size() - 1)));
				edge.machine = speed;
				loop.append(edge);
			}
			loopL.add(loop);
		}

		skeleton = new Skeleton(loopL, true);
		skeleton.skeleton();

		return skeleton;
	}

	/**
	 * Generate a Topographic-like isoline contour map of the shape. The "elevation"
	 * (or z values) of points is the euclidean distance between a point in the
	 * shape and the given high point.
	 * <p>
	 * Assigns each point feature a number equal to the distance between geometry's
	 * centroid and the point.
	 * 
	 * @param shape
	 * @param highPoint
	 * @param intervalSpacing
	 * @return
	 */
	public static PShape isolines(PShape shape, PVector highPoint, float intervalSpacing) {

		/**
		 * Also See:
		 * https://github.com/hageldave/JPlotter/blob/master/jplotter/src/main/java/hageldave/jplotter/misc/Contours.java
		 * https://blog.bruce-hill.com/meandering-triangles
		 * http://indiemaps.com/blog/2008/06/isolining-package-for-actionscript-3/
		 */

		// https://github.com/hageldave/JPlotter/blob/master/jplotter/src/main/java/hageldave/jplotter/misc/Contours.java

		Geometry g = fromPShape(shape);
		if (g.getCoordinates().length > 2000) {
			g = DouglasPeuckerSimplifier.simplify(g, 2);
		}
		int buffer = Math.max(10, Math.round(intervalSpacing) + 1);
		PreparedGeometry cache = PreparedGeometryFactory.prepare(g.buffer(10));

		List<Vertex> tinVertices = new ArrayList<Vertex>(200);
		double maxDist = 0;

		PoissonDistribution pd = new PoissonDistribution(0);
		Coordinate[] e = g.getEnvelope().getCoordinates(); // envelope/bounding box of shape

		/**
		 * Poisson a little faster, but isolines are more rough
		 */
//		ArrayList<PVector> randomPoints = pd.generate(e[0].x - buffer, e[0].y - buffer, e[3].x + buffer,
//				e[1].y + buffer, intervalSpacing, 6);
		ArrayList<PVector> randomPoints = generateGrid(e[0].x - buffer, e[0].y - buffer, e[3].x + buffer,
				e[1].y + buffer, intervalSpacing, intervalSpacing);

		for (PVector v : randomPoints) {
			/**
			 * Major bottleneck is isoline computation so reduce points to only those
			 * needed.
			 */
			if (cache.covers(PTS.pointFromPVector(v))) {
				double d = highPoint.dist(v);
				maxDist = Math.max(d, maxDist);
				tinVertices.add(new Vertex(v.x, v.y, d, 0));
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

		/**
		 * A null valuator tells the builder to just use the z values from the vertices
		 * rather than applying any adjustments to their values.
		 */
		final ContourBuilderForTin builder = new ContourBuilderForTin(tin, null, intervals, true);

		List<org.tinfour.contour.Contour> contours = builder.getContours();

		PShape parent = new PShape(PConstants.GROUP);
		parent.setKind(PConstants.GROUP);

		LineDissolver ld = new LineDissolver();
		int q = 0;
		for (org.tinfour.contour.Contour contour : contours) {
			Coordinate[] coords = new Coordinate[contour.getCoordinates().length / 2];
			for (int i = 0; i < contour.getCoordinates().length; i += 2) {
				float vx = (float) contour.getCoordinates()[i];
				float vy = (float) contour.getCoordinates()[i + 1];
				coords[i / 2] = new Coordinate(vx, vy);
			}
			ld.add(GEOM_FACTORY.createLineString(coords));
		}

		return toPShape(DouglasPeuckerSimplifier.simplify(ld.getResult(), 2).intersection(g)); // contains check
																								// instead?
	}

	/**
	 * 
	 * @param points               List of PVectors: the z coordinate for each
	 *                             PVector should define the contour height at that
	 *                             point
	 * @param intervalValueSpacing contour height distance represented by successive
	 *                             isolines
	 * @param isolineMin           minimum value represented by isolines
	 * @param isolineMax           maximum value represented by isolines
	 * @return
	 */
	public static PShape isolines(List<PVector> points, float intervalValueSpacing, float isolineMin,
			float isolineMax) {
		// lines = max-min/spacing
		final IncrementalTin tin = new IncrementalTin(10);
		points.forEach(point -> {
			tin.add(new Vertex(point.x, point.y, point.z));
		});

		double[] intervals = generateDoubleSequence(isolineMin, isolineMax, intervalValueSpacing);
		ContourBuilderForTin builder = null;
//		try {
		// catch org.tinfour.contour.PerimeterLink.addContourTip error if any vertex has
		// near-zero coordinate
		builder = new ContourBuilderForTin(tin, null, intervals, true);
//		} catch (Exception e) {
//			return new PShape();
//		}

		List<org.tinfour.contour.Contour> contours = builder.getContours();

		PShape parent = new PShape(PConstants.GROUP);
		parent.setKind(PConstants.GROUP);

		LineDissolver ld = new LineDissolver();
		for (org.tinfour.contour.Contour contour : contours) {
			Coordinate[] coords = new Coordinate[contour.getCoordinates().length / 2];
			for (int i = 0; i < contour.getCoordinates().length; i += 2) {
				float vx = (float) contour.getCoordinates()[i];
				float vy = (float) contour.getCoordinates()[i + 1];
				coords[i / 2] = new Coordinate(vx, vy);
			}
			ld.add(GEOM_FACTORY.createLineString(coords));
		}

		return toPShape(ld.getResult()); // contains check instead?
	}

	public static PShape isolinesJP(double[][] values, double isoValue) {
		PShape lines = new PShape(PShape.GEOMETRY);
		lines.setStroke(true);
		lines.setStrokeWeight(4);
		lines.setStroke(-1232222);
		lines.beginShape(LINES);
//		Contours.
		LineMerger m = new LineMerger();

		List<SegmentDetails> segments = Contours.computeContourLines(values, isoValue, RGB.composeclr(0, 255, 0, 255));
		segments.forEach(s -> {
			lines.vertex((float) s.p0.getX() * 2, (float) s.p0.getY() * 2);
			lines.vertex((float) s.p1.getX() * 2, (float) s.p1.getY() * 2);
//			m.add(GEOM_FACTORY.createLineString(new Coordinate[] { new Coordinate(s.p0.getX(), s.p0.getY()),
//					new Coordinate(s.p1.getX(), s.p1.getY()) }));
		});
//		m.getMergedLineStrings();
//		System.out.println(segments.size());
		lines.endShape();
		return lines;
	}

	/**
	 * 
	 * @param shape
	 * @param spacing spacing between offsets. should be >=1
	 * @param p
	 * @return
	 * @see #miteredOffset(PShape, double)
	 */
	public static PShape miteredOffset(PShape shape, double spacing) {
		Geometry g = fromPShape(shape);

		if (g.getCoordinates().length > 2000) {
			g = DouglasPeuckerSimplifier.simplify(g, 0.2);
		}

		final int joinStyle = BufferParameters.JOIN_MITRE; // TODO as input argument

		BufferParameters bufParams = new BufferParameters(4, BufferParameters.CAP_FLAT, joinStyle,
				BufferParameters.DEFAULT_MITRE_LIMIT);

//		bufParams.setSimplifyFactor(5); // can produce "poor" yet interesting results

		spacing = -Math.max(1, Math.abs(spacing)); // ensure negative and >=1

		/*
		 * Is is faster to buffer the shape inwards than to use the straight skeleton as
		 * the basis for mitered offsets. To this end a stack is used to handle times
		 * where the geometry breaks into multiple sub-geometries. The buffer is
		 * effectively applied in a recursive manner until all geometries are empty (in
		 * which case they are ignored for the next iteration).
		 */
		final ArrayDeque<Geometry> geometries = new ArrayDeque<>();
		geometries.push(g);

		PShape parent = new PShape(PConstants.GROUP);
		while (!geometries.isEmpty()) {
			g = geometries.pop();
			Coordinate[] coords = g.getCoordinates();
			if (coords.length > 1) {
				PShape lines = new PShape();
				lines.setFamily(PShape.PATH);
				lines.setStroke(true);
				lines.setStrokeWeight(2);
//				lines.setStroke(
//						RGB.composeclr((int) PApplet.map((float) g.getCoordinates()[0].x, 0, p.width, 0, 254), n * 10 % 255,
//								150,
//								255));
				lines.setStroke(-1232222);
				lines.beginShape();
				for (int i = 0; i < coords.length; i++) {
					Coordinate coord = coords[i];
					lines.vertex((float) coord.x, (float) coord.y);
				}
				lines.endShape();
				parent.addChild(lines);

				BufferOp b = new BufferOp(g, bufParams);
				g = b.getResultGeometry(spacing);

				for (int i = 0; i < g.getNumGeometries(); i++) {
					if (joinStyle == BufferParameters.JOIN_MITRE || coords.length < 20) {
						geometries.add(g.getGeometryN(i));
					} else {
						// because rounded miters produce LOADS of dense vertices
						geometries.add(DouglasPeuckerSimplifier.simplify(g.getGeometryN(i), 0.1));
					}
				}
			}
		}

		return parent;
	}

	/**
	 * A mitered offset that emanates from the shape
	 * 
	 * @param shape
	 * @param spacing
	 * @param curves  number of offset curves (including the original shape outline)
	 * @return A group PShape, where each child shape is a
	 * @see #miteredOffset(PShape, double)
	 */
	public static PShape miteredOffsetOutwards(PShape shape, double spacing, final int curves) {
		Geometry g = fromPShape(shape);

		if (g.getCoordinates().length > 2000) {
			g = DouglasPeuckerSimplifier.simplify(g, 0.5);
		}

		final int joinStyle = BufferParameters.JOIN_MITRE; // TODO as input argument

		BufferParameters bufParams = new BufferParameters(4, BufferParameters.CAP_FLAT, joinStyle,
				BufferParameters.DEFAULT_MITRE_LIMIT);

//		bufParams.setSimplifyFactor(5); // can produce "poor" yet interesting results

		spacing = Math.max(1, Math.abs(spacing)); // ensure positive and >=1

		final ArrayDeque<Geometry> geometries = new ArrayDeque<>();
		geometries.push(g);

		PShape parent = new PShape(PConstants.GROUP);
		int currentCurves = 0;
		while (!geometries.isEmpty() && currentCurves < curves) {
			g = geometries.poll();

			// need to iterate over rings individually
			for (LinearRing ring : new LinearRingIterator(g)) {
				Coordinate[] coords = ring.getCoordinates();
				if (coords.length > 1) {
					PShape lines = new PShape(PShape.PATH);
					lines.setStroke(true);
					lines.setStrokeWeight(2);
					lines.setStroke(-1232222);
					lines.beginShape();

					for (int i = 0; i < coords.length; i++) {
						Coordinate coord = coords[i];
						lines.vertex((float) coord.x, (float) coord.y);
					}
					lines.endShape();
					parent.addChild(lines);
				}
			}

			BufferOp b = new BufferOp(g, bufParams);
			g = b.getResultGeometry(spacing);
			for (int i = 0; i < g.getNumGeometries(); i++) {
				if (joinStyle == BufferParameters.JOIN_MITRE) {
					geometries.add(g.getGeometryN(i));
				} else {
					// because rounded miters produce LOADS of dense vertices
					geometries.add(DouglasPeuckerSimplifier.simplify(g.getGeometryN(i), 0.1));
				}
			}
			currentCurves++;
		}

		return parent;
	}

	/**
	 * Uniquely encodes two numbers (order-dependent) into a single natural number.
	 */
	private static double cantorPairing(double a, double b) {
		a = (a >= 0.0 ? 2.0 * a : (-2.0 * a) - 1.0); // enable negative input values
		b = (b >= 0.0 ? 2.0 * b : (-2.0 * b) - 1.0); // enable negative input values
		return (a + b) * (a + b + 1) / 2 + a;
	}

	/**
	 * Start inclusive; end exclusive
	 */
	public static double[] generateDoubleSequence(double start, double end, double step) {
		double[] sequence = new double[(int) Math.ceil((end - start) / step)];
		for (int i = 0; i < sequence.length; i++) {
			sequence[i] = start + i * step;
		}
		return sequence;
	}

	/**
	 * Generates a grid of points
	 * 
	 * @param minX
	 * @param minY
	 * @param maxX
	 * @param maxY
	 * @param spacingX
	 * @param spacingY
	 * @return
	 */
	private static ArrayList<PVector> generateGrid(double minX, double minY, double maxX, double maxY, double spacingX,
			double spacingY) {
		ArrayList<PVector> grid = new ArrayList<>();
		double[] y = generateDoubleSequence(minY, maxY, spacingY);
		double[] x = generateDoubleSequence(minX, maxX, spacingX);

		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < y.length; j++) {
				grid.add(new PVector((float) x[i], (float) y[j]));
			}
		}
		return grid;

	}

}
