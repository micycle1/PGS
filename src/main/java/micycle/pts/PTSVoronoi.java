package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import static micycle.pts.color.RGB.composeclr;
import static micycle.pts.PTS.GEOM_FACTORY;
import static processing.core.PConstants.LINES;
import static processing.core.PConstants.ROUND;

import java.util.ArrayList;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;

import de.alsclo.voronoi.Voronoi;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Rename to Mesh?
 * 
 * @author Michael Carleton
 *
 */
public class PTSVoronoi {

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
		return toPShape(out); // .intersection(g))
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

		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(ROUND);
		lines.setStroke(true);
		lines.setStrokeWeight(2);
		lines.setStroke(composeclr(100, 150, 200, 255));
		lines.beginShape(LINES);

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

	public static PShape voronoiDiagram(ArrayList<PVector> points, float tolerance) {
		ArrayList<Coordinate> coords = new ArrayList<>(points.size());
		for (PVector p : points) {
			coords.add(new Coordinate(p.x, p.y));
		}
		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setTolerance(tolerance);
		v.setSites(coords);
		Geometry out = v.getDiagram(GEOM_FACTORY);
		return toPShape(out);
	}

	/**
	 * Doesn't use JTS so that voronoi can be applied to (sub-divided) circles. Must
	 * manually insert boundary
	 */
	public static PShape voronoiDiagram2(PShape shape) {
		Geometry g = fromPShape(shape); // convert to geom to avoid "only works with PATH or GEOMETRY shapes"
		ArrayList<de.alsclo.voronoi.graph.Point> graphIn = new ArrayList<>();
		for (int i = 0; i < g.getCoordinates().length; i++) {
			Coordinate point = g.getCoordinates()[i];
			graphIn.add(new de.alsclo.voronoi.graph.Point(point.x, point.y));
		}

//		graphIn.add(new de.alsclo.voronoi.graph.Point(1000, 750));
//		graphIn.add(new de.alsclo.voronoi.graph.Point(1000, 0));
//		graphIn.add(new de.alsclo.voronoi.graph.Point(0, 0));
//		graphIn.add(new de.alsclo.voronoi.graph.Point(0, 750));
		// TODO ^

		Voronoi voronoi = new Voronoi(graphIn);
//		voronoi.applyBoundingBox(0, 0, 1000, 750);

		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(ROUND);
		lines.setStroke(true);
		lines.setStrokeWeight(2);
		lines.setStroke(-1232222);
		lines.beginShape(LINES);
		voronoi.getGraph().edgeStream().forEach(edge -> {
			if (edge.getA() != null && edge.getB() != null) {
				lines.vertex((float) edge.getA().getLocation().x, (float) edge.getA().getLocation().y);
				lines.vertex((float) edge.getB().getLocation().x, (float) edge.getB().getLocation().y);
			}
			// TODO doesn't draw edges outside bounds :Z

		});
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

		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(ROUND);
		lines.setStroke(true);
		lines.setStrokeWeight(2);
		lines.setStroke(-1232222);
		lines.beginShape(LINES);
		voronoi.getGraph().edgeStream().forEach(edge -> {
			if (edge.getA() != null && edge.getB() != null) {
				lines.vertex((float) edge.getA().getLocation().x, (float) edge.getA().getLocation().y);
				lines.vertex((float) edge.getB().getLocation().x, (float) edge.getB().getLocation().y);
			}

		});
		lines.endShape();
		return lines;

	}

}
