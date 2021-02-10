package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.fromPShapeVivid;
import static micycle.pts.Conversion.toPShape;
import static micycle.pts.color.RGB.composeclr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import org.apache.log4j.Level;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.geodelivery.jap.concavehull.SnapHull;
import org.geotools.coverage.CoverageFactoryFinder;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridCoverageFactory;
import org.geotools.data.DataUtilities;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.feature.SchemaException;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.geotools.process.vector.ContourProcess;
import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.dissolve.LineDissolver;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.geomgraph.Edge;
import org.locationtech.jts.geomgraph.GeometryGraph;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.shape.random.RandomPointsBuilder;
import org.locationtech.jts.shape.random.RandomPointsInGridBuilder;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.triangulate.ConformingDelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.IncrementalDelaunayTriangulator;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;
import org.locationtech.jts.util.GeometricShapeFactory;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.twak.camp.Corner;
import org.twak.camp.Machine;
import org.twak.camp.Skeleton;
import org.twak.utils.collections.Loop;

import de.alsclo.voronoi.Voronoi;
import de.incentergy.geometry.impl.GreedyPolygonSplitter;
import earcut4j.Earcut;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineSegment;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineString;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IPolygon;
import fr.ign.cogit.geoxygene.spatial.coordgeom.GM_Polygon;
import fr.ign.cogit.geoxygene.util.algo.JtsUtil;
import fr.ign.cogit.geoxygene.util.conversion.AdapterFactory;
import fr.ign.cogit.geoxygene.util.conversion.JtsGeOxygene;
import micycle.pts.color.Blending;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;
import uk.osgb.algorithm.concavehull.ConcaveHull;
import uk.osgb.algorithm.concavehull.TriCheckerChi;
import uk.osgb.algorithm.minkowski_sum.Minkowski_Sum;

/**
 * PTS | Processing Topology Suite
 * 
 * http://thecloudlab.org/processing/library.html
 * https://github.com/Rogach/jopenvoronoi https://github.com/IGNF/CartAGen
 * https://ignf.github.io/CartAGen/docs/algorithms/others/spinalize.html
 * https://github.com/twak/campskeleton
 * https://discourse.processing.org/t/straight-skeleton-or-how-to-draw-a-center-line-in-a-polygon-or-shape/17208/9
 * http://lastresortsoftware.blogspot.com/2010/12/smooth-as.html#note1
 * 
 * TODO https://ignf.github.io/CartAGen/docs/algorithms.html TODO take into
 * account strokweight (buffer by stroke amount?) TODO maintain pshape fill, etc
 * on output
 * 
 * @author MCarleton
 */
public class PTS implements PConstants {

	// TODO check for getCoordinates() in loops (and replace) (if lots of child
	// geometries)

	/**
	 * Calling Polygon#union repeatedly is one way to union several Polygons
	 * together. But here’s a trick that can be significantly faster (seconds rather
	 * than minutes) – add the Polygons to a GeometryCollection, then apply a buffer
	 * with zero distance
	 */

	protected static final int CURVE_SAMPLES = 20;

	protected static GeometryFactory geometryFactory = new GeometryFactory();

	static {
		// stop spina/skeleton console logging
		List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
		loggers.add(LogManager.getRootLogger());
		for (Logger logger : loggers) {
			logger.setLevel(Level.OFF);
		}
	}

	/**
	 * The Maximum Inscribed Circle is determined by a point in the interior of the
	 * area which has the farthest distance from the area boundary,along with a
	 * boundary point at that distance.
	 * 
	 * @param shape
	 */
	public static PShape maximumInscribedCircle(PShape shape, float tolerance) {
		MaximumInscribedCircle mic = new MaximumInscribedCircle(fromPShape(shape), tolerance);

		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(new Coordinate(mic.getCenter().getX(), mic.getCenter().getY()));
		shapeFactory.setWidth(mic.getRadiusLine().getLength() * 2); // r*2 for total width & height
		shapeFactory.setHeight(mic.getRadiusLine().getLength() * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());

	}

	/**
	 * Return the maximum circle (at a given centerpoint inside/outside the circle)
	 * 
	 * @param shape
	 * @param centerPoint
	 * @param tolerance
	 * @return
	 */
	public static PShape maximumInscribedCircle(PShape shape, PVector centerPoint) {
		Geometry g = fromPShape(shape);
//		Polygon poly = (Polygon) g;
		Point p = pointFromPVector(centerPoint);
		Coordinate center = DistanceOp.nearestPoints(g, p)[0];

		// or .distance(point)

		double radius = distance(geometryFactory.createPoint(center), p);
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(p.getCoordinate());
		shapeFactory.setWidth(radius * 2); // r*2 for total width & height
		shapeFactory.setHeight(radius * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * The MBC is the smallest circle which completely covers the input shape (this
	 * is also known as the Smallest Enclosing Circle).
	 */
	public static PShape minimumBoundingCircle(PShape shape) {
		MinimumBoundingCircle mbc = new MinimumBoundingCircle(fromPShape(shape));
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(new Coordinate(mbc.getCentre().getX(), mbc.getCentre().getY()));
		shapeFactory.setWidth(mbc.getRadius() * 2); // r*2 for total width & height
		shapeFactory.setHeight(mbc.getRadius() * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * 
	 * @param shape
	 * @param fit   0...1
	 * @return
	 */
	private static Polygon smooth(Polygon shape, float fit) {
		return (Polygon) JTS.smooth(shape, fit);
	}

	private static Geometry smooth(Geometry shape, float fit) {
		return JTS.smooth(shape, fit);
	}

	/**
	 * @return shape A - shape B
	 */
	public static PShape difference(PShape a, PShape b) {
		return toPShape(fromPShape(a).difference(fromPShape(b)));
	}

	/**
	 * Compute the parts that the shapes do not have in common.
	 * 
	 * @return A∪B - A∩B
	 */
	public static PShape symDifference(PShape a, PShape b) {
		return toPShape(fromPShape(a).symDifference(fromPShape(b)));
	}

	/**
	 * @return A∩B
	 */
	public static PShape intersection(PShape a, PShape b) {
		PShape out = toPShape(fromPShape(a).intersection(fromPShape(b)));
//		a.draw(p.getGraphics());
		out.setFill(Blending.screen(getPShapeFillColor(a), getPShapeFillColor(b))); // TODO
		return out;
	}

	private static final int getPShapeFillColor(final PShape sh) {
		try {
			final java.lang.reflect.Field f = PShape.class.getDeclaredField("fillColor");
			f.setAccessible(true);
			return f.getInt(sh);
		} catch (ReflectiveOperationException cause) {
			throw new RuntimeException(cause);
		}
	}

	/**
	 * @return A∪B
	 * @see #union(PShape...)
	 */
	public static PShape union(PShape a, PShape b) {
		return toPShape(fromPShape(a).union(fromPShape(b)));
	}

	public static PShape union(PShape... shapes) {
		ArrayList<Geometry> geoms = new ArrayList<>();
		for (int i = 0; i < shapes.length; i++) {
			geoms.add(fromPShape(shapes[i]));
		}
		return toPShape(UnaryUnionOp.union(geoms));
	}

	public static PShape buffer(PShape shape, float buffer) {
		// TODO read
		// https://locationtech.github.io/jts/javadoc/org/locationtech/jts/operation/buffer/BufferOp.html
		return toPShape(fromPShape(shape).buffer(buffer, 4));
	}

	/**
	 * Simplifies a PShape using the Douglas-Peucker algorithm.
	 * 
	 * @param shape
	 * @param distanceTolerance
	 * @return
	 */
	public static PShape simplify(PShape shape, float distanceTolerance) {
		return toPShape(DouglasPeuckerSimplifier.simplify(fromPShape(shape), distanceTolerance));
	}

	/**
	 * Densifies a Geometry by inserting extra vertices along the line segments
	 * contained in the geometry.
	 */
	public static PShape densify(PShape shape, float distanceTolerance) {
		Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(distanceTolerance);
		d.setValidate(false);
		return toPShape(d.getResultGeometry());
	}

	/**
	 * Smoothes a geometry
	 * 
	 * @param shape
	 * @param buffer 0...1
	 * @return
	 */
	public static PShape smooth(PShape shape, float fit) {
		return toPShape(smooth(fromPShape(shape), fit));
	}

	public static PShape flatten(PShape shape) {
		// TODO iterate over geometries shapes, then create group PShape
		Polygon poly = (Polygon) fromPShape(shape).union().getGeometryN(0);
		return toPShape(poly.getExteriorRing());
	}

	public static PShape convexHull(PShape... shapes) {
		Geometry g = fromPShape(shapes[0]);
		for (int i = 1; i < shapes.length; i++) {
			g = g.union(fromPShape(shapes[i]));
		}
		return toPShape(g.convexHull());
	}

	public static PShape concaveHull(ArrayList<PVector> points, float threshold) {

		// calls Ordnance Survey implementation

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

		Geometry g = geometryFactory.createPolygon(coords);
		ConcaveHull hull = new ConcaveHull(g);
		return toPShape(hull.getConcaveHullBFS(new TriCheckerChi(threshold), false, false).get(0));
	}

	/**
	 * Has a more "organic" structure compared to other concave method.
	 * 
	 * @param points
	 * @param threshold 0...1
	 * @return
	 */
	public static PShape concaveHull2(ArrayList<PVector> points, float threshold) {

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

		Geometry g = geometryFactory.createPolygon(coords);

		// TODO test AVG threshold heuristic
		org.geodelivery.jap.concavehull.ConcaveHull hull = new org.geodelivery.jap.concavehull.ConcaveHull(threshold);

		return toPShape(hull.transform(g));
	}

	/**
	 * Adjust segment factor to change between
	 * 
	 * @param shape
	 * @param segmentFactor
	 * @return
	 */
	public static PShape snapHull(PShape shape, float segmentFactor) {
		return toPShape(SnapHull.snapHull(fromPShape(shape), segmentFactor));
	}

	/**
	 * Minkowski addition a.k.a dilation
	 * 
	 * @return
	 */
	public static PShape minkSum(PShape source, PShape addition) {
//		Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
		Geometry sum = Minkowski_Sum.compMinkSum(fromPShape(source), fromPShape(addition), false, false);
		return toPShape(sum);
	}

	/**
	 * TODO check, a.k.a erosion
	 * 
	 * @param source
	 * @param addition
	 * @return
	 */
	public static PShape minkDifference(PShape source, PShape addition) {
//			Geometry sum = Minkowski_Sum.minkSum(fromPShape(source), fromPShape(addition));
		Geometry sum = Minkowski_Sum.compMinkDiff(fromPShape(source), fromPShape(addition), false, false);
		return toPShape(sum);
	}

	public static PShape earCutTriangulation(ArrayList<PVector> points) {
		double[] arrCoords = new double[points.size() * 2];

		for (int i = 0; i < points.size(); i++) {
			arrCoords[2 * i] = points.get(i).x;
			arrCoords[2 * i + 1] = points.get(i).y;
		}
//		arrCoords = new double[] { 0, 0, 100, 0, 100, 100, 0, 100, 20, 20, 80, 20, 80, 80, 20, 80 };

		List<Integer> triangles = Earcut.earcut(arrCoords, null, 2);

		PShape triangulation = new PShape();
		triangulation.setFamily(PShape.GEOMETRY);
//		triangulation.setStrokeCap(ROUND);
		triangulation.setStroke(true);
		triangulation.setStrokeWeight(2);
		triangulation.setStroke(-123222);
		triangulation.setFill(true);
		triangulation.setFill(micycle.pts.color.RGB.composeclr(0, 0, 0, 25));
		triangulation.beginShape(TRIANGLES);
		for (int i = 0; i < triangles.size(); i += 3) {
			int v1 = triangles.get(i);
			int v2 = triangles.get(i + 1);
			int v3 = triangles.get(i + 2);
//			triangulation.vertex((float) arrCoords[v1], (float) arrCoords[v1 + 1]);
			triangulation.vertex((float) arrCoords[v1], (float) arrCoords[v1 + 1]);
			triangulation.vertex((float) arrCoords[v2], (float) arrCoords[v2 + 1]);
			triangulation.vertex((float) arrCoords[v3], (float) arrCoords[v3 + 1]);
		}
		triangulation.endShape();
		return triangulation;
	}

	/**
	 * tolerance = 0
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape delaunayTriangulation(PShape shape) {
		Geometry g = fromPShape(shape);
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setSites(g);
		Geometry out = d.getTriangles(geometryFactory); // triangulates convex hull of points
		out = out.intersection(g); // get concave hull
		return toPShape(out);
	}

	protected static Coordinate coordFromPoint(Point p) {
		return new Coordinate(p.getX(), p.getY());
	}

	/**
	 * cirumcircle center must lie inside triangle
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	private static double smallestSide(Coordinate a, Coordinate b, Coordinate c) {
		double ab = Math.sqrt((b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x));
		double bc = Math.sqrt((c.y - b.y) * (c.y - b.y) + (c.x - b.x) * (c.x - b.x));
		double ca = Math.sqrt((a.y - c.y) * (a.y - c.y) + (a.x - c.x) * (a.x - c.x));
		return Math.min(Math.min(ab, bc), ca);
	}

	static Geometry refinedTriangulation(Geometry g, int nRefinements, double tolerance) {

		DelaunayTriangulationBuilder builder = new DelaunayTriangulationBuilder();
		builder.setSites(g); // set vertex sites
		builder.setTolerance(tolerance); // set tolerance for initial triangulation only

		Geometry triangulation = builder.getTriangles(geometryFactory); // initial triangulation

		HashSet<Coordinate> sites = new HashSet<>();
		for (int i = 0; i < triangulation.getCoordinates().length; i++) {
			sites.add(triangulation.getCoordinates()[i]);
		}

		for (int refinement = 0; refinement < nRefinements; refinement++) {
			for (int i = 0; i < triangulation.getNumGeometries(); i++) {
				Polygon triangle = (Polygon) triangulation.getGeometryN(i);

				if (triangle.getArea() > 100) { // skip small triangles
					sites.add(new Coordinate(triangle.getCentroid().getX(), triangle.getCentroid().getY()));
				}
			}
			builder = new DelaunayTriangulationBuilder();
			builder.setSites(sites);
			triangulation = builder.getTriangles(geometryFactory); // re-triangulate using new centroid sites
		}

		triangulation = triangulation.intersection(g); // restore concave hull and any holes
		return triangulation;
	}

	/**
	 * Calc from vertices of PShape
	 * 
	 * @param shape
	 * @param tolerance
	 * @return
	 */
	public static PShape delaunayTriangulation(PShape shape, float tolerance) {

		/**
		 * Assume each pair of points on exterior ring to be segments,
		 */
		// http://lin-ear-th-inking.blogspot.com/2011/04/polygon-triangulation-via-ear-clipping.html
//		Geometry g = fromPShape(shape);
//		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
//		d.setTolerance(tolerance);
//		d.setSites(g);
//		Geometry out = d.getTriangles(geometryFactory); // triangulates concave hull of points
//
////		Coordinate
//
////		var x = d.getSubdivision();
//		HashSet<Coordinate> coords = new HashSet<>();
//
//		// add from g
//		for (int i = 0; i < out.getCoordinates().length; i++) {
//			coords.add(out.getCoordinates()[i]);
//		}
//
//		for (int refinement = 0; refinement < 3; refinement++) {
//			// refinement: add new points (centroids)
//			for (int i = 0; i < out.getNumGeometries(); i++) {
//				Polygon triangle = (Polygon) out.getGeometryN(i);
//
//				if (triangle.getArea() > 50) { // skip small triangles
//					coords.add(coordFromPoint(triangle.getCentroid()));
//				}
//			}
//			d = new DelaunayTriangulationBuilder();
//			d.setSites(coords);
//			out = d.getTriangles(geometryFactory); // triangulates concave hull of points
//		}
//
//		// use d.getSubdivision() to get connected to a given segment to test
//		// encroachment
//
//		out = out.intersection(g); // get convex hull
		return toPShape(refinedTriangulation(fromPShape(shape), 3, 10));

//		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
//		b.setTolerance(tolerance);
//		b.setSites(g);
//		b.setConstraints(fromPShape(createSquircle(ThreadLocalRandom.current().nextInt(400, 450), 400, 200, 200)));

//		d.getSubdivision().getVoronoiCellPolygons(geometryFactory)
//		out = b.getTriangles(geometryFactory);
	}

	public static void incrementalDelaunay() {

		// TODO
		IncrementalDelaunayTriangulator i = new IncrementalDelaunayTriangulator(new QuadEdgeSubdivision(null, 5));
		i.insertSite(null);
	}

	/**
	 * Create a triangulation with that forces certain required segments into the
	 * triangulation from a 'constrain' shape.
	 * 
	 * Steiner point is a point that is not part of the input to a geometric
	 * optimization problem but is added during the solution of the problem
	 * 
	 * @param shape
	 * @param constraints points from this shape are used as 'steiner points',
	 *                    applied to the input
	 * @param tolerance
	 * @return
	 */
	public static PShape constrainedDelaunayTriangulation(PShape shape, PShape constraints, float tolerance) {
		// TODO point-set array constraints
		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
		Geometry g = fromPShape(shape);
		b.setSites(g);
		b.setTolerance(tolerance);
		b.setConstraints(fromPShape(constraints));
		Geometry out = b.getTriangles(geometryFactory); // triangulates concave hull of points
		out = out.intersection(g); // get convex hull
		return toPShape(out);
	}

	public static PShape constrainedDelaunayTriangulation(ArrayList<PVector> points, PShape constraints,
			float tolerance) {
		ConformingDelaunayTriangulationBuilder b = new ConformingDelaunayTriangulationBuilder();
		Coordinate[] coords = new Coordinate[points.size()];
		for (int j = 0; j < points.size(); j++) {
			coords[j] = new Coordinate(points.get(j).x, points.get(j).y);
		}
		b.setSites(geometryFactory.createPolygon(coords));
		b.setTolerance(tolerance);
//		b.setConstraints(fromPShape(constraints));
		Geometry out = b.getTriangles(geometryFactory); // triangulates concave hull of points
		out = out.intersection(fromPShape(constraints)); // get convex hull
		return toPShape(out);
	}

	/**
	 * 
	 * @param points
	 * @param tolerance default is 10, snaps pixels to the nearest N. higher values
	 *                  result in more coarse triangulation
	 * @return A PShape whose children are triangles making up the triangulation
	 */
	public static PShape delaunayTriangulation(ArrayList<PVector> points, float tolerance) {
		ArrayList<Coordinate> coords = new ArrayList<>(points.size());
		for (PVector p : points) {
			coords.add(new Coordinate(p.x, p.y));
		}
		/**
		 * Use ConformingDelau... for better result? (get constrained segments from edge
		 * ring) Code example: https://github.com/locationtech/jts/issues/320
		 */
		DelaunayTriangulationBuilder d = new DelaunayTriangulationBuilder();
		d.setTolerance(tolerance);
		d.setSites(coords);
		Geometry out = d.getTriangles(geometryFactory);
		return toPShape(out);
	}

	/**
	 * Contour interval on shape
	 */
	public static PShape contour(PShape shape) {

		Geometry g = fromPShape(shape);

		ContourProcess vcp = new ContourProcess();

		SimpleFeatureType TYPE = null;
		try {
			TYPE = DataUtilities.createType("", "the_geom:Point," + "name:String," + "number:Double");
		} catch (SchemaException e) {
			e.printStackTrace();
		}

		final SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(TYPE);
		final List<SimpleFeature> features = new ArrayList<>();

		/**
		 * Generate 100 random points within the geometry. Assign each point feature a
		 * number equal to the distance between geometry's centroid and the point.
		 */
		RandomPointsBuilder r = new RandomPointsBuilder();
		r.setExtent(g);
		r.setNumPoints(100);
		Geometry points = r.getGeometry();
		final Point shapeCentroid = g.getCentroid();

		for (Coordinate coord : points.getCoordinates()) {
			Point point = geometryFactory.createPoint(coord);
			featureBuilder.add(point);
			featureBuilder.add(String.valueOf(coord.hashCode()));
			featureBuilder.add(point.distance(shapeCentroid));
			SimpleFeature feature = featureBuilder.buildFeature(null);
			features.add(feature);
		}

		SimpleFeatureCollection collection = new ListFeatureCollection(TYPE, features);

		SimpleFeatureCollection results = vcp.execute(collection, "number", new double[] { 0, 25, 50, 100, 250 }, 50d,
				true, true, null);

		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(ROUND);
		lines.setStroke(true);
		lines.setStrokeWeight(10);
		lines.setStroke(-1232222);
		lines.beginShape(LINES);

		SimpleFeatureIterator contourIterator = results.features();

		while (contourIterator.hasNext()) {
			LineString l = (LineString) contourIterator.next().getDefaultGeometry();
			lines.vertex((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY());
			lines.vertex((float) l.getEndPoint().getX(), (float) l.getEndPoint().getY());
		}

		lines.endShape();
		return lines;
	}

	/**
	 * TODO set clip envelope?
	 * 
	 * @param shape
	 * @param tolerance snapping used in underlying triangulation algorithm
	 * @return
	 */
	public static PShape voronoiDiagram(PShape shape, float tolerance) {
		Geometry g = fromPShape(shape);
		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setTolerance(tolerance);
		v.setSites(g);
		v.setClipEnvelope(new Envelope(0, 1000, 0, 750)); // speeds up when lots of edges
//		v.setSites(new ArrayList<Coordinate>(Arrays.asList(g.getCoordinates())));
		Geometry out = v.getDiagram(geometryFactory);
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
		final Geometry g = fromPShape(shape);
		final PreparedGeometry cache = PreparedGeometryFactory.prepare(g);

		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setTolerance(tolerance);
		v.setSites(g);
//		v.setClipEnvelope(new Envelope(0, 1000, 0, 750)); // TODO

		final Geometry out = v.getDiagram(geometryFactory);

		final LineDissolver ld = new LineDissolver();
		ld.add(out);
		final Geometry d = ld.getResult();

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
			if (!cache.contains(l.getStartPoint()) && !cache.contains(l.getEndPoint())) {
				lines.vertex((float) l.getStartPoint().getX(), (float) l.getStartPoint().getY());
				lines.vertex((float) l.getEndPoint().getX(), (float) l.getEndPoint().getY());
//				}
			}
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
		Geometry out = v.getDiagram(geometryFactory);
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

	public static void prims(PShape shape) { // TODO
//		PrimMinimumSpanningTree<Vertex, Edge> prims = new PrimMinimumSpanningTree<>(null);
		GeometryGraph graph = new GeometryGraph(0, fromPShape(shape));
		for (Iterator<Edge> iterator = graph.getEdgeIterator(); iterator.hasNext();) {
			Edge type = iterator.next();
//		 System.out.println(type.getCoordinate().getX());
//		 fromPShape(shape).get
		}
//	Graph<Vertex, Edge> graph = new DefaultUndirectedGraph<Vertex, Edge>(Edge.class);
//	graph.
	}

	/**
	 * Set of points in space equidistant to 2 or more points on the surface. As
	 * density of boundary points goes to infinity, a voronoi diagram converges to a
	 * medial axis.
	 * 
	 * @param shape
	 * @param density          distance tolerance for boundary densification
	 *                         (smaller values more converge towards a more accurate
	 *                         axis, but slower), 5-20 is appropriate
	 * @param minimumCloseness the SQUARE of the
	 * @return
	 */
	public static PShape medialAxis(PShape shape, float density, float minimumCloseness) {
		final Geometry g = fromPShape(shape);
		final Densifier d = new Densifier(fromPShape(shape));
		d.setDistanceTolerance(density);
		d.setValidate(false); // don't perform validation processing (a little faster)
		final Geometry dense = d.getResultGeometry();

		VoronoiDiagramBuilder v = new VoronoiDiagramBuilder();
		v.setSites(dense);
		Geometry voronoi = v.getDiagram(geometryFactory);

		final Geometry small = dense.buffer(-minimumCloseness);

		PreparedGeometry cache = PreparedGeometryFactory.prepare(small); // provides MUCH faster contains() check

		// inline lines creation
		PShape lines = new PShape();
		lines.setFamily(PShape.GEOMETRY);
		lines.setStrokeCap(ROUND);
		lines.setStroke(true);
		lines.setStrokeWeight(2);
		lines.setStroke(-1232222);
		lines.beginShape(LINES);

//		LineMerger lm = new LineMerger();
		// TODO g.difference small
//		LineDissolver // TODO
		// apply LineMergeGraph to voronoi too?

		// TODO getCoordinates() call slow on cell too?

		// TODO compare both points at once?
		for (int i = 0; i < voronoi.getNumGeometries(); i++) {
			Polygon cell = (Polygon) voronoi.getGeometryN(i); // TODO .get(0) prevent occasional crash
			for (int j = 0; j < cell.getCoordinates().length - 1; j++) {
				Coordinate a = cell.getCoordinates()[j];
//				CoordinateSequence seq = geometryFactory.getCoordinateSequenceFactory().create(new Coordinate[] { a });
				if (cache.covers(geometryFactory.createPoint(a))) {
					Coordinate b = cell.getCoordinates()[j + 1];
//					seq = geometryFactory.getCoordinateSequenceFactory().create(new Coordinate[] { b });
					if (cache.covers(geometryFactory.createPoint(b))) {
						lines.vertex((float) a.x, (float) a.y);
						lines.vertex((float) b.x, (float) b.y);
//						lm.add(geometryFactory.createLineString(new Coordinate[] { a, b }));
					}
				}
			}
		}

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
	 * Not working: "Skeleton screw up"
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape spinalize(PShape shape) {
		// also see
		// https://github.com/IGNF/geoxygene/blob/master/geoxygene-spatial/src/main/java/fr/ign/cogit/geoxygene/util/algo/geometricAlgorithms/morphomaths/MorphologyTransform.java
//		GM_Polygon p = new GM_Polygon(fromPShape(null));
		try {
			GM_Polygon polygon = (GM_Polygon) AdapterFactory.toGM_Object(fromPShapeVivid(shape));

//			AdapterFactory.log
			List<IPolygon> l = new ArrayList<>();
			l.add(polygon);
			List<ILineString> lines = Spinalize.spinalize(l, 5, 100, false);
			System.out.println(lines.size());
			PShape skeleton = new PShape();
			skeleton.setFamily(PShape.GEOMETRY);
			skeleton.setStroke(true);
			skeleton.setStrokeWeight(3);
//			skeleton.setstrok
			skeleton.beginShape(PShape.LINES);
//			for (ILineString segment : lines) {
////				segment.get
////				System.out.println(segment.getStartPoint().getX());
//				skeleton.vertex((float) segment.getStartPoint().getX(), (float) segment.getStartPoint().getY());
//				skeleton.vertex((float) segment.getEndPoint().getX(), (float) segment.getEndPoint().getY());
//			}
			skeleton.endShape();
			return skeleton;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
//			return new ArrayList<ILineString>();
		}
//		JtsGeOxygene.makeGeOxygeneGeom(fromPShapeVivid(shape));
//		JTSGeomFactory f = new JTSGeomFactory();
//		f.createIPolygon((IRing) fromPShape(shape).getExteriorRing());
	}

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

		Polygon p = geometryFactory.createPolygon(coords); // reverse
		points.clear();

		for (Coordinate coordinate : p.getExteriorRing().getCoordinates()) {
			points.add(new PVector((float) coordinate.x, (float) coordinate.y));
		}
		points.remove(points.size() - 1); // remove closing point

		SolubSkeleton skeleton = new SolubSkeleton(points, 20);
		return skeleton;
	}

	/**
	 * 
	 * @param shape     a hull
	 * @param tolerance minimum closeness that skeleton "bone" is to nearest vertex
	 * @return
	 */
	public static SolubSkeleton solubSkeleton(PShape shape, float tolerance) {

		ArrayList<PVector> points = new ArrayList<>();

		Polygon p = (Polygon) fromPShape(shape);

		for (Coordinate coordinate : p.getExteriorRing().reverse().getCoordinates()) { // reverse
			points.add(new PVector((float) coordinate.x, (float) coordinate.y));
		}
		points.remove(0); // remove closing point
		points.remove(0); // remove closing point

		SolubSkeleton skeleton = new SolubSkeleton(points, tolerance);
		return skeleton;
	}

	/**
	 * Skeleton structures of shapes are line structures that represent the shape in
	 * some sense. Works, but sloooow
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape skeletonize(PShape shape) {

		// use Angle.interiorAngle(null, null, ) for reflex vertices
		try {
			GM_Polygon polygon = (GM_Polygon) AdapterFactory.toGM_Object(fromPShapeVivid(shape));
			Set<ILineSegment> lines = Skeletonize.skeletonizeStraightSkeleton(polygon);
//			System.out.println(lines.size());
			PShape skeleton = new PShape();
			skeleton.setFamily(PShape.GEOMETRY);
			skeleton.setStroke(true);
			skeleton.setStrokeWeight(3);
//			skeleton.setstrok
			skeleton.beginShape(PShape.LINES);
			for (ILineSegment segment : lines) {
//				System.out.println(segment.getStartPoint().getX());
				skeleton.vertex((float) segment.getStartPoint().getX(), (float) segment.getStartPoint().getY());
				skeleton.vertex((float) segment.getEndPoint().getX(), (float) segment.getEndPoint().getY());
			}
			skeleton.endShape();
			for (ILineString line : lines) {
				System.out.println(line.startPoint().getX());
			}
			return skeleton;
//			

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}

	public static PShape straightSkeleton(PShape shape) {
		Polygon p = (Polygon) fromPShape(shape);
		LineString l = p.getExteriorRing().reverse(); // ? | also, remove last coord?
//		l = (LineString) DouglasPeuckerSimplifier.simplify(l, 3);

		final Machine speed = new Machine(1);

		Loop<org.twak.camp.Edge> exteriorLoop = new Loop<>();
//		LoopL<Corner> vertices = new LoopL<>();

		Corner pCorner = new Corner(l.getCoordinateN(0).x, l.getCoordinateN(0).y);
		for (int i = 1; i < l.getCoordinates().length; i++) {
			Corner corner = new Corner(l.getCoordinateN(i).x, l.getCoordinateN(i).y);
			org.twak.camp.Edge e = new org.twak.camp.Edge(pCorner, corner);
			e.machine = speed;
			exteriorLoop.append(e);
			pCorner = corner;
			System.out.println(i);
		} // not closed loop
		int end = l.getCoordinates().length - 1;
		org.twak.camp.Edge e = new org.twak.camp.Edge(new Corner(l.getCoordinateN(end).x, l.getCoordinateN(end).y),
				new Corner(l.getCoordinateN(0).x, l.getCoordinateN(0).y));
		e.machine = speed;
		exteriorLoop.append(e);

		Skeleton skeleton = new Skeleton();
		skeleton.setupForEdges(exteriorLoop.singleton());
		skeleton.skeleton();

		return null;
	}

	/**
	 * TODO: outline of exterior ring only
	 * 
	 * @param shape
	 * @param points
	 * @param offsetDistance the distance the point is offset from the
	 *                       segment(positive is to the left, negative is to the
	 *                       right)
	 * @return
	 * @see #equidistantOutlineByDistance(PShape, float, float)
	 */
	public static PVector[] equidistantOutline(PShape shape, int points, float offsetDistance) {
		Polygon p = (Polygon) fromPShape(shape); // TODO CHECK CAST (ITERATE OVER GROUP)

		LengthIndexedLine l = new LengthIndexedLine(p.getExteriorRing());
		PVector[] outlinePoints = new PVector[points * (p.getNumInteriorRing() + 1)];
		// exterior ring
		for (int i = 0; i < points; i++) {
			Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
			outlinePoints[i] = new PVector((float) q.x, (float) q.y);
		}

		for (int j = 0; j < p.getNumInteriorRing(); j++) {
			l = new LengthIndexedLine(p.getInteriorRingN(j));
			for (int i = 0; i < points; i++) {
				Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
				outlinePoints[(j + 1) * points + i] = new PVector((float) q.x, (float) q.y);
			}
		}

		return outlinePoints;
	}

	/**
	 * 
	 * @param shape
	 * @param interPointDistance Distance between each point on outline
	 * @param offsetDistance
	 * @return nearest distance such that every distance is equal (maybe be
	 *         different due to rounding)
	 */
	public static PVector[] equidistantOutlineByDistance(PShape shape, float interPointDistance, float offsetDistance) {
		Polygon p = (Polygon) fromPShape(shape); // TODO CHECK CAST (ITERATE OVER GROUP)

		LengthIndexedLine l = new LengthIndexedLine(p.getExteriorRing());
		if (interPointDistance > l.getEndIndex()) {
			System.err.println("Interpoint greater than shape length");
			return null;
		}
		int points = (int) Math.round(l.getEndIndex() / interPointDistance);
		ArrayList<PVector> outlinePoints = new ArrayList<>();
		// exterior ring
		for (int i = 0; i < points; i++) {
			Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
			outlinePoints.add(new PVector((float) q.x, (float) q.y));
		}

		for (int j = 0; j < p.getNumInteriorRing(); j++) {
			l = new LengthIndexedLine(p.getInteriorRingN(j));
			points = (int) Math.round(l.getEndIndex() / interPointDistance);
			for (int i = 0; i < points; i++) {
				Coordinate q = l.extractPoint((i / (float) points) * l.getEndIndex(), offsetDistance);
				outlinePoints.add(new PVector((float) q.x, (float) q.y));
			}
		}

		final PVector[] out = new PVector[outlinePoints.size()];
		Arrays.setAll(out, outlinePoints::get);
		return out;
	}

	/**
	 * 
	 * @param shape
	 * @param distance 0...1 around shape perimeter; or -1...0 (other direction)
	 * @return
	 */
	public static PVector pointOnOutline(PShape shape, float distance, float offsetDistance) {
		distance %= 1;
		// TODO CHECK CAST (ITERATE OVER GROUP)
		LengthIndexedLine l = new LengthIndexedLine(((Polygon) fromPShape(shape)).getExteriorRing());
		Coordinate coord = l.extractPoint(distance * l.getEndIndex(), offsetDistance);
		return new PVector((float) coord.x, (float) coord.y);
	}

	public static PVector getCentroid(PShape shape) {
		Point point = fromPShape(shape).getCentroid();
		return new PVector((float) point.getX(), (float) point.getY());
	}

	/**
	 * Returns a hole-less version of the shape, or boundary of group of shapes
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape boundary(PShape shape) {
		return toPShape(fromPShape(shape).getBoundary());
	}

	/**
	 * Get N random points contained within the PShape region. Points are
	 * distributed randomly.
	 * 
	 * @param shape
	 * @param points
	 * @return
	 */
	public static ArrayList<PVector> generateRandomPoints(PShape shape, int points) {
		// or g.getInteriorPoint()?
		RandomPointsBuilder r = new RandomPointsBuilder();
		r.setExtent(fromPShape(shape));
		r.setNumPoints(points);

		ArrayList<PVector> vertices = new ArrayList<>();

		for (Coordinate coord : r.getGeometry().getCoordinates()) {
			vertices.add(new PVector((float) coord.x, (float) coord.y));
		}
		return vertices;
	}

	/**
	 * random points generated in a grid of cells (one point randomly located in
	 * each cell) from the envelope of the shape
	 * 
	 * @param shape
	 * @param points number of points, if this shape was its own envelope
	 * @return
	 */
	public static ArrayList<PVector> generateRandomGridPoints(PShape shape, int maxPoints) {
		Geometry g = fromPShape(shape);

		RandomPointsInGridBuilder r = new RandomPointsInGridBuilder();
		r.setConstrainedToCircle(true);
		r.setExtent(g.getEnvelopeInternal());
		r.setNumPoints(maxPoints);

		ArrayList<PVector> vertices = new ArrayList<>();

		for (Coordinate coord : r.getGeometry().getCoordinates()) {
			// manually prune envelope
			if (g.contains(geometryFactory.createPoint(coord))) {
				vertices.add(new PVector((float) coord.x, (float) coord.y));
			}
		}
		return vertices;
	}

	/**
	 * 
	 * @param n    number of vertices
	 * @param xMax
	 * @param yMax
	 * @return
	 */
	public static PShape randomPolygon(int n, float xMax, float yMax) {
		return pshapeFromPVector(RandomPolygon.generateRandomConvexPolygon(n, xMax, yMax));
	}

	/**
	 * Simple polygon from coordinate list
	 */
	private static PShape pshapeFromPVector(List<PVector> coords) {
		PShape shape = new PShape();
		shape.setFamily(PShape.GEOMETRY);
//		lines.setStrokeCap(ROUND);
		shape.setStroke(true);
		shape.setStrokeWeight(2);
		shape.setStroke(0);
		shape.setFill(-123712);
		shape.setFill(true);
		shape.beginShape();

		for (PVector v : coords) {
			shape.vertex(v.x, v.y);
		}

		shape.endShape(CLOSE);
		return shape;
	}

	/**
	 * TODO return list of points when shape is group
	 * 
	 * @param shape
	 * @param point
	 * @return
	 */
	public static PVector closestVertexToPoint(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		Coordinate coord = DistanceOp.nearestPoints(g, pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

	public static boolean contains(PShape outer, PShape inner) {
		return fromPShape(outer).contains(fromPShape(inner));
	}

	public static boolean containsPoint(PShape shape, PVector point) {
		return fromPShape(shape).contains(pointFromPVector(point));
	}

	public static float distance(PShape a, PShape b) {
		return (float) fromPShape(a).distance(fromPShape(b));
	}

	public static Point getPoint(float x, float y) {
		return geometryFactory.createPoint(new Coordinate(x, y));
	}

	public static float area(PShape shape) {
		return (float) fromPShape(shape).getArea();
	}

	/**
	 * * Calculates the Miller circularity index for a polygon. This index, between
	 * 0 and 1, is equal to 1 if the polygon is perfectly circular and tends towards
	 * 0 for a segment.
	 * 
	 * @param shape
	 * @return
	 */
	public static float getCircularity(PShape shape) {
		// TODO test
		return (float) JtsUtil.circularite((com.vividsolutions.jts.geom.Polygon) (fromPShapeVivid(shape)));
	}

	/**
	 * Same as getVertices for geometry PShapes; different (subdivides) for circles
	 * etc.
	 * 
	 * @param shape
	 * @return
	 */
	public static PVector[] vertices(PShape shape) {
		Coordinate[] coords = fromPShape(shape).getCoordinates();
		PVector[] vertices = new PVector[coords.length];
		for (int i = 0; i < coords.length; i++) {
			Coordinate coord = coords[i];
			vertices[i] = new PVector((float) coord.x, (float) coord.y);
		}
		return vertices;
	}

	/**
	 * aka envelope
	 * 
	 * @param shape
	 * @return X,Y,W,H
	 */
	public static float[] bound(PShape shape) {
		Envelope e = (Envelope) fromPShape(shape).getEnvelopeInternal();
		return new float[] { (float) e.getMinX(), (float) e.getMinY(), (float) e.getWidth(), (float) e.getHeight() };
	}

	/**
	 * aka envelope
	 * 
	 * @param shape
	 * @return X1, Y1, X2, Y2
	 */
	public static float[] boundCoords(PShape shape) {
		Envelope e = (Envelope) fromPShape(shape).getEnvelopeInternal();
		return new float[] { (float) e.getMinX(), (float) e.getMinY(), (float) e.getMaxX(), (float) e.getMaxY() };
	}

	/**
	 * 
	 * @param x      centre X
	 * @param y      center Y
	 * @param width  ?total? width
	 * @param height ?total? height
	 * @return
	 */
	public static PShape createSquircle(float x, float y, float width, float height) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(40);
		shapeFactory.setCentre(new Coordinate(x, y));
//		shapeFactory.setBase(new Coordinate(x, y));
//		shapeFactory.setRotation(1);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createSquircle());
	}

	public static PShape createArcPolygon(float x, float y, float width, float height) {
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(40);
		shapeFactory.setCentre(new Coordinate(x, y));
//		shapeFactory.setBase(new Coordinate(x, y));
//		shapeFactory.setRotation(1);
		shapeFactory.setWidth(width);
		shapeFactory.setHeight(height);
		return toPShape(shapeFactory.createArcPolygon(1, 2));
	}

	/**
	 * Splits a polygon into 4 equal quadrants
	 * 
	 * @param p
	 * @return
	 */
	public static ArrayList<PShape> split(PShape shape) {
		// https://stackoverflow.com/questions/64252638/how-to-split-a-jts-polygon
		Geometry p = fromPShape(shape);
		ArrayList<PShape> ret = new ArrayList<>();

		final Envelope envelope = p.getEnvelopeInternal();
		double minX = envelope.getMinX();
		double maxX = envelope.getMaxX();
		double midX = minX + (maxX - minX) / 2.0;
		double minY = envelope.getMinY();
		double maxY = envelope.getMaxY();
		double midY = minY + (maxY - minY) / 2.0;

		Envelope llEnv = new Envelope(minX, midX, minY, midY);
		Envelope lrEnv = new Envelope(midX, maxX, minY, midY);
		Envelope ulEnv = new Envelope(minX, midX, midY, maxY);
		Envelope urEnv = new Envelope(midX, maxX, midY, maxY);
		Geometry ll = JTS.toGeometry(llEnv).intersection(p);
		Geometry lr = JTS.toGeometry(lrEnv).intersection(p);
		Geometry ul = JTS.toGeometry(ulEnv).intersection(p);
		Geometry ur = JTS.toGeometry(urEnv).intersection(p);
		ret.add(toPShape(ll));
		ret.add(toPShape(lr));
		ret.add(toPShape(ul));
		ret.add(toPShape(ur));

		return ret;
	}

	/**
	 * Splits both convex and concave shapes into a number of parts with equal area,
	 * as long as they don't have any intersecting edges and are defined by a single
	 * exterior ring.
	 * 
	 * @param shape
	 * @param parts
	 * @return
	 * @deprecated not working!?
	 */
	public static ArrayList<PShape> nSplit(PShape shape, int parts) {
		Geometry g = fromPShape(shape);
		Polygon p;
		try {
			p = (Polygon) g; // try cast to polygon
		} catch (Exception e) {
			try {
				p = (Polygon) g.buffer(0); // try buffer and cast
			} catch (Exception e2) {
				return new ArrayList<>(); // geometry is not a single/mergeable polygon
			}
		}

		System.out.println("l " + p.getCoordinates().length);

		List<Polygon> splits = new GreedyPolygonSplitter().split(p, parts);

		ArrayList<PShape> out = new ArrayList<>();
		for (Polygon polygon : splits) {
			out.add(toPShape(polygon));
		}

		return out;
	}

	/**
	 * Slice a shape using a line given by its start and endpoints
	 * 
	 * @param shape
	 * @param p1    must be outside shape
	 * @param p2    must be outside shape
	 * @return
	 */
	public static ArrayList<PShape> slice(PShape shape, PVector p1, PVector p2) {
		// https://gis.stackexchange.com/questions/189976/
		Geometry poly = fromPShape(shape);
		LineString line = lineFromPVectors(p1, p2);
		Geometry nodedLinework = poly.getBoundary().union(line);
		Geometry polys = polygonize(nodedLinework);

		// Only keep polygons which are inside the input
		ArrayList<Polygon> slices = new ArrayList<>();
		for (int i = 0; i < polys.getNumGeometries(); i++) {
			Polygon candpoly = (Polygon) polys.getGeometryN(i);
			if (poly.contains(candpoly.getInteriorPoint())) {
				slices.add(candpoly);
			}
		}

		ArrayList<PShape> output = new ArrayList<>();

		for (Polygon polygon : slices) {
			output.add(toPShape(polygon));
		}
		return output;
	}

	/**
	 * Used by slice()
	 */
	@SuppressWarnings("unchecked")
	private static Geometry polygonize(Geometry geometry) {
		List<LineString> lines = LineStringExtracter.getLines(geometry);
		Polygonizer polygonizer = new Polygonizer();
		polygonizer.add(lines);
		Collection<Polygon> polys = polygonizer.getPolygons();
		Polygon[] polyArray = GeometryFactory.toPolygonArray(polys);
		return geometry.getFactory().createGeometryCollection(polyArray);
	}

	private static double distance(Coordinate a, Coordinate b) {
		double deltaX = a.y - b.y;
		double deltaY = a.x - b.x;
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	private static double distance(Point a, Point b) {
		double deltaX = a.getY() - b.getY();
		double deltaY = a.getX() - b.getX();
		return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
	}

	private static void removeCollinearVertices(Geometry g) {
		JTS.removeCollinearVertices(g);
	}

	private static Point pointFromPVector(PVector p) {
		return geometryFactory.createPoint(new Coordinate(p.x, p.y));
	}

	private static LineString lineFromPVectors(PVector a, PVector b) {
		return geometryFactory
				.createLineString(new Coordinate[] { new Coordinate(a.x, a.y), new Coordinate(b.x, b.y) });
	}

}
