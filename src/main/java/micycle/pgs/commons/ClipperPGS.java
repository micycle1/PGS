package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.util.GeometryExtracter;
import org.locationtech.jts.geom.util.LinearComponentExtracter;

import clipper2.Clipper;
import clipper2.core.ClipType;
import clipper2.core.FillRule;
import clipper2.core.Path64;
import clipper2.core.Paths64;
import clipper2.core.Point64;
import clipper2.engine.Clipper64;
import clipper2.engine.PolyPath64;
import clipper2.engine.PolyPathBase;
import clipper2.engine.PolyTree64;
import clipper2.offset.ClipperOffset;
import clipper2.offset.EndType;
import clipper2.offset.JoinType;

/**
 * Performs geometric boolean operations on JTS Geometries using Clipper2
 * subroutines, which are generally faster (on geometries with fewer vertices,
 * at least).
 * 
 * @author Michael Carleton
 *
 */
public class ClipperPGS {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));
	private static final double PRECISION = 10;

	private ClipperPGS() {
	}

	public static Geometry union(Geometry subject, Geometry clip) {
		return execute(subject, clip, ClipType.Union, FillRule.EvenOdd);
	}

	public static Geometry intersect(Geometry subject, Geometry clip) {
		return execute(subject, clip, ClipType.Intersection, FillRule.EvenOdd);
	}

	/**
	 * Intersects two polygons.
	 * <p>
	 * This method has less overhead than {@link #intersect(Geometry, Geometry)
	 * intersect()}, to be used when both geometries are known to be fairly simple
	 * polygons (at most one hole between them).
	 * 
	 * @param subject a polygon
	 * @param clip    a polygon
	 * @return
	 */
	public static Geometry intersectPolygons(Geometry subject, Geometry clip) {
		Clipper64 clipper = new Clipper64();
		Paths64 closedPaths = new Paths64();
		supplySubjects(subject, clipper);
		clipper.AddClip(toClipper(clip));
		clipper.Execute(ClipType.Intersection, FillRule.EvenOdd, closedPaths);
		return fromClipper(closedPaths, true);

	}

	public static Geometry difference(Geometry subject, Geometry clip) {
		return execute(subject, clip, ClipType.Difference, FillRule.EvenOdd);
	}

	public static Geometry symDifference(Geometry subject, Geometry clip) {
		return execute(subject, clip, ClipType.Xor, FillRule.EvenOdd);
	}

	public static Geometry execute(Geometry subject, Geometry clip, ClipType clipType, FillRule fillRule) {
		Clipper64 clipper = new Clipper64();

		/*
		 * If both geometries are polygonal and contain holes, output might be
		 * complicated (having nested polygons each with holes), so use polytree to
		 * preserve this structure faithfully during conversion to Geometry.
		 */
		int polysWithHoles = 0;
		countHoles: {
			for (Object p : GeometryExtracter.extract(subject, Geometry.TYPENAME_POLYGON)) {
				Polygon polygon = (Polygon) p;
				if (polygon.getNumInteriorRing() > 0) {
					polysWithHoles++;
					if (polysWithHoles > 1) {
						break countHoles;
					}
				}
			}
			if (clip != null) {
				for (Object p : GeometryExtracter.extract(clip, Geometry.TYPENAME_POLYGON)) {
					Polygon polygon = (Polygon) p;
					if (polygon.getNumInteriorRing() > 0) {
						polysWithHoles++;
					}
					if (polysWithHoles > 1) {
						break countHoles;
					}
				}
			}
		}

		boolean usePolyTree = polysWithHoles > 1; // use polytree when 2 or more polygons each have holes

		supplySubjects(subject, clipper);
		if (clip != null) {
			Paths64 clipContours = toClipper(clip);
			clipper.AddClip(clipContours);
		}

		Paths64 openPaths = new Paths64(); // unclosed lines (linework)

		if (usePolyTree) {
			PolyTree64 tree = new PolyTree64();
			clipper.Execute(clipType, fillRule, tree, openPaths);
			return fromPolyTree(tree, openPaths);
		} else {
			Paths64 closedPaths = new Paths64(); // closed paths -- outer polygon contours & holes
			clipper.Execute(clipType, fillRule, closedPaths, openPaths);
			return fromSolutionPaths(closedPaths, openPaths);
		}
	}

	static Geometry minkSum(Geometry pattern, Geometry path) {
		return fromClipper(Clipper.MinkowskiSum(linealToPath(pattern), linealToPath(path), true), true);
	}

	static Geometry minkDifference(Geometry pattern, Geometry path) {
		return fromClipper(Clipper.MinkowskiDiff(linealToPath(pattern), linealToPath(path), true), true);
	}

	static Geometry buffer(Geometry geometry, double delta, JoinType joinType, EndType endType) {
		Paths64 buffered = Clipper.InflatePaths(toClipper(geometry), delta * PRECISION, joinType, endType);
		return fromSolutionPaths(new Paths64(), new Paths64(buffered.subList(0, 1)));
	}

	static Geometry offset(Geometry geometry, double delta, JoinType joinType, EndType endType) {
		ClipperOffset offset = new ClipperOffset();
		offset.AddPaths(toClipper(geometry), joinType, endType);
		return fromSolutionPaths(offset.Execute(delta * PRECISION), new Paths64());
	}

	private static Geometry fromPolyTree(PolyTree64 tree, Paths64 openPaths) {
		List<Polygon> polygons = new ArrayList<>();
		for (PolyPathBase polyPath : tree) { // children at tree level, starting with polyon exteriors
			fromPolyPath((PolyPath64) polyPath, polygons);
		}

		if (polygons.isEmpty()) {
			if (openPaths.isEmpty()) {
				return GEOM_FACTORY.createEmpty(2);
			} else {
				return fromClipper(openPaths, false);
			}
		} else {
			Geometry polygonal;
			if (polygons.size() == 1) {
				polygonal = polygons.get(0);
			} else {
				polygonal = GEOM_FACTORY.createMultiPolygon(polygons.toArray(new Polygon[0]));
			}

			if (openPaths.isEmpty()) {
				return polygonal;
			} else {
				return GEOM_FACTORY.createGeometryCollection(new Geometry[] { polygonal, fromClipper(openPaths, false) });
			}
		}

	}

	/**
	 * 
	 * @param closedPaths closed paths solutions, which correspond to outer polygon
	 *                    contours and holes
	 * @param openPaths   open paths solutions
	 * @return
	 */
	private static Geometry fromSolutionPaths(Paths64 closedPaths, Paths64 openPaths) {
		if (closedPaths.isEmpty()) {
			if (openPaths.isEmpty()) {
				return GEOM_FACTORY.createEmpty(2);
			} else {
				return fromClipper(openPaths, false);
			}
		} else {
			if (openPaths.isEmpty()) {
				return fromClipper(closedPaths, true);
			} else {
				return GEOM_FACTORY
						.createGeometryCollection(new Geometry[] { fromClipper(closedPaths, true), fromClipper(openPaths, false) });
			}
		}
	}

	private static void fromPolyPath(PolyPath64 polyPath, List<Polygon> polygons) {
		LinearRing exterior = GEOM_FACTORY.createLinearRing(toCoord(polyPath.getPolygon(), true));
		LinearRing[] holes = new LinearRing[polyPath.getCount()];

		int holeIndex = 0;
		for (PolyPathBase hole : polyPath) {
			PolyPath64 hole64 = (PolyPath64) hole;
			holes[holeIndex++] = GEOM_FACTORY.createLinearRing(toCoord(hole64.getPolygon(), true));
			if (hole.getCount() > 0) {
				hole.forEach(subExterior -> fromPolyPath((PolyPath64) subExterior, polygons));
			}
		}

		polygons.add(GEOM_FACTORY.createPolygon(exterior, holes));
	}

	/**
	 * 
	 * @param solution
	 * @param closed   necessary because clipper outputs unclosed lines for
	 * @return
	 */
	private static Geometry fromClipper(List<Path64> solution, boolean closed) {
		if (closed) { // polygonal
			/*
			 * NOTE The subroutine below assumes the solution represents a single polygon
			 * exterior with any number holes. It may be the case that the solution actually
			 * corresponds to N polygon exteriors (each possibly having their own holes), in
			 * which case the shape structure of the output is wrong (a polytree is needed
			 * to preserve such structures during conversion).
			 */
			LinearRing exterior = GEOM_FACTORY.createLinearRing(toCoord(solution.get(0), true));
			LinearRing[] innerRings = new LinearRing[solution.size() - 1];
			for (int i = 0; i < innerRings.length; i++) {
				innerRings[i] = GEOM_FACTORY.createLinearRing(toCoord(solution.get(i + 1), true));
			}
			return GEOM_FACTORY.createPolygon(exterior, innerRings);
		} else { // lineal
			if (solution.size() > 1) {
				LineString[] lineStrings = new LineString[solution.size()];
				for (int i = 0; i < lineStrings.length; i++) {
					lineStrings[i] = GEOM_FACTORY.createLineString(toCoord(solution.get(i), false));
				}
				return GEOM_FACTORY.createMultiLineString(lineStrings);
			} else {
				return GEOM_FACTORY.createLineString(toCoord(solution.get(0), false));
			}
		}
	}

	@SuppressWarnings("unchecked")
	/**
	 * Adds (open) subject contours comprising the Geometry to the given clipper
	 * object, depending on whether they are open or not.
	 */
	private static void supplySubjects(Geometry g, Clipper64 clipper) {
		LinearComponentExtracter.getLines(g).forEach(l -> {
			LineString lineString = (LineString) l;
			Path64 path = linealToPath(lineString);
			if (lineString.isClosed()) {
				clipper.AddSubject(path);
			} else {
				clipper.AddOpenSubject(path);
			}
		});
	}

	@SuppressWarnings("unchecked")
	private static Paths64 toClipper(Geometry g) {
		Paths64 paths = new Paths64(g.getNumGeometries());
		LinearComponentExtracter.getLines(g).forEach(l -> {
			LineString lineString = (LineString) l;
			paths.add(linealToPath(lineString));
		});
		return paths;
	}

	private static Path64 linealToPath(Geometry lineal) {
		Path64 contour = new Path64(lineal.getNumPoints());
		for (Coordinate c : lineal.getCoordinates()) {
			contour.add(new Point64(c.x * PRECISION, c.y * PRECISION));
		}
		return contour;
	}

	private static Coordinate[] toCoord(List<Point64> path, boolean polygon) {
		Coordinate[] coords = new Coordinate[path.size() + (polygon ? 1 : 0)];
		for (int i = 0; i < path.size(); i++) {
			Point64 p = path.get(i);
			coords[i] = new Coordinate(p.x / PRECISION, p.y / PRECISION);
		}
		if (polygon) {
			coords[coords.length - 1] = coords[0]; // close polygon coords for JTS
		}
		return coords;
	}

}
