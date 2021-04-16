package micycle.pgs;

import static micycle.pgs.PGS.CURVE_SAMPLES;
import static micycle.pgs.PGS.GEOM_FACTORY;
import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.algorithm.construct.MaximumInscribedCircle;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.locationtech.jts.util.GeometricShapeFactory;

import processing.core.PShape;
import processing.core.PVector;

/**
 * Bounding volumes, inscribed areas, optimal distances, etc.
 * 
 * @author Michael Carleton
 *
 */
public class PGS_GeometricOptimisation {

	private PGS_GeometricOptimisation() {
	}

	/**
	 * The Maximum Inscribed Circle is determined by a point in the interior of the
	 * area which has the farthest distance from the area boundary, along with a
	 * boundary point at that distance.
	 * 
	 * @param shape
	 * @param tolerance the distance tolerance for computing the centre point
	 *                  (around 1)
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
	 * @return A circular PShape
	 */
	public static PShape maximumInscribedCircle(PShape shape, PVector centerPoint) {
		Geometry g = fromPShape(shape);
		Point p = PGS.pointFromPVector(centerPoint);
		Coordinate closestEdgePoint = DistanceOp.nearestPoints(g.getBoundary(), p)[0];

		double radius = PGS.distance(GEOM_FACTORY.createPoint(closestEdgePoint), p);
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(p.getCoordinate());
		shapeFactory.setWidth(radius * 2); // r*2 for total width & height
		shapeFactory.setHeight(radius * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * Computes the Minimum Bounding Circle (MBC) for the points in a Geometry. The
	 * MBC is the smallest circle which covers all the vertices of the input shape
	 * (this is also known as the Smallest Enclosing Circle). This is equivalent to
	 * computing the Maximum Diameter of the input vertex set.
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
	 * Computes the minimum bounding rectangle that encloses a shape. Unlike the
	 * envelope for a shape, the rectangle returned by this method can have any
	 * orientation.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumBoundingRectangle(PShape shape) {
		Polygon md = (Polygon) MinimumDiameter.getMinimumRectangle(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Computes the minimum diameter of a shape.
	 * <p>
	 * The minimum diameter is defined to be the width of the smallest band that
	 * contains the shape, where a band is a strip of the plane defined by two
	 * parallel lines. This can be thought of as the smallest hole that the geometry
	 * can bemoved through, with a single rotation.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumDiameter(PShape shape) {
		LineString md = (LineString) MinimumDiameter.getMinimumDiameter(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * Returns the nearest point of the shape to the given point. If the shape is
	 * has multiple children/geometries, the single closest point is returned.
	 * 
	 * @param shape
	 * @param point
	 * @return
	 * @see #closestPoints(PShape, PVector)
	 */
	public static PVector closestPoint(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		Coordinate coord = DistanceOp.nearestPoints(g, PGS.pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

	/**
	 * Returns the nearest point for each "island" in the input shape.
	 * 
	 * @param shape
	 * @param point
	 * @return list of closest points for each child shape. Output is identical to
	 *         {@link #closestPoint(PShape, PVector)} if the input shape
	 * @see #closestPoint(PShape, PVector)
	 */
	public static List<PVector> closestPoints(PShape shape, PVector point) {
		Geometry g = fromPShape(shape);
		ArrayList<PVector> points = new ArrayList<>();
		for (int i = 0; i < g.getNumGeometries(); i++) {
			Coordinate coord = DistanceOp.nearestPoints(g.getGeometryN(i), PGS.pointFromPVector(point))[0];
			points.add(new PVector((float) coord.x, (float) coord.y));
		}
		return points;
	}

}
