package micycle.pts;

import static micycle.pts.Conversion.fromPShape;
import static micycle.pts.Conversion.toPShape;
import static micycle.pts.PTS.CURVE_SAMPLES;
import static micycle.pts.PTS.GEOM_FACTORY;

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
 * Bounding Volumes, inscribed areas, Optimal Distances, etc.
 * 
 * @author MCarleton
 *
 */
public class PTSGeometricOptimisation {

	/**
	 * The Maximum Inscribed Circle is determined by a point in the interior of the
	 * area which has the farthest distance from the area boundary, along with a
	 * boundary point at that distance.
	 * 
	 * @param shape
	 * @param tolerance the distance tolerance for computing the centre point
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
	 * @return
	 */
	public static PShape maximumInscribedCircle(PShape shape, PVector centerPoint) {
		Geometry g = fromPShape(shape);
		Point p = PTS.pointFromPVector(centerPoint);
		Coordinate closestEdgePoint = DistanceOp.nearestPoints(g.getBoundary(), p)[0];

		double radius = PTS.distance(GEOM_FACTORY.createPoint(closestEdgePoint), p);
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory();
		shapeFactory.setNumPoints(CURVE_SAMPLES * 4); // TODO magic constant
		shapeFactory.setCentre(p.getCoordinate());
		shapeFactory.setWidth(radius * 2); // r*2 for total width & height
		shapeFactory.setHeight(radius * 2); // r*2 for total width & height
		return toPShape(shapeFactory.createEllipse());
	}

	/**
	 * Computes the Minimum Bounding Circle (MBC)for the points in a Geometry. The
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
	 * Gets the minimum rectangle enclosing a shape.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumBoundingRectangle(PShape shape) {
		Polygon md = (Polygon) MinimumDiameter.getMinimumRectangle(fromPShape(shape));
		return toPShape(md);
	}

	/**
	 * 
	 * 
	 * Computes the minimum diameter of a shape. The minimum diameter is defined to
	 * be the width of the smallest band that contains the shape, where a band is a
	 * strip of the plane defined by two parallel lines. This can be thought of as
	 * the smallest hole that the geometry can bemoved through, with a single
	 * rotation.
	 * 
	 * @param shape
	 * @return
	 */
	public static PShape minimumDiameter(PShape shape) {
		LineString md = (LineString) MinimumDiameter.getMinimumDiameter(fromPShape(shape));
		return toPShape(md);
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
		Coordinate coord = DistanceOp.nearestPoints(g, PTS.pointFromPVector(point))[0];
		return new PVector((float) coord.x, (float) coord.y);
	}

}
