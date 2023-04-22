package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;

/**
 * Adapts {@link org.locationtech.jts.algorithm.construct.LargestEmptyCircle
 * LargestEmptyCircle}, allowing for repeated calls to find the N largest empty
 * circles in an optimised manner.
 * <p>
 * In this adaptation circle circumferences are constrained to lie within the
 * boundary (originally only circle center points must lie within the boundary).
 * This adaption also supports polygonal obstacles; if a boundary is given,
 * circles will not lie within the polygonal obstacles.
 * 
 * @author Martin Davis
 * @author Michael Carleton
 */
public class LargestEmptyCircles {

	private final Geometry obstacles;
	private final Geometry boundary;
	private final double tolerance;

	private final GeometryFactory factory;
	private IndexedPointInAreaLocator obstaclesPointLocator; // when obstacles are polygonal
	private IndexedPointInAreaLocator boundsPointLocator;
	private IndexedFacetDistance obstacleDistance;
	private IndexedFacetDistance boundaryDistance;
	private Envelope gridEnv;
	private Cell farthestCell;

	// Priority queue of cells, ordered by decreasing distance from constraints
	private LinkedList<Cell> cellQueue = new LinkedList<>();
	private List<Cell> nextIterCells = new ArrayList<>();
	private List<double[]> circles = new ArrayList<>();

	/**
	 * Constructs a new Largest Empty Circles (LEC) instance, ensuring that the
	 * circles are interior-disjoint to a set of obstacle geometries and (optional)
	 * contained within a polygonal boundary.
	 * <p>
	 * <li>If the provided boundary is null and obstacles are linear/pointal, the
	 * convex hull of the obstacles is used as the boundary.</li>
	 * <li>If the provided boundary is null and obstacles are polygonal, the
	 * obstacles themselves form the boundary (effectively function as a "Maximum
	 * Inscribed Circles" algorithm.).</li>
	 * 
	 * @param obstacles geometry representing the obstacles; if null, the boundary
	 *                  is used instead
	 * @param boundary  a polygonal geometry (may be null)
	 * @param tolerance a positive distance tolerance for computing the circle
	 *                  center point
	 * @throws IllegalArgumentException if the obstacles geometry or the tolerance
	 *                                  is non-positive
	 */
	public LargestEmptyCircles(Geometry obstacles, Geometry boundary, double tolerance) {
		if (obstacles == null || obstacles.isEmpty()) {
			throw new IllegalArgumentException("Obstacles geometry is null or empty.");
		}
		if (boundary != null && !(boundary instanceof Polygonal)) {
			throw new IllegalArgumentException("A non-null boundary must be polygonal.");
		}
		if (tolerance <= 0) {
			throw new IllegalArgumentException("Accuracy tolerance is non-positive: " + tolerance);
		}
		this.tolerance = tolerance;

		if (obstacles instanceof Polygonal && boundary != null) {
			// used to exclude circle from lying within the obstacles
			obstaclesPointLocator = new IndexedPointInAreaLocator(obstacles);
		}
		/*
		 * If no boundary given, use convex hull of obstacles as boundary.
		 */
		if (boundary == null || boundary.isEmpty()) {
			if (obstacles instanceof Polygonal) {
				boundary = obstacles;
			} else {
				boundary = obstacles.convexHull();
			}
		}

		this.obstacles = obstacles;
		this.boundary = boundary;
		this.factory = obstacles.getFactory();

		/*
		 * Circle *radii* will always be bounded by the boundary (in the JTS
		 * implementation, only circle center points must lie within the boundary).
		 */
		final Geometry distGeom = obstacles.getFactory().createGeometryCollection(new Geometry[] { obstacles, boundary });
		obstacleDistance = new IndexedFacetDistance(distGeom);

	}

	private void initBoundary() {
		gridEnv = boundary.getEnvelopeInternal();
		// if bounds does not enclose an area cannot create a boundsPointLocator
		if (boundary.getDimension() >= 2) {
			boundsPointLocator = new IndexedPointInAreaLocator(boundary);
			boundaryDistance = new IndexedFacetDistance(boundary);
		}

		createInitialGrid(gridEnv, cellQueue);
	}

	/**
	 * Computes the signed distance from a point to the constraints (obstacles and
	 * boundary). Points outside the boundary polygon are assigned a negative
	 * distance. Their containing cells will be last in the priority queue (but will
	 * still end up being tested since they may be refined).
	 * 
	 * @param p the point to compute the distance for
	 * @return the signed distance to the constraints (negative indicates outside
	 *         the boundary)
	 */
	private double distanceToConstraints(Point p) {
		Coordinate c = p.getCoordinate();
		boolean isOutsideBounds = boundsPointLocator.locate(c) == Location.EXTERIOR;
		if (isOutsideBounds) {
			double boundaryDist = boundaryDistance.distance(p);
			return -boundaryDist;
		}
		double dist = obstacleDistance.distance(p);

		/*
		 * If obstacles are polygonal, ensure circles do not lie within their interior.
		 * Only applies when the given boundary is not null.
		 */
		if (obstaclesPointLocator != null && (obstaclesPointLocator.locate(c) == Location.INTERIOR)) {
			dist = -dist;
		}
		return dist;
	}

	private double distanceToConstraints(double x, double y) {
		Coordinate coord = new Coordinate(x, y);
		Point pt = factory.createPoint(coord);
		return distanceToConstraints(pt);
	}

	/**
	 * Computes the (next) N largest empty circles.
	 * 
	 * @param n number of circles
	 * @return array of circles; each circle is represented as: [x,y,r]
	 */
	public double[][] findLECs(int n) {
		double[][] lecs = new double[n][3];
		for (int i = 0; i < n; i++) {
			lecs[i] = findNextLEC();
		}
		return lecs;
	}

	/**
	 * Computes the next largest empty circle.
	 * 
	 * @return an array representing the circle: [x,y,r]
	 */
	public double[] findNextLEC() {

		double farthestD;
		if (gridEnv == null) { // first iteration
			initBoundary();
			// use the area centroid as the initial candidate center point
			farthestCell = createCentroidCell(obstacles);
			farthestD = farthestCell.getDistance();
		} else {
			nextIterCells.forEach(c -> c.updateDistance(circles.get(circles.size() - 1)));
			cellQueue = new LinkedList<>(nextIterCells);
			nextIterCells.clear();
			farthestD = Double.MIN_VALUE;
		}

		/*
		 * Carry out the branch-and-bound search of the cell space.
		 */
		while (!cellQueue.isEmpty()) {
			// pick the cell with greatest distance from the queue
			Cell cell = cellQueue.removeLast();

			// update the center cell if the candidate is further from the constraints
			if (cell.getDistance() > farthestD) {
				farthestCell = cell;
				farthestD = farthestCell.getDistance();
			}

			/*
			 * If this cell may contain a better approximation to the center of the empty
			 * circle, then refine it (partition into subcells which are added into the
			 * queue for further processing). Otherwise the cell is pruned (not investigated
			 * further), since no point in it can be further than the current farthest
			 * distance.
			 */
			if (!cell.isFullyOutside()) {
				/*
				 * The cell is outside, but overlaps the boundary so it may contain a point
				 * which should be checked. This is only the case if the potential overlap
				 * distance is larger than the tolerance.
				 */
				if (cell.isOutside()) {
					boolean isOverlapSignificant = cell.getMaxDistance() > tolerance;
					if (isOverlapSignificant) {
						enqueueChildren(cell);
					}
				} else {
					/*
					 * Cell is inside the boundary. It may contain the center if the maximum
					 * possible distance is greater than the current distance (up to tolerance).
					 */
					double potentialIncrease = cell.getMaxDistance() - farthestCell.getDistance();
					if (potentialIncrease > tolerance) {
						enqueueChildren(cell);
					} else {
						nextIterCells.add(cell);
					}
				}
			}
		}
		// the farthest cell is the best approximation to the LEC center
		final Cell lecCell = farthestCell;
		final double r = lecCell.distance;
		final double[] circle = new double[] { lecCell.getX(), lecCell.getY(), r };
		circles.add(circle);

		return circle;
	}

	// split the cell into four sub-cells
	private void enqueueChildren(final Cell cell) {
		final double h2 = cell.getHSide() / 2;

		cellQueue.add(createCell(cell.x - h2, cell.y - h2, h2));
		cellQueue.add(createCell(cell.x + h2, cell.y - h2, h2));
		cellQueue.add(createCell(cell.x - h2, cell.y + h2, h2));
		cellQueue.add(createCell(cell.x + h2, cell.y + h2, h2));
	}

	/**
	 * Initializes the queue with a grid of cells covering the extent of the area.
	 * 
	 * @param env       the area extent to cover
	 * @param cellQueue the queue to initialize
	 */
	private void createInitialGrid(Envelope env, Collection<Cell> cellQueue) {
		double minX = env.getMinX();
		double maxX = env.getMaxX();
		double minY = env.getMinY();
		double maxY = env.getMaxY();
		double width = env.getWidth();
		double height = env.getHeight();
		double cellSize = Math.min(width, height);
		double hSize = cellSize / 2.0;

		// compute initial grid of cells to cover area
		for (double x = minX; x < maxX; x += cellSize) {
			for (double y = minY; y < maxY; y += cellSize) {
				cellQueue.add(createCell(x + hSize, y + hSize, hSize));
			}
		}
	}

	private Cell createCell(final double x, final double y, final double h) {
		Cell c = new Cell(x, y, h, distanceToConstraints(x, y));
		c.updateDistance(circles);
		return c;
	}

	// create a cell centered on area centroid
	private Cell createCentroidCell(Geometry geom) {
		Point p = geom.getCentroid();
		return new Cell(p.getX(), p.getY(), 0, distanceToConstraints(p));
	}

	/**
	 * A square grid cell centered on a given point with a given side half-length,
	 * and having a given distance from the center point to the constraints. The
	 * maximum possible distance from any point in the cell to the constraints can
	 * be computed. This is used as the ordering and upper-bound function in the
	 * branch-and-bound algorithm.
	 */
	private static class Cell implements Comparable<Cell> {

		private static final double SQRT2 = 1.4142135623730951;

		private double x;
		private double y;
		private double hSide;
		private double distance;
		private double maxDist;

		Cell(double x, double y, double hSide, double distanceToConstraints) {
			this.x = x; // cell center x
			this.y = y; // cell center y
			this.hSide = hSide; // half the cell size

			// the distance from cell center to constraints
			distance = distanceToConstraints;

			/*
			 * The maximum possible distance to the constraints for points in this cell is
			 * the center distance plus the radius (half the diagonal length).
			 */
			this.maxDist = distance + hSide * SQRT2;
		}

		public void updateDistance(double[] c) {
			double deltaX = x - c[0];
			double deltaY = y - c[1];
			// for now d is circle-cell-center dist
			double d = Math.sqrt(deltaX * deltaX + deltaY * deltaY);
			double r = c[2];

			if (d < r) {
				d = -(r - d); // negative (inside the circle)
			} else {
				// d is distance from cell center to circle boundary
				d -= r;
			}

			if (d < distance) {
				distance = d;
				maxDist = distance + hSide * SQRT2;
			}
		}

		/**
		 * Updates the distance of this cell based on circle constraints. Distance
		 * between circle boundary is used.
		 */
		public void updateDistance(List<double[]> circles) {
			double minCircleDist = Double.MAX_VALUE;
			for (double[] c : circles) {
				double deltaX = x - c[0];
				double deltaY = y - c[1];
				// for now d is circle-cell-center dist
				double d = Math.sqrt(deltaX * deltaX + deltaY * deltaY);
				double r = c[2];

				if (d < r) {
					d = -(r - d); // negative (inside the circle)
				} else {
					// d is distance from cell center to circle boundary
					d -= r;
				}
				minCircleDist = Math.min(minCircleDist, d);
			}

			if (minCircleDist < distance) {
				distance = minCircleDist;
				maxDist = distance + hSide * SQRT2;
			}
		}

		public boolean isFullyOutside() {
			return getMaxDistance() < 0;
		}

		public boolean isOutside() {
			return distance < 0;
		}

		public double getMaxDistance() {
			return maxDist;
		}

		public double getDistance() {
			return distance;
		}

		public double getHSide() {
			return hSide;
		}

		public double getX() {
			return x;
		}

		public double getY() {
			return y;
		}

		/**
		 * A cell is greater if its maximum distance is larger.
		 */
		@Override
		public int compareTo(Cell o) {
			return (int) (o.maxDist - this.maxDist);
		}
	}

}
