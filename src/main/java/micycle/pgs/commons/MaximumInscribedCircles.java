package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;
import org.tinspin.index.PointDistanceFunction;
import org.tinspin.index.PointEntryDist;
import org.tinspin.index.PointIndex;
import org.tinspin.index.covertree.CoverTree;

import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;
import processing.core.PVector;

/**
 * An bespoke version of
 * {@link org.locationtech.jts.algorithm.construct.MaximumInscribedCircle
 * MaximumInscribedCircle} to find N largest maximum inscribed circles in an
 * optimised manner.
 * 
 * @author Michael Carleton
 *
 */
public class MaximumInscribedCircles {

	private static final boolean USE_TREE = false;

	private final Geometry inputGeom;
	private final double tolerance;

	private GeometryFactory factory;
	private IndexedPointInAreaLocator ptLocater;
	private IndexedFacetDistance indexedDistance;
	private final Point centroid;
	private final PointIndex<Coordinate> tree = CoverTree.create(3, 2, circleDistanceMetric);
	private final ArrayList<PVector> circles = new ArrayList<>();

	private double maxCircleRadius = 0;

	private final Map<Coordinate, Double> distanceCache;

	/**
	 * Creates a new instance of a Maximum Inscribed Circles computation.
	 * 
	 * @param polygonal an areal geometry
	 * @param tolerance the distance tolerance for computing centre points (must be
	 *                  positive)
	 * @throws IllegalArgumentException if the tolerance is non-positive, or the
	 *                                  input geometry is non-polygonal or empty.
	 */
	public MaximumInscribedCircles(Geometry polygonal, double tolerance) {
		if (tolerance <= 0) {
			throw new IllegalArgumentException("Tolerance must be positive");
		}
		if (!(polygonal instanceof Polygon || polygonal instanceof MultiPolygon)) {
			throw new IllegalArgumentException("Input geometry must be a Polygon or MultiPolygon");
		}
		if (polygonal.isEmpty()) {
			throw new IllegalArgumentException("Empty input geometry is not supported");
		}

		this.inputGeom = polygonal;
		this.centroid = polygonal.getCentroid();
		this.factory = polygonal.getFactory();
		this.tolerance = tolerance;
		ptLocater = new IndexedPointInAreaLocator(polygonal);
		indexedDistance = new IndexedFacetDistance(polygonal.getBoundary());
		distanceCache = new HashMap<>();
	}

	/**
	 * Get the next largest MIC.
	 * 
	 * @return
	 */
	public double[] getNextLargestCircle() {
		return compute();
	}

	private double[] compute() {
		// Priority queue of cells, ordered by maximum distance from boundary
		ObjectHeapPriorityQueue<Cell> cellQueue = new ObjectHeapPriorityQueue<>();
		/*
		 * TODO optimise this -- rather than recreating grid each time, prune, update
		 * then use grid cells from previous iteration. Prunes if they lie inside the
		 * last circle placed.
		 */
		createInitialGrid(inputGeom.getEnvelopeInternal(), cellQueue);
		// use the area centroid as the initial candidate center point
		Cell farthestCell = createCentroidCell();

		/*
		 * Carry out the branch-and-bound search of the cell space
		 */
		while (!cellQueue.isEmpty()) {
			// pick the most promising cell from the queue
			final Cell cell = cellQueue.dequeue();

			// update the center cell if the candidate is further from the boundary
			if (cell.getDistance() > farthestCell.getDistance()) {
				farthestCell = cell;
			}
			/*
			 * Refine this cell if the potential distance improvement is greater than the
			 * required tolerance. Otherwise the cell is pruned (not investigated further),
			 * since no point in it is further than the current farthest distance.
			 */
			double potentialIncrease = cell.getMaxDistance() - farthestCell.getDistance();
			if (potentialIncrease > tolerance) {
				// split the cell into four sub-cells
				double h2 = cell.getHSide() / 2;
				cellQueue.enqueue(createCell(cell.getX() - h2, cell.getY() - h2, h2));
				cellQueue.enqueue(createCell(cell.getX() + h2, cell.getY() - h2, h2));
				cellQueue.enqueue(createCell(cell.getX() - h2, cell.getY() + h2, h2));
				cellQueue.enqueue(createCell(cell.getX() + h2, cell.getY() + h2, h2));
			}
		}
		// the farthest cell is the best approximation to the MIC center
		Cell centerCell = farthestCell;
		maxCircleRadius = Math.max(maxCircleRadius, farthestCell.getDistance());
		if (USE_TREE) {
			tree.insert(new double[] { centerCell.getX(), centerCell.getY(), farthestCell.getDistance() },
					new Coordinate(centerCell.getX(), centerCell.getY()));
		} else {
			circles.add(new PVector((float) centerCell.getX(), (float) centerCell.getY(), (float) farthestCell.getDistance()));
		}

		return new double[] { centerCell.getX(), centerCell.getY(), farthestCell.getDistance() };
	}

	/**
	 * Initializes the queue with a grid of cells covering the extent of the area.
	 * 
	 * @param env       the area extent to cover
	 * @param cellQueue the queue to initialize
	 */
	private void createInitialGrid(Envelope env, ObjectHeapPriorityQueue<Cell> cellQueue) {
		double minX = env.getMinX();
		double maxX = env.getMaxX();
		double minY = env.getMinY();
		double maxY = env.getMaxY();
		double width = env.getWidth();
		double height = env.getHeight();
		double cellSize = Math.min(width, height);

		// Check for flat collapsed input and if so short-circuit
		// Result will just be centroid
		if (cellSize == 0) {
			return;
		}

		final double hSide = cellSize / 2.0;
		// compute initial grid of cells to cover area
		for (double x = minX; x < maxX; x += cellSize) {
			for (double y = minY; y < maxY; y += cellSize) {
				cellQueue.enqueue(createCell(x + hSide, y + hSide, hSide));
			}
		}
	}

	private double distanceToBoundary(double x, double y) {
		Coordinate coord = new Coordinate(x, y);
		return distanceToBoundary(coord);
	}

	/**
	 * Computes the signed distance from a point to the area boundary or the
	 * perimeter of the nearest existing MIC (whichever is smallest).
	 * <p>
	 * Points outside the polygon are assigned a negative distance. Their containing
	 * cells will be last in the priority queue (but may still end up being tested
	 * since they may need to be refined).
	 * 
	 * @param coord the point to compute the distance for
	 * @return the signed distance to the area boundary (negative indicates outside
	 *         the area)
	 */
	private double distanceToBoundary(Coordinate coord) {
		double boundaryDist = distanceCache.computeIfAbsent(coord, c -> {
			final double d = indexedDistance.distance(factory.createPoint(coord));
			final boolean isOutide = Location.EXTERIOR == ptLocater.locate(coord);
			return isOutide ? -d : d;
		});
		if (boundaryDist < 0) {
			return boundaryDist;
		}

		/*
		 * Brute-force (above) appears to be faster than cover tree query (probably
		 * because circle count is fairly low.
		 */
		if (USE_TREE) {
			final PointEntryDist<Coordinate> query = tree.query1NN(new double[] { coord.getX(), coord.getY(), maxCircleRadius });
			if (tree.getNodeCount() == 0) {
				return boundaryDist;
			}

			final double circleRadius = query.point()[2];
			final double distance = query.value().distance(coord);
			if (distance < circleRadius) {
				return -(circleRadius - distance);
			}

			return Math.min(boundaryDist, query.dist() - maxCircleRadius);
		} else {
			final PVector q = new PVector((float) coord.getX(), (float) coord.getY()); // query
			double dist = boundaryDist;
			for (PVector c : circles) {
				q.z = c.z; // so dist isnt affected
				if (q.dist(c) < c.z) { // inside
					return -(c.z - q.dist(c)); // dist to circle edge
				}
				dist = Math.min(dist, q.dist(c) - c.z);
			}
			return dist;
		}
	}

	private Cell createCell(double x, double y, double hSide) {
		return new Cell(x, y, hSide, distanceToBoundary(x, y));
	}

	// create a cell centered on area centroid
	private Cell createCentroidCell() {
		return new Cell(centroid.getX(), centroid.getY(), 0, distanceToBoundary(centroid.getCoordinate()));
	}

	/**
	 * A square grid cell centered on a given point, with a given half-side size,
	 * and having a given distance to the area boundary. The maximum possible
	 * distance from any point in the cell to the boundary can be computed, and is
	 * used as the ordering and upper-bound function in the branch-and-bound
	 * algorithm.
	 *
	 */
	private static class Cell implements Comparable<Cell> {

		private static final double SQRT2 = 1.4142135623730951;

		private double x;
		private double y;
		private double hSide;
		private double distance;
		private double maxDist;

		Cell(double x, double y, double hSide, double distanceToBoundary) {
			this.x = x; // cell center x
			this.y = y; // cell center y
			this.hSide = hSide; // half the cell size

			// the distance from cell center to area boundary
			distance = distanceToBoundary;

			// the maximum possible distance to area boundary for points in this cell
			this.maxDist = distance + hSide * SQRT2;
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
		 * A cell is greater if its maximum possible distance is larger.
		 */
		@Override
		public int compareTo(Cell o) {
			return Double.compare(o.maxDist, this.maxDist);
		}

	}

	private static final PointDistanceFunction circleDistanceMetric = (p1, p2) -> {
		final double dx = p1[0] - p2[0];
		final double dy = p1[1] - p2[1];
		return Math.sqrt(dx * dx + dy * dy) + Math.abs(p1[2] - p2[2]);
	};

}
