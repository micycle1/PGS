package micycle.pgs.commons;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;

import com.github.micycle1.geoblitz.YStripesPointInAreaLocator;

/**
 * Adapts {@link org.locationtech.jts.algorithm.construct.LargestEmptyCircle
 * LargestEmptyCircle}, allowing for repeated calls to find the N largest empty
 * circles in an optimised manner.
 * <p>
 * In this adaptation circle circumferences are constrained to lie within the
 * boundary (originally only circle center points must lie within the boundary).
 * This adaption also supports polygonal obstacles; if a boundary is provided,
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
	private PointOnGeometryLocator obstaclesPointLocator; // when obstacles are polygonal
	private PointOnGeometryLocator boundsPointLocator;
	private IndexedFacetDistance obstacleDistance;
	private IndexedFacetDistance boundaryDistance;
	private Envelope gridEnv;
	private Cell farthestCell;

	private ArrayDeque<Cell> cellStack = new ArrayDeque<>();
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
			obstaclesPointLocator = new YStripesPointInAreaLocator(obstacles);
		}
		if (boundary == null || boundary.isEmpty()) {
			if (obstacles instanceof Polygonal) {
				boundary = obstacles;
			} else {
				/*
				 * If no boundary given, use convex hull of obstacles as boundary.
				 */
				boundary = obstacles.convexHull();
			}
		}

		this.obstacles = obstacles;
		this.boundary = boundary;
		this.factory = obstacles.getFactory();

		/*
		 * Combine, in case the nearest obstacle is farther away than the nearest
		 * boundary. (it's faster to make one call on a larger index than 2 separate
		 * calls to each index).
		 */
		final Geometry distGeom = obstacles.getFactory().createGeometryCollection(new Geometry[] { obstacles, boundary });
		obstacleDistance = new IndexedFacetDistance(distGeom);
//	    obstacleDistance = new IndexedFacetDistance(obstacles);
	}

	private void initBoundary() {
		gridEnv = boundary.getEnvelopeInternal();
		if (boundary.getDimension() >= 2) {
			boundsPointLocator = new YStripesPointInAreaLocator(boundary);
			boundaryDistance = new IndexedFacetDistance(boundary);
		}
		createInitialGrid(gridEnv, cellStack);
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
		/*
		 * If obstacles are polygonal, ensure circles do not lie within their interior.
		 * Only applies when the given boundary is not null.
		 */
		double dist = obstacleDistance.distance(p);
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

	public double[][] findLECs(int n) {
		double[][] lecs = new double[n][3];
		for (int i = 0; i < n; i++) {
			lecs[i] = findNextLEC();
		}
		return lecs;
	}

	public double[] findNextLEC() {

		double farthestD;
		if (gridEnv == null) { // first iteration
			initBoundary();
			// pick best seed from initial grid instead of only centroid
			farthestCell = createCentroidCell(obstacles);
			farthestD = farthestCell.getDistance();
			for (Cell c : cellStack) {
				if (c.getDistance() > farthestD) {
					farthestCell = c;
					farthestD = c.getDistance();
				}
			}
		} else {
			// Update remaining candidates with the newly-placed circle
			nextIterCells.forEach(c -> c.updateDistance(circles.get(circles.size() - 1)));
			cellStack = new ArrayDeque<>(nextIterCells);
			nextIterCells.clear();

			// Seed best for this iteration from what we already have
			farthestD = Double.NEGATIVE_INFINITY;
			for (Cell c : cellStack) {
				if (c.getDistance() > farthestD) {
					farthestCell = c;
					farthestD = c.getDistance();
				}
			}
		}

		// Branch-and-bound
		while (!cellStack.isEmpty()) {
			// LIFO pop for DFS-like behavior
			Cell cell = cellStack.removeLast();

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
					double potentialIncrease = cell.getMaxDistance() - farthestD;
					if (potentialIncrease > tolerance) {
						enqueueChildren(cell);
					} else {
						nextIterCells.add(cell);
					}
				}
			}
		}

		final Cell lecCell = farthestCell;
		final double r = lecCell.distance;
		final double[] circle = new double[] { lecCell.getX(), lecCell.getY(), r };
		circles.add(circle);

		return circle;
	}

	private void enqueueChildren(final Cell cell) {
		final double h2 = cell.getHSide() / 2;
		final double parentDist = cell.getDistance();
		final double farthestD = (farthestCell != null) ? farthestCell.getDistance() : Double.NEGATIVE_INFINITY;

		// The max potential increase from parent's center to any point in a child cell
		// is dist(parent_center, child_corner) = sqrt((h2)^2 + (h2)^2) = h2*sqrt(2)
		// The max distance in a child cell is at its corner, which is h2*sqrt(2) from
		// its center.
		// So, an upper bound on a child's maxDist is parentDist + 2 * h2 * SQRT2
		double maxChildPotential = parentDist + 2 * h2 * Math.sqrt(2);

		if (maxChildPotential <= farthestD + tolerance) {
			// Even the most optimistic estimate for any child of this cell
			// won't beat the current best. So we don't need to subdivide.
			// We might still need to keep this cell for the next iteration.
			nextIterCells.add(cell);
			return;
		}

		Cell c1 = createCellIfPromising(cell.x - h2, cell.y - h2, h2, farthestD);
		Cell c2 = createCellIfPromising(cell.x + h2, cell.y - h2, h2, farthestD);
		Cell c3 = createCellIfPromising(cell.x - h2, cell.y + h2, h2, farthestD);
		Cell c4 = createCellIfPromising(cell.x + h2, cell.y + h2, h2, farthestD);

		Cell[] kids = new Cell[] { c1, c2, c3, c4 };
		Arrays.sort(kids, (a, b) -> {
			if (a == null && b == null)
				return 0;
			if (a == null)
				return -1; // nulls go first
			if (b == null)
				return 1;
			return Double.compare(a.getMaxDistance(), b.getMaxDistance());
		});

		for (Cell k : kids) {
			if (k != null) {
				cellStack.addLast(k);
			}
		}
	}

	// Helper method to create a cell only if it's worth investigating
	private Cell createCellIfPromising(final double x, final double y, final double h, double farthestD) {
		// We can't use the Lipschitz bound here because we don't know the parent's
		// distance,
		// but we can check the cell after creation before adding it.
		// The main pruning is the check in enqueueChildren. This is a secondary check.
		Cell c = createCell(x, y, h);
		if (c.getMaxDistance() > farthestD + tolerance) {
			return c;
		}
		// If not promising, but might be useful for the next LEC search, add it there.
		if (!c.isFullyOutside()) {
			nextIterCells.add(c);
		}
		return null;
	}

	private void createInitialGrid(Envelope env, Collection<Cell> target) {
		double minX = env.getMinX();
		double maxX = env.getMaxX();
		double minY = env.getMinY();
		double maxY = env.getMaxY();
		double width = env.getWidth();
		double height = env.getHeight();
		double cellSize = Math.min(width, height);
		double hSize = cellSize / 2.0;

		for (double x = minX; x < maxX; x += cellSize) {
			for (double y = minY; y < maxY; y += cellSize) {
				target.add(createCell(x + hSize, y + hSize, hSize));
			}
		}
	}

	private Cell createCell(final double x, final double y, final double h) {
		Cell c = new Cell(x, y, h, distanceToConstraints(x, y));
		c.updateDistance(circles);
		return c;
	}

	private Cell createCentroidCell(Geometry geom) {
		Point p = geom.getCentroid();
		return new Cell(p.getX(), p.getY(), 0, distanceToConstraints(p));
	}

	private static class Cell implements Comparable<Cell> {

		private static final double SQRT2 = 1.4142135623730951;

		private double x;
		private double y;
		private double hSide;
		private double distance;
		private double maxDist;

		Cell(double x, double y, double hSide, double distanceToConstraints) {
			this.x = x;
			this.y = y;
			this.hSide = hSide;
			distance = distanceToConstraints;
			this.maxDist = distance + hSide * SQRT2;
		}

		// CHANGED: sqrt-free early-rejection for circle; only take sqrt if it might
		// improve
		public void updateDistance(double[] c) {
			final double dx = x - c[0];
			final double dy = y - c[1];
			final double dsq = dx * dx + dy * dy;
			final double r = c[2];

			double D = distance;
			double t = D + r; // improvement only possible if sqrt(dsq) < D + r
			if (t > 0) {
				double tsq = t * t;
				if (dsq < tsq) {
					double d = Math.sqrt(dsq) - r; // signed (negative when inside)
					if (d < D) {
						distance = d;
						maxDist = distance + hSide * SQRT2;
					}
				}
			}
		}

		// CHANGED: sqrt-free early-rejection loop over all circles
		public void updateDistance(List<double[]> circles) {
			double D = distance;
			for (double[] c : circles) {
				final double r = c[2];
				double t = D + r;
				if (t <= 0) {
					// No circle can improve when D <= -r for this circle
					continue;
				}
				final double dx = x - c[0];
				final double dy = y - c[1];
				final double dsq = dx * dx + dy * dy;
				final double tsq = t * t;

				if (dsq < tsq) {
					double d = Math.sqrt(dsq) - r; // signed (negative when inside)
					if (d < D) {
						D = d;
					}
				}
			}
			if (D < distance) {
				distance = D;
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

		@Override
		public int compareTo(Cell o) {
			return (int) (o.maxDist - this.maxDist);
		}
	}

}
