package micycle.pgs.commons;

import java.util.Arrays;

/**
 * A spatial index for circles that supports fast incremental insertion and
 * nearest-circle queries.
 * <p>
 * Each entry is stored as a circle {@code (x, y, r)} plus an associated value.
 * The structure is designed for workloads where circles are added over time and
 * queried frequently, such as stochastic circle packing or geometric sampling.
 * <p>
 * Internally this uses a cover tree with the metric
 * {@code d(c1, c2) = hypot(dx, dy) + |dr|}. Point queries are evaluated using
 * an embedding based on the current {@linkplain #maxRadius() maximum inserted
 * radius}, which makes nearest-neighbor search correspond to minimizing circle
 * clearance {@code hypot(q - center) - r}.
 * <p>
 * This implementation is a fusion of the circle-distance metric described in
 * <a href="https://stackoverflow.com/a/21975136/">this Stack Overflow
 * answer</a> and the TinSpin utilities {@code CoverTree} class. The metric is
 * hard-wired into a fast, specialized reimplementation of TinSpin's cover tree,
 * rather than being supplied as a generic distance function over
 * {@code (x, y, r)} tuples.
 *
 * @author Michael Carleton
 */
public final class CircleCoverTree<T> {

	// tree state
	private Node<T> root;
	private Node<T> anyLeaf; // for quick detach during rare root-growth
	private int size;
	private double maxR;

	// scratch (reused, no allocations during queries)
	private Node<T>[] stack = null;
	private int[] order = null;
	private double[] keys = null;

	/**
	 * Returns the number of circles currently stored in the tree.
	 *
	 * @return the number of indexed circles
	 */
	public int size() {
		return size;
	}

	/**
	 * Returns the largest radius inserted so far.
	 *
	 * @return the maximum stored circle radius
	 */
	public double maxRadius() {
		return maxR;
	}

	// circle-circle metric for the cover tree structure
	private static double dCircleCircle(Node<?> a, Node<?> b) {
		double dx = a.x - b.x;
		double dy = a.y - b.y;
		double dr = a.r - b.r;
		return Math.sqrt(dx * dx + dy * dy) + Math.abs(dr);
	}

	// Query embedding distance with R>=r: sqrt(d2) + (R - r)
	// We try to avoid sqrt most of the time by using squared pruning.
	private static double queryOffset(double R, Node<?> n) {
		return (R - n.r); // >=0
	}

	/**
	 * Inserts a circle and its associated value into the index.
	 *
	 * @param x     circle center x-coordinate
	 * @param y     circle center y-coordinate
	 * @param r     circle radius
	 * @param value value associated with the circle
	 */
	public void insert(double x, double y, double r, T value) {
		maxR = Math.max(maxR, r);

		if (root == null) {
			root = new Node<>(x, y, r, value, 0);
			anyLeaf = root;
			size = 1;
			return;
		}

		Node<T> xn = new Node<>(x, y, r, value, -1);
		anyLeaf = xn; // new node is a leaf at insertion time

		if (!root.hasChildren()) {
			double dist = dCircleCircle(root, xn);
			int level = levelFromDist(dist);
			root.level = level + 1;
			xn.level = level;
			root.addChild(xn, dist);
			root.adjustMaxDist(dist);
			size++;
			return;
		}

		insert(root, xn);
		size++;
	}

	private void insert(Node<T> p, Node<T> x) {
		double distPX = dCircleCircle(p, x);
		double covP = covDist(p.level);

		if (distPX > covP) {
			// Grow tree upward until x fits.
			while (distPX > covP) { // BASE=2 => (BASE-1)*covP == covP
				Node<T> q = detachAnyLeaf(); // O(1) typical
				q.level = p.level + 1;
				double dqp = dCircleCircle(q, p);
				q.addChild(p, dqp);
				q.adjustMaxDist(dqp + p.maxDist);
				p = q;
				distPX = dCircleCircle(p, x);
				covP = covDist(p.level);
			}

			x.level = p.level + 1;
			Node<T> newRoot = x;
			newRoot.addChild(p, distPX);
			newRoot.adjustMaxDist(distPX + p.maxDist);
			root = newRoot;
			return;
		}

		insertNearestAncestor(p, x, distPX);
	}

	private void insertNearestAncestor(Node<T> p, Node<T> x, double distPX) {
		double covChild = covDist(p.level - 1);

		Node<T> best = null;
		double bestDist = Double.POSITIVE_INFINITY;

		// pick closest child that can cover x
		for (int i = 0; i < p.nChildren; i++) {
			Node<T> q = p.children[i];
			double dqx = dCircleCircle(q, x);
			if (dqx <= covChild && dqx < bestDist) {
				bestDist = dqx;
				best = q;
			}
		}

		if (best != null) {
			insertNearestAncestor(best, x, bestDist);
			p.adjustMaxDist(distPX);
			return;
		}

		x.level = p.level - 1;
		p.addChild(x, distPX);
		p.adjustMaxDist(distPX);
	}

	private Node<T> detachAnyLeaf() {
		Node<T> leaf = anyLeaf;

		// If leaf is root (tree size 1), root-growth case won't happen anyway.
		if (leaf == null || leaf.parent == null) {
			// fallback (rare)
			leaf = root;
			while (leaf.hasChildren()) {
				leaf = leaf.children[0];
			}
		}

		Node<T> parent = leaf.parent;
		parent.removeChildAt(leaf.parentIndex);
		leaf.parent = null;
		leaf.parentIndex = -1;

		// pick a new "anyLeaf" cheaply
		Node<T> cur = parent;
		while (cur != null && cur.hasChildren()) {
			cur = cur.children[0];
		}
		anyLeaf = (cur != null) ? cur : root;

		// NOTE: We are not tightening maxDist after removal (conservative). Correct,
		// slightly less pruning.
		return leaf;
	}

	/**
	 * Tests whether the query point lies inside any indexed circle.
	 *
	 * @param qx query point x-coordinate
	 * @param qy query point y-coordinate
	 * @return {@code true} if the point is contained in at least one stored circle;
	 *         {@code false} otherwise
	 */
	public boolean existsInside(double qx, double qy) {
		if (root == null) {
			return false;
		}
		final double R = maxR;
		return existsWithinMetric(qx, qy, R, R);
	}

	/**
	 * Returns whether any stored circle is within a given metric distance of the
	 * query.
	 * <p>
	 * This is the lower-level query primitive used by the more user-friendly
	 * helpers. For example, with {@code R = maxRadius}, {@code limitMetric = R}
	 * tests whether the point lies inside any circle, and
	 * {@code limitMetric = R + t} tests whether some circle has clearance at most
	 * {@code t}.
	 *
	 * @param qx          query point x-coordinate
	 * @param qy          query point y-coordinate
	 * @param R           query embedding radius, typically {@link #maxRadius()}
	 * @param limitMetric inclusive metric-distance threshold
	 * @return {@code true} if any indexed circle satisfies the threshold;
	 *         {@code false} otherwise
	 */
	public boolean existsWithinMetric(double qx, double qy, double R, double limitMetric) {
		if (root == null) {
			return false;
		}

		ensureStack();

		int sp = 0;
		stack[sp++] = root;

		while (sp != 0) {
			Node<T> p = stack[--sp];

			// Check p itself within limit without sqrt:
			// D = sqrt(d2) + (R - r) <= limit <=> sqrt(d2) <= limit - (R - r)
			double rhs = limitMetric - queryOffset(R, p);
			if (rhs >= 0) {
				double dx = qx - p.x, dy = qy - p.y;
				double d2 = dx * dx + dy * dy;
				if (d2 <= rhs * rhs) {
					return true;
				}
			}

			if (!p.hasChildren()) {
				continue;
			}

			// Push children that cannot be pruned (also sqrt-free):
			// Need: (sqrt(d2) + (R-rChild)) - maxDistChild <= limit
			// <=> sqrt(d2) <= limit + maxDistChild - (R - rChild)
			for (int i = 0; i < p.nChildren; i++) {
				Node<T> c = p.children[i];
				double bound = limitMetric + c.maxDist - queryOffset(R, c);
				if (bound <= 0) {
					continue; // impossible
				}
				double dx = qx - c.x, dy = qy - c.y;
				double d2 = dx * dx + dy * dy;
				if (d2 <= bound * bound) {
					stack[sp++] = c;
					if (sp == stack.length) {
						growStack();
					}
				}
			}
		}
		return false;
	}

	/**
	 * Finds the stored circle with minimum clearance to the query point.
	 * <p>
	 * Clearance is {@code hypot(q - center) - radius}, so negative values mean the
	 * point lies inside the returned circle.
	 *
	 * @param qx query point x-coordinate
	 * @param qy query point y-coordinate
	 * @return the nearest-circle result, or {@code null} if the tree is empty
	 */
	public Nearest<T> nearest(double qx, double qy) {
		if (root == null) {
			return null;
		}
		final double R = maxR;

		ensureStack();

		Nearest<T> best = new Nearest<>();
		// compute exact root distance once
		double dx0 = qx - root.x, dy0 = qy - root.y;
		double dRoot = Math.sqrt(dx0 * dx0 + dy0 * dy0) + queryOffset(R, root);
		best.set(root, dRoot, R);

		int sp = 0;
		stack[sp++] = root;

		while (sp != 0) {
			Node<T> p = stack[--sp];

			if (!p.hasChildren()) {
				continue;
			}

			// Optional near-first: gather candidates and visit smaller d2 first (no alloc).
			int nCand = 0;
			ensureOrderAndKeys(p.nChildren);

			for (int i = 0; i < p.nChildren; i++) {
				Node<T> c = p.children[i];

				// sqrt-free prune test derived from:
				// explore if (dist(query,c) - c.maxDist) < best.metricDist
				//
				// dist(query,c) = sqrt(d2) + (R - c.r)
				// => sqrt(d2) < best.metricDist + c.maxDist - (R - c.r)
				double t = best.metricDist + c.maxDist - queryOffset(R, c);
				if (t <= 0) {
					continue;
				}

				double dx = qx - c.x, dy = qy - c.y;
				double d2 = dx * dx + dy * dy;
				if (d2 >= t * t) {
					continue;
				}

				// Keep candidate; key by d2 (approx near-first)
				order[nCand] = i;
				keys[nCand] = d2;
				nCand++;
			}

			// insertion sort (fast for small degrees; cover trees typically have small
			// branching)
			for (int i = 1; i < nCand; i++) {
				int oi = order[i];
				double ki = keys[i];
				int j = i - 1;
				while (j >= 0 && keys[j] > ki) {
					keys[j + 1] = keys[j];
					order[j + 1] = order[j];
					j--;
				}
				keys[j + 1] = ki;
				order[j + 1] = oi;
			}

			// push in reverse so nearest visited first (stack LIFO)
			for (int idx = nCand - 1; idx >= 0; idx--) {
				Node<T> c = p.children[order[idx]];

				// Now compute exact dist (sqrt) only for survivors
				double dx = qx - c.x, dy = qy - c.y;
				double d2 = dx * dx + dy * dy;
				double dist = Math.sqrt(d2) + queryOffset(R, c);

				if (dist < best.metricDist) {
					best.set(c, dist, R);
				}

				// We already know it's worth exploring based on squared prune,
				// but best may have improved; optional re-check:
				if (dist - c.maxDist < best.metricDist) {
					stack[sp++] = c;
					if (sp == stack.length) {
						growStack();
					}
				}
			}
		}

		return best;
	}

	@SuppressWarnings("unchecked")
	private void ensureStack() {
		if (stack == null) {
			stack = new Node[64];
		}
	}

	private void growStack() {
		stack = Arrays.copyOf(stack, stack.length * 2);
	}

	private void ensureOrderAndKeys(int needed) {
		if (order == null || order.length < needed) {
			int cap = 1;
			while (cap < needed) {
				cap <<= 1;
			}
			order = new int[cap];
			keys = new double[cap];
		}
	}

	private static double covDist(int level) {
		return Math.scalb(1.0, level); // 2^level
	}

	private static int levelFromDist(double dist) {
		return Math.getExponent(dist); // floor(log2(dist)) for dist>0
	}

	private static final class Node<T> {
		double x, y, r;
		int level;

		// metric bookkeeping
		double distToParent; // circle-circle metric distance to parent
		double maxDist; // max metric distance from this node to any descendant

		// parent pointers (speeds leaf detaches)
		Node<T> parent;
		int parentIndex = -1;

		T value;

		Node<T>[] children;
		int nChildren;

		@SuppressWarnings("unchecked")
		Node(double x, double y, double r, T value, int level) {
			this.x = x;
			this.y = y;
			this.r = r;
			this.value = value;
			this.level = level;
			this.children = new Node[0];
		}

		boolean hasChildren() {
			return nChildren != 0;
		}

		void ensureChildCap(int cap) {
			if (children.length >= cap) {
				return;
			}
			int newCap = Math.max(4, children.length * 2);
			if (newCap < cap) {
				newCap = cap;
			}
			children = Arrays.copyOf(children, newCap);
		}

		void addChild(Node<T> c, double distThisToChild) {
			ensureChildCap(nChildren + 1);
			int idx = nChildren++;
			children[idx] = c;

			c.parent = this;
			c.parentIndex = idx;
			c.distToParent = distThisToChild;
		}

		void removeChildAt(int idx) {
			int last = --nChildren;
			if (idx != last) {
				Node<T> moved = children[last];
				children[idx] = moved;
				moved.parentIndex = idx;
			}
			children[last] = null;
		}

		void adjustMaxDist(double maybe) {
			if (maybe > maxDist) {
				maxDist = maybe;
			}
		}
	}

	public static final class Nearest<T> {
		/** Value associated with the nearest circle. */
		public T value;

		/** X-coordinate of the nearest circle center. */
		public double cx;

		/** Y-coordinate of the nearest circle center. */
		public double cy;

		/** Radius of the nearest circle. */
		public double cr;

		/** Metric distance {@code D((qx, qy, R=maxR), (cx, cy, r))}. */
		public double metricDist;

		/**
		 * Signed Euclidean distance from the query point to the nearest circle
		 * boundary: {@code hypot(q - c) - r = metricDist - R}.
		 * <p>
		 * Positive values mean the query lies outside the circle, zero means it lies on
		 * the boundary, and negative values mean it lies inside the circle.
		 */
		public double clearance;

		void set(Node<T> n, double metricDist, double R) {
			this.value = n.value;
			this.cx = n.x;
			this.cy = n.y;
			this.cr = n.r;
			this.metricDist = metricDist;
			this.clearance = metricDist - R;
		}
	}
}
