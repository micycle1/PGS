package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;

import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import processing.core.PVector;

/**
 * Circle packing of rectangle boundaries using the front-chain packing
 * algorithm from 'Visualization of Large Hierarchical Data by Circle Packing'.
 * <p>
 * The algorithm initialises circles in the boundary center and builds up the
 * packing in a spiral pattern from this center until it reaches the rectangle
 * boundary.
 * 
 * @author Mike Bostock
 * @author Java port and modifications by Michael Carleton
 *
 */
public class FrontChainPacker {

	// https://observablehq.com/@mbostock/packing-circles-inside-a-rectangle
	// 'Visualization of large hierarchical data by circle packing'

	private final float width, height;
	private final float offsetX, offsetY;
	private final float radiusMin, radiusMax;
	private final XoRoShiRo128PlusRandom rand;

	/**
	 * The square of the max euclidean distance between a circle center (of
	 * radiusMax) and the center of the boundary. The packing terminates when the
	 * center point of the next circle to place falls outside this distance.
	 */
	private final float maxDistSq;

	private final List<PVector> circles;

	/**
	 * Creates a FrontChainPacker instance. Circles are packed upon initialisation.
	 * <p>
	 * Each circle in the output packing is prescribed a random radius between the
	 * range given.
	 * 
	 * @param width     width of rectangle boundary to pack
	 * @param height    height of rectangle boundary to pack
	 * @param radiusMin minimum radius of circles in the packing
	 * @param radiusMax maximum radius of circles in the packing#
	 * @see #FrontChainPacker(float, float, float, float, float, float)
	 */
	public FrontChainPacker(float width, float height, float radiusMin, float radiusMax) {
		this(width, height, radiusMin, radiusMax, 0, 0, System.nanoTime());
	}

	/**
	 * Creates a FrontChainPacker instance. Circles are packed upon initialisation.
	 * Each circle in the output packing is prescribed a random radius between a
	 * radius range given by its minimum and maximum values.
	 * 
	 * @param width     width of rectangle boundary to pack
	 * @param height    height of rectangle boundary to pack
	 * @param radiusMin minimum radius of circles in the packing
	 * @param radiusMax maximum radius of circles in the packing
	 * @see #FrontChainPacker(float, float, float, float)
	 */
	public FrontChainPacker(float width, float height, float radiusMin, float radiusMax, float offsetX, float offsetY, long seed) {
		this.width = width;
		this.height = height;
		this.radiusMin = Math.max(1f, Math.min(radiusMin, radiusMax)); // choose min and constrain
		this.radiusMax = Math.max(1f, Math.max(radiusMin, radiusMax)); // choose max and constrain
		this.offsetX = offsetX;
		this.offsetY = offsetY;
		this.maxDistSq = this.width * this.height / 2 + this.radiusMax * this.radiusMax;
		rand = new XoRoShiRo128PlusRandom(seed);

		this.circles = pack(new ArrayList<>());
	}

	public List<PVector> getCircles() {
		return circles;
	}

	private List<PVector> pack(List<PVector> circles) {
		// init first chain of 3
		circles.add(new PVector(0, 0, randomRadius()));
		circles.add(new PVector(0, 0, randomRadius()));
		circles.add(new PVector(0, 0, randomRadius()));

		PVector A, B, C;
		A = circles.get(0); // Place the first circle.
		B = circles.get(1); // Place the second circle.
		A.x = -B.z;
		B.x = A.z;
		B.y = 0;

		place(B, A, C = circles.get(2));

		Node a, b, c;

		// Initialize the front-chain using the first three circles a, b and c.
		a = new Node(A);
		b = new Node(B);
		c = new Node(C);
		a.next = c.previous = b;
		b.next = a.previous = c;
		c.next = b.previous = a;

		circles.add(new PVector(0, 0, randomRadius()));

		int iter = 0;
		pack: while (place(a.c, b.c, C = circles.get(circles.size() - 1))) {
			c = new Node(C);

			// Find the closest intersecting circle on the front-chain, if any.
			// “Closeness” is determined by linear distance along the front-chain.
			// “Ahead” or “behind” is likewise determined by linear distance.
			Node j = b.next;
			Node k = a.previous;
			float sj = b.c.z;
			float sk = a.c.z;
			do {
				if (++iter > 10) {
					break; // prevent infinite loop when large difference between the given min & max radii
				}
				if (sj <= sk) {
					if (intersects(j.c, c.c)) {
						b = j;
						a.next = b;
						b.previous = a;
						continue pack;
					}
					sj += j.c.z;
					j = j.next;
				} else {
					if (intersects(k.c, c.c)) {
						a = k;
						a.next = b;
						b.previous = a;
						continue pack;
					}
					sk += k.c.z;
					k = k.previous;
				}
			} while (j != k.next);
			iter = 0;

			// Success! Insert the new circle c between a and b.
			c.previous = a;
			c.next = b;
			a.next = b.previous = b = c;

			// Compute the new closest circle pair to the centroid.
			float aa = score(a);
			float ca;
			while ((c = c.next) != b) {
				if ((ca = score(c)) < aa) {
					a = c;
					aa = ca;
				}
			}
			b = a.next;
			circles.add(new PVector(0, 0, randomRadius())); // last was within bounds; add new circle
		}

		circles.forEach(cl -> { // translate to corner, then offset (if applicable)
			cl.x += (width / 2f) + offsetX;
			cl.y += (height / 2f) + offsetY;
		});

		return circles;
	}

	/**
	 * Compute the position of a new PVector c, given two other PVectors in the
	 * chain same chain, a, b.
	 * 
	 * @return true if
	 */
	private boolean place(PVector b, PVector a, PVector c) {
		final float dx = b.x - a.x;
		final float dy = b.y - a.y;
		final float d2 = dx * dx + dy * dy;
		final float cx, cy; // coordinates to be computed for c
		if (d2 != 0) {
			final float a2 = (a.z + c.z) * (a.z + c.z);
			final float b2 = (b.z + c.z) * (b.z + c.z);
			if (a2 > b2) {
				final float x = (d2 + b2 - a2) / (2 * d2);
				final float y = (float) Math.sqrt(Math.max(0f, b2 / d2 - x * x));
				cx = b.x - x * dx - y * dy;
				cy = b.y - x * dy + y * dx;
			} else {
				final float x = (d2 + a2 - b2) / (2 * d2);
				final float y = (float) Math.sqrt(Math.max(0f, a2 / d2 - x * x));
				cx = a.x + x * dx - y * dy;
				cy = a.y + x * dy + y * dx;
			}
		} else {
			cx = a.x + c.z;
			cy = a.y;
		}

		if (withinBounds(c)) {
			c.x = cx;
			c.y = cy;
			return true;
		}
		return false;
	}

	/**
	 * Determines whether the candidate circle (represented by PVector) cannot cover
	 * the packing region.
	 * 
	 * @return true if circle can cover region
	 */
	private boolean withinBounds(PVector v) {
		return v.x * v.x + v.y * v.y < maxDistSq;
	}

	private static boolean intersects(PVector a, PVector b) {
		final float dr = a.z + b.z - 1e-6f;
		final float dx = b.x - a.x;
		final float dy = b.y - a.y;
		return dr > 0 && dr * dr > dx * dx + dy * dy;
	}

	private float score(Node node) {
		final PVector a = node.c;
		final PVector b = node.next.c;
		final float ab = a.z + b.z;
		final float cx = (a.x * b.z + b.x * a.z) / ab;
		final float cy = (a.y * b.z + b.y * a.z) / ab;
		return Math.max(Math.abs(cx * height), Math.abs(cy * width));
	}

	private float randomRadius() {
		return radiusMin == radiusMax ? radiusMin : rand.nextFloat(radiusMin, radiusMax);
	}

	private static class Node {

		final PVector c;
		Node next, previous;

		Node(PVector circle) {
			this.c = circle;
			this.next = null;
			this.previous = null;
		}
	}

}
