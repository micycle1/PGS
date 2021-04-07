package micycle.pgs.utility;

import static processing.core.PApplet.round;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Stream;

import processing.core.PVector;

/**
 * Algorithm:
 * <p>
 * 1. Calculate the bisecting vector at each vertex of a chosen polygon. Make it
 * a ray.<br>
 * 2. Find the intersection points with the previous and the next neighboring
 * ray.<br>
 * 3. Find all pairs of rays that intersect sooner than either intersects its
 * other neighbor.<br>
 * 4. Sort them by distance to their intersection point.<br>
 * 5. Take the pair with the shortest distance and replace it with a single new
 * ray originating from their intersection point. That child ray must be aiming
 * towards the average direction of its neighboring polygon edges.<br>
 * 6. Repeat from step 2 until there are only 2 rays left.<br>
 * 7. Connect the originating points of the last 2 rays to complete the
 * skeleton. *
 * 
 * @author Solub
 * @author Java port from Python by Michael Carleton
 */
public class SolubSkeleton {

	/**
	 * https://discourse.processing.org/t/17208/7
	 */

	private final float tol;

	public ArrayList<PVector> vertices;
	public ArrayList<Ray> rays;
	/**
	 * Branches connect the inner bones to the vertices of the polygon.
	 */
	public ArrayList<Bone> branches;
	/**
	 * Bones are the inner straight skeleton.
	 */
	public ArrayList<Bone> bones;

	public ArrayList<Edge> edges;

	/**
	 * 
	 * @param vertices unclosed list of vertices going clockwise
	 * @param tol
	 */
	public SolubSkeleton(List<PVector> vertices, float tol) {
		this.vertices = new ArrayList<>(vertices);
		this.tol = tol;
		rays = new ArrayList<>();
		bones = new ArrayList<>();
		branches = new ArrayList<>();
		edges = new ArrayList<>();
	}

	public void run() {
		init();
		if (vertices.size() > 3) {
			reduce();
		}
	}

	/**
	 * Starts by computing the rays originating from the corners of the polygon.
	 * Each ray object also stores the PVectors of its neighboring polygon edges
	 * (previous and next).
	 */
	private void init() {
		int lv = vertices.size();

		for (int i = 0; i < vertices.size(); i++) {
			PVector pv = vertices.get(Math.floorMod(i - 1, lv)); // previous vertex
			PVector cv = vertices.get(i); // current vertex
			PVector nv = vertices.get((i + 1) % lv); // next vertex

			PVector b = bisector(pv, cv, cv, nv);
			PVector ept = b.setMag(100000).add(cv); // end-point of the ray (far from start-point)

			Edge in = new Edge(pv, cv);
			Edge out = new Edge(cv, nv);
			edges.add(in);
			edges.add(out);

			Ray r = new Ray(cv, ept, in, out);
			rays.add(r);
		}

	}

	private void reduce() {

		/*
		 * Check for ray-ray collision
		 */
		while (rays.size() > 2) {

			ArrayList<Intersection> intersections = new ArrayList<>(); // Nearest intersection points from parent
			int lr = rays.size();

			for (int i = 0; i < lr; i++) {
				Ray pRay = rays.get(Math.floorMod(i - 1, lr)); // previous ray
				Ray cr = rays.get(i); // current ray
				Ray nRay = rays.get((i + 1) % lr); // next ray

				float minD = 100000; // min-distance
				PVector X = null; // intersection point (optional)

				// check current ray for intersection with prev and next

				PVector x = intersect(cr, pRay); // check with prev
				if (x != null) { // rays intersected
					float d = x.dist(cr.startPoint);
					if (d < minD) {
						minD = d;
						X = x;
					}
				}

				x = intersect(cr, nRay); // check with next
				if (x != null) { // rays intersected
					float d = x.dist(cr.startPoint);
					if (d < minD) {
						minD = d;
						X = x;
					}
				}

				intersections.add(new Intersection(X, minD));
			}

			ArrayList<Candidate> candidates = new ArrayList<>();

			for (int i = 0; i < intersections.size(); i++) {

				/**
				 * zip(intersections, intersections[1:]+intersections[:1]) is equivalent to
				 * iterating over consecutive pairs of Iterations
				 */

				Intersection i1 = intersections.get(i);
				Intersection i2 = intersections.get((i + 1) % intersections.size()); // wrap around to 0 for last pair
				if ((i1.X() != null && i2.X() != null)) { // if D not null and equivalent

					if (i1.X().equals(i2.X())) {
						// sum of distances from selected rays (s-pt) to their intersection point
						float dSum = i1.minD + i2.minD;
						candidates.add(new Candidate(i, i1.X(), dSum));
					}
				}
			}

			if (candidates.isEmpty()) {
				// making sure that NO consecutive rays intersect each other before either
				// intersects its other neighbor
				if (new HashSet<>(intersections).size() == intersections.size()) {
					// similar loop to above
					for (int id = 0; id < intersections.size(); id++) {
						// iterate over consecutive pairs of Iterations
						Intersection i1 = intersections.get(id);
						Intersection i2 = intersections.get((id + 1) % intersections.size()); // wrap around to 0
						if (i1.X() != null && i2.X() != null) {
							// sending to candidates anyway
							float dSum = i1.minD + i2.minD; // sum of distances from selected rays (s-pt) to their
															// intersection point
							candidates.add(new Candidate(id, i1.X(), dSum));
						}
					}

				}
			}

			// sort rays by distance from their s-pt to their intersection pt
			ArrayList<Candidate> srays = new ArrayList<>(candidates);
			Collections.sort(srays);

			ArrayList<Candidate> selectedRays = new ArrayList<>();
			for (Candidate r : srays) {
				if (r.minD == srays.get(0).minD) {
					selectedRays.add(r);
				}
			}

			if (!selectedRays.isEmpty()) {
				int offset = 0;

				for (Candidate ray : selectedRays) {
					int id = ray.id;
					PVector X = ray.intersectionPoint;
					Ray r1 = rays.get(id - offset);
					Ray r2 = rays.get((id + 1 - offset) % lr);

					// stores bones (segments from parent rays to intersection point 'X')
					Stream.of(r1, r2).forEach(r -> {
						if (!vertices.contains(r.startPoint)) {
							float d1 = X.dist(r1.prevEdge.p2);
							float d2 = X.dist(r2.nextEdge.p1);
							if ((d1 + d2) / 2f > tol) {
								bones.add(new Bone(r.startPoint, X));
							} else {
								branches.add(new Bone(r.startPoint, X));
							}
						} else {
							branches.add(new Bone(r.startPoint, X));
						}
					});

					// compute direction of the new ray
					PVector b = bisector(r1, r2);
					PVector ept = X.copy().add(b.setMag(100000));

					// check if ray points in the right direction
					int si = vertices.indexOf(r2.nextEdge.p1); // #start-pt index
					int ei = vertices.indexOf(r1.prevEdge.p2); // #end-pt index

					ArrayList<PVector> slice = new ArrayList<>();
					slice.add(X); // moved from '[X] + slice' to here
					if (ei > si) {
						for (int i = si; i < ei; i++) {
							slice.add(vertices.get(i)); // vertices[si:ei]
						}
					} else {
						for (int i = si; i < vertices.size(); i++) {
							slice.add(vertices.get(i)); // vertices[si:]
						}
						for (int i = 0; i < ei; i++) {
							slice.add(vertices.get(i)); // self.vertices[:ei]
						}
					}

					if (!contains(slice, X.copy().add(b.setMag(1)))) {
						ept.mult(-1);
					}

					// delete parents rays from array list and insert their 'child' instead
					rays.set(id - offset, new Ray(X, ept, r1.prevEdge, r2.nextEdge));
					rays.remove((id + 1 - offset) % lr);
					offset++;
					lr = rays.size();

				}
			} else {
				System.err.println(
						"Error: no rays have been found for reduction. A shared intersection point is probably missing.");
				break;
			}
		}
		bones.add(new Bone(rays.get(0).startPoint, rays.get(1).startPoint)); // connect start-points of the last 2 rays
	}

	/**
	 * Computes the bisector of a corner between edge p1-p2 and edge p3-p4. p2 & p3
	 * are usually the same.
	 */
	private static PVector bisector(PVector p1, PVector p2, PVector p3, PVector p4) {

		PVector dir1 = PVector.sub(p2, p1).normalize(); // direction between 1st point and 2nd point of edge 1
		PVector dir2 = PVector.sub(p4, p3).normalize(); // direction between 1st point and 2nd point of edge 2

		PVector dsum = PVector.add(dir1, dir2);

		return new PVector(dsum.y, -dsum.x).normalize();
	}

	private static PVector bisector(Ray r1, Ray r2) {
		return bisector(r1.prevEdge.p1, r1.prevEdge.p2, r2.nextEdge.p1, r2.nextEdge.p2);
	}

	/**
	 * Checks if 2 lines are intersecting. Optional: returns location of
	 * intersection point.
	 */
	private static PVector intersect(PVector s1, PVector e1, PVector s2, PVector e2) {

		float uA = ((e2.x - s2.x) * (s1.y - s2.y) - (e2.y - s2.y) * (s1.x - s2.x))
				/ ((e2.y - s2.y) * (e1.x - s1.x) - (e2.x - s2.x) * (e1.y - s1.y));
		float uB = ((e1.x - s1.x) * (s1.y - s2.y) - (e1.y - s1.y) * (s1.x - s2.x))
				/ ((e2.y - s2.y) * (e1.x - s1.x) - (e2.x - s2.x) * (e1.y - s1.y));

		float secX, secY;
		if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {

			secX = s1.x + (uA * (e1.x - s1.x));
			secY = s1.y + (uA * (e1.y - s1.y));
			return new PVector(round(secX), round(secY)); // TODO round?
		}
		return null; // no intersection
	}

	private static PVector intersect(Ray r1, Ray r2) {
		return intersect(r1.startPoint, r1.endPoint, r2.startPoint, r2.endPoint);
	}

	private static boolean contains(ArrayList<PVector> verts, PVector pt) {

		int N = verts.size();
		boolean isInside = false;

		for (int i = 0; i < N; i++) {
			PVector v1 = verts.get((i + 1) % N);
			PVector v2 = verts.get(i);

			if ((v2.y > pt.y) != (v1.y > pt.y)) {
				if (pt.x < (v1.x - v2.x) * (pt.y - v2.y) / (v1.y - v2.y) + v2.x) {
					isInside = !isInside;
				}
			}
		}

		return isInside;
	}

	/**
	 * Simple container.
	 * 
	 * Rays originate from the corners of the polygon, moving inwards.
	 */
	public class Ray {

		public PVector startPoint; // the vertex ray originated from
		public PVector endPoint; // live end point of ray
		public Edge prevEdge; // edge from previous vertex to startpoint of this ray
		public Edge nextEdge; // from to next vertex

		public Ray(PVector startPoint, PVector endPoint, Edge e1, Edge e2) {
			this.startPoint = startPoint;
			this.endPoint = endPoint;
			this.prevEdge = e1;
			this.nextEdge = e2;
		}
	}

	/**
	 * Simple container
	 */
	public class Edge {

		public PVector p1; // vertex from
		public PVector p2; // vertex to

		public Edge(PVector p1, PVector p2) {
			this.p1 = p1;
			this.p2 = p2;
		}
	}

	public class Bone {
	
		public PVector sp1; // startPoint of ray 1
		public PVector sp2; // startPoint of ray 2
	
		public Bone(PVector sp1, PVector sp2) {
			this.sp1 = sp1;
			this.sp2 = sp2;
		}
	}

	/**
	 * Simple container
	 */
	private class Intersection {

		float minD;
		PVector intersectionPoint;

		public Intersection(PVector intersectionPoint, float minD) {
			this.minD = minD;
			this.intersectionPoint = intersectionPoint;
		}

		/**
		 * Where the rays cross
		 */
		PVector X() {
			return intersectionPoint;
		}
	}

	/**
	 * Simple container
	 */
	private class Candidate implements Comparable<Candidate> {

		int id;
		PVector intersectionPoint;
		float minD;

		public Candidate(int id, PVector intersectionPoint, float minD) {
			this.id = id;
			this.intersectionPoint = intersectionPoint;
			this.minD = minD;
		}

		@Override
		public int compareTo(Candidate o) {
			if (minD > o.minD) { // this bigger
				return 1;
			}
			if (minD < o.minD) { // this smaller
				return -1;
			}
			return 0; // equal
		}
	}

}
