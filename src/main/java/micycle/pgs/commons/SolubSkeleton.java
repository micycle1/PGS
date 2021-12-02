package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Stream;

import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.geom.Coordinate;

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

	private final double tol;

	public List<Coordinate> vertices;
	private List<Ray> rays;
	/**
	 * Branches connect the inner bones to the vertices of the polygon.
	 */
	public List<Bone> branches;
	/**
	 * Bones are the inner straight skeleton.
	 */
	public List<Bone> bones;

	public List<Edge> edges;

	/**
	 * 
	 * @param vertices unclosed list of vertices going clockwise
	 * @param tol
	 */
	public SolubSkeleton(List<Coordinate> vertices, double tol) {
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
	 * Each ray object also stores the Coordinates of its neighboring polygon edges
	 * (previous and next).
	 */
	private void init() {
		int lv = vertices.size();

		for (int i = 0; i < vertices.size(); i++) {
			Coordinate pv = vertices.get(Math.floorMod(i - 1, lv)); // previous vertex
			Coordinate cv = vertices.get(i); // current vertex
			Coordinate nv = vertices.get((i + 1) % lv); // next vertex

			Coordinate b = bisector(pv, cv, cv, nv);
			Coordinate ept = add(cv, setMag(b, 100000)); // end-point of the ray (far from start-point)

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

				double minD = 100000; // min-distance
				Coordinate X = null; // intersection point (optional)

				// check current ray for intersection with prev and next

				Coordinate x = intersect(cr, pRay); // check with prev
				if (x != null) { // rays intersected
					double d = x.distance(cr.startPoint);
					if (d < minD) {
						minD = d;
						X = x;
					}
				}

				x = intersect(cr, nRay); // check with next
				if (x != null) { // rays intersected
					double d = x.distance(cr.startPoint);
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
						double dSum = i1.minD + i2.minD;
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
							double dSum = i1.minD + i2.minD; // sum of distances from selected rays (s-pt) to their
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
					Coordinate X = ray.intersectionPoint;
					Ray r1 = rays.get(id - offset);
					Ray r2 = rays.get((id + 1 - offset) % lr);

					// stores bones (segments from parent rays to intersection point 'X')
					Stream.of(r1, r2).forEach(r -> {
						if (!vertices.contains(r.startPoint)) {
							double d1 = X.distance(r1.prevEdge.p2);
							double d2 = X.distance(r2.nextEdge.p1);
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
					Coordinate b = bisector(r1, r2);
					Coordinate ept = add(X.copy(), setMag(b, 100000));

					// check if ray points in the right direction
					int si = vertices.indexOf(r2.nextEdge.p1); // #start-pt index
					int ei = vertices.indexOf(r1.prevEdge.p2); // #end-pt index

					ArrayList<Coordinate> slice = new ArrayList<>();
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

					if (!contains(slice, add(X.copy(), setMag(b, 1)))) {
						ept = mult(ept, -1);
					}

					// delete parents rays from array list and insert their 'child' instead
					rays.set(id - offset, new Ray(X, ept, r1.prevEdge, r2.nextEdge));
					rays.remove((id + 1 - offset) % lr);
					offset++;
					lr = rays.size();

				}
			} else {
				System.err.println("Error: no rays have been found for reduction. A shared intersection point is probably missing.");
				break;
			}
		}
		bones.add(new Bone(rays.get(0).startPoint, rays.get(1).startPoint)); // connect start-points of the last 2 rays
	}

	/**
	 * Computes the bisector of a corner between edge p1-p2 and edge p3-p4. p2 & p3
	 * are usually the same.
	 */
	private static Coordinate bisector(Coordinate p1, Coordinate p2, Coordinate p3, Coordinate p4) {

		Coordinate dir1 = normalize(sub(p2, p1)); // direction between 1st point and 2nd point of edge 1
		Coordinate dir2 = normalize(sub(p4, p3)); // direction between 1st point and 2nd point of edge 2

		Coordinate dsum = add(dir1, dir2);

		return normalize(new Coordinate(dsum.y, -dsum.x));
	}

	private static Coordinate sub(Coordinate a, Coordinate b) {
		return new Coordinate(a.x - b.x, a.y - b.y);
	}

	private static Coordinate add(Coordinate a, Coordinate b) {
		return new Coordinate(a.x + b.x, a.y + b.y);
	}

	private static Coordinate normalize(Coordinate a) {
		double mag = Math.sqrt(a.x * a.x + a.y * a.y);
		return new Coordinate(a.x / mag, a.y / mag);
	}

	private static Coordinate mult(Coordinate a, double mult) {
		return new Coordinate(a.x * mult, a.y * mult);
	}

	private static Coordinate setMag(Coordinate a, double mag) {
		return mult(normalize(a), mag);
	}

	private static Coordinate bisector(Ray r1, Ray r2) {
		return bisector(r1.prevEdge.p1, r1.prevEdge.p2, r2.nextEdge.p1, r2.nextEdge.p2);
	}

	/**
	 * Checks if 2 lines are intersecting. Optional: returns location of
	 * intersection point.
	 */
	private static Coordinate intersect(Coordinate s1, Coordinate e1, Coordinate s2, Coordinate e2) {
		LineIntersector lineIntersector = new RobustLineIntersector();
		lineIntersector.computeIntersection(s1, e1, s2, e2);

		if (lineIntersector.hasIntersection()) {
			Coordinate x = lineIntersector.getIntersection(0);
			return new Coordinate(Math.round(x.x), Math.round(x.y));
		}
		return null;
	}

	private static Coordinate intersect(Ray r1, Ray r2) {
		return intersect(r1.startPoint, r1.endPoint, r2.startPoint, r2.endPoint);
	}

	private static boolean contains(ArrayList<Coordinate> verts, Coordinate pt) {

		int N = verts.size();
		boolean isInside = false;

		for (int i = 0; i < N; i++) {
			Coordinate v1 = verts.get((i + 1) % N);
			Coordinate v2 = verts.get(i);

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
	private class Ray {

		Coordinate startPoint; // the vertex ray originated from
		Coordinate endPoint; // live end point of ray
		Edge prevEdge; // edge from previous vertex to startpoint of this ray
		Edge nextEdge; // from to next vertex

		private Ray(Coordinate startPoint, Coordinate endPoint, Edge e1, Edge e2) {
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

		public Coordinate p1; // vertex from
		public Coordinate p2; // vertex to

		public Edge(Coordinate p1, Coordinate p2) {
			this.p1 = p1;
			this.p2 = p2;
		}
	}

	public class Bone {

		public Coordinate sp1; // startPoint of ray 1
		public Coordinate sp2; // startPoint of ray 2

		public Bone(Coordinate sp1, Coordinate sp2) {
			this.sp1 = sp1;
			this.sp2 = sp2;
		}
	}

	/**
	 * Simple container
	 */
	private class Intersection {

		double minD;
		Coordinate intersectionPoint;

		public Intersection(Coordinate intersectionPoint, double minD) {
			this.minD = minD;
			this.intersectionPoint = intersectionPoint;
		}

		/**
		 * Where the rays cross
		 */
		Coordinate X() {
			return intersectionPoint;
		}

//		@Override
//		public int hashCode() {
//			return intersectionPoint.hashCode()^(int)Double.doubleToRawLongBits(minD);
//		}
//
//		@Override
//		public boolean equals(Object obj) {
//			if (obj instanceof Intersection) {
//				Intersection o = (Intersection) obj;
//				return o.intersectionPoint.equals2D(intersectionPoint) && o.minD == minD;
//			}
//			return false;
//		}
	}

	/**
	 * Simple container
	 */
	private class Candidate implements Comparable<Candidate> {

		int id;
		Coordinate intersectionPoint;
		double minD;

		public Candidate(int id, Coordinate intersectionPoint, double minD) {
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

		@Override
		public boolean equals(Object obj) {
			if (obj instanceof Candidate) {
				Candidate other = (Candidate) obj;
				return other.minD == minD && other.intersectionPoint.equals(intersectionPoint);
			}
			return false;
		}

		@Override
		public int hashCode() {
			return 1327 * id + intersectionPoint.hashCode();
		}
	}

}
