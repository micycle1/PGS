package micycle.pgs.utility;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import processing.core.PVector;

/**
 * Quadrangulations of Planar Point Sets via the Spiraling Rotating Calipers
 * Algorithm.
 * <p>
 * Implements 'Characterizing and efficiently computing quadrangulations of
 * planar point sets' by Prosenjit Bose and Godfried T. Toussaint.
 * <p>
 * This (Processing) implementation is derived from a highly coupled old Java
 * applet written by Martin Blais.
 * 
 * @author Martin Blais
 * @author Michael Carleton
 *
 */
public class SpiralQuadrangulation {

	// http://www-cgrl.cs.mcgill.ca/~godfried/teaching/cg-projects/97/Blais/quadrang/algorithm.html

	private List<PVector> vertices;
	private PVector vMinx;
	private List<PEdge> spiral;
	private List<PEdge> closedChain;
	private PVector vPh;
	private List<PEdge> innerDiagonals; // From inside (center) towards outside
	private PVector vY;
	private List<PEdge> chainDiagonals; // From outside towards inside (center)
	private List<PEdge> quadrangulation;

	/**
	 * Instantiates a Spiral Quadrangulation. The quadrangulation is computed upon
	 * instantiation.
	 * 
	 * @param points no duplicates
	 */
	public SpiralQuadrangulation(List<PVector> points) {
		vertices = new ArrayList<>(points);
		spiral = new ArrayList<>();
		closedChain = new ArrayList<>();
		innerDiagonals = new ArrayList<>();
		chainDiagonals = new ArrayList<>();
		quadrangulation = new ArrayList<>();

		vertices.sort((a, b) -> a.x > b.x ? 1 : 0);
		for (int ii = 0; ii < vertices.size() - 1; ii++) {
			PVector v1 = vertices.get(ii);
			PVector v2 = vertices.get(ii + 1);
			if (v1.equals(v2)) {
				vertices.remove(ii);
				ii--;
			}
		}

		computeConvexSpiral();
		computeCloseChain();
		computeInnerDiagonals();
		computeChainDiagonals();
		computeRemoveDiagonals();
		computeSteiner();
		computeFinalQuadrang();
	}

	public List<PEdge> getQuadrangulationEdges() {
		return quadrangulation;
	}

	/**
	 * Computes a spiral polygonal chain of the point set by using the rotating
	 * calipers method (the convex spiral of the set S) in O(n^2).
	 */
	private void computeConvexSpiral() {
		if (vertices.size() <= 2) {
			return;
		}

		// Find min x vertex O(n)
		vMinx = new PVector(Float.MAX_VALUE, Float.MAX_VALUE);
		vertices.forEach(v -> {
			if (v.x < vMinx.x || (v.x == vMinx.x && v.y < vMinx.y)) {
				vMinx = v;
			}
		});

		// Clone the vertex array, in order to be able to mark the vertices
		List<PVector> rgVertices = new ArrayList<>(vertices);
		List<PVector> rgSpiral = new ArrayList<>();

		// Start spiral with minx vertex
		rgSpiral.add(vMinx);
		double anglecur;

		while (!rgVertices.isEmpty()) {
			if (rgSpiral.size() < 2) {
				anglecur = Math.PI / 2;
			} else {
				// Compute angle of last two vertices in spiral chain
				PVector last = rgSpiral.get(rgSpiral.size() - 1);
				PVector last2 = rgSpiral.get(rgSpiral.size() - 2);

				anglecur = Math.atan(-(last.y - last2.y) / (last.x - last2.x));
				if (last.x < last2.x) {
					anglecur += Math.PI;
					// anglecur is in [-pi/2,+3pi/2]
				}
			}

			// Jarvis march step
			//
			// For each of the other vertices, wrap angle around to find next
			// CH vertex O(n)
			PVector vcur = rgSpiral.get(rgSpiral.size() - 1);
			double anglemin = anglecur + 3 * Math.PI / 2 + 0.0001; // REVIEW
			PVector vnext = null;

			for (PVector v : rgVertices) {
				// Compute wrap angle
				double angle = Math.atan(-(v.y - vcur.y) / (v.x - vcur.x));
				if (v.x < vcur.x) {
					angle += Math.PI;
				}

				angle = anglecur - angle;
				// Make sure difference is positive
				if (angle < 0) {
					angle += 2 * Math.PI;
				}

				if (angle < anglemin) {
					anglemin = angle;
					vnext = v;
				}
			}
			rgVertices.remove(vnext);
			rgSpiral.add(vnext);
		}

		// Add vertices to spiral
		Iterator<PVector> iterator = rgSpiral.iterator();
		PVector v1 = iterator.next();
		PVector v2;
		while (iterator.hasNext()) {
			v2 = iterator.next();
			spiral.add(new PEdge(v1, v2));
			v1 = v2;
		}
	}

	/**
	 * Close the polygonal chain between its start and end.
	 */
	private void computeCloseChain() {
		if (vertices.size() <= 2) {
			return;
		}

		// Start spiral with minx vertex
		double anglecur = Math.PI / 2;

		// Jarvis march step
		// For each of the other vertices, wrap angle around to find vertex in O(n)
		double anglemax = anglecur - Math.PI / 2 - 0.0001;
		vPh = null;

		for (PVector v : vertices) {
			// Compute wrap angle
			double angle = Math.atan(-(v.y - vMinx.y) / (v.x - vMinx.x));
			if (v.x < vMinx.x) {
				angle += Math.PI;
			}

			angle = anglecur - angle;
			// Make sure difference is positive
			if (angle < 0) {
				angle += 2 * Math.PI;
			}

			if (angle > anglemax) {
				anglemax = angle;
				vPh = v;
			}
		}

		closedChain.add(new PEdge(vMinx, vPh));
	}

	/**
	 * Split the spiral into inner/outer regions by extending the line through the
	 * polygonal chain. The inner region is star-shaped, and therefore we
	 * triangulate it easily.
	 */
	private void computeInnerDiagonals() {
		if (vertices.size() <= 2) {
			return;
		}

		// Create long long long long long final line segment (extend)
		final int FARAWAYMULT = 100;
		PEdge elast = spiral.get(spiral.size() - 1);
		if (elast != null) {
			elast = elast.clone();
		} else {
			return;
		}

		PVector vfaraway = PVector.add(PVector.mult(elast.a, 1 - FARAWAYMULT), PVector.mult(elast.b, FARAWAYMULT));
		elast.b.set(vfaraway);

		// Find the line segment which intersects with the extended line segment
		int ii;
		for (ii = spiral.size() - 3; ii >= 0; ii--) {
			PEdge ecomp = spiral.get(ii);

			if (intersectSegments(elast, ecomp)) {
				break;
			}
		}

		elast = spiral.get(spiral.size() - 1);
		if (ii < 0) {
			// Vertices are all on convex hull
			PEdge efirst = spiral.get(0);
			vY = efirst.b;
		} else {
			// Some vertices are inside the convex hull, normal case
			// Search for last vertex parallel to last segment (using rotating calipers)

			double prevdist = -1;
			PEdge edg = null;
			for (double dist = 0; dist > prevdist; ii--, prevdist = dist) {
				edg = spiral.get(ii);
				dist = distance(elast.a, elast.b, edg.a);
			}
			// Add the vertex to close the inner star-shaped region.
			vY = edg.a;
		}

		// Triangulate the star-shaped polygon
		int iii = spiral.size() - 2;
		PEdge edg;
		do {
			edg = spiral.get(iii--);
			innerDiagonals.add(new PEdge(elast.b, edg.a));
		} while (edg.a != vY);
	}

	/**
	 * Triangulate the outer region of the spiral, starting at the inner end.
	 */
	private void computeChainDiagonals() {
		if ((vertices.size() <= 2) || (vY == null)) {
			return;
		}

		// Create polygonal annulusses
		List<PVector> rgAnnulusOuter = new ArrayList<>();
		List<PVector> rgAnnulusInner = new ArrayList<>();

		///////////////////////////////////////////////////////////////////////
		// Create the outer annulus [p_1, ... , Y]

		Iterator<PEdge> e = spiral.iterator();
		PEdge edg;
		do {
			edg = e.next();
			rgAnnulusOuter.add(edg.a);
		} while (edg.a != vY && e.hasNext());

		///////////////////////////////////////////////////////////////////////
		// Create the inner annulus [p_h, ... , p_n]

		// Skip vertices before Ph
		e = spiral.iterator();
		PVector vv;
		do {
			edg = e.next();
			vv = edg.b;
		} while (vv != vPh && e.hasNext());

		rgAnnulusInner.add(vPh);
		// Add all the remaining edges to the form the inner annulus
		while (e.hasNext()) {
			PEdge edg1 = e.next();
			rgAnnulusInner.add(edg1.b);
		}

		///////////////////////////////////////////////////////////////////////
		// Run Spiraling Rotating Calipers algorithm on polygonal annulusses
		//
		PVector vouter = rgAnnulusOuter.get(0);
		PVector vinner = rgAnnulusInner.get(0);
		double angleSRC = angleBetween(vinner, vouter);
		double angleouter;
		double angleinner;
		int iouter = 0;
		int iinner = 0;
		PVector vnextouter = null;
		PVector vnextinner = null;
		// Compute angle of next two vertices in annulusses
		while (iinner != rgAnnulusInner.size() - 1 || iouter != rgAnnulusOuter.size() - 1) {
			///////////////////////////////////////////////////////////////////
			// Special case 1: Simply complete advance on outer chain

			if (iinner == rgAnnulusInner.size() - 1) {
				vnextouter = rgAnnulusOuter.get(iouter + 1);
				chainDiagonals.add(new PEdge(vinner, vnextouter));
				iouter++;
				vouter = vnextouter;
				continue;
			}

			///////////////////////////////////////////////////////////////////
			// Special case 1: Simply complete advance on inner chain

			if (iouter == rgAnnulusOuter.size() - 1) {
				vnextinner = rgAnnulusInner.get(iinner + 1);
				chainDiagonals.add(new PEdge(vnextinner, vouter));
				iinner++;
				vinner = vnextinner;
				continue;
			}

			///////////////////////////////////////////////////////////////////
			// Normal case

			// Compute angles for next segments on polygonal chains
			vnextinner = rgAnnulusInner.get(iinner + 1);
			vnextouter = rgAnnulusOuter.get(iouter + 1);

			angleinner = angleBetween(vinner, vnextinner);
			angleouter = angleBetween(vouter, vnextouter);

			// Compute differences between angles and current SRC angle
			double diffouter = angleSRC - angleouter;
			if (diffouter < 0) {
				diffouter += 2 * Math.PI;
			}
			double diffinner = angleSRC - angleinner;
			if (diffinner < 0) {
				diffinner += 2 * Math.PI;
			}

			// Choose the caliper side with the smallest difference
			if (diffinner < diffouter) {
				// Inner was hit by caliper first
				chainDiagonals.add(new PEdge(vnextinner, vouter));
				iinner++;
				vinner = vnextinner;
				angleSRC = angleinner;
			} else {
				// Outer was hit by caliper first
				chainDiagonals.add(new PEdge(vinner, vnextouter));
				iouter++;
				vouter = vnextouter;
				angleSRC = angleouter;
			}
		}

		// Remove last element, it was already inserted in the inner
		// star-shaped chain
		chainDiagonals.remove(chainDiagonals.size() - 1);
	}

	/**
	 * Remove one out of 2 diagonals of the whole spiral triangulation, starting at
	 * the center. This way we obtain quadrilaterals.
	 */
	private void computeRemoveDiagonals() {
		if (vertices.size() <= 2) {
			return;
		}

		int correct = 0;
		if (innerDiagonals.size() % 2 == 1) {
			correct = 1;
		}

		///////////////////////////////////////////////////////////////////////
		// Start removing one out of two diagonals from inner star-shaped region
		//
		List<PEdge> rgNewInnerDiagonals = new ArrayList<>();
		for (int ii = 0; ii < innerDiagonals.size() - 1; ii += 2) {
			rgNewInnerDiagonals.add(innerDiagonals.get(ii + 1));
		}
		innerDiagonals = rgNewInnerDiagonals;

		///////////////////////////////////////////////////////////////////////
		// Start removing one out of two diagonals from outer star-shaped
		// poylgonal annulus
		//

		List<PEdge> rgNewChainDiagonals = new ArrayList<>();
		for (int ii = chainDiagonals.size() - 1 + correct; ii >= 1; ii -= 2) {
			rgNewChainDiagonals.add(chainDiagonals.get(ii - 1));
		}
		chainDiagonals = rgNewChainDiagonals;
	}

	/**
	 * Add a Steiner point to complete the last polygon of the chain, which may
	 * otherwise be a triangle.
	 */
	private void computeSteiner() {
		if ((vertices.size() <= 2) || (vPh == null)) {
			return;
		}

		// Count the number of vertices on the convex hull [P1, ..., Ph] incl.
		Iterator<PEdge> e = spiral.iterator();
		PEdge edg;
		int inbverts = 0;
		do {
			edg = e.next();
			inbverts++;
		} while (edg.b != vPh);
		inbverts++;

		if (inbverts % 2 == 1) {
			// The number of vertices on the CH is odd, add a Steiner point
			edg = closedChain.get(0);
			PVector vSteiner = new PVector((edg.a.x + edg.b.x) / 2, (edg.a.y + edg.b.y) / 2);
			vertices.add(vSteiner);

			// Change the closing chain to include the Steiner point
			closedChain.clear();
			closedChain.add(new PEdge(edg.a, vSteiner));
			closedChain.add(new PEdge(vSteiner, edg.b));
		}
	}

	private void computeFinalQuadrang() {
		if (vertices.size() <= 2) {
			return;
		}
		quadrangulation.clear();

		// Simply put all the edges indiscriminately in a quadrangulation
		spiral.forEach(quadrangulation::add);
		closedChain.forEach(quadrangulation::add);
		innerDiagonals.forEach(quadrangulation::add);
		chainDiagonals.forEach(quadrangulation::add);
	}

	/**
	 * Determines whether two segments intersect.
	 */
	private static boolean intersectSegments(PEdge e1, PEdge e2) {
		boolean f11 = (signedArea(e1.a, e1.b, e2.a) >= 0);
		boolean f12 = (signedArea(e1.a, e1.b, e2.b) >= 0);
		boolean f21 = (signedArea(e2.a, e2.b, e1.a) >= 0);
		boolean f22 = (signedArea(e2.a, e2.b, e1.b) >= 0);
		return f11 != f12 && f21 != f22;
	}

	/**
	 * Signed area of a triangle.
	 */
	private static double signedArea(PVector v1, PVector v2, PVector v3) {
		return (v1.x * (v2.y - v3.y) + v2.x * (v3.y - v1.y) + v3.x * (v1.y - v2.y)) / 2.0;
	}

	/**
	 * Perpendicular distance from an infinite line to a point.
	 * <p>
	 * Note: Distance is positive if the point is to the left of the line. This
	 * convention has been chosen so that it is convenient to think of it as the
	 * right-hand rule.
	 */
	private static double distance(PVector vLine1, PVector vLine2, PVector vPoint) {
		// b*h/2 = area... h = area*2/b
		final PVector vline = PVector.sub(vLine2, vLine1);
		return signedArea(vLine1, vLine2, vPoint) / vline.mag();
	}

	/**
	 * Finds the angle determined by points v1 and v2.
	 * 
	 * @param v1
	 * @param v2
	 * @return an angle between [-pi/2, +3pi/2]
	 */
	private static double angleBetween(PVector v1, PVector v2) {
		double angle = Math.atan(-(v2.y - v1.y) / (v2.x - v1.x));
		if (v2.x < v1.x) {
			angle += Math.PI;
		}
		return angle;
	}

}
