package org.locationtech.jts.operation.polygonize;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;

/**
 * A quicker version of
 * {@link org.locationtech.jts.operation.polygonize.Polygonizer Polygonizer}
 * that forgoes sorting the shells list (which is not needed when polygonizer
 * processes a mesh-like input).
 * 
 * @author Michael Carleton
 *
 */
@SuppressWarnings("unchecked")
public class QuickPolygonizer extends Polygonizer {

	private boolean isCheckingRingsValid = true;
	private boolean extractOnlyPolygonal;

	/**
	 * Creates a polygonizer, specifying whether a valid polygonal geometry must be
	 * created. If the argument is <code>true</code> then areas may be discarded in
	 * order to ensure that the extracted geometry is a valid polygonal geometry.
	 * 
	 * @param extractOnlyPolygonal true if a valid polygonal geometry should be
	 *                             extracted
	 */
	public QuickPolygonizer(boolean extractOnlyPolygonal) {
		this.extractOnlyPolygonal = extractOnlyPolygonal;
	}

	@Override
	public Collection<Geometry> getPolygons() {
		quickPolygonize();
		return polyList;
	}

	/**
	 * Performs the polygonization, if it has not already been carried out.
	 */

	private void quickPolygonize() {
		// check if already computed
		if (polyList != null) {
			return;
		}
		polyList = new ArrayList<Geometry>();

		// if no geometries were supplied it's possible that graph is null
		if (graph == null) {
			return;
		}

		dangles = graph.deleteDangles();
		cutEdges = graph.deleteCutEdges();
		List<EdgeRing> edgeRingList = graph.getEdgeRings();

		// Debug.printTime("Build Edge Rings");

		List<EdgeRing> validEdgeRingList = new ArrayList<EdgeRing>();
		invalidRingLines = new ArrayList<LineString>();
		if (isCheckingRingsValid) {
			findValidRings(edgeRingList, validEdgeRingList, invalidRingLines);
		} else {
			validEdgeRingList = edgeRingList;
		}
		// Debug.printTime("Validate Rings");

		findShellsAndHoles(validEdgeRingList);
		HoleAssigner.assignHolesToShells(holeList, shellList);

		// order the shells to make any subsequent processing deterministic
//		Collections.sort(shellList, new EdgeRing.EnvelopeComparator()); // NOTE skip sort

		// Debug.printTime("Assign Holes");

		boolean includeAll = true;
		if (extractOnlyPolygonal) {
			findDisjointShells(shellList);
			includeAll = false;
		}
		polyList = extractPolygons(shellList, includeAll);
	}

	private void findValidRings(List<EdgeRing> edgeRingList, List<EdgeRing> validEdgeRingList, List<LineString> invalidRingList) {
		for (Iterator<EdgeRing> i = edgeRingList.iterator(); i.hasNext();) {
			EdgeRing er = i.next();
			if (er.isValid()) {
				validEdgeRingList.add(er);
			} else {
				invalidRingList.add(er.getLineString());
			}
		}
	}

	private void findShellsAndHoles(List<EdgeRing> edgeRingList) {
		holeList = new ArrayList<EdgeRing>();
		shellList = new ArrayList<EdgeRing>();
		for (Iterator<EdgeRing> i = edgeRingList.iterator(); i.hasNext();) {
			EdgeRing er = i.next();
			er.computeHole();
			if (er.isHole()) {
				holeList.add(er);
			} else {
				shellList.add(er);
			}
		}
	}

	private static void findDisjointShells(List<EdgeRing> shellList) {
		findOuterShells(shellList);

		boolean isMoreToScan;
		do {
			isMoreToScan = false;
			for (Iterator<EdgeRing> i = shellList.iterator(); i.hasNext();) {
				EdgeRing er = i.next();
				if (er.isIncludedSet()) {
					continue;
				}
				er.updateIncluded();
				if (!er.isIncludedSet()) {
					isMoreToScan = true;
				}
			}
		} while (isMoreToScan);
	}

	/**
	 * For each outer hole finds and includes a single outer shell. This seeds the
	 * traversal algorithm for finding only polygonal shells.
	 * 
	 * @param shellList the list of shell EdgeRings
	 */
	private static void findOuterShells(List<EdgeRing> shellList) {

		for (Iterator<EdgeRing> i = shellList.iterator(); i.hasNext();) {
			EdgeRing er = i.next();
			EdgeRing outerHoleER = er.getOuterHole();
			if (outerHoleER != null && !outerHoleER.isProcessed()) {
				er.setIncluded(true);
				outerHoleER.setProcessed(true);
			}
		}
	}

	private static List<Geometry> extractPolygons(List<EdgeRing> shellList, boolean includeAll) {
		List<Geometry> polyList = new ArrayList<Geometry>();
		for (Iterator<EdgeRing> i = shellList.iterator(); i.hasNext();) {
			EdgeRing er = i.next();
			if (includeAll || er.isIncluded()) {
				polyList.add(er.getPolygon());
			}
		}
		return polyList;
	}

}
