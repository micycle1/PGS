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
		if (polyList != null) { // check if already computed
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
		List<EdgeRing> validEdgeRingList = new ArrayList<>();
		invalidRingLines = new ArrayList<LineString>();
		if (isCheckingRingsValid) {
			findValidRingz(edgeRingList, validEdgeRingList, invalidRingLines);
		} else {
			validEdgeRingList = edgeRingList;
		}

		findShellsAndHolez(validEdgeRingList);
		HoleAssigner.assignHolesToShells(holeList, shellList);

		// order the shells to make any subsequent processing deterministic
//		Collections.sort(shellList, new EdgeRing.EnvelopeComparator()); // NOTE skip sort

		boolean includeAll = true;
		if (extractOnlyPolygonal) {
			findDisjointShells(shellList);
			includeAll = false;
		}
		polyList = extractPolygons(shellList, includeAll);
	}

	private void findValidRingz(List<EdgeRing> edgeRingList, List<EdgeRing> validEdgeRingList, List<LineString> invalidRingList) {
		for (Iterator<EdgeRing> i = edgeRingList.iterator(); i.hasNext();) {
			EdgeRing er = i.next();
			if (er.isValid()) {
				validEdgeRingList.add(er);
			} else {
				invalidRingList.add(er.getLineString());
			}
		}
	}

	private void findShellsAndHolez(List<EdgeRing> edgeRingList) {
		holeList = new ArrayList<>();
		shellList = new ArrayList<>();
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
		shellList.forEach(er -> {
			EdgeRing outerHoleER = er.getOuterHole();
			if (outerHoleER != null && !outerHoleER.isProcessed()) {
				er.setIncluded(true);
				outerHoleER.setProcessed(true);
			}
		});
	}

	private static List<Geometry> extractPolygons(List<EdgeRing> shellList, boolean includeAll) {
		List<Geometry> polyList = new ArrayList<>();
		shellList.forEach(ring -> {
			if (includeAll || ring.isIncluded()) {
				polyList.add(ring.getPolygon());
			}
		});
		return polyList;
	}

}
