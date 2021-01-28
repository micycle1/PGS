package micycle.pts;

/*******************************************************************************
 * This software is released under the licence CeCILL
 *
 * see Licence_CeCILL-C_fr.html see Licence_CeCILL-C_en.html
 *
 * see <a href="http://www.cecill.info/">http://www.cecill.info/a>
 *
 * @copyright IGN
 * @copyright twak
 ******************************************************************************/

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.twak.camp.Corner;
import org.twak.camp.Edge;
import org.twak.camp.Machine;
import org.twak.camp.Output.SharedEdge;
import org.twak.camp.Skeleton;
import org.twak.utils.collections.Loop;
import org.twak.utils.collections.LoopL;

import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IDirectPosition;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IDirectPositionList;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineSegment;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineString;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IPolygon;
import fr.ign.cogit.geoxygene.api.spatial.geomaggr.IMultiPoint;
import fr.ign.cogit.geoxygene.api.spatial.geomprim.IPoint;
import fr.ign.cogit.geoxygene.api.spatial.geomprim.IRing;
import fr.ign.cogit.geoxygene.generalisation.Filtering;
import fr.ign.cogit.geoxygene.spatial.coordgeom.DirectPosition;
import fr.ign.cogit.geoxygene.spatial.coordgeom.GM_LineSegment;
import fr.ign.cogit.geoxygene.spatial.coordgeom.GM_LineString;
import fr.ign.cogit.geoxygene.util.algo.geometricAlgorithms.CommonAlgorithmsFromCartAGen;
import fr.ign.cogit.geoxygene.util.algo.geomstructure.Vector2D;

/**
 * Extract the skeleton of polygons by different methods
 *
 * @author Guillaume
 *
 */
public class Skeletonize {

	class MultiSkeleton {
		private Set<ILineSegment> segments = new HashSet<ILineSegment>();
		private List<ILineSegment> segFace = new ArrayList<ILineSegment>();
		private Set<IDirectPosition> isolatedPts = new HashSet<IDirectPosition>();

		public Set<ILineSegment> getSegments() {
			return segments;
		}

		public void setSegments(Set<ILineSegment> segments) {
			this.segments = segments;
		}

		public List<ILineSegment> getSegFace() {
			return segFace;
		}

		public void setSegFace(List<ILineSegment> segFace) {
			this.segFace = segFace;
		}

		public Set<IDirectPosition> getIsolatedPts() {
			return isolatedPts;
		}

		public void setIsolatedPts(Set<IDirectPosition> isolatedPts) {
			this.isolatedPts = isolatedPts;
		}

	}

	/**
	 * Raw Straight Skeleton algorithm using the campskeleton project
	 * implementation. Returns all the skeleton segments that do not touch the
	 * polygon outline. Be careful, the computation time may be extremely long with
	 * very large polygons. Does not work with holes.
	 *
	 * @param polygon
	 * @return
	 */
	public static Set<ILineSegment> skeletonizeStraightSkeleton(IPolygon polygon) {
		Set<ILineSegment> skeletonSegments = new HashSet<ILineSegment>();
		var directionMachine = new Machine();

		// when the geometry is too big, it needs to be simplified first
		var geom = polygon;
		if (polygon.numPoints() > 500) {
			geom = (IPolygon) Filtering.DouglasPeucker(polygon, 15.0);
		}
		if (polygon.numPoints() > 1000) {
			geom = (IPolygon) Filtering.DouglasPeucker(polygon, 30.0);
		}

		var p = geom; // original impl has .reverse() here
		var input = new LoopL<Edge>();

		var rExt = p.getExterior();

		var loop = new Loop<Edge>();
		var lEExt = fromDPLToEdges(rExt.coord());
		for (Edge e : lEExt) {
			loop.append(e);
		}
		for (Edge e : lEExt) {
			e.machine = directionMachine;
		}

		input.add(loop);

		for (IRing rInt : p.getInterior()) {

			var loopIn = new Loop<Edge>();
			input.add(loopIn);
			var lInt = fromDPLToEdges(rInt.coord());
			for (Edge e : lInt) {
				loop.append(e);
			}

			for (Edge e : lInt) {
				e.machine = directionMachine;
			}
		}

		var ske = new Skeleton(input, true);
		ske.skeleton();

		for (SharedEdge edge : ske.output.edges.map.keySet()) {
			ILineSegment segment = new GM_LineSegment(
					new DirectPosition(edge.getStart(edge.left).x, edge.getStart(edge.left).y),
					new DirectPosition(edge.getEnd(edge.left).x, edge.getEnd(edge.left).y));
			if (segment.intersects(polygon.getExterior().getPrimitive())) {
				continue; // TODO optional
			}
			for (IRing hole : polygon.getInterior()) {
				if (segment.intersects(hole.getPrimitive())) {
					continue;
				}
			}
			if (polygon.disjoint(segment)) {
				continue;
			}
			skeletonSegments.add(segment);
		}

		return skeletonSegments;
	}

	/**
	 * Convertit un positon en corner
	 *
	 * @param dp
	 * @return
	 */
	private static Corner fromPositionToCorner(IDirectPosition dp) {
		return new Corner(dp.getX(), dp.getY(), 0);
	}

	/**
	 * Convertit une liste de sommets formant un cycle en arr�tes
	 *
	 * @param dpl
	 * @return
	 */
	private static List<Edge> fromDPLToEdges(IDirectPositionList dpl) {

		var nbPoints = dpl.size();
		List<Edge> lEOut = new ArrayList<Edge>();
		List<Corner> lC = new ArrayList<Corner>();

		for (var i = 0; i < nbPoints - 1; i++) {
			lC.add(fromPositionToCorner(dpl.get(i)));
		}

		lC.add(lC.get(0));

		for (var i = 0; i < nbPoints - 1; i++) {
			lEOut.add(new Edge(lC.get(i), lC.get(i + 1)));
		}

		return lEOut;
	}

	/**
	 * Connect a skeleton computed by any method to the nearest edges of the
	 * skeletonized polygon.
	 *
	 * @param skeleton
	 * @param polygon
	 */
	public static Set<ILineString> connectSkeletonToPolygon(Set<ILineSegment> skeleton, IPolygon polygon) {
		Set<ILineString> extendedSkeleton = new HashSet<ILineString>();
		for (ILineSegment skeSeg : skeleton) {
			ILineString line = new GM_LineString(skeSeg.coord());
			// first check if the start node has to be extended
			var extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg)) {
					continue;
				}
				if (other.coord().contains(skeSeg.coord().get(0))) {
					extend = false;
					break;
				}
			}
			// extend it if necessary
			if (extend) {
				var vect = new Vector2D(skeSeg.coord().get(1), skeSeg.coord().get(0));
				var firstPt = CommonAlgorithmsFromCartAGen.projection(skeSeg.coord().get(0),
						polygon.exteriorLineString(), vect);
				line.addControlPoint(0, firstPt);
			}

			extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg)) {
					continue;
				}
				if (other.coord().contains(skeSeg.coord().get(1))) {
					extend = false;
					break;
				}
			}
			// then, extend the last node
			if (extend) {
				var vect = new Vector2D(skeSeg.coord().get(0), skeSeg.coord().get(1));
				var lastPt = CommonAlgorithmsFromCartAGen.projection(skeSeg.coord().get(1),
						polygon.exteriorLineString(), vect);
				line.addControlPoint(lastPt);
			}
			extendedSkeleton.add(line);
		}
		return extendedSkeleton;
	}

	/**
	 * Connect a skeleton computed by any method to the nearest edges of the given
	 * linear network. If no network edge can be found at a skeleton extremity, it's
	 * extended to polygon's outline.
	 *
	 * @param skeleton
	 * @param polygon
	 */
	public static Set<ILineString> connectSkeletonToNetwork(Set<ILineSegment> skeleton, Set<ILineString> network,
			IPolygon polygon) {
		Set<ILineString> extendedSkeleton = new HashSet<ILineString>();
		// first, find the intersection between the network and the polygon
		Set<IPoint> intersections = new HashSet<IPoint>();
		for (ILineString line : network) {
			if (line.intersects(polygon)) {
				if (line.intersection(polygon) instanceof IPoint) {
					intersections.add((IPoint) line.intersection(polygon));
				} else if (line.intersection(polygon) instanceof IMultiPoint) {
					var inter = (IMultiPoint) line.intersection(polygon);
					for (var i = 0; i < inter.getList().size(); i++) {
						intersections.add(inter.get(i));
					}
				}
			}
		}
		// then get the skeleton extremities
		Map<IDirectPosition, ILineString> skeIni = new HashMap<IDirectPosition, ILineString>();
		Map<IDirectPosition, ILineString> skeFin = new HashMap<IDirectPosition, ILineString>();
		// loop on the skeleton edges to find possible extensions
		for (ILineSegment skeSeg : skeleton) {
			ILineString line = new GM_LineString(skeSeg.coord());
			// first check if the start node has to be extended
			var extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg)) {
					continue;
				}
				if (other.coord().contains(skeSeg.coord().get(0))) {
					extend = false;
					break;
				}
			}
			// extend it if necessary
			if (extend) {
				skeIni.put(skeSeg.coord().get(0), line);
			}

			extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg)) {
					continue;
				}
				if (other.coord().contains(skeSeg.coord().get(1))) {
					extend = false;
					break;
				}
			}
			// then, extend the last node
			if (extend) {
				skeFin.put(skeSeg.coord().get(1), line);
			}

			extendedSkeleton.add(line);
		}

		// now, loop on the intersection points to get the nearest skeleton
		// extremity
		for (IPoint pt : intersections) {
			IDirectPosition nearest = null;
			var maxDist = polygon.perimeter();
			for (IDirectPosition ptIni : skeIni.keySet()) {
				if (pt.getPosition().distance2D(ptIni) < maxDist) {
					nearest = ptIni;
					maxDist = pt.getPosition().distance2D(ptIni);
				}
			}
			for (IDirectPosition ptIni : skeFin.keySet()) {
				if (pt.getPosition().distance2D(ptIni) < maxDist) {
					nearest = ptIni;
					maxDist = pt.getPosition().distance2D(ptIni);
				}
			}

			// now modify the line that has to be extended
			if (skeIni.containsKey(nearest)) {
				skeIni.get(nearest).addControlPoint(0, pt.getPosition());
			} else {
				skeFin.get(nearest).addControlPoint(pt.getPosition());
			}
		}

		return extendedSkeleton;
	}

}