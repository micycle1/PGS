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
package micycle.pts;

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

	// https://github.com/IGNF/CartAGen/blob/master/cartagen-core/src/main/java/fr/ign/cogit/cartagen/algorithms/polygon/Skeletonize.java

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
		Machine directionMachine = new Machine();

		// when the geometry is too big, it needs to be simplified first
		IPolygon geom = polygon;
		if (polygon.numPoints() > 500) {
			geom = (IPolygon) Filtering.DouglasPeucker(polygon, 15.0);
		}
		if (polygon.numPoints() > 1000) {
			geom = (IPolygon) Filtering.DouglasPeucker(polygon, 30.0);
		}

		IPolygon p = (IPolygon) geom.reverse();
		LoopL<Edge> input = new LoopL<Edge>();

		IRing rExt = p.getExterior();

		Loop<Edge> loop = new Loop<Edge>();
		List<Edge> lEExt = fromDPLToEdges(rExt.coord());
		for (Edge e : lEExt)
			loop.append(e);
		for (Edge e : lEExt)
			e.machine = directionMachine;

		input.add(loop);

		for (IRing rInt : p.getInterior()) {

			Loop<Edge> loopIn = new Loop<Edge>();
			input.add(loopIn);
			List<Edge> lInt = fromDPLToEdges(rInt.coord());
			for (Edge e : lInt)
				loop.append(e);

			for (Edge e : lInt)
				e.machine = directionMachine;
		}

		Skeleton ske = new Skeleton(input, true);
		ske.skeleton();

		for (SharedEdge edge : ske.output.edges.map.keySet()) {
			ILineSegment segment = new GM_LineSegment(
					new DirectPosition(edge.getStart(edge.left).x, edge.getStart(edge.left).y),
					new DirectPosition(edge.getEnd(edge.left).x, edge.getEnd(edge.left).y));
			if (segment.intersects(polygon.getExterior().getPrimitive()))
				continue;
			for (IRing hole : polygon.getInterior()) {
				if (segment.intersects(hole.getPrimitive()))
					continue;
			}
			if (polygon.disjoint(segment))
				continue;
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
	 * Convertit une liste de sommets formant un cycle en arrêtes
	 * 
	 * @param dpl
	 * @return
	 */
	private static List<Edge> fromDPLToEdges(IDirectPositionList dpl) {

		int nbPoints = dpl.size();
		List<Edge> lEOut = new ArrayList<Edge>();
		List<Corner> lC = new ArrayList<Corner>();

		for (int i = 0; i < nbPoints - 1; i++)
			lC.add(fromPositionToCorner(dpl.get(i)));

		lC.add(lC.get(0));

		for (int i = 0; i < nbPoints - 1; i++)
			lEOut.add(new Edge(lC.get(i), lC.get(i + 1)));

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
			boolean extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg))
					continue;
				if (other.coord().contains(skeSeg.coord().get(0))) {
					extend = false;
					break;
				}
			}
			// extend it if necessary
			if (extend) {
				Vector2D vect = new Vector2D(skeSeg.coord().get(1), skeSeg.coord().get(0));
				IDirectPosition firstPt = CommonAlgorithmsFromCartAGen.projection(skeSeg.coord().get(0),
						polygon.exteriorLineString(), vect);
				line.addControlPoint(0, firstPt);
			}

			extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg))
					continue;
				if (other.coord().contains(skeSeg.coord().get(1))) {
					extend = false;
					break;
				}
			}
			// then, extend the last node
			if (extend) {
				Vector2D vect = new Vector2D(skeSeg.coord().get(0), skeSeg.coord().get(1));
				IDirectPosition lastPt = CommonAlgorithmsFromCartAGen.projection(skeSeg.coord().get(1),
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
				if (line.intersection(polygon) instanceof IPoint)
					intersections.add((IPoint) line.intersection(polygon));
				else if (line.intersection(polygon) instanceof IMultiPoint) {
					IMultiPoint inter = (IMultiPoint) line.intersection(polygon);
					for (int i = 0; i < inter.getList().size(); i++) {
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
			boolean extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg))
					continue;
				if (other.coord().contains(skeSeg.coord().get(0))) {
					extend = false;
					break;
				}
			}
			// extend it if necessary
			if (extend)
				skeIni.put(skeSeg.coord().get(0), line);

			extend = true;
			for (ILineSegment other : skeleton) {
				if (other.equals(skeSeg))
					continue;
				if (other.coord().contains(skeSeg.coord().get(1))) {
					extend = false;
					break;
				}
			}
			// then, extend the last node
			if (extend)
				skeFin.put(skeSeg.coord().get(1), line);

			extendedSkeleton.add(line);
		}

		// now, loop on the intersection points to get the nearest skeleton
		// extremity
		for (IPoint pt : intersections) {
			IDirectPosition nearest = null;
			double maxDist = polygon.perimeter();
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

//	/**
//	 * Skeletonize a polygon and returns all segments of the skeleton as a
//	 * tree-graph.
//	 * 
//	 * @param polygon
//	 * @param densStep densifies polygon with a vertex every densStep. Put -1 to
//	 *                 avoid densification.
//	 * @return
//	 */
//	public static TreeGraph skeletonizeTINGraph(IPolygon polygon, double densStep) {
//		// initialisation of useful collections
//		List<TriangulationPoint> points = new ArrayList<TriangulationPoint>();
//		TreeGraph graph = new TreeGraph("skeleton");
//
//		// compute the outline of the polygon to extract the triangulation
//		// vertices and constraining segments
//		ILineString contour = null;
//		IPolygon densePol = null;
//		if (densStep == -1) {
//			contour = new GM_LineString(polygon.coord());
//			densePol = polygon;
//		} else {
//			contour = LineDensification.densification(new GM_LineString(polygon.coord()), densStep);
//			densePol = GeometryEngine.getFactory().createIPolygon(contour.coord());
//		}
//
//		// the list contains all the points arround the face
//		for (int i = 0; i < contour.numPoints(); i++) {
//			TriangulationPointImpl point = new TriangulationPointImpl(contour.getControlPoint(i));
//			points.add(point);
//		}
//
//		// the triangulation is generated
//		Triangulation tri = new Triangulation(points, new TriangulationSegmentFactoryImpl(),
//				new TriangulationTriangleFactoryImpl());
//		tri.compute(true);
//		Collection<TriangulationTriangle> triangles = tri.getTriangles();
//
//		if (triangles.size() != 2) {
//			// There are three cases of triangles
//			for (TriangulationTriangle triangle : triangles) {
//
//				// System.out.println(triangle.getGeom().centroid());
//				// System.out.println(polygon.relate(triangle.getGeom()));
//				if (JTSAlgorithms.coversPredicate(densePol, triangle.getGeom())) {
//
//					TriangulationPointImpl point1 = new TriangulationPointImpl(triangle.getPoint1().getPosition());
//					TriangulationPointImpl point2 = new TriangulationPointImpl(triangle.getPoint2().getPosition());
//					TriangulationPointImpl point3 = new TriangulationPointImpl(triangle.getPoint3().getPosition());
//					TriangulationSegmentImpl seg1 = new TriangulationSegmentImpl(point1, point2);
//					TriangulationSegmentImpl seg2 = new TriangulationSegmentImpl(point2, point3);
//					TriangulationSegmentImpl seg3 = new TriangulationSegmentImpl(point3, point1);
//					int nbSeg = 0, contourSeg1 = 0, contourSeg2 = 0, contourSeg3 = 0;
//
//					if (contour.buffer(0.1).contains(seg1.getGeom()) == true) {
//						contourSeg1 = 1;
//					}
//					if (contour.buffer(0.1).contains(seg2.getGeom()) == true) {
//						contourSeg2 = 1;
//					}
//					if (contour.buffer(0.1).contains(seg3.getGeom()) == true) {
//						contourSeg3 = 1;
//					}
//
//					nbSeg = contourSeg1 + contourSeg2 + contourSeg3;
//
//					// The end triangles (the end line starts and finishes where these
//					// triangles are) => We don't make anything now
//
//					// The most frequent triangles (1 side in common with the separator
//					// outline)
//					if (nbSeg == 1) {
//						IDirectPosition pos1 = new DirectPosition();
//						IDirectPosition pos2 = new DirectPosition();
//
//						if (contourSeg1 == 1) {
//							pos1 = seg2.getGeom().centroid();
//							pos2 = seg3.getGeom().centroid();
//							INode node1 = graph.getNodeAt(pos1);
//							if (node1 == null) {
//								node1 = new Node(pos1.toGM_Point());
//								node1.setGraph(graph);
//								graph.getNodes().add(node1);
//							}
//							INode node2 = graph.getNodeAt(pos2);
//							if (node2 == null) {
//								node2 = new Node(pos2.toGM_Point());
//								node2.setGraph(graph);
//								graph.getNodes().add(node2);
//							}
//							graph.addEdge(node1, node2);
//						} else if (contourSeg2 == 1) {
//							pos1 = seg1.getGeom().centroid();
//							pos2 = seg3.getGeom().centroid();
//							INode node1 = graph.getNodeAt(pos1);
//							if (node1 == null) {
//								node1 = new Node(pos1.toGM_Point());
//								node1.setGraph(graph);
//								graph.getNodes().add(node1);
//							}
//							INode node2 = graph.getNodeAt(pos2);
//							if (node2 == null) {
//								node2 = new Node(pos2.toGM_Point());
//								node2.setGraph(graph);
//								graph.getNodes().add(node2);
//							}
//							graph.addEdge(node1, node2);
//						} else {
//							pos1 = seg1.getGeom().centroid();
//							pos2 = seg2.getGeom().centroid();
//							INode node1 = graph.getNodeAt(pos1);
//							if (node1 == null) {
//								node1 = new Node(pos1.toGM_Point());
//								node1.setGraph(graph);
//								graph.getNodes().add(node1);
//							}
//							INode node2 = graph.getNodeAt(pos2);
//							if (node2 == null) {
//								node2 = new Node(pos2.toGM_Point());
//								node2.setGraph(graph);
//								graph.getNodes().add(node2);
//							}
//							graph.addEdge(node1, node2);
//						}
//					}
//
//					// Third case : a fork is created (no side in common with the
//					// separator outline)
//					if (nbSeg == 0) {
//						IDirectPosition pos1 = new DirectPosition();
//						IDirectPosition pos2 = new DirectPosition();
//						IDirectPosition pos3 = new DirectPosition();
//						IDirectPosition barycentre = new DirectPosition();
//						pos1 = seg1.getGeom().centroid();
//						pos2 = seg2.getGeom().centroid();
//						pos3 = seg3.getGeom().centroid();
//						barycentre = triangle.getGeom().centroid();
//
//						Node node = new Node(barycentre.toGM_Point());
//						node.setGraph(graph);
//						graph.getNodes().add(node);
//
//						INode node1 = graph.getNodeAt(pos1);
//						if (node1 == null) {
//							node1 = new Node(pos1.toGM_Point());
//							node1.setGraph(graph);
//							graph.getNodes().add(node1);
//						}
//						graph.addEdge(node, node1);
//
//						INode node2 = graph.getNodeAt(pos2);
//						if (node2 == null) {
//							node2 = new Node(pos2.toGM_Point());
//							node2.setGraph(graph);
//							graph.getNodes().add(node2);
//						}
//						graph.addEdge(node, node2);
//
//						INode node3 = graph.getNodeAt(pos3);
//						if (node3 == null) {
//							node3 = new Node(pos3.toGM_Point());
//							node3.setGraph(graph);
//							graph.getNodes().add(node3);
//						}
//						graph.addEdge(node, node3);
//					}
//
//				}
//			}
//		}
//		return graph;
//	}
//
//	/**
//	 * @param polygon
//	 * @param multiSke
//	 * @param densStep
//	 */
//	private static void computeTINSkeSegments(IPolygon polygon, MultiSkeleton multiSke, double densStep) {
//		// initialisation of useful collections
//		List<TriangulationPoint> points = new ArrayList<TriangulationPoint>();
//		Set<ILineSegment> hashSegFork = new HashSet<ILineSegment>();
//
//		// compute the outline of the polygon to extract the triangulation
//		// vertices and constraining segments
//		ILineString contour = LineDensification.densification(new GM_LineString(polygon.coord()), densStep);
//
//		// the list contains all the points arround the face
//		for (int i = 0; i < contour.numPoints(); i++) {
//			TriangulationPointImpl point = new TriangulationPointImpl(contour.getControlPoint(i));
//			points.add(point);
//		}
//
//		// the triangulation is generated
//		Triangulation tri = new Triangulation(points, new TriangulationSegmentFactoryImpl(),
//				new TriangulationTriangleFactoryImpl());
//		tri.compute(true);
//		Collection<TriangulationTriangle> triangles = tri.getTriangles();
//
//		// Cas particulier où la triangulation de la face contient deux triangles
//		if (triangles.size() == 2) {
//			for (TriangulationSegment seg : tri.getSegments()) {
//				if ((polygon.contains(seg.getGeom()))
//						&& (contour.buffer(0.05).intersection(seg.getGeom()).equals(seg.getGeom()) == false)) {
//					multiSke.getIsolatedPts().add(seg.getGeom().centroid());
//				}
//			}
//		}
//
//		else {
//			// There are three case of triangles
//			for (TriangulationTriangle triangle : triangles) {
//
//				if (polygon.contains(triangle.getGeom())) {
//
//					TriangulationPointImpl point1 = new TriangulationPointImpl(triangle.getPoint1().getPosition());
//					TriangulationPointImpl point2 = new TriangulationPointImpl(triangle.getPoint2().getPosition());
//					TriangulationPointImpl point3 = new TriangulationPointImpl(triangle.getPoint3().getPosition());
//					TriangulationSegmentImpl seg1 = new TriangulationSegmentImpl(point1, point2);
//					TriangulationSegmentImpl seg2 = new TriangulationSegmentImpl(point2, point3);
//					TriangulationSegmentImpl seg3 = new TriangulationSegmentImpl(point3, point1);
//					int nbSeg = 0, contourSeg1 = 0, contourSeg2 = 0, contourSeg3 = 0;
//
//					if (contour.buffer(0.1).contains(seg1.getGeom()) == true) {
//						contourSeg1 = 1;
//					}
//					if (contour.buffer(0.1).contains(seg2.getGeom()) == true) {
//						contourSeg2 = 1;
//					}
//					if (contour.buffer(0.1).contains(seg3.getGeom()) == true) {
//						contourSeg3 = 1;
//					}
//
//					nbSeg = contourSeg1 + contourSeg2 + contourSeg3;
//
//					// The end triangles (the end line starts and finishes where these
//					// triangles are => We don't make anything now
//
//					// The most frequent triangles (1 side in common with the separator
//					// outline)
//					if (nbSeg == 1) {
//						IDirectPosition pos1 = new DirectPosition();
//						IDirectPosition pos2 = new DirectPosition();
//
//						if (contourSeg1 == 1) {
//							pos1 = seg2.getGeom().centroid();
//							pos2 = seg3.getGeom().centroid();
//							ILineSegment segment = new GM_LineSegment(pos1, pos2);
//							multiSke.segFace.add(segment);
//						} else if (contourSeg2 == 1) {
//							pos1 = seg1.getGeom().centroid();
//							pos2 = seg3.getGeom().centroid();
//							ILineSegment segment = new GM_LineSegment(pos1, pos2);
//							multiSke.segFace.add(segment);
//						} else {
//							pos1 = seg1.getGeom().centroid();
//							pos2 = seg2.getGeom().centroid();
//							ILineSegment segment = new GM_LineSegment(pos1, pos2);
//							multiSke.segFace.add(segment);
//						}
//					}
//
//					// Third case : a fork is created (no side in common with the
//					// separator outline)
//					if (nbSeg == 0) {
//						IDirectPosition pos1 = new DirectPosition();
//						IDirectPosition pos2 = new DirectPosition();
//						IDirectPosition pos3 = new DirectPosition();
//						IDirectPosition barycentre = new DirectPosition();
//						pos1 = seg1.getGeom().centroid();
//						pos2 = seg2.getGeom().centroid();
//						pos3 = seg3.getGeom().centroid();
//						barycentre = triangle.getGeom().centroid();
//
//						ILineSegment segSuivant1 = new GM_LineSegment(barycentre, pos1);
//						ILineSegment segSuivant2 = new GM_LineSegment(barycentre, pos2);
//						ILineSegment segSuivant3 = new GM_LineSegment(barycentre, pos3);
//						hashSegFork.add(segSuivant1);
//						hashSegFork.add(segSuivant2);
//						hashSegFork.add(segSuivant3);
//					}
//
//				}
//			}
//		}
//
//		// Removal of the useless segments created with the forks by filtering
//		for (ILineSegment segmentFork : hashSegFork) {
//			IDirectPosition debutFork = segmentFork.getStartPoint();
//			IDirectPosition finFork = segmentFork.getEndPoint();
//			boolean boolFork = false;
//
//			for (ILineSegment segment : multiSke.segFace) {
//				IDirectPosition debut = segment.getStartPoint();
//				IDirectPosition fin = segment.getEndPoint();
//				if (debutFork.equals(debut) || debutFork.equals(fin) || finFork.equals(debut) || finFork.equals(fin)) {
//					boolFork = true;
//				}
//			}
//
//			if (boolFork == true) {
//				multiSke.segments.add(segmentFork);
//			}
//		}
//
//		multiSke.segments.addAll(multiSke.segFace);
//		multiSke.segFace.clear();
//	}

	/**
	 * Connect a skeleton computed by any method to the nearest edges of the
	 * skeletonized polygon.
	 * 
	 * @param skeleton
	 * @param polygon
	 */
	/**
	 * @param skeleton
	 * @param polygon
	 */
	private static void connectSkeletonToPolygon(ILineString skeleton, IPolygon polygon) {
		// first extend the start node
		Vector2D vect = new Vector2D(skeleton.coord().get(1), skeleton.coord().get(0));
		IDirectPosition firstPt = CommonAlgorithmsFromCartAGen.projection(skeleton.coord().get(0),
				polygon.exteriorLineString(), vect);
		skeleton.addControlPoint(0, firstPt);

		// then, extend the last node
		vect = new Vector2D(skeleton.coord().get(skeleton.numPoints() - 2),
				skeleton.coord().get(skeleton.numPoints() - 1));
		IDirectPosition lastPt = CommonAlgorithmsFromCartAGen.projection(skeleton.coord().get(skeleton.numPoints() - 1),
				polygon.exteriorLineString(), vect);
		skeleton.addControlPoint(lastPt);
	}

}