/** This library provides method to compute Concave Hull of a set of points
 * 
 * Currently three criteria are provided:
 * 
 * Alpha shape in Edelsbrunner, Herbert; Kirkpatrick, David G.; Seidel, Raimund (1983), "On the shape of a set of points in the plane", IEEE Transactions on Information Theory, 29 (4): 551�559, doi:10.1109/TIT.1983.1056714 
 * 
 * Chi criterion in Matt Duckham et al 2008 Efficient generation of simple polygons for characterizing the shape of a set of points in the plane
 * 
 * Park edge ratio criterion in JIN-SEO PARK AND SE-JONG OH (2012) A New Concave Hull Algorithm and Concaveness Measure for n-dimensional Datasets, JOURNAL OF INFORMATION SCIENCE AND ENGINEERING 28, 587-600 
 * 
 * Alpha shape is effectively to use a disc of radius R = 1/alpha to remove space among points without enclosing any point in disc's interior
 * Therefore, if a triangle has a circumcircle with a radius large than the given R, the triangle may be removed from the initial hull. For an infinite R (alpha = 0), it becomes the convex hull of the point set 
 * 
 * Chi-criterion is the edge length threshold. If an edge of a triangle is longer than the threshold, the triangle (in fact, often two) may be removed from the initial hull.
 * 
 * Park edge ratio is the ratio between the length of the "outter" edge and the shorter inner edge. Unlike the above two criteria, it is scale-independent.    
 *
 * In this library we provide several different implementations.
 * 
 * Criteria for triangle removal is parametrised and implemented as the TriangleChecker interface.
 * 
 * getConcaveHullDFS supports all three criteria. It starts the "digging" from the first qualified boundary edge and follows a depth-first strategy. The "dig" stops only if no more dig-able edges left inside a "cave", then it will start digging again with another qualified boundary edge. 
 * 
 * getConcaveHullBFS support all three criteria with an option of whether allow multiple parts to be generated. It follows a breadth-first strategy and at each step it will "dig" the longest qualified edge.    
 * 
 * getConcaveHullWithHolesAlpha supports alpha shape criterion and will multiple parts as well as generate holes if applicable.
 * 
 * getConcaveHullWithHolesChi supports Chi criterion and also generate multiple parts and holes if applicable. 
 * 
 * 
 * This library is built entirely on top of JTS geometry and Delaunay triangulation libraries. Therefore, it will have the same numerical robustness issues as has JTS.  
 *  
 * 
 * Author: Sheng Zhou (Sheng.Zhou@os.uk)
 * 
 * version 0.5
 * 
 * Date: 2019-04-24
 * 
 * Copyright (C) 2019 Ordnance Survey
 *
 * Licensed under the Open Government Licence v3.0 (the "License");
 * 
 * you may not use this file except in compliance with the License.
 * 
 * You may obtain a copy of the License at
 *
 *     http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
*/

package uk.osgb.algorithm.concavehull;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.quadedge.QuadEdge;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;
import org.locationtech.jts.triangulate.quadedge.Vertex;

import org.locationtech.jts.operation.union.UnaryUnionOp;

public class ConcaveHull {

	/**
	 * A coordinate on the convex hull of the dataset, used as a seed coordinate for
	 * various operations
	 */
	Coordinate hullCoord = null;
	/**
	 * Initial hull of the dataset (closed, first==last), built from JTS DT. It is
	 * normally (but not always) the convex hull (due to JTS DT's use of a bounding
	 * super triangle)
	 */
	LinkedList<Coordinate> hullDT = null;
	QuadEdgeSubdivision sd = null;
	GeometryFactory gf = null;
	Geometry convexHull = null;

	/**
	 * constructor that takes a JTS geometry as input
	 * 
	 * @param geom
	 */
	public ConcaveHull(Geometry geom) {
		Coordinate[] coordArray = geom.getCoordinates();
		initArray(coordArray, geom.getFactory());
	}

	/**
	 * constructor that takes a JTS geometry as input, with optional densification
	 * 
	 * @param geom       input geometry
	 * @param densifyTol tolerance for densification of linear geometry or boundary
	 *                   of areal geometry. 0.0 for none
	 */
	public ConcaveHull(Geometry geom, double densifyTol) {
		if (densifyTol > 0.0) {
			geom = Densifier.densify(geom, densifyTol);
		}
		Coordinate[] coordArray = geom.getCoordinates();
		initArray(coordArray, geom.getFactory());
	}

	/**
	 * constructor that takes a coordinate collection as input
	 * 
	 * @param coordCol input JTS Coordinate collection
	 * @param gf       GeometryFactory to be used for geometry construction
	 */
	public ConcaveHull(Collection<Coordinate> coordCol, GeometryFactory gf) {
		initCol(coordCol, gf);
	}

	/**
	 * constructor that takes a coordinate array as input
	 * 
	 * @param coordArray
	 * @param gf
	 */
	public ConcaveHull(Coordinate[] coordArray, GeometryFactory gf) {
		initArray(coordArray, gf);
	}

	/**
	 * constructor that takes a Geometry collection as input, with densification
	 * 
	 * @param geomCol
	 * @param densifyTol tolerance for densification of linear geometry or boundary
	 *                   of areal geometry. 0.0 for none
	 */
	public ConcaveHull(Collection<Geometry> geomCol, double densifyTol) {
		if (geomCol != null && !geomCol.isEmpty()) {
			Collection<Geometry> gc = null;
			if (densifyTol > 0.0) {
				gc = new ArrayList<>(geomCol.size());
				for (Geometry geom : geomCol) {
					gc.add(Densifier.densify(geom, densifyTol));
				}
			} else {
				gc = geomCol;
			}
			Collection<Coordinate> coordCol = geomCol2Coordinate(gc);
			initCol(coordCol, geomCol.iterator().next().getFactory());
		}
	}

	/**
	 * dig depth-first, faster but may generate weird results
	 * 
	 * @param triChecker Triangle Checker to be used
	 * @return concave hull of input data
	 */
	public Geometry getConcaveHullDFS(TriangleChecker triChecker) {
		if (triChecker != null) {
			TreeSet<Coordinate> nodeSet = new TreeSet<>(hullDT);
			LinkedList<Coordinate> hullList = new LinkedList<>(hullDT);
			boolean modified = true;
			while (modified) {
				modified = false;
				ListIterator<Coordinate> iter = hullList.listIterator(0);
				while (iter.hasNext()) {
					Coordinate nodeS = iter.next(); // iter after nodeS
					if (iter.hasNext()) {
						Coordinate nodeE = iter.next(); // iter afer nodeE
						QuadEdge qe = sd.locate(nodeS, nodeE);
						if (qe == null) {
							System.out.println("fail to find edge: " + nodeS.toString() + " -" + nodeE.toString());
							return null;
						}
						QuadEdge lnext = qe.lNext();
						Vertex verO = lnext.dest();
						Coordinate nodeO = verO.getCoordinate();
						if (nodeSet.contains(nodeO)) {
							iter.previous();
							continue; // will break areal connectivity so omit?
						}
						if (triChecker.removeable(nodeS, nodeE, nodeO)) {
							modified = true;
							nodeSet.add(nodeO);
							iter.previous(); // move iter before nodeE
							iter.add(nodeO); // insert nodeO between nodeS and nodeE (and before iter)
							iter.previous(); // iter before nodeO
							iter.previous(); // iter before nodeS
						} else {
							iter.previous(); // iter before nodeE
						}
					}
				}
			}
			//
			Coordinate[] newCoords = new Coordinate[hullList.size()];
			hullList.toArray(newCoords);
			return gf.createPolygon(newCoords);
		} else {
			return null;
		}
	}

	/**
	 * breadth first digging, supports multi parts
	 * 
	 * @param triChecker      triangle checker to be used
	 * @param allowMultiParts if multiple parts are to be generated
	 * @param keepLineSeg     if degenerated line segments should be kept
	 * @return a collection of geometry that form the concave hull of the input data
	 *         (may contains linestring as degenerated segments)
	 */
	public List<Geometry> getConcaveHullBFS(TriangleChecker triChecker, boolean allowMultiParts, boolean keepLineSeg) {
		if (triChecker != null) {
			LinkedList<DLCirList<Coordinate>> hulls = new LinkedList<>(); //
			ArrayList<DLCirList<Coordinate>> rltHulls = new ArrayList<>(); // finished hulls

			ArrayList<LineString> rltLS = new ArrayList<>();
			ArrayList<Point> rltPt = new ArrayList<>();

			TreeSet<HullEdgeCir> edgeIdx = new TreeSet<>();

			TreeMap<Coordinate, DLNode<Coordinate>> coordNodeMap = new TreeMap<>();
			HashMap<DLNode<Coordinate>, HullEdgeCir> nodeEdgeMap = new HashMap<>();

			DLCirList<Coordinate> hull = generateHullEdgeRep(hullDT);

			hulls.add(hull);
			while (!hulls.isEmpty()) {
				hull = hulls.pollFirst();
				generateHullIndices(hull, coordNodeMap, edgeIdx, nodeEdgeMap);
				boolean addHull = true;
				while (!edgeIdx.isEmpty()) {
					HullEdgeCir edge = edgeIdx.pollLast(); // the longest

					DLNode<Coordinate> sn = edge.sn;
					DLNode<Coordinate> en = sn.getNext();
					Coordinate nodeS = sn.getObj();
					Coordinate nodeE = en.getObj();
					QuadEdge qe = sd.locate(nodeS, nodeE);
					if (qe == null) {
						System.out.println("fail to find edge: " + nodeS.toString() + " -" + nodeE.toString());
						return new ArrayList<Geometry>();
					}
					QuadEdge lnext = qe.lNext();
					Vertex verO = lnext.dest(); // "opposite" vertex
					Coordinate nodeO = verO.getCoordinate();
					if (triChecker.removeable(nodeS, nodeE, nodeO)) {
						if (hull.size() > 3) {
							if (coordNodeMap.containsKey(nodeO)) { // nodeO is a boundary node, split takes place
								if (allowMultiParts) {
									DLNode<Coordinate> on = coordNodeMap.get(nodeO);
									// check if a tri-corner
									if (on == en.getNext()) {// trim en
										HullEdgeCir ee = nodeEdgeMap.get(en);
										edgeIdx.remove(ee); // se popped out already
										coordNodeMap.remove(nodeE);
										nodeEdgeMap.remove(en);
										nodeEdgeMap.remove(sn);
										hull.remove(en); // sn now connected to on
										//
										HullEdgeCir seNew = new HullEdgeCir(sn);
										nodeEdgeMap.put(sn, seNew);
										edgeIdx.add(seNew);
										if (keepLineSeg) {
											Coordinate[] lsCoords = new Coordinate[2];
											lsCoords[0] = nodeE;
											lsCoords[1] = nodeO;
											rltLS.add(gf.createLineString(lsCoords));
										}
										rltPt.add(gf.createPoint(en.getObj()));
									} else if (on == sn.getPrev()) { // trim sn
										HullEdgeCir oe = nodeEdgeMap.get(on);
										edgeIdx.remove(oe); // se popped out already
										coordNodeMap.remove(nodeS);
										nodeEdgeMap.remove(sn);
										nodeEdgeMap.remove(on);
										hull.remove(sn); // on now connected to en
										//
										HullEdgeCir oeNew = new HullEdgeCir(on);
										nodeEdgeMap.put(on, oeNew);
										edgeIdx.add(oeNew);
										if (keepLineSeg) {
											Coordinate[] lsCoords = new Coordinate[2];
											lsCoords[0] = nodeO;
											lsCoords[1] = nodeS;
											rltLS.add(gf.createLineString(lsCoords));
										}
										rltPt.add(gf.createPoint(sn.getObj()));
									} else {// split
										addHull = false;
										// split
										DLCirList<Coordinate> hull2 = hull.split(sn, en, on);
										hulls.add(hull);
										hulls.add(hull2);
										//
										break;
									}
								} else {// multiparts not allowed
										// don't do anything
								}
							} else { // normal digging
								nodeEdgeMap.remove(sn);
								DLNode<Coordinate> on = new DLNode<>(nodeO);
								coordNodeMap.put(nodeO, on);
								hull.addAfter(sn, on);
								//
								HullEdgeCir seNew = new HullEdgeCir(sn);
								HullEdgeCir oe = new HullEdgeCir(on);
								edgeIdx.add(seNew);
								edgeIdx.add(oe);
								nodeEdgeMap.put(sn, seNew);
								nodeEdgeMap.put(on, oe);
							}
						} else {// 3 vertices only, clapse to a line segement, may be saved separately if
								// needed?
							addHull = false;
							break;
						}
					}
				}
				
				if (addHull) {
					// add hull to result
					rltHulls.add(hull);
				}
				// clear indices
				coordNodeMap.clear();
				nodeEdgeMap.clear();
				edgeIdx.clear();
			}
			ArrayList<Geometry> rtn = new ArrayList<>(hulls.size() + rltLS.size() + rltPt.size());
			for (DLCirList<Coordinate> h : rltHulls) {
				Coordinate[] coords = new Coordinate[h.size() + 1];
				int cnt = 0;
				DLNode<Coordinate> node = h.getNode();
				do {
					Coordinate coord = node.getObj();
					coords[cnt++] = coord;
					node = node.getNext();
				} while (node != h.getNode());
				coords[cnt] = new Coordinate(coords[0]);
				Geometry geom = gf.createPolygon(coords);
				rtn.add(geom);
			}
			rtn.addAll(rltLS);
			rtn.addAll(rltPt);

			return rtn;
		} else {
			return new ArrayList<>();
		}
	}

	/**
	 * Experimental Concave hull construction with an extra metric to control the
	 * order of boundary edge "digging"
	 * 
	 * @param triChecker
	 * @param m               TriMetricLength for controlling the order of digging
	 *                        (default is the length of the boundary edge)
	 * @param allowMultiParts
	 * @param keepLineSeg
	 * @return
	 */
	public Collection<Geometry> getConcaveHullMetric(TriangleChecker triChecker, TriMetricLength m, boolean allowMultiParts,
			boolean keepLineSeg) {
		if (triChecker != null) {
			LinkedList<DLCirList<Coordinate>> hulls = new LinkedList<>(); //
			ArrayList<DLCirList<Coordinate>> rltHulls = new ArrayList<>(); // finished hulls

			ArrayList<LineString> rltLS = new ArrayList<LineString>();

			TreeSet<HullEdgeCir> edgeIdx = new TreeSet<HullEdgeCir>();

			TreeMap<Coordinate, DLNode<Coordinate>> coordNodeMap = new TreeMap<>();
			HashMap<DLNode<Coordinate>, HullEdgeCir> nodeEdgeMap = new HashMap<>();

			DLCirList<Coordinate> hull = generateHullEdgeRep(hullDT);

			hulls.add(hull);
			while (!hulls.isEmpty()) {
				hull = hulls.pollFirst();
				generateHullIndices(hull, m, coordNodeMap, edgeIdx, nodeEdgeMap);
				boolean addHull = true;
				while (!edgeIdx.isEmpty()) {
					HullEdgeCir edge = edgeIdx.pollLast(); // the longest

					DLNode<Coordinate> sn = edge.sn;
					DLNode<Coordinate> en = sn.getNext();
					Coordinate nodeS = sn.getObj();
					Coordinate nodeE = en.getObj();
					QuadEdge qe = sd.locate(nodeS, nodeE);
					if (qe == null) {
						System.out.println("fail to find edge: " + nodeS.toString() + " -" + nodeE.toString());
						return null;
					}
					QuadEdge lnext = qe.lNext();
					Vertex verO = lnext.dest();
					Coordinate nodeO = verO.getCoordinate();
					if (triChecker.removeable(nodeS, nodeE, nodeO)) {
						if (hull.size() > 3) {
							if (coordNodeMap.containsKey(nodeO)) { // nodeO is a boundary node, split takes place
								DLNode<Coordinate> on = coordNodeMap.get(nodeO);
								// check if a tri-corner
								if (on == en.getNext()) {// trim en
									HullEdgeCir ee = nodeEdgeMap.get(en);
									edgeIdx.remove(ee); // se popped out already
									coordNodeMap.remove(nodeE);
									nodeEdgeMap.remove(en);
									nodeEdgeMap.remove(sn);
									hull.remove(en); // sn now connected to on
									//
									HullEdgeCir seNew = new HullEdgeCir(sn, m);
									nodeEdgeMap.put(sn, seNew);
									edgeIdx.add(seNew);
									if (keepLineSeg) {
										Coordinate[] lsCoords = new Coordinate[2];
										lsCoords[0] = nodeE;
										lsCoords[1] = nodeO;
										rltLS.add(gf.createLineString(lsCoords));
									}
								} else if (on == sn.getPrev()) { // trim sn
									HullEdgeCir oe = nodeEdgeMap.get(on);
									edgeIdx.remove(oe); // se popped out already
									coordNodeMap.remove(nodeS);
									nodeEdgeMap.remove(sn);
									nodeEdgeMap.remove(on);
									hull.remove(sn); // on now connected to en
									//
									HullEdgeCir oeNew = new HullEdgeCir(on, m);
									nodeEdgeMap.put(on, oeNew);
									edgeIdx.add(oeNew);
									if (keepLineSeg) {
										Coordinate[] lsCoords = new Coordinate[2];
										lsCoords[0] = nodeO;
										lsCoords[1] = nodeS;
										rltLS.add(gf.createLineString(lsCoords));
									}
								} else {// split
									if (allowMultiParts) {
										addHull = false;
										// split
										DLCirList<Coordinate> hull2 = hull.split(sn, en, on);
										hulls.add(hull);
										hulls.add(hull2);
										//
										break;
									}
								}
								//
							} else { // normal digging
								//
								nodeEdgeMap.remove(sn);
								DLNode<Coordinate> on = new DLNode<>(nodeO);
								coordNodeMap.put(nodeO, on);
								hull.addAfter(sn, on);
								//
								HullEdgeCir seNew = new HullEdgeCir(sn, m);
								HullEdgeCir oe = new HullEdgeCir(on, m);
								edgeIdx.add(seNew);
								edgeIdx.add(oe);
								nodeEdgeMap.put(sn, seNew);
								nodeEdgeMap.put(on, oe);
							}
						} else {// 3 vertices only, clapse to a line segement, may be saved separately if
								// needed?
							addHull = false;
							break;
						}
					}
				}
				if (addHull) {
					// add hull to result
					rltHulls.add(hull);
				}
				// clear indices
				coordNodeMap.clear();
				nodeEdgeMap.clear();
				edgeIdx.clear();
			}
			ArrayList<Geometry> rtn = new ArrayList<>(hulls.size() + rltLS.size());
			for (DLCirList<Coordinate> h : rltHulls) {
				Coordinate[] coords = new Coordinate[h.size() + 1];
				int cnt = 0;
				DLNode<Coordinate> node = h.getNode();
				do {
					Coordinate coord = node.getObj();
					coords[cnt++] = coord;
					node = node.getNext();
				} while (node != h.getNode());
				coords[cnt] = new Coordinate(coords[0]);
				Geometry geom = gf.createPolygon(coords);
				rtn.add(geom);
			}
			rtn.addAll(rltLS);
			return rtn;
		} else {
			return null;
		}
	}

	/**
	 * Concave hull with multiple parts and holes (going through all triangles and
	 * use union operators to generate results)
	 * 
	 * @param tcAlpha triangle checker using Alpha shape criterion
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public Geometry getConcaveHullWithHolesAlpha(TriCheckerAlpha tcAlpha) {
		Coordinate[] hullCoords = new Coordinate[hullDT.size()];
		hullDT.toArray(hullCoords);
		Geometry shell = gf.createPolygon(hullCoords);

		List<Coordinate[]> triCoords = sd.getTriangleCoordinates(false);
		for (Coordinate[] coords : triCoords) {
			boolean rmvable = tcAlpha.removeable(coords[0], coords[1], coords[2]);
			if (rmvable) {
				Geometry hole = gf.createPolygon(coords);
				shell = shell.difference(hole);
			}
		}
		return shell;
	}

	//
	/**
	 * Concave hull with multiple parts and holes (going through all triangles and
	 * use union operators to generate results)
	 * 
	 * @param tcChi triangle checker using Chi criterion (edge length)
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public Geometry getConcaveHullWithHolesChi(TriCheckerChi tcChi) {
		Coordinate[] hullCoords = new Coordinate[hullDT.size()];
		hullDT.toArray(hullCoords);
		Geometry shell = gf.createPolygon(hullCoords);

		List<Coordinate[]> triCoords = sd.getTriangleCoordinates(false);
		Collection<Geometry> holeCol = new ArrayList<>();
		for (Coordinate[] coords : triCoords) {
			boolean rmvable0 = tcChi.removeable(coords[0], coords[1], coords[2]);
			boolean rmvable1 = tcChi.removeable(coords[1], coords[2], coords[0]);
			boolean rmvable2 = tcChi.removeable(coords[2], coords[0], coords[1]);
			if (rmvable0 || rmvable1 || rmvable2) {
				Geometry hole = gf.createPolygon(coords);
				holeCol.add(hole);
			}
		}
		try {
			Geometry holes = UnaryUnionOp.union(holeCol);
			if (holes != null && holes.getDimension() > 1) {
				shell = shell.difference(holes);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return shell;
	}

	public Geometry exportTriangles(GeometryFactory gf) {
		if (sd != null) {
			return sd.getTriangles(gf);
		} else {
			return null;
		}
	}

	//
	private void initArray(Coordinate[] coordArray, GeometryFactory gf) {
		this.gf = gf;
		init(Arrays.asList(coordArray), coordArray);
	}

	private void initCol(Collection<Coordinate> coordCol, GeometryFactory gf) {
		this.gf = gf;
		Coordinate[] coordArray = new Coordinate[coordCol.size()];
		coordCol.toArray(coordArray);
		init(coordCol, coordArray);
	}

	private void init(Collection<Coordinate> coordCol, Coordinate[] coordArray) {
		ConvexHull ch = new ConvexHull(coordArray, gf);
		Geometry chGeom = ch.getConvexHull();
		convexHull = chGeom;
		coordArray = null;
		hullCoord = chGeom.getCoordinate();
		DelaunayTriangulationBuilder builder = new DelaunayTriangulationBuilder();
		builder.setSites(coordCol);
		sd = builder.getSubdivision();
		hullDT = getHull(sd, hullCoord);
	}

	/**
	 * Flattens the coordinates of geometries in a geometry collection.
	 */
	private static Collection<Coordinate> geomCol2Coordinate(Collection<Geometry> geomCol) {
		ArrayList<Coordinate> coordCol = new ArrayList<>();
		if (geomCol != null) {
			for (Geometry geom : geomCol) {
				coordCol.addAll(Arrays.asList(geom.getCoordinates()));
			}
		}
		return coordCol;
	}

	/**
	 * knowing coord is in DT, locate the quad-edge whose origin is at coord
	 * 
	 * @param sd
	 * @param coord
	 * @return
	 */
	private static QuadEdge locateVertexInDT(QuadEdgeSubdivision sd, Coordinate coord) {
		QuadEdge qe = sd.locate(coord);
		if (qe == null) { // something wrong
			return null;
		}
		if (qe.orig().getCoordinate().equals2D(coord)) {
			return qe;
		} else if (qe.dest().getCoordinate().equals2D(coord)) {
			return qe.sym();
		} else {
			return null;
		}
	}

	/**
	 * get the hull from current DT, which is the CH of data (or almost the CH of
	 * data)
	 * 
	 * @param sd
	 * @param baseCoord
	 * @return
	 */
	private static LinkedList<Coordinate> getHull(QuadEdgeSubdivision sd, Coordinate baseCoord) {
		QuadEdge qe = locateVertexInDT(sd, baseCoord);
		if (qe == null) {
			return new LinkedList<>();
		}
		Vertex baseVer = qe.orig();
		// find a frame vertex
		while (!sd.isFrameVertex(qe.dest())) {
			qe = qe.oNext();
		}
		// turn to a hull edge
		while (sd.isFrameVertex(qe.dest())) {
			qe = qe.oNext();
		}
		LinkedList<Coordinate> hull = new LinkedList<>();
		do {
			Vertex ver = qe.orig();
			hull.add(ver.getCoordinate());
			qe = qe.rPrev();
			while (sd.isFrameVertex(qe.dest())) {
				qe = qe.oNext();
			}
		} while (qe.orig() != baseVer);
		hull.add(hull.getFirst()); // close the ring, now first == last
		return hull;
	}

	//
	/**
	 * generate a circular list edge representation of a hull in the form of a
	 * CLOSED coordinate list. the order of the coordinates should be CCW
	 * 
	 * @param hullCoords
	 * @return
	 */
	private DLCirList<Coordinate> generateHullEdgeRep(LinkedList<Coordinate> hullCoords) {
		DLCirList<Coordinate> rtn = new DLCirList<>();
		ListIterator<Coordinate> iter = hullCoords.listIterator(1); // 2nd in list (should have at least 3 coords)
		while (iter.hasNext()) {
			Coordinate coord = iter.next();
			rtn.add(coord);
		}
		return rtn;
	}

	//
	/**
	 * @param hullCL
	 * @param coordNodeMap
	 * @param edgeIdx
	 * @param nodeEdgeMap
	 */
	private void generateHullIndices(DLCirList<Coordinate> hullCL, Map<Coordinate, DLNode<Coordinate>> coordNodeMap,
			Set<HullEdgeCir> edgeIdx, Map<DLNode<Coordinate>, HullEdgeCir> nodeEdgeMap) {
		DLNode<Coordinate> node = hullCL.getNode();
		do {
			Coordinate coord = node.getObj();
			coordNodeMap.put(coord, node);
			HullEdgeCir e = new HullEdgeCir(node);
			edgeIdx.add(e);
			nodeEdgeMap.put(node, e);
			node = node.getNext();
		} while (node != hullCL.getNode());
	}

	private void generateHullIndices(DLCirList<Coordinate> hullCL, TriMetricLength m, Map<Coordinate, DLNode<Coordinate>> coordNodeMap,
			Set<HullEdgeCir> edgeIdx, Map<DLNode<Coordinate>, HullEdgeCir> nodeEdgeMap) {
		DLNode<Coordinate> node = hullCL.getNode();
		do {
			Coordinate coord = node.getObj();
			coordNodeMap.put(coord, node);
			HullEdgeCir e = new HullEdgeCir(node, m);
			edgeIdx.add(e);
			nodeEdgeMap.put(node, e);
			node = node.getNext();
		} while (node != hullCL.getNode());

	}

	private class HullEdgeCir implements Comparable<Object> {
		DLNode<Coordinate> sn;
		double metric;

		public HullEdgeCir(DLNode<Coordinate> sn) {
			this.sn = sn;
			metric = sn.getObj().distance(sn.getNext().getObj());
		}

		public HullEdgeCir(DLNode<Coordinate> sn, TriMetricLength m) {
			this.sn = sn;
			DLNode<Coordinate> en = sn.getNext();
			Coordinate s = sn.getObj();
			Coordinate e = en.getObj();
			QuadEdge qe = sd.locate(s, e);
			if (qe == null) {
				System.out.println("fail to find edge: " + s.toString() + " -" + e.toString());
			} else {
				QuadEdge lnext = qe.lNext();
				Vertex verO = lnext.dest();
				Coordinate o = verO.getCoordinate();

				metric = m.compMetric(s, e, o);
			}
		}

		public int compareTo(Object o) {
			HullEdgeCir other = (HullEdgeCir) o;
			if (metric < other.metric) {
				return -1;
			} else if (metric > other.metric) {
				return 1;
			} else {
				int rlt = sn.getObj().compareTo(other.sn.getObj());
				if (rlt == 0) {
					return (sn.getNext().getObj().compareTo(other.sn.getNext().getObj()));
				} else {
					return rlt;
				}
			}
		}
	}

	private class DLNode<T> {

		DLNode<T> prev, next;
		T obj;

		public DLNode(T o) {
			prev = next = null;
			obj = o;
		}

		public DLNode<T> getPrev() {
			return prev;
		}

		public void setPrev(DLNode<T> prev) {
			this.prev = prev;
		}

		public DLNode<T> getNext() {
			return next;
		}

		public T getObj() {
			return obj;
		}

		/**
		 * @param node
		 */
		public void insertAfter(DLNode<T> node) {
			if (node != null) {
				if (this.next != null) {
					this.next.setPrev(node);
				}
				node.prev = this;
				node.next = this.next;
				this.next = node;
			}
		}

		public DLNode<T> remove() {
			if (this.prev != null) {
				this.prev.next = this.next;
			}
			if (this.next != null) {
				this.next.prev = this.prev;
			}
			DLNode<T> rtn = null;
			if (this.prev != null) {
				rtn = this.prev;
			} else if (this.next != null) {
				rtn = this.next;
			}
			this.prev = this.next = null;
			return rtn;
		}

	}

	private class DLCirList<T> {

		DLNode<T> anchor = null;
		int size = 0;

		public DLCirList() {

		}

		public DLCirList(DLNode<T> node) {
			anchor = node;
			if (anchor.next == null && anchor.prev == null) {
				anchor.next = anchor.prev = anchor;
			}
			calculateSize();
		}

		public DLNode<T> getNode() {
			return anchor;
		}

		private void init(T o) {
			anchor = new DLNode<>(o);
			anchor.next = anchor.prev = anchor;
			size++;
		}

		public void add(T o) {
			if (anchor == null) {
				init(o);
			} else {
				DLNode<T> newNode = new DLNode<>(o);
				anchor.insertAfter(newNode);
				anchor = newNode;
				size++;
			}
		}

		public void addAfter(DLNode<T> node, DLNode<T> newNode) {
			node.insertAfter(newNode);
			size++;
		}

		public void remove(DLNode<T> node) {
			DLNode<T> rtn = node.remove();
			if (node == anchor) {
				anchor = rtn;
			}
			size--;
		}

		public int calculateSize() {
			if (anchor == null) {
				size = 0;
				return size;
			}
			DLNode<T> node = anchor;
			int cnt = 1;
			while (node.next != null && node.next != anchor) {
				node = node.next;
				cnt++;
			}
			size = cnt;
			return size;
		}

		public int size() {
			return size;
		}

		/**
		 * split, remove ns-ne and connect ns-no and no-ne, new circular list containing
		 * no-ne is returned
		 * 
		 * @param ns
		 * @param ne
		 * @param no
		 * @return
		 */
		public DLCirList<T> split(DLNode<T> ns, DLNode<T> ne, DLNode<T> no) {
			if (ns.next != ne) {
				return null;
			}
			if (size < 3) {
				calculateSize();
				if (size < 3) {
					return null;
				}
			}
			ns.next = no;
			DLNode<T> no2 = new DLNode<>(no.obj);
			no2.prev = no.prev;
			no.prev.next = no2;
			no2.next = ne;
			ne.prev = no2;
			no.prev = ns;
			anchor = ns;
			calculateSize();
			return new DLCirList<T>(ne);
		}
	}
}
