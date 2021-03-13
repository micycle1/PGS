package micycle.pts.utility;
//package micycle.pts;
//
///*******************************************************************************
// * This software is released under the licence CeCILL
// * 
// * see Licence_CeCILL-C_fr.html see Licence_CeCILL-C_en.html
// * 
// * see <a href="http://www.cecill.info/">http://www.cecill.info/a>
// * 
// * @copyright IGN
// * @copyright twak
// ******************************************************************************/
//
//import java.util.ArrayList;
//import java.util.HashMap;
//import java.util.HashSet;
//import java.util.Iterator;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//
//import fr.ign.cogit.geoxygene.api.feature.IFeature;
//import fr.ign.cogit.geoxygene.api.feature.IFeatureCollection;
//import fr.ign.cogit.geoxygene.api.feature.IPopulation;
//import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IDirectPosition;
//import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IDirectPositionList;
//import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineSegment;
//import fr.ign.cogit.geoxygene.api.spatial.coordgeom.ILineString;
//import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IPolygon;
//import fr.ign.cogit.geoxygene.api.spatial.geomaggr.IMultiCurve;
//import fr.ign.cogit.geoxygene.api.spatial.geomaggr.IMultiPoint;
//import fr.ign.cogit.geoxygene.api.spatial.geomprim.IPoint;
//import fr.ign.cogit.geoxygene.api.spatial.geomroot.IGeometry;
//import fr.ign.cogit.geoxygene.contrib.cartetopo.Arc;
//import fr.ign.cogit.geoxygene.contrib.cartetopo.CarteTopo;
//import fr.ign.cogit.geoxygene.contrib.cartetopo.CarteTopoFactory;
//import fr.ign.cogit.geoxygene.contrib.cartetopo.Groupe;
//import fr.ign.cogit.geoxygene.contrib.cartetopo.Noeud;
//import fr.ign.cogit.geoxygene.contrib.delaunay.TriangulationJTS;
//import fr.ign.cogit.geoxygene.contrib.geometrie.Operateurs;
//import fr.ign.cogit.geoxygene.feature.DefaultFeature;
//import fr.ign.cogit.geoxygene.feature.FT_FeatureCollection;
//import fr.ign.cogit.geoxygene.spatial.coordgeom.DirectPositionList;
//import fr.ign.cogit.geoxygene.spatial.coordgeom.GM_LineString;
//import fr.ign.cogit.geoxygene.spatial.geomprim.GM_Ring;
//import fr.ign.cogit.geoxygene.util.algo.geometricAlgorithms.LineDensification;
//
///**
// * Extract the spinal column of polygons
// * 
// * @author JFG
// * 
// */
//public class Spinalize {
//
//	// Triangulated irregular network
//
//	/**
//	 * Beastie method to skeletonize (only the spinal column) a list of polygons
//	 * using JTS. The method takes into account adjacent borders between polygons,
//	 * in order to generate a connected network of linestrings.
//	 * 
//	 * @param a          list of polygon
//	 * @param lengthMin, overSample, removeHoles
//	 * @return a list of linestring
//	 */
//	public static List<ILineString> spinalize(List<IPolygon> listPoly, double lengthMin, double overSample,
//			boolean removeHoles) {
//
//		List<ILineString> skeletonLines = new ArrayList<ILineString>();
//
//		// Eliminates holes in the polygons
//		if (removeHoles == true) {
//			for (IPolygon poly : listPoly) {
//				if (poly.getInterior().size() > 0) {
//					for (int it = 0; it < poly.getInterior().size(); it++) {
//						poly.removeInterior(it);
//					}
//				}
//			}
//		}
//
//		IFeatureCollection<IFeature> ftColVoronoiTotal = new FT_FeatureCollection<IFeature>();
//		IFeatureCollection<IFeature> ftColCheminsTotal = new FT_FeatureCollection<IFeature>();
//
//		// 1 - Identification of adjacent borders
//		// Detect adjacent polygons and get the middle of common borders
//		Map<IPolygon, IDirectPositionList> mapPtMedParPoly = new HashMap<IPolygon, IDirectPositionList>();
//		for (IPolygon polyRef : listPoly) {
//			IDirectPositionList listPtsMed = new DirectPositionList();
//			mapPtMedParPoly.put(polyRef, listPtsMed);
//			for (IPolygon polyComp : listPoly) {
//				if (polyRef.equals(polyComp)) {
//					continue;
//				}
//				if (polyRef.intersects(polyComp)) {
//					if (polyRef.intersection(polyComp).isMultiCurve()) {
//						@SuppressWarnings("unchecked")
//						IMultiCurve<ILineString> multiLs = (IMultiCurve<ILineString>) polyRef.intersection(polyComp);
//						IFeatureCollection<IFeature> ftColLs = new FT_FeatureCollection<IFeature>();
//						for (ILineString ls : multiLs.getList()) {
//							ftColLs.add(new DefaultFeature(ls));
//						}
//						CarteTopo carteTopoLs = CarteTopoFactory.newCarteTopo(ftColLs);
//						carteTopoLs.filtreNoeudsSimples();
//						carteTopoLs.filtreNoeudsIsoles();
//						for (Arc arc : carteTopoLs.getPopArcs()) {
//							IDirectPosition dpMilieu = Operateurs.milieu(arc.getGeometrie());
//							listPtsMed.add(dpMilieu);
//						}
//					}
//					if (polyRef.intersection(polyComp).isLineString()) {
//						ILineString ls = (ILineString) polyRef.intersection(polyComp);
//						listPtsMed.add(ls.centroid());
//					}
//				} else {
//					continue;
//				}
//			}
//		}
//
//		// 2 - Computation of the skeleton
//		// Compute the Voronoi diagram for each polygon and the path between
//		// extremities and middle of adjacent borders;
//		DirectPositionList dplPtsMedProj = new DirectPositionList();
//		IFeatureCollection<IFeature> ftColArcsProjIsoles = new FT_FeatureCollection<IFeature>();
//		Iterator<IPolygon> itPoly;
//		itPoly = mapPtMedParPoly.keySet().iterator();
//		while (itPoly.hasNext()) {
//			IPolygon poly = itPoly.next();
//			IDirectPositionList dplist = mapPtMedParPoly.get(poly);
//			if (!(poly.getExterior().coord().size() < 5)) {
//				IFeatureCollection<IFeature> ftColVoronoi = Spinalize.computeVoronoiDiagram(poly, overSample);
//				// Create a segment between the middle of adjacent borders and the
//				// closest point of the Voronoi diagram
//				for (IDirectPosition dp : dplist) {
//					IDirectPositionList dplVoronoi = new DirectPositionList();
//					for (IFeature ft : ftColVoronoi) {
//						IDirectPosition dp0 = ft.getGeom().coord().get(0);
//						IDirectPosition dp1 = ft.getGeom().coord().get(1);
//						dplVoronoi.add(dp0);
//						dplVoronoi.add(dp1);
//					}
//					double distanceMin = Double.MAX_VALUE;
//					IDirectPosition dpPtProche = null;
//					for (IDirectPosition dpVoronoi : dplVoronoi) {
//						double distance = dp.distance2D(dpVoronoi);
//						if (distance < distanceMin) {
//							distanceMin = distance;
//							dpPtProche = dpVoronoi;
//						}
//					}
//					dplPtsMedProj.add(dpPtProche);
//					ILineString ls = new GM_LineString(new DirectPositionList(dp, dpPtProche));
//					ftColCheminsTotal.add(new DefaultFeature(ls));
//					ftColArcsProjIsoles.add(new DefaultFeature(ls));
//					ftColVoronoi.add(new DefaultFeature(ls));
//				}
//
//				// For isolated polygons
//				if (dplist.size() == 0) {
//					// Compute the path between extremities (e.g. longest euclidian
//					// distance)
//					double distanceMax = 0.0;
//					IDirectPosition dpStart = null;
//					IDirectPosition dpEnd = null;
//					for (IDirectPosition dp1 : poly.coord()) {
//						for (IDirectPosition dp2 : poly.coord()) {
//							double distance = dp1.distance2D(dp2);
//							if (distanceMax < distance) {
//								distanceMax = distance;
//								dpStart = dp1;
//								dpEnd = dp2;
//							}
//						}
//					}
//					ftColCheminsTotal.addAll(Spinalize.computeShortestPath(dpStart, dpEnd, ftColVoronoi));
//					ftColCheminsTotal.addAll(Spinalize.computeShortestPath(dpEnd, dpStart, ftColVoronoi));
//
//					// For connected polygons
//				} else if (dplist.size() > 0) {
//					// Compute the path between extremities...
//					Double distanceMax = 0.0;
//					IDirectPosition dpStart = null;
//					IDirectPosition dpEnd = null;
//					for (IDirectPosition dp1 : poly.coord()) {
//						for (IDirectPosition dp2 : poly.coord()) {
//							double distance = dp1.distance2D(dp2);
//							if (distanceMax < distance) {
//								distanceMax = distance;
//								dpStart = dp1;
//								dpEnd = dp2;
//							}
//						}
//					}
//					ftColCheminsTotal.addAll(Spinalize.computeShortestPath(dpStart, dpEnd, ftColVoronoi));
//					ftColCheminsTotal.addAll(Spinalize.computeShortestPath(dpEnd, dpStart, ftColVoronoi));
//
//					// ... and from the middle of each border with the most far away point
//					// of the polygon geometry
//					for (IDirectPosition dp : dplist) {
//						Double distanceMedMax = 0.0;
//						IDirectPosition dpMedStart = dp;
//						IDirectPosition dpMedEnd = null;
//						for (IDirectPosition dpPoly : poly.coord()) {
//							double distance = dpMedStart.distance2D(dpPoly);
//							if (distanceMedMax < distance) {
//								distanceMedMax = distance;
//								dpMedEnd = dpPoly;
//							}
//						}
//						ftColCheminsTotal.addAll(Spinalize.computeShortestPath(dpMedStart, dpMedEnd, ftColVoronoi));
//						ftColCheminsTotal.addAll(Spinalize.computeShortestPath(dpMedEnd, dpMedStart, ftColVoronoi));
//					}
//				}
//				ftColVoronoiTotal.addAll(ftColVoronoi);
//			}
//		}
//
//		// 3 - Cleaning of the skeleton
//		// Get all computed paths and generate the topology
//		CarteTopo carteTopoChemins = CarteTopoFactory.newCarteTopo(ftColCheminsTotal);
//		carteTopoChemins.filtreNoeudsSimples();
//		carteTopoChemins.filtreNoeudsIsoles();
//
//		// Pruning of the skeleton using the lengthMin parameter
//		IFeatureCollection<IFeature> ftColCheminsSelect = new FT_FeatureCollection<IFeature>();
//		for (Arc arc : carteTopoChemins.getPopArcs()) {
//			int degreStart = arc.getNoeudIni().arcs().size();
//			int degreEnd = arc.getNoeudFin().arcs().size();
//
//			// Isolated arc (e.g. isolated polygon) then it is conserved
//			if (degreStart == 1 && degreEnd == 1) {
//				ftColCheminsSelect.add(arc);
//			}
//			// Connected arc on both sides, then it is conserved
//			else if (degreStart > 1 && degreEnd > 1) {
//				ftColCheminsSelect.add(arc);
//			}
//			// Dead-end, only conserved if it is longer than the lengthMin parameter
//			else {
//				if (degreStart == 1 && arc.getGeom().length() > lengthMin) {
//					ftColCheminsSelect.add(arc);
//				}
//				if (degreEnd == 1 && arc.getGeom().length() > lengthMin) {
//					ftColCheminsSelect.add(arc);
//				}
//			}
//		}
//
//		// Generate the topology of the selected arcs of the skeleton
//		CarteTopo carteTopoCheminsClean = CarteTopoFactory.newCarteTopo(ftColCheminsSelect);
//		carteTopoCheminsClean.filtreNoeudsSimples();
//		carteTopoCheminsClean.filtreNoeudsIsoles();
//
//		// Get segments of connexion between borders and each Voronoi diagram...
//		CarteTopo carteTopoArcsIsoles = CarteTopoFactory.newCarteTopo(ftColArcsProjIsoles);
//		carteTopoArcsIsoles.filtreNoeudsSimples();
//		carteTopoArcsIsoles.filtreNoeudsIsoles();
//
//		// ... and remove isolated arcs of the skeleton.
//		Iterator<Arc> itArc;
//		itArc = carteTopoCheminsClean.getPopArcs().iterator();
//		while (itArc.hasNext()) {
//			Arc arc = itArc.next();
//			for (Arc arcProj : carteTopoArcsIsoles.getPopArcs()) {
//				if (arc.getGeom().equals(arcProj.getGeom())) {
//					itArc.remove();
//				}
//			}
//		}
//
//		IFeatureCollection<IFeature> ftColSkeletons = new FT_FeatureCollection<IFeature>();
//		// Gaussian smoothing of selected arcs of the skeleton
//		for (Arc arc : carteTopoCheminsClean.getPopArcs()) {
//			// ILineString lsCheminClean = Filtering.DouglasPeuckerLineString(
//			// GaussianFilter.gaussianFilter(new GM_LineString(arc.getCoord()), 30,
//			// 1), 20);
//			// skeletonLines.add(lsCheminClean);
//			ftColSkeletons.add(arc);
//		}
//
//		for (IPolygon poly : listPoly) {
//			IGeometry intersection = poly.intersection(ftColSkeletons.getGeomAggregate());
//			if (intersection.isLineString()) {
//				ILineString ls = (ILineString) intersection;
//				// ILineString lsCheminClean = Filtering.DouglasPeuckerLineString(
//				// GaussianFilter.gaussianFilter(ls, 30, 1), 20);
//				// skeletonLines.add(lsCheminClean);
//				skeletonLines.add(ls);
//			}
//			if (intersection.isMultiCurve()) {
//				IMultiCurve<ILineString> multiLs = (IMultiCurve<ILineString>) intersection;
//				for (ILineString ls : multiLs.getList()) {
//					// ILineString lsCheminClean = Filtering.DouglasPeuckerLineString(
//					// GaussianFilter.gaussianFilter(ls, 30, 1), 20);
//					// skeletonLines.add(lsCheminClean);
//					skeletonLines.add(ls);
//				}
//			}
//		}
//
//		return skeletonLines;
//	}
//
//	/**
//	 * Compute the Voronoi diagram of a polygon using JTS triangulation library
//	 * 
//	 * @param polygon
//	 * @param threshold (pas de sur�chantillonnage)
//	 * @return
//	 */
//	private static IFeatureCollection<IFeature> computeVoronoiDiagram(IPolygon polygon, double overSample) {
//
//		// Densification of the polygon contour
//		IFeatureCollection<IFeature> ftColArcsVoronoi = new FT_FeatureCollection<IFeature>();
//		if (polygon.getExterior().coord().size() < 5) {
//			return ftColArcsVoronoi;
//		}
//		ILineString contourDense = LineDensification.densification(new GM_LineString(polygon.getExterior().coord()),
//				overSample);
//		polygon.setExterior(new GM_Ring(contourDense));
//
//		// Triangulation of the densified geometry using JTS
//		IFeatureCollection<IFeature> ftcolPoints = new FT_FeatureCollection<IFeature>();
//		for (IDirectPosition dp : polygon.coord()) {
//			ftcolPoints.add(new DefaultFeature(dp.toGM_Point()));
//		}
//		TriangulationJTS triangule = new TriangulationJTS("TriangulationJTS");
//		triangule.importAsNodes(ftcolPoints);
//		try {
//			triangule.triangule("v");
//		} catch (Exception e1) {
//			e1.printStackTrace();
//		}
//
//		// Compute Voronoi diagram
//		IPopulation<Arc> popArcVoronoi = triangule.getPopVoronoiEdges();
//		triangule.getVoronoiDiagram().filtreNoeudsSimples();
//		for (Arc arc : popArcVoronoi) {
//			if (arc.getGeom().within(polygon)) {
//				ftColArcsVoronoi.add(arc);
//			}
//		}
//		return ftColArcsVoronoi;
//	}
//
//	/**
//	 * Compute the shortest path between two points using the Voronoi diagram
//	 * 
//	 * @param dpStart,     dpEnd
//	 * @param ftcolVoronoi
//	 * @return
//	 * 
//	 **/
//	private static IFeatureCollection<IFeature> computeShortestPath(IDirectPosition dpStart, IDirectPosition dpEnd,
//			IFeatureCollection<IFeature> ftColVoronoi) {
//
//		IFeatureCollection<IFeature> ftColVoronoiClean = new FT_FeatureCollection<IFeature>();
//		for (IFeature ft : ftColVoronoi) {
//			if (!(ft.getGeom().coord().size() == 1)) {
//				ftColVoronoiClean.add(ft);
//			}
//		}
//
//		// Generate the topology of the Voronoi diagram
//		CarteTopo carteTopo = CarteTopoFactory.newCarteTopo(ftColVoronoiClean);
//		carteTopo.filtreNoeudsSimples();
//		carteTopo.filtreNoeudsIsoles();
//
//		// Identify closest nodes of start and end points of the path in the Voronoi
//		// diagram
//		Noeud noeudStart = null, noeudEnd = null;
//		double distanceDpStartMin = Double.MAX_VALUE;
//		double distanceDpEndMin = Double.MAX_VALUE;
//		for (Noeud noeud : carteTopo.getPopNoeuds()) {
//			double distance1 = dpStart.distance2D(noeud.getCoord());
//			double distance2 = dpEnd.distance2D(noeud.getCoord());
//			if (distance1 < distanceDpStartMin) {
//				distanceDpStartMin = distance1;
//				noeudStart = noeud;
//			}
//			if (distance2 < distanceDpEndMin) {
//				distanceDpEndMin = distance2;
//				noeudEnd = noeud;
//			}
//		}
//		if (noeudStart == null) {
//			return new FT_FeatureCollection<IFeature>();
//		}
//
//		// Compute the shortest path between the two points (if possible)
//		Groupe groupePCC = noeudStart.plusCourtChemin(noeudEnd, 0);
//		IFeatureCollection<IFeature> listChemins = new FT_FeatureCollection<IFeature>();
//		if (groupePCC == null) {
//			System.out.println("Foirage du squelette...");
//			return listChemins;
//		}
//		List<Arc> listArcs = groupePCC.getListeArcs();
//		for (Arc arc : listArcs) {
//			listChemins.add(new DefaultFeature(arc.getGeom()));
//		}
//		return listChemins;
//
//	}
//
//	/**
//	 * Connect a linear network to the nearest node of a skeleton.
//	 * 
//	 * @param skeleton
//	 * @param network
//	 * @param polygon
//	 */
//	public static Set<ILineString> connectSkeletonToNetwork(Set<ILineSegment> skeleton, Set<ILineString> network,
//			IPolygon polygon) {
//		Set<ILineString> extendedSkeleton = new HashSet<ILineString>();
//		// first, find the intersection between the network and the polygon
//		Set<IPoint> intersections = new HashSet<IPoint>();
//		for (ILineString line : network) {
//			if (line.intersects(polygon)) {
//				if (line.intersection(polygon) instanceof IPoint) {
//					intersections.add((IPoint) line.intersection(polygon));
//				} else if (line.intersection(polygon) instanceof IMultiPoint) {
//					IMultiPoint inter = (IMultiPoint) line.intersection(polygon);
//					for (int i = 0; i < inter.getList().size(); i++) {
//						intersections.add(inter.get(i));
//					}
//				}
//			}
//		}
//
//		for (IPoint point : intersections) {
//			double distanceMin = Double.MAX_VALUE;
//			IDirectPosition dpMin = null;
//			for (ILineSegment segment : skeleton) {
//				extendedSkeleton.add(segment);
//				double distanceMinSegment;
//				IDirectPosition dpMinSegment = null;
//				IDirectPosition dpStart = segment.startPoint();
//				IDirectPosition dpEnd = segment.endPoint();
//				double distanceStart = point.getPosition().distance2D(dpStart);
//				double distanceEnd = point.getPosition().distance2D(dpEnd);
//				if (distanceStart > distanceEnd) {
//					distanceMinSegment = distanceEnd;
//					dpMinSegment = dpEnd;
//				} else {
//					distanceMinSegment = distanceStart;
//					dpMinSegment = dpStart;
//				}
//				if (distanceMinSegment < distanceMin) {
//					distanceMin = distanceMinSegment;
//					dpMin = dpMinSegment;
//				}
//			}
//			ILineString lsConnect = new GM_LineString(point.getPosition(), dpMin);
//			extendedSkeleton.add(lsConnect);
//		}
//
//		return extendedSkeleton;
//	}
//
//}