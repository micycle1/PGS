/*
 *
 * version 0.3
 * 
 * date: 2019-01-19
 * 
 * author: Sheng Zhou (Sheng.Zhou@os.uk)
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
//
package uk.osgb.algorithm.minkowski_sum;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.operation.union.CascadedPolygonUnion;

/**
 * This is an experimental implementation of Minkowski sum and difference based
 * on JTS geometric functionality.
 * <p>
 * Current implementation supports Minkowski sums between a "source" geometry
 * (polygon/linestring/multipolygon/multilinestring/GeometryCollection) and a
 * "reference" geometry (polygon/linestring), and Minkowski difference between a
 * "source" geometry (polygon/multipolygon/GeometryCollection) and a "reference"
 * geometry (polygon/linestring)
 * <p>
 * Polygons may be concave. The "source" polygon may contain holes.
 * <p>
 * Any holes in "reference" polygon are ignored (in most cases it doesn't make
 * practical sense anyway).
 *
 * @author Sheng Zhou
 * @author Small improvements by Michael Carleton
 */
public class Minkowski_Sum {

	private static GeometryFactory gf = new GeometryFactory();

	/**
	 * Minkowski sum of two general geometries, computed without any presumption
	 * 
	 * @param src source geometry which might be
	 *            polygon/multipolygon/linestring/multilinestring
	 * @param ref reference geometry, which might be a polygon or a linestring
	 * @return minkowski sum of src and ref (or an empty geometry)
	 */
	public static Geometry minkSum(Geometry src, Geometry ref) {
		return minkSum(src, ref, false, false);
	}

	/*************************************************************
	 * 
	 * Minkowski Sum of a Polygon/MultiPolygon/LineString/MultiLineString and a
	 * Polygon/LineString
	 * 
	 *************************************************************/
	/**
	 * Minkowski sum of two general geometries, dispatched by manually testing types
	 * 
	 * @param src          source geometry which might be
	 *                     polygon/multipolygon/linestring/multilinestring
	 * @param ref          reference geometry, which might be a polygon or a
	 *                     linestring
	 * @param doReflection If reflection over origin is performed first on reference
	 *                     polygon (e.g. for collision detection purpose).
	 * @return the sum or an empty polygon (for un-supported types)
	 */
	public static Geometry minkSum(Geometry src, Geometry ref, boolean doReflection, boolean isRefConvex) {
		if (src == null || ref == null) {
			return gf.createPolygon();
		}
		if (ref instanceof Polygon) {
			if (src instanceof Polygon) {
				return compMinkSumPlgPlg((Polygon) src, (Polygon) ref, doReflection, isRefConvex);
			} else if (src instanceof LineString) {
				return compMinkSumLSPlg((LineString) src, (Polygon) ref, doReflection, isRefConvex);
			} else if (src instanceof MultiPolygon) {
				return compMinkSumMultiPlgPlg((MultiPolygon) src, (Polygon) ref, doReflection, isRefConvex);
			} else if (src instanceof MultiLineString) {
				return compMinkSumMultiLSPlg((MultiLineString) src, (Polygon) ref, doReflection, isRefConvex);
			} else if (src instanceof GeometryCollection) {
				return compMinkSumGeometryCollection((GeometryCollection) src, ref, doReflection, isRefConvex);
			} else if (src instanceof Point) {
				return compMinkSumPoint((Point) src, ref, doReflection, isRefConvex);
			}
		} else if (ref instanceof LineString || ref instanceof LinearRing) {
			if (src instanceof Polygon) {
				return compMinkSumPlgLS((Polygon) src, (LineString) ref, doReflection);
			} else if (src instanceof LineString) {
				return compMinkSumLSLS((LineString) src, (LineString) ref, doReflection);
			} else if (src instanceof MultiPolygon) {
				return compMinkSumMultiPlgLS((MultiPolygon) src, (LineString) ref, doReflection);
			} else if (src instanceof MultiLineString) {
				return compMinkSumMultiLSLS((MultiLineString) src, (LineString) ref, doReflection);
			} else if (src instanceof GeometryCollection) {
				return compMinkSumGeometryCollection((GeometryCollection) src, ref, doReflection, isRefConvex);
			} else if (src instanceof Point) {
				return compMinkSumPoint((Point) src, ref, doReflection, isRefConvex);
			}
		}
		return src.getFactory().createPolygon();
	}

	//
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkSumPoint(Point src, Geometry ref, boolean doReflection, boolean isRefConvex) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric on coordinate origin, the reflected geometry is the
							// same as original
			if (ref instanceof Polygon) {
				ref = reflectionOrgPlg((Polygon) ref);
			} else if (ref instanceof LineString || ref instanceof LinearRing) {
				ref = reflectionOrgLS((LineString) ref);
			}
		}
		AffineTransformation af = new AffineTransformation();
		af.translate(src.getX(), src.getY());
		return af.transform(ref);
	}

	/**
	 * Minkowski sum of two polygons. Any holes in reference polygon are ignored
	 * 
	 * @param src          Source polygon, may contain holes
	 * @param ref          Reference polygon, holes are ignored
	 * @param doReflection If reflection is performed first on reference polygon
	 *                     (e.g. for collision detection purpose). For repeated
	 *                     computation, reflection may be performed first before
	 *                     calling this method.
	 * @param isRefConvex  whether the reference polygon is convex (if unknown,
	 *                     please use false)
	 * @return
	 */
	private static Geometry compMinkSumPlgPlg(Polygon src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric on coordinate origin, the reflected geometry is the
							// same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkSumPlg(src, refCoords, isRefConvex);
	}

	/**
	 * Minkowski sum of a multipolygon and a polygon
	 * 
	 * @param src          source multipolygon, may contain holes
	 * @param ref          reference polygon, holes are ignored
	 * @param doReflection if reflection is performed first on reference polygon
	 *                     (for repeated computation, reflection may be performed
	 *                     first before calling this method)
	 * @param isRefConvex  whether the reference polygon is convex (if unknown,
	 *                     please use false)
	 * @return
	 */
	private static Geometry compMinkSumMultiPlgPlg(MultiPolygon src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the
							// same as original
			ref = reflectionOrgPlg(ref);
		}
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry sum1 = compMinkSumPlgPlg(plg, ref, false, isRefConvex);
			try {
				rlt = rlt.union(sum1);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}

	/**
	 * Minkowski sum between a polygon and a reference linestring
	 * 
	 * @param src          Soource polygon, may contain holes
	 * @param ref          reference linestring
	 * @param doReflection if reflection is performed first on reference linestring
	 *                     (e.g. for collision detection purpose). For repeated
	 *                     computation, reflection may be performed first before
	 *                     calling this method.
	 * @return
	 */
	private static Geometry compMinkSumPlgLS(Polygon src, LineString ref, boolean doReflection) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the
							// same as original
			ref = reflectionOrgLS(ref);
		}

		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkSumPlg(src, refCoords, false); // always treat LS ref as non-convex
	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	private static Geometry compMinkSumMultiPlgLS(MultiPolygon src, LineString ref, boolean doReflection) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric over coordinate origin, the geometry is the same as
							// original
			ref = reflectionOrgLS(ref);
		}
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry sum1 = compMinkSumPlgLS(plg, ref, false);
			try {
				rlt = rlt.union(sum1);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkSumLSPlg(LineString src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the
							// same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkSumLS(src, refCoords, isRefConvex);
	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	private static Geometry compMinkSumLSLS(LineString src, LineString ref, boolean doReflection) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the
							// same as original
			ref = reflectionOrgLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkSumLS(src, refCoords, false);

	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	private static Geometry compMinkSumMultiLSLS(MultiLineString src, LineString ref, boolean doReflection) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the
							// same as original
			ref = reflectionOrgLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			LineString ls = (LineString) src.getGeometryN(i);
			Geometry sum1 = compMinkSumLS(ls, refCoords, false);
			try {
				rlt = rlt.union(sum1);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkSumMultiLSPlg(MultiLineString src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if (doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the
							// same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			LineString ls = (LineString) src.getGeometryN(i);
			Geometry sum1 = compMinkSumLS(ls, refCoords, isRefConvex);
			try {
				rlt = rlt.union(sum1);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}

	/**
	 * Minkowskie sum of an arbitrary GeometryCollection and a polygon
	 * 
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkSumGeometryCollection(GeometryCollection src, Geometry ref, boolean doReflection, boolean isRefConvex) {
		if (src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			Geometry part = src.getGeometryN(i);
			Geometry sum1 = minkSum(part, ref, doReflection, isRefConvex);
			try {
				rlt = rlt.union(sum1);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}

	/**
	 * Minkowski sum between a linestring and a reference array of coordinates,
	 * assuming src is not null or emtpy (handled elsewehere)
	 * 
	 * @param src         source linestring
	 * @param refCoords   coordinate sequence, which forms either an open LineString
	 *                    or a closed LinearRing
	 * @param isRefConvex if refCoords is a closed LinearRing, whether it is convex
	 *                    or not
	 * @return
	 */
	private static Geometry compMinkSumLS(LineString src, Coordinate[] refCoords, boolean isRefConvex) {
		Coordinate[] coords = src.getCoordinates();
		// expanding the shell of original geometry (on two directions)
		if (refCoords[0].equals2D(refCoords[refCoords.length - 1])) { // reference is a linear ring
			return coordArrayVectorAddition(coords, refCoords, isRefConvex, src.getFactory());
		} else {
			return coordArrayVectorAddition(coords, refCoords, src.getFactory());
		}
	}

	/**
	 * @param src
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkSumPlg(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if (src == null || src.isEmpty()) {
			return src;
		}

		final LineString extRing = src.getExteriorRing();
		final Coordinate[] extCoords = extRing.getCoordinates();
		Geometry extSum = expansionShell(gf.createPolygon(extCoords), refCoords, isRefConvex);
		if (extSum == null) {
			return gf.createEmpty(2); // rather than return null
		}
		final int numHoles = src.getNumInteriorRing();
		if (numHoles > 0) {
			for (int k = 0; k < numHoles; ++k) {
				LineString holeRing = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeRing.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords);
				Geometry intDiff = shrinkHole(holePlg, refCoords, isRefConvex);
				extSum = extSum.difference(intDiff);
			}
		}
		return extSum;
	}

	/**
	 * Compute the "vector addition" of refCoords and the exterior ring of src and
	 * return the exterior ring of the result (including holes outside original src
	 * geometry)
	 * 
	 * @param src         Polygon
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry expansionShell(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if (src == null || src.isEmpty()) {
			return src;
		}
		boolean isRefRing = false;
		if (refCoords[0].equals2D(refCoords[refCoords.length - 1])) { // reference is a linear ring
			isRefRing = true;
		}

		Geometry sum = null;

		LineString extRing = src.getExteriorRing();
		Coordinate[] extCoords = extRing.getCoordinates();
		// expanding the shell of original geometry (on two directions)
		if (isRefRing) {
			sum = coordArrayVectorAddition(extCoords, refCoords, isRefConvex, src.getFactory());
		} else {
			sum = coordArrayVectorAddition(extCoords, refCoords, src.getFactory());
		}
		// if holes are outside original geometry, they are generated by boundary
		// expansion and should be kept as holes in final results
		Vector<LinearRing> ringVec = null;

		LinearRing[] shellHoles = getHoles((Polygon) sum);
		if (shellHoles != null) {
			Polygon geomNoHole = gf.createPolygon(extCoords);// polygon formed by the exterior ring of src geometry
			for (int i = 0; i < shellHoles.length; ++i) {
				// Point ep = shellHoles[i].getEndPoint(); // a point on the hole boundary
				Point ep = gf.createPolygon(shellHoles[i]).getInteriorPoint();
				if (!geomNoHole.contains(ep)) {// holes outside original geometry, generated by expansion, should be
												// part of the result as a hole
					// if(!geomNoHole.intersects(ep)){// holes outside original geometry, generated
					// by expansion, should be part of the result as a hole
					if (ringVec == null) {
						ringVec = new Vector<LinearRing>();
					}
					ringVec.add(shellHoles[i]);
				}
			}
		}
		// assemble the result
		LinearRing shell = getShellRing((Polygon) sum);
		if (ringVec != null) {
			// holes of sumHoles
			LinearRing[] holes = new LinearRing[ringVec.size()];
			ringVec.toArray(holes);
			return gf.createPolygon(shell, holes);
		} else {// no holes as result of expansion from original holes
			return gf.createPolygon(shell.getCoordinates());
		}
	}

	/******************************************************************
	 * 
	 * Minkowski Difference between Polygon/MultiPolygon and Polyon/LineString
	 * 
	 * 
	 *****************************************************************/

	public static Geometry minkDiff(Geometry src, Geometry ref) {
		return minkDiff(src, ref, false, false);
	}

	/**
	 * Minkowski difference of two geometries.
	 * 
	 * @param src          Source geometry, may be polygon or multipolygon
	 * @param ref          Reference geometry, may be polygon or
	 *                     linestring/linearring
	 * @param refSymmetric true if ref is symmetric over origin
	 * @param isRefConvex
	 * @return null (not implemented yet)
	 */
	public static Geometry minkDiff(Geometry src, Geometry ref, boolean isRefSymmetric, boolean isRefConvex) {
		if (src == null) {
			return gf.createPolygon();
		} else if (ref == null) {
			return src;
		}
		if (ref instanceof Polygon) {
			if (src instanceof Polygon) {
				return compMinkDiffPlgPlg((Polygon) src, (Polygon) ref, isRefSymmetric, isRefConvex);
			} else if (src instanceof MultiPolygon) {
				return compMinkDiffMultiPlgPlg((MultiPolygon) src, (Polygon) ref, isRefSymmetric, isRefConvex);
			} else if (src instanceof GeometryCollection) {
				return compMinkDiffGeometryCollection((GeometryCollection) src, ref, isRefSymmetric);
			} else {
				return src.getFactory().createPolygon();
			}
		} else if (ref instanceof LineString || ref instanceof LinearRing) {
			if (src instanceof Polygon) {
				return compMinkDiffPlgLS((Polygon) src, (LineString) ref, isRefSymmetric);
			} else if (src instanceof MultiPolygon) {
				return compMinkDiffMultiPlgLS((MultiPolygon) src, (LineString) ref, isRefSymmetric);
			} else if (src instanceof GeometryCollection) {
				return compMinkDiffGeometryCollection((GeometryCollection) src, ref, isRefSymmetric);
			} else {
				return src.getFactory().createPolygon(); // unsupported source type
			}
		} else {
			return src.getFactory().createPolygon();
		}
	}

	/**
	 * Minkowski difference of polygon (may contain holes) and reference geometry,
	 * using difference between minkDiff(exterior ring) and minkSum(holes)
	 * 
	 * @param src
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkDiffPlg(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if (src.isEmpty()) {
			return src;
		}

		// bufferring exterior ring of source geometry
		LineString extRing = src.getExteriorRing();
		Geometry extDiff = shrinkHole(gf.createPolygon(extRing.getCoordinates()), refCoords, isRefConvex);
		int numHoles = src.getNumInteriorRing();
		if (numHoles > 0) {
			for (int k = 0; k < numHoles; ++k) {
				LineString holeSrc = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeSrc.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords); // hole polygon
				Geometry intSum = expansionShell(holePlg, refCoords, isRefConvex);
				try {
					extDiff = extDiff.difference(intSum);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		return extDiff;
	}

	/**
	 * @param src         polygon (holes, if any, are ignored)
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry shrinkHole(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if (src.isEmpty()) {
			return src;
		}
		GeometryFactory gf = src.getFactory();
		boolean isRefRing = false;
		if (refCoords[0].equals2D(refCoords[refCoords.length - 1])) { // reference is a linear ring so isRefConvex is
																		// considered
			isRefRing = true;
		}
		Geometry sum = null;
		LineString extRing = src.getExteriorRing();
		Coordinate[] extCoords = extRing.getCoordinates();
		if (isRefRing) {
			sum = coordArrayVectorAddition(extCoords, refCoords, isRefConvex, src.getFactory());
		} else {
			sum = coordArrayVectorAddition(extCoords, refCoords, src.getFactory());
		}
		// holes in buffering result
		LinearRing[] shellHoles = getHoles((Polygon) sum);
		if (shellHoles == null) {// no holes, shrinked completely
			return src.getFactory().createPolygon();
		}
		//
		// need to check if the holes are inside the shell of the original geometry
		// holes could be generated from outwards expansion as well
		Vector<LinearRing> ringVec = new Vector<LinearRing>(shellHoles.length);
		Polygon geomNoHole = gf.createPolygon(extCoords); // original geom
		for (int i = 0; i < shellHoles.length; ++i) {
			// Point ep = shellHoles[i].getEndPoint();
			Point ep = gf.createPolygon(shellHoles[i]).getInteriorPoint();
			if (geomNoHole.contains(ep)) { // the hole is inside the original source geometry, not a hole generated by
											// outwards expansion
				// if(geomNoHole.intersects(ep)){ // the hole is inside the original source
				// geometry, not a hole generated by outwards expansion
				ringVec.add(shellHoles[i]);
			}
		}
		//
		if (ringVec.size() == 0) {
			return gf.createPolygon();
		} else {
			shellHoles = new LinearRing[ringVec.size()];
			ringVec.toArray(shellHoles);
		}
		// shell of result
		Geometry rltShell = null;
		if (shellHoles.length > 1) {// multipolygon
			Polygon[] plgs = new Polygon[shellHoles.length];
			for (int i = 0; i < plgs.length; ++i) {
				plgs[i] = gf.createPolygon(shellHoles[i]);
			}
			rltShell = gf.createMultiPolygon(plgs);
		} else {// polygon
			rltShell = gf.createPolygon(shellHoles[0]);
		}
		// holes
		if (rltShell != null) {
			return rltShell;
		} else {
			return gf.createPolygon();
		}
	}

	/**
	 * Minkowski difference of two polygons (src - ref). Holes in ref are ignored.
	 * 
	 * @param src          Source polygon
	 * @param ref          reference polgyon (without hole)
	 * @param doReflection whether ref should be reflected (for repeated calls, ref
	 *                     may be reflected before calling this method). If ref is
	 *                     symmetric with respect to the coordinate origin, there is
	 *                     no need to perform reflection either.
	 * @PARAM isRefConvex If you know the ref polygon is convex, set this to true
	 *        will improve performance and robustness.
	 * @return
	 */
	private static Geometry compMinkDiffPlgPlg(Polygon src, Polygon ref, boolean isRefSymmetric, boolean isRefConvex) {
		if (!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the
								// same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkDiffPlg(src, refCoords, isRefConvex);
	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	private static Geometry compMinkDiffPlgLS(Polygon src, LineString ref, boolean isRefSymmetric) {
		if (!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the
								// same as original
			ref = reflectionOrgLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkDiffPlg(src, refCoords, false);
	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkDiffMultiPlgPlg(MultiPolygon src, Polygon ref, boolean isRefSymmetric, boolean isRefConvex) {
		if (!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the
								// same as original
			ref = reflectionOrgPlg(ref);
		}
		Geometry rlt = null;
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry diff1 = compMinkDiffPlgPlg(plg, ref, false, isRefConvex);
			if (rlt == null) {
				rlt = diff1;
			} else {
				try {
					rlt = rlt.union(diff1);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		return rlt;
	}

	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	private static Geometry compMinkDiffMultiPlgLS(MultiPolygon src, LineString ref, boolean isRefSymmetric) {
		if (!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the
								// same as original
			ref = reflectionOrgLS(ref);
		}
		Geometry rlt = null;
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry diff1 = compMinkDiffPlgLS(plg, ref, false);
			if (rlt == null) {
				rlt = diff1;
			} else {
				try {
					rlt = rlt.union(diff1);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		return rlt;
	}

	//
	private static Geometry compMinkDiffGeometryCollection(GeometryCollection src, Geometry ref, boolean isRefSymmetric) {
		Geometry rlt = null;
		int numParts = src.getNumGeometries();
		for (int i = 0; i < numParts; ++i) {
			Geometry geom = src.getGeometryN(i);
			Geometry diff1 = minkDiff(geom, ref, isRefSymmetric, false);
			if (rlt == null) {
				rlt = diff1;
			} else {
				try {
					rlt = rlt.union(diff1);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		return rlt;
	}

	/**
	 * @param plg
	 * @return
	 */
	private static LinearRing getShellRing(Polygon plg) {
		return plg.getFactory().createLinearRing(plg.getExteriorRing().getCoordinates());
	}

	/**
	 * @param plg
	 * @return
	 */
	private static LinearRing[] getHoles(Polygon plg) {
		int numHoles = plg.getNumInteriorRing();
		if (numHoles > 0) {
			LinearRing[] holes = new LinearRing[numHoles];
			for (int i = 0; i < numHoles; ++i) {
				holes[i] = plg.getFactory().createLinearRing(plg.getInteriorRingN(i).getCoordinates());
			}
			return holes;
		} else {
			return null;
		}
	}

	/**
	 * Symmetry of JTS polygon with respect to coordinate origin
	 * 
	 * @param plg
	 * @return
	 */
	private static Polygon reflectionOrgPlg(Polygon plg) {
		Coordinate[] extNeg = reflectionOrgCoordArray(plg.getExteriorRing().getCoordinates());
		GeometryFactory gf = plg.getFactory();
		int numHoles = plg.getNumInteriorRing();
		if (numHoles > 0) {
			LinearRing[] holes = new LinearRing[numHoles];
			for (int i = 0; i < numHoles; ++i) {
				Coordinate[] holeNeg = reflectionOrgCoordArray(plg.getInteriorRingN(i).getCoordinates());
				holes[i] = gf.createLinearRing(holeNeg);
			}
			return gf.createPolygon(gf.createLinearRing(extNeg), holes);
		} else {
			return gf.createPolygon(extNeg);
		}
	}

	private static LineString reflectionOrgLS(LineString ls) {
		Coordinate[] coordsNeg = reflectionOrgCoordArray(ls.getCoordinates());
		return ls.getFactory().createLineString(coordsNeg);
	}

	/**
	 * symmetry of a Coordinate array with respect to coordinate origin
	 * 
	 * @param coords
	 * @return
	 */
	private static Coordinate[] reflectionOrgCoordArray(Coordinate[] coords) {
		if (coords != null && coords.length > 0) {
			Coordinate[] coordsNeg = new Coordinate[coords.length];
			for (int i = 0; i < coords.length; ++i) {
				Coordinate coord = coords[i];
				coordsNeg[i] = new Coordinate(-coord.x, -coord.y);
			}
			return coordsNeg;
		} else {
			return null;
		}
	}

	/**
	 * Given a line segment segSp-segEp, and coordinate array of a polygon
	 * refCoords, return the point set union of: vector sum of polygon refCoords and
	 * segSp; vector sum of polygon refCoords and segEp; all points covered by the
	 * polygon when it moves from segSp to segEp
	 * 
	 * @param refCoords   Coordinate array of the reference polygon (without holes)
	 * @param segSp       starting point of the line segment
	 * @param segEp       ending point of the line segment
	 * @param isRefConvex if the reference polygon refCoords is convex. If reference
	 *                    polygon is convex, the computation is simplier and robust
	 * @return
	 */
	private static Polygon segPlgAddition(Coordinate[] refCoords, Coordinate segSp, Coordinate segEp, boolean isRefConvex,
			GeometryFactory gf) {
		if (isRefConvex) {
			int numCoord = refCoords.length;
			Coordinate[] coords = new Coordinate[numCoord * 2];
			for (int i = 0; i < refCoords.length; ++i) {
				Coordinate ref = refCoords[i];
				coords[i] = new Coordinate(ref.x + segSp.x, ref.y + segSp.y);
				coords[i + numCoord] = new Coordinate(ref.x + segEp.x, ref.y + segEp.y);
			}
			ConvexHull ch = new ConvexHull(coords, gf);
			return (Polygon) ch.getConvexHull();
		} else {
			int numCoord = refCoords.length;
			Coordinate[] coords = new Coordinate[numCoord * 2];
			Coordinate[] coords1 = new Coordinate[numCoord];
			Coordinate[] coords2 = new Coordinate[numCoord];
			for (int i = 0; i < refCoords.length; ++i) {
				Coordinate ref = refCoords[i];
				coords[i] = new Coordinate(ref.x + segSp.x, ref.y + segSp.y);
				coords[i + numCoord] = new Coordinate(ref.x + segEp.x, ref.y + segEp.y);
				coords1[i] = (Coordinate) coords[i].clone();
				coords2[i] = (Coordinate) coords[i + numCoord].clone();
			}
			ConvexHull chAll = new ConvexHull(coords, gf);
			ConvexHull ch1 = new ConvexHull(coords1, gf);
			ConvexHull ch2 = new ConvexHull(coords2, gf);
			Geometry core = chAll.getConvexHull().difference(ch1.getConvexHull()).difference(ch2.getConvexHull());
			ConvexHull coreHull = new ConvexHull(core);
			Geometry coreNew = coreHull.getConvexHull();
			Polygon plg1 = gf.createPolygon(coords1);
			Polygon plg2 = gf.createPolygon(coords2);
			return (Polygon) coreNew.union(plg1).union(plg2);
		}
	}

	/**
	 * vector addition of a segment and a linestring (open, although it will handle
	 * linearing as well)
	 * 
	 * @param refCoords Coordinate array of the reference linestring (whether it is
	 *                  a ring is checked outside this method)
	 * @param segSp
	 * @param segEp
	 * @return the sum (may be an emtpy polygon)
	 */
	private static Polygon segLSAddition(Coordinate segSp, Coordinate segEp, Coordinate[] refCoords, GeometryFactory gf) {
		int numPts = refCoords.length;
		List<Polygon> parts = new ArrayList<>();
		for (int i = 0; i < numPts - 1; ++i) {
			Coordinate refSp = refCoords[i];
			Coordinate refEp = refCoords[i + 1];
			Polygon sum1 = segVectorAddition(segSp, segEp, refSp, refEp, gf);
			if (sum1 != null) {
				parts.add(sum1);
			}
		}
		return (Polygon) CascadedPolygonUnion.union(parts);
	}

	/**
	 * vector addition of a segment on reference linestring and a segment on the
	 * source geometry. This is the basic operation upon which all other methods are
	 * built Convexhull is used to construct result to improve robustness (avoiding
	 * self-intersection, hopefully).
	 * 
	 * @param segSp
	 * @param segEp
	 * @param refSp
	 * @param refEp
	 * @return a polygon (empty polygon if the sum is 0D or 1D)
	 */
	private static Polygon segVectorAddition(Coordinate segSp, Coordinate segEp, Coordinate refSp, Coordinate refEp, GeometryFactory gf) {
		Coordinate[] coords = new Coordinate[4];
		coords[0] = new Coordinate(refSp.x + segSp.x, refSp.y + segSp.y);
		coords[1] = new Coordinate(refEp.x + segSp.x, refEp.y + segSp.y);
		coords[2] = new Coordinate(refEp.x + segEp.x, refEp.y + segEp.y);
		coords[3] = new Coordinate(refSp.x + segEp.x, refSp.y + segEp.y);
		ConvexHull ch = new ConvexHull(coords, gf);
		Polygon hull = (Polygon) ch.getConvexHull();
		if (hull.getDimension() > 1) {
			return hull; // a polygon
		} else {
			return gf.createPolygon(); // return an empty polygon
		}
	}

	/**
	 * computing the vector sum of a linestring/polygon (without holes) and a
	 * polygon in form of coordinate array.
	 * 
	 * @param geomCoords   Coordinate array of the source linestring/polygon
	 *                     (without holes)
	 * @param refPlgCoords the reference polygon
	 * @param isRefConvex  whether the reference polygon is convex
	 * @return
	 */
	private static Geometry coordArrayVectorAddition(Coordinate[] geomCoords, Coordinate[] refPlgCoords, boolean isRefConvex,
			GeometryFactory gf) {
		List<Polygon> parts = new ArrayList<>();
		for (int j = 0; j < geomCoords.length - 1; ++j) {
			Coordinate segSp = geomCoords[j];
			Coordinate segEp = geomCoords[j + 1];
			if (!segSp.equals2D(segEp)) {
				Polygon sum1 = segPlgAddition(refPlgCoords, segSp, segEp, isRefConvex, gf);
				if (sum1 != null) {
					parts.add(sum1);
				}
			}
		}
		return CascadedPolygonUnion.union(parts);
	}

	/**
	 * computing the vector sum of a linestring or polygon (without holes) and a
	 * linestring (maybe closed) in form of coordinate array.
	 * 
	 * @param geomCoords Coordinate array of source geometry (linestring, or polygon
	 *                   without holes)
	 * @param refCoords  Coordinate array of reference geometry, an open linestring
	 *                   or a closed lienarring
	 * @return
	 */
	private static Geometry coordArrayVectorAddition(Coordinate[] geomCoords, Coordinate[] refCoords, GeometryFactory gf) {
		int numPts = refCoords.length;
		if (refCoords[0].equals2D(refCoords[numPts - 1])) { // linear ring
			return coordArrayVectorAddition(geomCoords, refCoords, false, gf);
		}
		List<Polygon> parts = new ArrayList<>();
		for (int j = 0; j < geomCoords.length - 1; ++j) {
			Coordinate segSp = geomCoords[j];
			Coordinate segEp = geomCoords[j + 1];
			if (!segSp.equals2D(segEp)) {
				Polygon sum1 = segLSAddition(segSp, segEp, refCoords, gf);
				if (sum1 != null) {
					parts.add(sum1);
				}
			}
		}
		return CascadedPolygonUnion.union(parts);
	}

	public static void setGeometryFactory(GeometryFactory geomFactory) {
		gf = geomFactory;
	}

}