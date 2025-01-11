package micycle.pgs.commons;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Polygon;

/**
 * 
 * The ProcrustesAlignment class provides methods for performing
 * ProcrustesAlignment analysis, which is a technique for aligning and comparing
 * geometric shapes. ProcrustesAlignment analysis aims to find the best
 * transformation (translation, rotation, and scaling) that minimizes the
 * differences between corresponding points of two shapes.
 * <p>
 * This class is particularly useful in shape matching and analysis tasks where
 * the only permitted deformation modes are uniform scaling, rotation, and
 * translation. It is commonly used in various fields such as computer vision,
 * image processing, pattern recognition, and bioinformatics.
 * 
 * @author Michael Carleton
 */
public class ProcrustesAlignment {

	private ProcrustesAlignment() {
	}

	/**
	 * Performs <i>ProcrustesAlignment Analysis</i> to align two polygons.
	 * <p>
	 * Finds the optimal scaling, translation and rotation to best align
	 * <code>transformPolygon</code> with respect to <code>referencePolygon</code>.
	 * <p>
	 * Note: the polygons should have the same number of vertices.
	 *
	 * @param referencePolygon the first polygon
	 * @param transformPolygon the polygon to transform/align
	 * @return an array having 4 values: the optimal translation (x, y), scale, and
	 *         rotation angle (radians, clockwise) to best align the transform
	 *         polygon to the reference polygon.
	 */
	public static double[] transform(final Polygon referencePolygon, final Polygon transformPolygon) {
		final Coordinate[] coordsA = referencePolygon.getExteriorRing().getCoordinates();
		final Coordinate[] coordsB = transformPolygon.getExteriorRing().getCoordinates();

		if (coordsA.length != coordsB.length) {
			throw new IllegalArgumentException(
					"Polygon exterior rings are different lengths (" + coordsA.length + ", " + coordsB.length + ")!");
		}

		// Find optimal translation
		Coordinate t = findTranslation(referencePolygon, transformPolygon);

		// Shift to origin (required for scaling & rotation)
		Coordinate ca = referencePolygon.getCentroid().getCoordinate();
		Coordinate cb = transformPolygon.getCentroid().getCoordinate();
		translate(coordsA, -ca.x, -ca.y);
		translate(coordsB, -cb.x, -cb.y);

		// Find optimal scale
		double scale = findScale(coordsA, coordsB);

		// Find optimal rotation
		double theta = findRotation(coordsA, coordsB);

		// Return transformation parameters
		return new double[] { t.x, t.y, scale, theta };
	}

	/**
	 * Translates an array of coordinates by dx and dy.
	 * 
	 * @param coords the array of coordinates to translate
	 * @param dx     the translation in x
	 * @param dy     the translation in y
	 */
	private static void translate(final Coordinate[] coords, final double dx, final double dy) {
		for (Coordinate c : coords) {
			c.x += dx;
			c.y += dy;
		}
	}

	/**
	 * Finds the optimal translation vector between two polygons.
	 *
	 * @param a the first polygon
	 * @param b the second polygon
	 * @return the translation vector as a Coordinate
	 */
	private static Coordinate findTranslation(Polygon a, Polygon b) {
		Coordinate ca = a.getCentroid().getCoordinate();
		Coordinate cb = b.getCentroid().getCoordinate();

		return new Coordinate(ca.x - cb.x, ca.y - cb.y);

	}

	/**
	 * Finds the optimal scale factor to match the size of two coordinate arrays.
	 * 
	 * @param aCoords the first coordinate array
	 * @param bCoords the second coordinate array
	 * @return the scale factor
	 */
	private static double findScale(final Coordinate[] aCoords, final Coordinate[] bCoords) {
		double rmsd1 = getRMSD(aCoords);
		double rmsd2 = getRMSD(bCoords);

		return rmsd1 / rmsd2;
	}

	/**
	 * Calculates the RMSD (root mean square deviation) of a coordinate array.
	 *
	 * @param coords the array of coordinates
	 * @return the RMSD
	 */
	private static double getRMSD(final Coordinate[] coords) {
		double sum = 0;

		for (Coordinate c : coords) {
			sum += c.x * c.x + c.y * c.y;
		}

		return Math.sqrt(sum / coords.length);
	}

	/**
	 * Finds the optimal rotation angle to align two coordinate arrays.
	 * 
	 * @param c1 the first coordinate array
	 * @param c2 the second coordinate array
	 * @return the rotation angle in radians
	 */
	private static double findRotation(final Coordinate[] c1, final Coordinate[] c2) {
		RealMatrix m = computeOptimalRotationMatrix(c1, c2);
		double cosTheta = m.getEntry(0, 0);
		double sinTheta = m.getEntry(1, 0);

		double theta = Math.atan2(sinTheta, cosTheta);
		return -theta; // NOTE negate for CW rotation
	}

	/**
	 * Computes the optimal rotation matrix to align two coordinate arrays using
	 * Kabsch algorithm.
	 * <p>
	 * Each array must be translated first, such that its centroid coincides with
	 * (0, 0).
	 *
	 * @param c1 the first coordinate array
	 * @param c2 the second coordinate array
	 * @return the 2x2 rotation matrix
	 */
	private static RealMatrix computeOptimalRotationMatrix(final Coordinate[] c1, final Coordinate[] c2) {
		final double[][] covarianceMatrix = new double[2][2];
		for (int i = 0; i < c1.length; i++) {
			final double[] translatedSourcePoint = { c1[i].x, c1[i].y };
			final double[] translatedTargetPoint = { c2[i].x, c2[i].y };

			covarianceMatrix[0][0] += translatedSourcePoint[0] * translatedTargetPoint[0];
			covarianceMatrix[0][1] += translatedSourcePoint[0] * translatedTargetPoint[1];
			covarianceMatrix[1][0] += translatedSourcePoint[1] * translatedTargetPoint[0];
			covarianceMatrix[1][1] += translatedSourcePoint[1] * translatedTargetPoint[1];
		}

		final RealMatrix covarianceMatrixMatrix = MatrixUtils.createRealMatrix(covarianceMatrix);

		final SingularValueDecomposition svd = new SingularValueDecomposition(covarianceMatrixMatrix);
		final RealMatrix rotationMatrix = svd.getV().multiply(svd.getUT());
		return rotationMatrix;
	}

}