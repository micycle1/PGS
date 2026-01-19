package micycle.pgs.commons;

import java.util.Objects;

import org.locationtech.jts.algorithm.distance.DiscreteHausdorffDistance;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;

import com.github.micycle1.geoblitz.ProHausdorffDistance;

/**
 * Interpolates ("morphs") between two planar shapes using a Hausdorff-distance
 * construction.
 * <p>
 * <b>Intuition:</b> Treat {@code A} and {@code B} as two ink blots. To obtain
 * an in-between shape, inflate {@code A} a little and inflate {@code B} a
 * little (by complementary amounts), then keep only the region where these two
 * inflated blots overlap. When {@code α=0} the overlap is essentially
 * {@code A}; when {@code α=1} it is essentially {@code B}; intermediate
 * {@code α} values yield intermediate shapes.
 * <p>
 * More formally, computes {@code S_α(A,B) = (A ⊕ D_{α d}) ∩ (B ⊕ D_{(1-α) d})},
 * where {@code d} is the (undirected) Hausdorff distance between {@code A} and
 * {@code B}, {@code ⊕} is the Minkowski sum, and {@code D_r} is a disk of
 * radius {@code r}. When {@code d} is the true Hausdorff distance, this is the
 * "maximal" Hausdorff morph from: <br>
 * Marc van Kreveld, Tillmann Miltzow, Tim Ophelders, Willem Sonke, Jordi L.
 * Vermeulen, <i>Between Shapes, Using the Hausdorff Distance</i>.
 * 
 * @author Michael Carleton
 */
public final class HausdorffInterpolator {

	private HausdorffInterpolator() {
	}

	/**
	 * Interpolates between two shapes using the Hausdorff morph: {@code S_α(A,B) =
	 * (A ⊕ D_{α d}) ∩ (B ⊕ D_{(1-α) d})}.
	 * <p>
	 * If {@code d} is the true Hausdorff distance, then the resulting shape is the
	 * maximal set that is at Hausdorff distance {@code α d} from {@code A} and
	 * {@code (1-α) d} from {@code B}.
	 *
	 * @param a                 input geometry {@code A} (typically an area
	 *                          geometry)
	 * @param b                 input geometry {@code B} (typically an area
	 *                          geometry)
	 * @param alpha             interpolation parameter in {@code [0,1]}
	 *                          ({@code 0 -> A}, {@code 1 -> B})
	 * @param hausdorffDistance {@code d} in coordinate units; if {@code 0}, returns
	 *                          {@code A} (fixed)
	 * @param quadSegs          number of quadrant segments used to approximate
	 *                          round buffers (higher is smoother/slower)
	 * @return the interpolated geometry {@code S_α(A,B)} (fixed for robustness)
	 * @throws NullPointerException     if {@code a} or {@code b} is {@code null}
	 * @throws IllegalArgumentException if {@code alpha} is not in {@code [0,1]} or
	 *                                  {@code hausdorffDistance < 0}
	 */
	public static Geometry interpolate(Geometry a, Geometry b, double alpha, double hausdorffDistance, int quadSegs) {
		Objects.requireNonNull(a, "a");
		Objects.requireNonNull(b, "b");
		if (alpha < 0.0 || alpha > 1.0) {
			throw new IllegalArgumentException("alpha must be in [0,1]");
		}
		if (hausdorffDistance < 0.0) {
			throw new IllegalArgumentException("hausdorffDistance must be >= 0");
		}

		// Normalise trivial cases
		if (alpha == 0.0) {
			return a;
		}
		if (alpha == 1.0) {
			return b;
		}
		if (hausdorffDistance == 0.0) {
			return a; // A and B coincide (or distance not meaningful)
		}
		Geometry aFix = GeometryFixer.fix(a);
		Geometry bFix = GeometryFixer.fix(b);

		double rA = alpha * hausdorffDistance;
		double rB = (1.0 - alpha) * hausdorffDistance;

		BufferParameters bp = new BufferParameters();
		bp.setQuadrantSegments(quadSegs);
		bp.setEndCapStyle(BufferParameters.CAP_ROUND);
		bp.setJoinStyle(BufferParameters.JOIN_ROUND);

		Geometry aOff = BufferOp.bufferOp(aFix, rA, bp);
		Geometry bOff = BufferOp.bufferOp(bFix, rB, bp);

		Geometry s = aOff.intersection(bOff);

		return s;
	}

	/**
	 * Convenience overload that first estimates {@code d = d_H(A,B)}
	 * (approximately) using {@link DiscreteHausdorffDistance}, then calls
	 * {@link #interpolate(Geometry, Geometry, double, double, int)}.
	 * <p>
	 * The estimate is sampling-based: smaller {@code densifyFraction} increases
	 * sampling density along edges (slower, typically closer to the true Hausdorff
	 * distance).
	 *
	 * @param a                input geometry {@code A}
	 * @param b                input geometry {@code B}
	 * @param alpha            interpolation parameter in {@code [0,1]}
	 * @param maxSegmentLength if > 0, densifies geometry segments so consecutive
	 *                         sample points are at most this far apart (in geometry
	 *                         units). If <= 0, uses only the geometry's existing
	 *                         vertices.
	 * @param quadSegs         buffer roundness accuracy (quadrant segments)
	 * @return the interpolated geometry using the estimated Hausdorff distance
	 * @throws IllegalArgumentException if {@code densifyFraction} is not in
	 *                                  {@code (0,1]}
	 */
	public static Geometry interpolateUsingEstimatedHausdorff(Geometry a, Geometry b, double alpha, double maxSegmentLength, int quadSegs) {
		double d = estimateHausdorffDistance(a, b, maxSegmentLength);
		return interpolate(a, b, alpha, d, quadSegs);
	}

	/**
	 * Estimates the (undirected) Hausdorff distance between {@code a} and {@code b}
	 * using JTS {@link DiscreteHausdorffDistance}.
	 *
	 * @param a                input geometry {@code A}
	 * @param b                input geometry {@code B}
	 * @param maxSegmentLength if > 0, densifies geometry segments so consecutive
	 *                         sample points are at most this far apart (in geometry
	 *                         units). If <= 0, uses only the geometry's existing
	 *                         vertices.
	 * @return an approximation of {@code d_H(A,B)} in coordinate units
	 */
	public static double estimateHausdorffDistance(Geometry a, Geometry b, double maxSegmentLength) {
		return ProHausdorffDistance.distance(a, b, 0.1, maxSegmentLength); // undirected discrete Hausdorff approximation
	}
}