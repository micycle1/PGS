package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;

/**
 * <p>
 * Morphs between two or more planar rings by evaluating a bivariate
 * Newton&ndash;Thiele interpolation surface <code>R(t, y)</code> over the
 * rings&rsquo; boundary coordinates.
 * </p>
 *
 * <p>
 * <strong>Overview:</strong>
 * </p>
 * <ol>
 * <li>Normalize and resample each input ring to <code>nIn</code> unique points
 * (CCW, uniform arc-length).</li>
 * <li>Align cyclic rotation of all rings to the first ring to fix the
 * <code>y</code>-parameter seam.</li>
 * <li>For each boundary index <code>j</code>, compute Newton divided
 * differences over time nodes <code>t</code> in <code>[0,1]</code>.</li>
 * <li>For each Newton level <code>p</code>, build a Thiele (rational)
 * interpolant over <code>y</code> and evaluate it at the output <code>y</code>
 * grid.</li>
 * <li>Evaluate the resulting Newton polynomial at time <code>t</code> for each
 * output <code>y</code> to form the intermediate ring.</li>
 * </ol>
 *
 * <p>
 * <strong>Notes:</strong>
 * </p>
 * <ul>
 * <li>Works best with 3+ keyframes; for 2 keyframes the time interpolation is
 * linear.</li>
 * <li><code>y ∈ [0,1)</code> excludes the closure point to avoid duplicating
 * the first vertex.</li>
 * <li>Output can be retessellated to any vertex count <code>nOut</code>.</li>
 * </ul>
 *
 * <p>
 * <strong>Reference:</strong> <cite>B. Chen, J. Chen, “Algorithm of Shape
 * Morphing Based on Bivariate Non-linear Interpolation,” ACM CSAE 2020</cite>.
 * <a href=
 * "https://doi.org/10.1145/3424978.3424991">https://doi.org/10.1145/3424978.3424991</a>
 * </p>
 *
 * @author Michael Carleton
 */
public final class NewtonThieleRingMorpher {

	/*-
	 * TODO:
	 * Cut shape with holes (hairline gap) to ensure each geometry has genus=1.
	 * Cuts are discussed in 'Polygon Vertex Set Matching Algorithm for Shapefile Tweening'
	 * Break-up single polygon to interpolate with a multipolygon.
	 * See https://github.com/veltman/openvis/blob/master/README.md
	 * See 'Guaranteed intersection-free polygon morphing'
	 * See CGAl Shape Deformation: https://doc.cgal.org/latest/Barycentric_coordinates_2/index.html#title10
	 * https://homepages.inf.ed.ac.uk/tkomura/cav/presentation14_2018.pdf
	 * RAP C++ : https://github.com/catherinetaylor2/Shape_Interpolation/blob/master/rigid_interp.cpp
	 * and https://github.com/deliagander/ARAPShapeInterpolation
	 */

	private final GeometryFactory gf;

	private final int m; // number of keyframes
	private final int nIn; // input resampled vertices (unique, not closed)
	private final int nOut; // output vertices (unique, not closed)

	private final double[] tNodes; // size m
	private final double[] yIn; // size nIn, in [0,1)
	private final double[] yOut; // size nOut, in [0,1)

	// Precomputed Newton coefficients evaluated at yOut:
	// coeffXAtOut[p][k] = coefficient of Newton level p at output vertex k
	// (parameter yOut[k])
	private final double[][] coeffXAtOut; // [m][nOut]
	private final double[][] coeffYAtOut; // [m][nOut]

	/**
	 * Convenience constructor for N rings.
	 * 
	 * @param rings the rings to interpolate (can be more than two)
	 */
	public NewtonThieleRingMorpher(LinearRing... rings) {
		this(Arrays.asList(rings), maxUniqueVertexCount(rings), // resampleVertices
				maxUniqueVertexCount(rings), // outputVertices
				false // useThieleInY (can be false since nOut==nIn)
		);
	}

	/**
	 * General constructor for multiple keyframes.
	 *
	 * @param keyframes        rings in temporal order
	 * @param resampleVertices vertices to resample each ring to (unique points, no
	 *                         closure point included)
	 * @param outputVertices   output vertex count (unique)
	 * @param useThieleInY     if false and outputVertices==resampleVertices, skips
	 *                         Thiele and uses lattice values
	 */
	public NewtonThieleRingMorpher(List<LinearRing> keyframes, int resampleVertices, int outputVertices, boolean useThieleInY) {
		this.gf = keyframes.get(0).getFactory();
		Objects.requireNonNull(keyframes, "keyframes");
		if (keyframes.size() < 2) {
			throw new IllegalArgumentException("Need at least 2 keyframes");
		}
		if (resampleVertices < 3) {
			throw new IllegalArgumentException("resampleVertices must be >= 3");
		}
		if (outputVertices < 3) {
			throw new IllegalArgumentException("outputVertices must be >= 3");
		}

		this.m = keyframes.size();
		this.nIn = resampleVertices;
		this.nOut = outputVertices;

		this.tNodes = linspace01(m); // [0..1], inclusive endpoints
		this.yIn = linspace01Open(nIn); // [0..1), no 1.0
		this.yOut = linspace01Open(nOut); // [0..1), no 1.0

		// 1) Normalize CCW + resample all rings to nIn vertices (unique)
		List<Coordinate[]> rings = new ArrayList<>(m);
		for (LinearRing ring : keyframes) {
			LinearRing ccw = ensureCCW(ring);
			Coordinate[] pts = resampleClosedRingUniform(ccw, nIn);
			rings.add(pts);
		}

		// 2) Align cyclic rotation of each ring to the first ring (keeps y-parameter
		// seam consistent)
		Coordinate[] ref = rings.get(0);
		for (int i = 1; i < m; i++) {
			Coordinate[] aligned = rotateToBestMatch(ref, rings.get(i));
			rings.set(i, aligned);
		}

		// 3) Build X[i][j], Y[i][j] grid
		double[][] X = new double[m][nIn];
		double[][] Y = new double[m][nIn];
		for (int i = 0; i < m; i++) {
			Coordinate[] pts = rings.get(i);
			for (int j = 0; j < nIn; j++) {
				X[i][j] = pts[j].x;
				Y[i][j] = pts[j].y;
			}
		}

		// 4) Newton divided differences along t for each boundary sample j
		// coeffX[p][j] = Newton coefficient at level p for vertex j (same for Y)
		double[][] coeffX = new double[m][nIn];
		double[][] coeffY = new double[m][nIn];

		double[] tmp = new double[m];
		for (int j = 0; j < nIn; j++) {
			for (int i = 0; i < m; i++) {
				tmp[i] = X[i][j];
			}
			double[] nx = newtonCoefficients(tNodes, tmp);
			for (int p = 0; p < m; p++) {
				coeffX[p][j] = nx[p];
			}

			for (int i = 0; i < m; i++) {
				tmp[i] = Y[i][j];
			}
			double[] ny = newtonCoefficients(tNodes, tmp);
			for (int p = 0; p < m; p++) {
				coeffY[p][j] = ny[p];
			}
		}

		// 5) Thiele recursion along y (boundary parameter) for each Newton level p,
		// then precompute those coefficients on yOut.
		this.coeffXAtOut = new double[m][nOut];
		this.coeffYAtOut = new double[m][nOut];

		boolean canSkipThiele = !useThieleInY && (nOut == nIn);

		if (canSkipThiele) {
			// lattice mode: coeff at yOut[k] is just coeff[*][k]
			for (int p = 0; p < m; p++) {
				System.arraycopy(coeffX[p], 0, coeffXAtOut[p], 0, nOut);
				System.arraycopy(coeffY[p], 0, coeffYAtOut[p], 0, nOut);
			}
		} else {
			// full mode: build Thiele interpolators for each level p and evaluate at yOut
			for (int p = 0; p < m; p++) {
				Thiele1D thX = Thiele1D.fromSamples(yIn, coeffX[p]);
				Thiele1D thY = Thiele1D.fromSamples(yIn, coeffY[p]);
				for (int k = 0; k < nOut; k++) {
					double y = yOut[k];
					coeffXAtOut[p][k] = thX.eval(y);
					coeffYAtOut[p][k] = thY.eval(y);
				}
			}
		}
	}

	/**
	 * Build the intermediate ring at time t in [0,1]. Output is a valid closed
	 * LinearRing (first point repeated at end).
	 */
	public LinearRing interpolate(double t) {
		if (Double.isNaN(t) || Double.isInfinite(t)) {
			throw new IllegalArgumentException("t must be finite");
		}
		t = clamp01(t);

		Coordinate[] out = new Coordinate[nOut + 1];
		for (int k = 0; k < nOut; k++) {
			double x = evalNewtonAtT(t, tNodes, coeffXAtOut, k);
			double y = evalNewtonAtT(t, tNodes, coeffYAtOut, k);
			out[k] = new Coordinate(x, y);
		}
		out[nOut] = new Coordinate(out[0]); // close ring
		return gf.createLinearRing(out);
	}

	/**
	 * Compute Newton interpolation coefficients a[0..m-1] from samples f[0..m-1] at
	 * nodes x[0..m-1]. a[] are the divided differences: P(x)=a0 + a1(x-x0) +
	 * a2(x-x0)(x-x1)+...
	 */
	private static double[] newtonCoefficients(double[] x, double[] f) {
		int n = x.length;
		double[] a = Arrays.copyOf(f, n); // will be overwritten into divided differences
		for (int k = 1; k < n; k++) {
			for (int i = n - 1; i >= k; i--) {
				double den = x[i] - x[i - k];
				if (den == 0.0) {
					throw new IllegalArgumentException("Duplicate x nodes in Newton interpolation");
				}
				a[i] = (a[i] - a[i - 1]) / den;
			}
		}
		// now a[k] is divided difference f[x0..xk]
		return a;
	}

	/**
	 * Evaluate Newton polynomial at time t for a fixed output vertex index k.
	 * coeff[p][k] are the Newton coefficients at level p for this vertex.
	 */
	private static double evalNewtonAtT(double t, double[] tNodes, double[][] coeff, int k) {
		int n = tNodes.length;
		double v = coeff[n - 1][k];
		for (int p = n - 2; p >= 0; p--) {
			v = coeff[p][k] + (t - tNodes[p]) * v;
		}
		return v;
	}

	/**
	 * Thiele 1D rational interpolation (continued fraction). Builds
	 * reciprocal-difference coefficients c[] and evaluates: R(x)=c0 + (x-x0)/(c1 +
	 * (x-x1)/(c2 + ...))
	 */
	private static final class Thiele1D {
		private static final double EPS = 1e-12;

		private final double[] x; // nodes
		private final double[] c; // continued fraction coefficients

		private Thiele1D(double[] x, double[] c) {
			this.x = x;
			this.c = c;
		}

		static Thiele1D fromSamples(double[] x, double[] f) {
			if (x.length != f.length) {
				throw new IllegalArgumentException("x and f length mismatch");
			}
			int n = x.length;
			if (n < 2) {
				throw new IllegalArgumentException("Need at least 2 points for Thiele");
			}

			double[] coeff = new double[n];
			double[] r = Arrays.copyOf(f, n); // r[i] holds reciprocal differences in-place

			coeff[0] = r[0];
			for (int k = 1; k < n; k++) {
				// update r[i] for i>=k using previous stage values
				for (int i = n - 1; i >= k; i--) {
					double denom = r[i] - r[k - 1];
					if (Math.abs(denom) < EPS) {
						// guard singularity; preserve sign
						denom = (denom >= 0.0) ? EPS : -EPS;
					}
					r[i] = (x[i] - x[k - 1]) / denom;
				}
				coeff[k] = r[k];
			}
			return new Thiele1D(Arrays.copyOf(x, n), coeff);
		}

		double eval(double xq) {
			int n = c.length;
			double v = c[n - 1];
			for (int k = n - 2; k >= 0; k--) {
				double denom = v;
				if (Math.abs(denom) < EPS) {
					denom = (denom >= 0.0) ? EPS : -EPS;
				}
				v = c[k] + (xq - x[k]) / denom;
			}
			return v;
		}
	}

	// ------------------------- Geometry helpers -------------------------

	private static LinearRing ensureCCW(LinearRing ring) {
		Coordinate[] c = ring.getCoordinates();
		if (!Orientation.isCCW(c)) {
			// reverse() returns a LineString; for a LinearRing it will still be a
			// LinearRing in practice,
			// but JTS types return LineString so we rebuild below if needed.
			Coordinate[] rev = reverseCoords(c);
			// ensure it's still a proper ring coordinate array:
			if (!rev[0].equals2D(rev[rev.length - 1])) {
				rev = Arrays.copyOf(rev, rev.length + 1);
				rev[rev.length - 1] = new Coordinate(rev[0]);
			}
			return ring.getFactory().createLinearRing(rev);
		}
		return ring;
	}

	/**
	 * Resample a closed ring uniformly by arc length into n unique points (no
	 * closure point).
	 */
	private static Coordinate[] resampleClosedRingUniform(LinearRing ring, int n) {
		Coordinate[] coords = ring.getCoordinates();
		if (coords.length < 4) {
			throw new IllegalArgumentException("Ring must have at least 4 coordinates (including closure)");
		}

		// Drop closure coordinate if present
		int last = coords.length - 1;
		boolean closed = coords[0].equals2D(coords[last]);
		Coordinate[] unique = closed ? Arrays.copyOf(coords, coords.length - 1) : Arrays.copyOf(coords, coords.length);

		unique = removeConsecutiveDuplicates(unique);
		if (unique.length < 3) {
			throw new IllegalArgumentException("Ring has too few unique points after cleanup");
		}

		// Compute perimeter
		double perimeter = 0.0;
		for (int i = 0; i < unique.length; i++) {
			Coordinate a = unique[i];
			Coordinate b = unique[(i + 1) % unique.length];
			perimeter += a.distance(b);
		}
		if (perimeter == 0.0) {
			throw new IllegalArgumentException("Degenerate ring (zero perimeter)");
		}

		// Sample n points at distances k*perimeter/n
		Coordinate[] out = new Coordinate[n];

		int seg = 0;
		double segStartDist = 0.0;
		double segLen = unique[0].distance(unique[1 % unique.length]);

		for (int k = 0; k < n; k++) {
			double target = (k * perimeter) / n;

			while (segStartDist + segLen < target) {
				segStartDist += segLen;
				seg++;
				Coordinate a = unique[seg % unique.length];
				Coordinate b = unique[(seg + 1) % unique.length];
				segLen = a.distance(b);
				if (segLen == 0.0) {
					// skip zero-length segments
					continue;
				}
			}

			Coordinate a = unique[seg % unique.length];
			Coordinate b = unique[(seg + 1) % unique.length];
			double u = (segLen == 0.0) ? 0.0 : (target - segStartDist) / segLen;
			out[k] = lerp(a, b, u);
		}

		return out;
	}

	private static Coordinate[] rotateToBestMatch(Coordinate[] ref, Coordinate[] candidate) {
		if (ref.length != candidate.length) {
			throw new IllegalArgumentException("Length mismatch for rotation alignment");
		}
		int n = ref.length;

		int bestShift = 0;
		double best = Double.POSITIVE_INFINITY;

		for (int shift = 0; shift < n; shift++) {
			double sse = 0.0;
			for (int i = 0; i < n; i++) {
				Coordinate a = ref[i];
				Coordinate b = candidate[(i + shift) % n];
				double dx = a.x - b.x;
				double dy = a.y - b.y;
				sse += dx * dx + dy * dy;
			}
			if (sse < best) {
				best = sse;
				bestShift = shift;
			}
		}

		if (bestShift == 0) {
			return candidate;
		}

		Coordinate[] out = new Coordinate[n];
		for (int i = 0; i < n; i++) {
			out[i] = candidate[(i + bestShift) % n];
		}
		return out;
	}

	private static Coordinate lerp(Coordinate a, Coordinate b, double t) {
		return new Coordinate(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
	}

	private static Coordinate[] removeConsecutiveDuplicates(Coordinate[] pts) {
		ArrayList<Coordinate> out = new ArrayList<>(pts.length);
		Coordinate prev = null;
		for (Coordinate c : pts) {
			if (prev == null || !c.equals2D(prev)) {
				out.add(c);
				prev = c;
			}
		}
		// also avoid last==first in unique list
		if (out.size() >= 2 && out.get(0).equals2D(out.get(out.size() - 1))) {
			out.remove(out.size() - 1);
		}
		return out.toArray(new Coordinate[0]);
	}

	private static Coordinate[] reverseCoords(Coordinate[] c) {
		Coordinate[] r = new Coordinate[c.length];
		for (int i = 0; i < c.length; i++) {
			r[i] = c[c.length - 1 - i];
		}
		return r;
	}

	private static double clamp01(double t) {
		if (t <= 0.0) {
			return 0.0;
		}
		if (t >= 1.0) {
			return 1.0;
		}
		return t;
	}

	/**
	 * tNodes in [0,1], inclusive endpoints: size m.
	 */
	private static double[] linspace01(int m) {
		if (m == 1) {
			return new double[] { 0.0 };
		}
		double[] x = new double[m];
		for (int i = 0; i < m; i++) {
			x[i] = i / (double) (m - 1);
		}
		return x;
	}

	/**
	 * y nodes in [0,1), excludes 1.0 to avoid duplicating ring closure point.
	 */
	private static double[] linspace01Open(int n) {
		double[] x = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = i / (double) n; // last is (n-1)/n < 1
		}
		return x;
	}

	private static int maxUniqueVertexCount(LinearRing... rings) {
		int max = 3;
		if (rings == null || rings.length == 0) {
			return max;
		}
		for (LinearRing r : rings) {
			if (r == null) {
				continue;
			}
			int c = uniqueVertexCount(r);
			if (c > max) {
				max = c;
			}
		}
		return max;
	}

	private static int uniqueVertexCount(LinearRing ring) {
		Coordinate[] c = ring.getCoordinates();
		if (c.length == 0) {
			return 0;
		}

		// drop closure if present
		int len = c.length;
		if (len >= 2 && c[0].equals2D(c[len - 1])) {
			len--;
		}

		// remove consecutive duplicates (same logic as resampling pre-clean)
		int count = 0;
		Coordinate prev = null;
		for (int i = 0; i < len; i++) {
			Coordinate cur = c[i];
			if (prev == null || !cur.equals2D(prev)) {
				count++;
				prev = cur;
			}
		}

		// if last equals first after dedupe, drop it
		if (count >= 2) {
			// cheap check using original endpoints; good enough for sizing N
			if (c[0].equals2D(c[len - 1])) {
				count--;
			}
		}

		// ensure minimum sensible ring size
		return Math.max(count, 3);
	}
}