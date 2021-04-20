package micycle.pgs.utility;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Minimum Volume Enclosing Ellipsoid
 *
 * @author Peter (pborissow)
 *
 */
public class MinimumBoundingEllipse {

	// https://stackoverflow.com/questions/1768197/bounding-ellipse

	private double[] center;
	private double[][] A;
	private double l1;
	private double l2;
	private double thu;

	/**
	 * @param P An array of points. Each entry in the array contains an x,y
	 *          coordinate. The points must form a closed loop. Consider
	 *          pre-processing a geometry to find its convex hull, and supply thos
	 *          coordinates to compute the MBE.
	 */
	public MinimumBoundingEllipse(double[][] P, double tolerance) {

		// Dimension of the points
		double d = 2;

		// Number of points
		int N = P.length;

		// Rotate the array of points
		P = transpose(P);

		// Add a row of 1s to the 2xN matrix P - so Q is 3xN now.
		// Q = [P;ones(1,N)]
		double[][] Q = merge(P, ones(1, N));

		// Initialize
		int count = 1;
		double err = 1;

		// u is an Nx1 vector where each element is 1/N
		// u = (1/N) * ones(N,1)
		double[] u = new double[N];
		for (int i = 0; i < u.length; i++) {
			u[i] = 1D / N;
		}

		// Khachiyan Algorithm
		while (err > tolerance) {

			// Matrix multiplication:
			// diag(u) : if u is a vector, places the elements of u
			// in the diagonal of an NxN matrix of zeros
			// X = Q*diag(u)*Q'; // Q' - transpose of Q
			double[][] X = multiply(multiply(Q, diag(u)), transpose(Q));

			// inv(X) returns the matrix inverse of X
			// diag(M) when M is a matrix returns the diagonal vector of M
			// M = diag(Q' * inv(X) * Q); // Q' - transpose of Q
			double[] M = diag(multiply(multiply(transpose(Q), inv(X)), Q));

			// Find the value and location of the maximum element in the vector M
			double maximum = max(M);
			int j = find_maximum_value_location(M, maximum);

			// Calculate the step size for the ascent
			double step_size = (maximum - d - 1) / ((d + 1) * (maximum - 1));

			// Calculate the new_u:
			// Take the vector u, and multiply all the elements in it by (1-step_size)
			double[] new_u = multiply(1 - step_size, u);

			// Increment the jth element of new_u by step_size
			new_u[j] = new_u[j] + step_size;

			// Calculate error by taking finding the square root of the SSD
			// between new_u and u
			err = Math.sqrt(ssd(new_u, u));

			// Increment count and replace u
			count = count + 1;
			u = new_u;
		}

		// Compute center point
		// c = P * u
		double[][] c = multiply(P, u);
		center = transpose(c)[0];

		// Put the elements of the vector u into the diagonal of a matrix
		// U with the rest of the elements as 0
		double[][] U = diag(u);

		// Compute the A-matrix
		// A = (1/d) * inv(P * U * P' - (P * u)*(P*u)' );
		double[][] pup = multiply(multiply(P, U), transpose(P));
		double[][] pupu = multiply(multiply(P, u), transpose(multiply(P, u)));
		double[][] pup_pupu = minus(pup, pupu);
		A = multiply(1 / d, inv(pup_pupu));

		// Compute Eigen vectors and values
		// A=inv(A);
		// [Ve,De]=eig(A);
		RealMatrix m = new Array2DRowRealMatrix(inv(A)); // TODO
		EigenDecomposition eig = new EigenDecomposition(m);
//		Eigen eig = new Eigen(inv(A));
		double[][] Ve = eig.getV().getData(); // eigenvalues
		double[][] De = eig.getD().getData(); // right eigenvectors
		reorderEigenVectors(De);
		reorderEigenValues(Ve);

		// v=sqrt(diag(De));
		double[] v = sqrt(diag(De));

		// [l1,Ie] = max(v);
		l1 = max(v);
		int Ie = find_maximum_value_location(v, l1); // off by one from MatLab but I think it's ok here

		// veig=Ve(:,Ie);
		double[] veig = new double[Ve.length];
		for (int i = 0; i < veig.length; i++) {
			veig[i] = Ve[Ie][i];
		}

		// thu=atan2(veig(2),veig(1));
		thu = Math.atan2(veig[1], veig[0]);

		// l2=v(setdiff([1 2],Ie));
		l2 = v[setdiff(new int[] { 0, 1 }, Ie)];
	}

	// **************************************************************************
	// ** getCenter
	// **************************************************************************
	/**
	 * Returns the center point of the ellipse
	 */
	public double[] getCenter() {
		double[] pt = new double[2];
		pt[0] = center[0];
		pt[1] = center[1];
		return pt;
	}

	// **************************************************************************
	// ** getMatrix
	// **************************************************************************
	/**
	 * Returns a matrix containing all the information regarding the shape of the
	 * ellipsoid. To get the radii and orientation of the ellipsoid take the
	 * Singular Value Decomposition of the matrix.
	 */
	public double[][] getMatrix() {
		return A;
	}

	// **************************************************************************
	// ** getBoundingCoordinates
	// **************************************************************************
	/**
	 * Returns a list of coordinates that can be used to render the ellipse.
	 *
	 * @param numPoints The number of points used to represent the ellipse. The
	 *                  higher the number the more dense the ellipse outline, the
	 *                  more accurate the shape.
	 */
	public double[][] getBoundingCoordinates(int numPoints) {

		// tq=linspace(-pi,pi,50);
		double[] tq = linspace(-Math.PI, Math.PI, numPoints);

		// U=[cos(thu) -sin(thu);sin(thu) cos(thu)]*[l1*cos(tq);l2*sin(tq)];
		double[][] U = multiply(
				new double[][] { createVector(Math.cos(thu), -Math.sin(thu)),
						createVector(Math.sin(thu), Math.cos(thu)) },
				new double[][] { multiply(l1, cos(tq)), multiply(l2, sin(tq)) });
		// System.out.println(toString(transpose(U)));

		double[][] coords = transpose(U);
		for (int i = 0; i < coords.length; i++) {
			double x = coords[i][0] + center[0];
			double y = coords[i][1] + center[1];

			coords[i][0] = x;
			coords[i][1] = y;
		}

		return coords;
	}

	// **************************************************************************
	// ** reorderEigenVectors
	// **************************************************************************
	/**
	 * Eigen values generated from Apache Common Math and JAMA are different than
	 * MatLab. The vectors are in the reverse order than expected. This function
	 * will update the array to what we expect to see in MatLab.
	 */
	private void reorderEigenVectors(double[][] De) {
		rotateMatrix(De);
		rotateMatrix(De);
	}

	// **************************************************************************
	// ** reorderEigenValues
	// **************************************************************************
	/**
	 * Eigen values generated from Apache Common Math and JAMA are different than
	 * MatLab. The vectors are in reverse order than expected and with an opposite
	 * sign. This function will update the array to what we expect to see in MatLab.
	 */
	private void reorderEigenValues(double[][] Ve) {
		rotateMatrix(Ve);
		for (int i = 0; i < Ve.length; i++) {
			for (int j = 0; j < Ve[i].length; j++) {
				Ve[i][j] = -Ve[i][j];
			}
		}
	}

	// **************************************************************************
	// ** linspace
	// **************************************************************************
	private double[] linspace(double min, double max, int points) {
		double[] d = new double[points];
		for (int i = 0; i < points; i++) {
			d[i] = min + i * (max - min) / (points - 1);
		}
		return d;
	}

	// **************************************************************************
	// ** ssd
	// **************************************************************************
	/**
	 * Returns the sum-of-square-differences between tow arrays. Takes two vectors
	 * of the same size, creates a new vector by finding the difference between
	 * corresponding elements, squaring each difference and adding them all
	 * together. So if the vectors were: a = [1 2 3] and b = [5 4 6], then: SSD =
	 * (1-5)^2 + (2-4)^2 + (3-6)^2;
	 */
	private double ssd(double[] a, double[] b) {
		double ssd = 0;
		for (int i = 0; i < a.length; i++) {
			final double n = a[i] - b[i];
			ssd += n * n;
		}
		return ssd;
	}

	// **************************************************************************
	// ** ones
	// **************************************************************************
	/**
	 * Creates an array of all ones. For example, ones(2,3) returns a 2-by-3 array
	 * of ones.
	 *
	 * <pre>
	     1 1 1
	     1 1 1
	 * </pre>
	 *
	 * Reference: https://www.mathworks.com/help/matlab/ref/ones.html
	 */
	private double[][] ones(int rows, int cols) {
		double[][] arr = new double[rows][];
		for (int i = 0; i < arr.length; i++) {
			double[] row = new double[cols];
			for (int j = 0; j < row.length; j++) {
				row[j] = 1;
			}
			arr[i] = row;
		}
		return arr;
	}

	// **************************************************************************
	// ** merge
	// **************************************************************************
	/**
	 * Used to combine two arrays into one
	 */
	private double[][] merge(double[][] m1, double[][] m2) {
		int x = 0;
		double[][] out = new double[m1.length + m2.length][];
		for (int i = 0; i < m1.length; i++) {
			out[x] = m1[i];
			x++;
		}
		for (int i = 0; i < m2.length; i++) {
			out[x] = m2[i];
			x++;
		}
		return out;
	}

	// **************************************************************************
	// ** multiply
	// **************************************************************************
	/**
	 * Used to multiply all the values in the vector (arr) by n. This is called
	 * scalar multiplication.
	 */
	private double[] multiply(double n, double[] arr) {
		double[] out = new double[arr.length];
		for (int i = 0; i < arr.length; i++) {
			out[i] = arr[i] * n;
		}
		return out;
	}

	// **************************************************************************
	// ** multiply
	// **************************************************************************
	/**
	 * Used to multiply all the values in the matrix (arr) by n
	 */
	private double[][] multiply(double n, double[][] arr) {
		double[][] out = new double[arr.length][];
		for (int i = 0; i < arr.length; i++) {
			double[] row = arr[i];
			double[] r = new double[row.length];
			for (int j = 0; j < row.length; j++) {
				r[j] = row[j] * n;
			}
			out[i] = r;
		}
		return out;
	}

	// **************************************************************************
	// ** multiply
	// **************************************************************************
	/**
	 * Multiply a matrix with a vector by converting the vector to a matrix
	 */
	private double[][] multiply(double[][] P, double[] u) {
		double[][] m2 = new double[u.length][];
		for (int i = 0; i < m2.length; i++) {
			double[] row = new double[1];
			row[0] = u[i];
			m2[i] = row;
		}
		return multiply(P, m2);
	}

	// **************************************************************************
	// ** multiply
	// **************************************************************************
	/**
	 * Used to multiply two matrices. Credit: https://stackoverflow.com/a/23817780
	 */
	private double[][] multiply(double[][] m1, double[][] m2) {
		int m1ColLength = m1[0].length; // m1 columns length
		int m2RowLength = m2.length; // m2 rows length
		if (m1ColLength != m2RowLength) {
			return null; // matrix multiplication is not possible
		}
		int mRRowLength = m1.length; // m result rows length
		int mRColLength = m2[0].length; // m result columns length
		double[][] mResult = new double[mRRowLength][mRColLength];
		for (int i = 0; i < mRRowLength; i++) { // rows from m1
			for (int j = 0; j < mRColLength; j++) { // columns from m2
				for (int k = 0; k < m1ColLength; k++) { // columns from m1
					mResult[i][j] += m1[i][k] * m2[k][j];
				}
			}
		}
		return mResult;
	}

	// **************************************************************************
	// ** diag
	// **************************************************************************
	/**
	 * Returns a matrix for a given vector. The values in the vector will appear
	 * diagonally in the output. Reference:
	 * https://www.mathworks.com/help/matlab/ref/diag.html
	 */
	private double[][] diag(double[] arr) {
		double[][] out = new double[arr.length][];
		for (int i = 0; i < arr.length; i++) {
			double[] row = new double[arr.length];
			for (int j = 0; j < row.length; j++) {
				if (j == i) {
					row[j] = arr[i];
				} else {
					row[j] = 0;
				}
			}
			out[i] = row;
		}
		return out;
	}

	// **************************************************************************
	// ** diag
	// **************************************************************************
	/**
	 * Returns a vector representing values that appear diagonally in the given
	 * matrix. Reference: https://www.mathworks.com/help/matlab/ref/diag.html
	 */
	private double[] diag(double[][] arr) {
		double[] out = new double[arr.length];
		for (int i = 0; i < arr.length; i++) {
			out[i] = arr[i][i];
		}
		return out;
	}

	// **************************************************************************
	// ** transpose
	// **************************************************************************
	/**
	 * Interchanges the row and column index for each element Reference:
	 * https://www.mathworks.com/help/matlab/ref/transpose.html
	 */
	private double[][] transpose(double[][] arr) {
		int rows = arr.length;
		int cols = arr[0].length;

		double[][] out = new double[cols][rows];
		for (int x = 0; x < cols; x++) {
			for (int y = 0; y < rows; y++) {
				out[x][y] = arr[y][x];
			}
		}

		return out;
	}

	// **************************************************************************
	// ** inv
	// **************************************************************************
	/**
	 * Returns the inverse of a matrix. Relies on 2 different implementations. The
	 * first implementation is more accurate (passes inverse check) but has the
	 * potential to fail. If so, falls back to second method that relies on
	 * partial-pivoting Gaussian elimination. Reference:
	 * https://www.mathworks.com/help/matlab/ref/inv.html
	 */
	private double[][] inv(double[][] matrix) {
		try {
			return inv1(matrix);
		} catch (Exception e) {
			try {
				return inv2(matrix);
			} catch (Exception ex) {
				throw new RuntimeException(ex);
			}
		}
	}

	// **************************************************************************
	// ** inv1
	// **************************************************************************
	/**
	 * Returns the inverse of a matrix. This implementation passes inverse check so
	 * I think it's valid but it has a tendency to fail. For example, the following
	 * matrix fails with a ArrayIndexOutOfBoundsException in the determinant method.
	 *
	 * <pre>
	     1171.18     658.33
	      658.33    1039.55
	 * </pre>
	 *
	 * Credit: https://github.com/rchen8/Algorithms/blob/master/Matrix.java
	 */
	private double[][] inv1(double[][] matrix) {
		double[][] inverse = new double[matrix.length][matrix.length];

		// minors and cofactors
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				inverse[i][j] = Math.pow(-1, i + j) * determinant(minor(matrix, i, j));
			}
		}

		// adjugate and determinant
		double det = 1.0 / determinant(matrix);
		for (int i = 0; i < inverse.length; i++) {
			for (int j = 0; j <= i; j++) {
				double temp = inverse[i][j];
				inverse[i][j] = inverse[j][i] * det;
				inverse[j][i] = temp * det;
			}
		}

		return inverse;
	}

	private static double determinant(double[][] matrix) {
		if (matrix.length != matrix[0].length) {
			throw new IllegalStateException("invalid dimensions");
		}

		if (matrix.length == 2) {
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		}

		double det = 0;
		for (int i = 0; i < matrix[0].length; i++) {
			det += Math.pow(-1, i) * matrix[0][i] * determinant(minor(matrix, 0, i));
		}
		return det;
	}

	private static double[][] minor(double[][] matrix, int row, int column) {
		double[][] minor = new double[matrix.length - 1][matrix.length - 1];

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; i != row && j < matrix[i].length; j++) {
				if (j != column) {
					minor[i < row ? i : i - 1][j < column ? j : j - 1] = matrix[i][j];
				}
			}
		}
		return minor;
	}

	// **************************************************************************
	// ** inv2
	// **************************************************************************
	/**
	 * Returns the inverse of a matrix. This implementation successfully executes
	 * but does not pass the inverse check. Credit:
	 * https://www.sanfoundry.com/java-program-find-inverse-matrix/
	 */
	public static double[][] inv2(double a[][]) {

		int n = a.length;
		double x[][] = new double[n][n];
		double b[][] = new double[n][n];
		int index[] = new int[n];

		for (int i = 0; i < n; ++i) {
			b[i][i] = 1;
		}

		// Transform the matrix into an upper triangle
		gaussian(a, index);

		// Update the matrix b[i][j] with the ratios stored
		for (int i = 0; i < n - 1; ++i) {
			for (int j = i + 1; j < n; ++j) {
				for (int k = 0; k < n; ++k) {
					b[index[j]][k] -= a[index[j]][i] * b[index[i]][k];
				}
			}
		}

		// Perform backward substitutions
		for (int i = 0; i < n; ++i) {
			x[n - 1][i] = b[index[n - 1]][i] / a[index[n - 1]][n - 1];
			for (int j = n - 2; j >= 0; --j) {
				x[j][i] = b[index[j]][i];
				for (int k = j + 1; k < n; ++k) {
					x[j][i] -= a[index[j]][k] * x[k][i];
				}
				x[j][i] /= a[index[j]][j];
			}
		}

		return x;
	}

	// Method to carry out the partial-pivoting Gaussian
	// elimination. Here index[] stores pivoting order.
	public static void gaussian(double a[][], int index[]) {

		int n = index.length;
		double c[] = new double[n];

		// Initialize the index
		for (int i = 0; i < n; ++i) {
			index[i] = i;
		}

		// Find the rescaling factors, one from each row
		for (int i = 0; i < n; ++i) {
			double c1 = 0;
			for (int j = 0; j < n; ++j) {
				double c0 = Math.abs(a[i][j]);
				if (c0 > c1) {
					c1 = c0;
				}
			}
			c[i] = c1;
		}

		// Search the pivoting element from each column
		int k = 0;
		for (int j = 0; j < n - 1; ++j) {
			double pi1 = 0;
			for (int i = j; i < n; ++i) {
				double pi0 = Math.abs(a[index[i]][j]);
				pi0 /= c[index[i]];
				if (pi0 > pi1) {
					pi1 = pi0;
					k = i;
				}
			}

			// Interchange rows according to the pivoting order
			int itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;
			for (int i = j + 1; i < n; ++i) {
				double pj = a[index[i]][j] / a[index[j]][j];

				// Record pivoting ratios below the diagonal
				a[index[i]][j] = pj;

				// Modify other elements accordingly
				for (int l = j + 1; l < n; ++l) {
					a[index[i]][l] -= pj * a[index[j]][l];
				}
			}
		}
	}

	// **************************************************************************
	// ** max
	// **************************************************************************
	/**
	 * Returns the max value in a vector
	 */
	private double max(double[] arr) {
		double max = arr[0];
		for (double d : arr) {
			max = Math.max(d, max);
		}
		return max;
	}

	// **************************************************************************
	// ** find_maximum_value_location
	// **************************************************************************
	/**
	 * Returns the index of the max value in a vector
	 */
	private int find_maximum_value_location(double[] arr, double max) {
		for (int i = 0; i < arr.length; i++) {
			if (arr[i] == max) {
				return i;
			}
		}
		return 0;
	}

	// **************************************************************************
	// ** minus
	// **************************************************************************
	/**
	 * Used to subtract array B from array A and returns the result Reference:
	 * https://www.mathworks.com/help/matlab/ref/minus.html
	 */
	private double[][] minus(double[][] a, double[][] b) {
		double[][] out = new double[a.length][];
		for (int i = 0; i < out.length; i++) {
			double[] row = new double[a[i].length];
			for (int j = 0; j < row.length; j++) {
				row[j] = a[i][j] - b[i][j];
			}
			out[i] = row;
		}
		return out;
	}

	// **************************************************************************
	// ** sqrt
	// **************************************************************************
	/**
	 * Returns the square root of each element in a vector Reference:
	 * https://www.mathworks.com/help/matlab/ref/sqrt.html
	 */
	private double[] sqrt(double[] arr) {
		double[] out = new double[arr.length];
		for (int i = 0; i < out.length; i++) {
			out[i] = Math.sqrt(arr[i]);
		}
		return out;
	}

	private double[] cos(double[] arr) {
		double[] out = new double[arr.length];
		for (int i = 0; i < out.length; i++) {
			out[i] = Math.cos(arr[i]);
		}
		return out;
	}

	private double[] sin(double[] arr) {
		double[] out = new double[arr.length];
		for (int i = 0; i < out.length; i++) {
			out[i] = Math.sin(arr[i]);
		}
		return out;
	}

	// **************************************************************************
	// ** setdiff
	// **************************************************************************
	/**
	 * Partial implementation of setdiff
	 */
	private int setdiff(int[] arr, int x) {
		for (int i : arr) {
			if (i != x) {
				return i;
			}
		}
		return 0; // ?
	}

	// **************************************************************************
	// ** rotateMatrix
	// **************************************************************************
	private void rotateMatrix(double mat[][]) {
		int N = mat[0].length;

		// Consider all squares one by one
		for (int x = 0; x < N / 2; x++) {
			// Consider elements in group of 4 in
			// current square
			for (int y = x; y < N - x - 1; y++) {
				// store current cell in temp variable
				double temp = mat[x][y];

				// move values from right to top
				mat[x][y] = mat[y][N - 1 - x];

				// move values from bottom to right
				mat[y][N - 1 - x] = mat[N - 1 - x][N - 1 - y];

				// move values from left to bottom
				mat[N - 1 - x][N - 1 - y] = mat[N - 1 - y][x];

				// assign temp to left
				mat[N - 1 - y][x] = temp;
			}
		}
	}

	// **************************************************************************
	// ** createVector
	// **************************************************************************
	/**
	 * Used to generate a vector for testing purposes
	 */
	private double[] createVector(double... d) {
		double[] arr = new double[d.length];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = d[i];
		}
		return arr;
	}
}