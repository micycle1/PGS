package micycle.pgs.commons;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;

import net.jafama.FastMath;

/**
 * Convexity Rule for Shape Decomposition Based on Discrete Contour Evolution.
 * 
 * @author Diego Catalano
 * @author Michael Carleton
 */
public class DiscreteCurveEvolution {

	// from github.com/DiegoCatalano/Catalano-Framework/

	private int vertices;

	public DiscreteCurveEvolution() {
		this(20);
	}

	/**
	 * 
	 * @param vertices the number of vertices to preserve
	 */
	public DiscreteCurveEvolution(int vertices) {
		this.vertices = vertices;
	}

	public Coordinate[] process(LineString lineString) {
		return process(lineString.getCoordinates());
	}

	public Coordinate[] process(Coordinate[] coords) {
		if (vertices > coords.length) {
			throw new IllegalArgumentException("Number of points left must be higher than number of the shape.");
		}

		List<Complex> complex = new ArrayList<>();
		for (int i = 0; i < coords.length; i++) {
			complex.add(new Complex(coords[i].x, coords[i].y));
		}

		for (int i = 0; i < coords.length - vertices; i++) {
			double[] angleMeasure = angle(complex);
			int index = minIndex(angleMeasure);
			complex.remove(index);
		}

		Coordinate[] newShape = new Coordinate[complex.size()];

		for (int i = 0; i < complex.size(); i++) {
			newShape[i] = new Coordinate(complex.get(i).getReal(), complex.get(i).getImaginary());
		}

		return newShape;
	}

	private double[] angle(List<Complex> z) {
		int n = z.size();
		double max = -Double.MAX_VALUE;

		double[] his = new double[n];
		for (int j = 1; j < n - 1; j++) {
			Complex c = z.get(j - 1).subtract(z.get(j + 1));
			double lm = magnitude(c);

			c = z.get(j).subtract(z.get(j + 1));
			double lr = magnitude(c);

			c = z.get(j - 1).subtract(z.get(j));
			double ll = magnitude(c);

			double alpha = FastMath.acos((lr * lr + ll * ll - lm * lm) / (2 * lr * ll));

			// turning angle (0-180)
			double a = 180 - alpha * 180 / Math.PI;

			// the original relevance measure
			his[j] = a * lr * ll / (lr + ll);

			if (his[j] > max) {
				max = his[j];
			}
		}

		his[0] = max;
		his[n - 1] = max;

		return his;
	}

	private static double magnitude(Complex c) {
		return Math.sqrt(c.getReal() * c.getReal() + c.getImaginary() * c.getImaginary());
	}

	private static int minIndex(double[] matrix) {
		int index = 0;
		double min = Double.MAX_VALUE;
		for (int i = 0; i < matrix.length; i++) {
			double currentValue = Math.min(min, matrix[i]);
			if (currentValue < min) {
				min = currentValue;
				index = i;
			}
		}
		return index;
	}

}