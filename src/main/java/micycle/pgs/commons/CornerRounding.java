package micycle.pgs.commons;

import static java.lang.Math.PI;

import java.util.ArrayList;
import java.util.List;

import micycle.pgs.color.Colors;
import net.jafama.FastMath;
import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * A utility class for rounding the corners of a polygon represented as a
 * {@link PShape}.
 * <p>
 * The implementation is based on the algorithm described in the ObservableHQ
 * notebook: <a href=
 * "https://observablehq.com/@daformat/rounding-polygon-corners">Rounding
 * Polygon Corners</a>. It calculates the necessary arc points for each corner
 * and constructs a new {@link PShape} with rounded corners.
 * <p>
 * The class supports different rounding styles, such as strictly geometric
 * rounding, natural rounding, and freehand-style rounding, through the
 * {@link RoundingStyle} enum.
 *
 * @author Mathieu Jouhet
 * @author Michael Carleton
 */
public class CornerRounding {

	/**
	 * An enum representing the rounding style for the corners of the polygon.
	 * <ul>
	 * <li>{@link #CIRCLE}: Strictly geometric rounding, using perfect circular
	 * arcs.</li>
	 * <li>{@link #APPROX}: Natural rounding, approximating circular arcs with
	 * Bézier curves.</li>
	 * <li>{@link #HAND}: Freehand-style rounding, providing a more organic
	 * look.</li>
	 * </ul>
	 */
	public enum RoundingStyle {
		CIRCLE, // Strictly geometric rounding
		APPROX, // Natural rounding
		HAND // Freehand-style rounding
	}

	/**
	 * Rounds the corners of a given {@link PShape} using the specified radius and
	 * rounding style. The method generates circular arcs for each corner and
	 * constructs a new {@link PShape} with the rounded corners.
	 *
	 * @param original The original {@link PShape} whose corners are to be rounded.
	 * @param radius   The radius of the circular arc used to round each corner.
	 *                 This determines how much a circle of the given radius "cuts
	 *                 into" the corner. The effective radius is bounded by the
	 *                 lengths of the edges forming the corner: If the radius is
	 *                 larger than half the length of either edge, it is clamped to
	 *                 the smaller of the two half-lengths to prevent overlapping or
	 *                 invalid geometry.
	 * @param style    The rounding style to apply, as defined in the
	 *                 {@link RoundingStyle} enum.
	 * @return A new {@link PShape} with rounded corners. If the input shape is null
	 *         or the radius is zero, the original shape is returned unchanged.
	 */
	public static PShape roundCorners(PShape original, double radius, RoundingStyle style) {
		if (original == null || radius == 0) {
			return original;
		}

		List<PVector> vertices = extractAndFilterVertices(original);

		int numVertices = vertices.size();
		if (numVertices == 0) {
			return original;
		}

		List<CornerData> corners = new ArrayList<>();

		for (int i = 0; i < numVertices; i++) {
			PVector c1 = vertices.get(i);
			PVector c2 = vertices.get((i + 1) % numVertices);
			PVector c3 = vertices.get((i + 2) % numVertices);

			PVector vC1c2 = PVector.sub(c1, c2);
			PVector vC3c2 = PVector.sub(c3, c2);

			float dx1 = vC1c2.x;
			float dy1 = vC1c2.y;
			float mag1 = vC1c2.mag();
			if (mag1 == 0) {
				continue;
			}
			float unitX1 = dx1 / mag1;
			float unitY1 = dy1 / mag1;

			float dx3 = vC3c2.x;
			float dy3 = vC3c2.y;
			float mag3 = vC3c2.mag();
			if (mag3 == 0) {
				continue;
			}
			float unitX3 = dx3 / mag3;
			float unitY3 = dy3 / mag3;

			float cross = dx1 * dy3 - dy1 * dx3;
			float dot = dx1 * dx3 + dy1 * dy3;
			final float angle = abs(atan2(cross, dot)); // == angleBetween(vC1c2, vC3c2)

			float cornerLength = min((float) radius, mag1 / 2, mag3 / 2);
			if (cornerLength <= 0) {
				continue;
			}

			float a = cornerLength * tan(angle / 2);

			float idealCPDistance = 0;
			idealCPDistance = switch (style) {
				case CIRCLE -> {
					double numPointsCircle = (2 * PI) / (PI - angle);
					yield (4.0f / 3) * tan(PI / (2 * numPointsCircle)) * a;
				}
				case APPROX -> {
					float multiplier = angle < PI / 2 ? 1 + cos(angle) : 2 - sin(angle);
					yield (4.0f / 3) * tan(PI / (2 * ((2 * PI) / angle))) * cornerLength * multiplier;
				}
				case HAND -> (4.0f / 3) * tan(PI / (2 * ((2 * PI) / angle))) * cornerLength * (2 + sin(angle));
				default -> (4.0f / 3) * cornerLength;
			};

			float cpDistance = cornerLength - idealCPDistance;

			PVector c1CurvePoint = new PVector(c2.x + unitX1 * cornerLength, c2.y + unitY1 * cornerLength);
			PVector c1CP = new PVector(c2.x + unitX1 * cpDistance, c2.y + unitY1 * cpDistance);
			PVector c3CurvePoint = new PVector(c2.x + unitX3 * cornerLength, c2.y + unitY3 * cornerLength);
			PVector c3CP = new PVector(c2.x + unitX3 * cpDistance, c2.y + unitY3 * cpDistance);

			corners.add(new CornerData(c1CurvePoint, c1CP, c3CP, c3CurvePoint));
		}

		if (corners.isEmpty()) {
			return original;
		}

		PVector startPoint = corners.get(corners.size() - 1).c3CurvePoint;
		PShape rounded = new PShape(PShape.PATH);
		rounded.setStroke(Colors.PINK);
		rounded.setStroke(true);
		rounded.setFill(true);
		rounded.setFill(Colors.WHITE);
		rounded.beginShape();
		rounded.vertex(startPoint.x, startPoint.y);

		for (CornerData corner : corners) {
			rounded.vertex(corner.c1CurvePoint.x, corner.c1CurvePoint.y);
			rounded.bezierVertex(corner.c1CP.x, corner.c1CP.y, corner.c3CP.x, corner.c3CP.y, corner.c3CurvePoint.x, corner.c3CurvePoint.y);
		}

		rounded.endShape(PConstants.CLOSE);
		return rounded;
	}

	private static List<PVector> extractAndFilterVertices(PShape shape) {
		List<PVector> vertices = new ArrayList<>();
		int vertexCount = shape.getVertexCount();
		if (vertexCount > 0) {
			PVector firstVertex = new PVector(shape.getVertexX(0), shape.getVertexY(0));
			vertices.add(firstVertex);
			PVector prevVertex = firstVertex;

			for (int i = 1; i < vertexCount; i++) {
				PVector currentVertex = new PVector(shape.getVertexX(i), shape.getVertexY(i));
				if (currentVertex.x != prevVertex.x || currentVertex.y != prevVertex.y) {
					vertices.add(currentVertex);
					prevVertex = currentVertex;
				}
			}

			// Check if last vertex is the same as the first (for closed shapes)
			if (vertexCount > 1 && firstVertex.x == prevVertex.x && firstVertex.y == prevVertex.y) {
				vertices.remove(vertices.size() - 1);
			}
		}
		return vertices;
	}

	private record CornerData(PVector c1CurvePoint, PVector c1CP, PVector c3CP, PVector c3CurvePoint) {
	}

	private static float min(float a, float b, float c) {
		return Math.min(a, Math.min(b, c));
	}

	private static float abs(float val) {
		return Math.abs(val);
	}

	private static float atan2(float y, float x) {
		return (float) FastMath.atan2(y, x);
	}

	private static float tan(double angle) {
		return (float) FastMath.tan(angle);
	}

	private static float cos(double angle) {
		return (float) FastMath.cos(angle);
	}

	private static float sin(double angle) {
		return (float) FastMath.sin(angle);
	}
}