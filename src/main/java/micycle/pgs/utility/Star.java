package micycle.pgs.utility;

/*
 * www.javagl.de - Geom - Geometry utilities
 *
 * Copyright (c) 2013-2016 Marco Hutter - http://www.javagl.de
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import java.util.ArrayList;
import java.util.List;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Methods to create star shapes
 * 
 * @author Marco Hutter
 * @author Processing integration by Michael Carleton
 */
public class Star {

	// from
	// https://github.com/javagl/Geom/blob/master/src/main/java/de/javagl/geom/Stars.java

	/**
	 * Epsilon for floating point computations
	 */
	private static final float EPSILON = 1e-8f;

	/**
	 * Create a star shape from the given parameters.
	 * 
	 * @param centerX     The x coordinate of the center
	 * @param centerY     The y coordinate of the center
	 * @param innerRadius The inner radius of the star
	 * @param outerRadius The outer radius of the star
	 * @param numRays     The number of rays that the star should have
	 * @return The star shape
	 */
	public static PShape createStarShape(double centerX, double centerY, double innerRadius, double outerRadius,
			int numRays, double roundness) {
		return createStarShape(centerX, centerY, innerRadius, outerRadius, numRays, 0.5 * Math.PI / numRays, roundness,
				roundness);
	}

	/**
	 * Create a star shape from the given parameters.
	 * 
	 * @param centerX        The x coordinate of the center
	 * @param centerY        The y coordinate of the center
	 * @param innerRadius    The inner radius of the star
	 * @param outerRadius    The outer radius of the star
	 * @param numRays        The number of rays that the star should have
	 * @param startAngleRad  The angle, in radians, where the first ray should start
	 *                       (referring to the positive x-axis)
	 * @param innerRoundness A roundness value between 0.0 and 1.0, for the inner
	 *                       corners of the star.
	 * @param outerRoundness A roundness value between 0.0 and 1.0, for the outer
	 *                       corners (ray tips) of the star.
	 * @return The star shape
	 */
	private static PShape createStarShape(double centerX, double centerY, double innerRadius, double outerRadius,
			int numRays, double startAngleRad, double innerRoundness, double outerRoundness) {
		if (numRays < 2) {
			throw new IllegalArgumentException("The number of rays must be at least 2, but is " + numRays);
		}

		// Create the list containing the inner and outer tip points
		List<PVector> points = new ArrayList<>();
		double deltaAngleRad = Math.PI / numRays;
		for (int i = 0; i < numRays * 2; i++) {
			double angleRad = startAngleRad + i * deltaAngleRad;
			double ca = Math.cos(angleRad);
			double sa = Math.sin(angleRad);
			double relX = ca;
			double relY = sa;
			if ((i & 1) == 0) {
				relX *= outerRadius;
				relY *= outerRadius;
			} else {
				relX *= innerRadius;
				relY *= innerRadius;
			}
			PVector p = new PVector((float) (centerX + relX), (float) (centerY + relY));
			points.add(p);
		}

		// Create a path based on the inner and outer tip points,
		// based on the roundness values
		PVector prevCenter;
		PVector nextCenter;
		PVector step0;
		PVector step1;

		final PShape path2 = new PShape(PShape.PATH);
		path2.beginShape();

		for (int i = 0; i < numRays * 2; i++) {
			// For the current point, compute the previous and next point
			int iPrev = (i - 1 + points.size()) % points.size();
			int iCurr = i;
			int iNext = (i + 1) % points.size();
			PVector pPrev = points.get(iPrev);
			PVector pCurr = points.get(iCurr);
			PVector pNext = points.get(iNext);

			// Compute the center between the previous and the current
			// point, and between the current and the next
			prevCenter = PVector.lerp(pPrev, pCurr, 0.5f);
			nextCenter = PVector.lerp(pCurr, pNext, 0.5f);

			if (i == 0) {
				path2.vertex(prevCenter.x, prevCenter.y);
			}

			// Pick the roundness, depending on whether the current
			// point is an inner or an outer point
			float roundness = (float) innerRoundness;
			if ((i & 1) == 0) {
				roundness = (float) outerRoundness;
			}

			if (Math.abs(roundness) < EPSILON) {
				// For a roundness of 0.0, just walk from the previous center
				// to the current point, and then to the next center
				path2.vertex(pCurr.x, pCurr.y);
				path2.vertex(nextCenter.x, nextCenter.y);

			} else if (Math.abs(roundness - 1.0) < EPSILON) {
				// For a roundness of 1.0, just draw a quad from the
				// previous center to the next, using the current
				// point as the control point
				path2.quadraticVertex(pCurr.x, pCurr.y, nextCenter.x, nextCenter.y);
			} else {
				// Compute interpolated points on the segment between
				// the previous center and the current point, and the
				// current point and the next center, based on the
				// roundness
				step0 = PVector.lerp(prevCenter, pCurr, 1 - roundness);
				step1 = PVector.lerp(pCurr, nextCenter, roundness);

				// Connect the interpolated points using a quad
				path2.vertex(step0.x, step0.y);
				path2.quadraticVertex(pCurr.x, pCurr.y, step1.x, step1.y);
				path2.vertex(nextCenter.x, nextCenter.y);
			}
		}
		path2.endShape(PConstants.CLOSE);

		return path2;
	}

	/**
	 * Private constructor to prevent instantiation
	 */
	private Star() {
		// Private constructor to prevent instantiation
	}

}