package org.geotools.process.vector;

/*
 *    GeoTools - The Open Source Java GIS Toolkit
 *    http://geotools.org
 *
 *    (C) 2020, Open Source Geospatial Foundation (OSGeo)
 *
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation;
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 */

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.process.ProcessException;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.util.ProgressListener;


public class ContourProcess implements VectorProcess {

	/**
	 * 
	 * @param features         Input point feature collection
	 * @param propertyName     PropertyName to be contoured
	 * @param levels           Values of levels at which to generate contours
	 * @param interval         Interval between contour values (ignored if levels
	 *                         parameter is supplied)
	 * @param simplify         Indicates whether contour lines are simplified
	 * @param smooth           Indicates whether contour lines are smoothed using
	 *                         Bezier smoothing
	 * @param progressListener
	 * @return The contours of the input features
	 * @throws ProcessException
	 */
	public SimpleFeatureCollection execute(SimpleFeatureCollection features, String propertyName, double[] levels,
			Double interval, Boolean simplify, Boolean smooth, ProgressListener progressListener)
			throws ProcessException {

		Contours contours = new Contours();
		if (smooth != null) {
			contours.setSmooth(smooth.booleanValue());
		}
		if (simplify != null) {
			contours.setSimplify(simplify.booleanValue());
		}
		if (progressListener != null) {
			contours.setProgressListener(progressListener);
		}
		if (levels.length > 0) {
			contours.setLevels(levels);
		} else {
			double min = Double.POSITIVE_INFINITY;
			double max = Double.NEGATIVE_INFINITY;
			try (SimpleFeatureIterator it = (SimpleFeatureIterator) features.features()) {
				while (it.hasNext()) {
					SimpleFeature feature = (SimpleFeature) it.next();
					double value = (Double) feature.getAttribute(propertyName);
					min = Math.min(min, value);
					max = Math.max(max, value);
				}
			}
			int nSteps = (int) Math.ceil((max - min) / interval);
			double[] l = new double[nSteps];
			for (int i = 0; i < nSteps; i++) {
				l[i] = i * interval;
			}
			contours.setLevels(l);
		}

		SimpleFeatureCollection result = contours.contour(features, propertyName);

		return result;
	}
}