/* For concave hull construction based on the criterion described in:
 * JIN-SEO PARK AND SE-JONG OH (2012) 
 * A New Concave Hull Algorithm and Concaveness Measure for n-dimensional Datasets, 
 * JOURNAL OF INFORMATION SCIENCE AND ENGINEERING 28, 587-600
 * 
 * Author: Sheng Zhou (Sheng.Zhou@os.uk)
 * 
 * version 0.4
 * 
 * Date: 2019-01-31
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
package uk.osgb.algorithm.concavehull;

import org.locationtech.jts.geom.Coordinate;

public class TriCheckerPark implements TriangleChecker {
	double n, cm;

	public TriCheckerPark(double N, double CM) {
		n = N;
		cm = CM;
	}

	public boolean removeable(Coordinate coordS, Coordinate coordE, Coordinate coordO) {
		double bl = coordS.distance(coordE);
		double ds = coordS.distance(coordO);
		double de = coordE.distance(coordO);
		double d = ds < de ? ds : de;
		double ratio = bl / d;
		return ratio > n;
	}

}
