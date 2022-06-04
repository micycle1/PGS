[![](https://jitpack.io/v/micycle1/PGS.svg)](https://jitpack.io/#micycle1/PGS) [![Lines of Code](https://sonarcloud.io/api/project_badges/measure?project=micycle1_PTS&metric=ncloc)](https://sonarcloud.io/dashboard?id=micycle1_PTS)

# Processing Geometry Suite

*Processing Geometry Suite* is a software project that provides easy access to 2D geometric algorithms in the form of a [Processing](https://processing.org/) library.

The focus of the library is on visualisation rather than providing underlying data structures. To this end all methods in the library are static and most of them take in and return [`PShapes`](https://processing.org/reference/PShape.html) or [`PVectors`](https://processing.org/reference/PVector.html).

Docs are hosted via *GitHub Pages* [here](https://micycle1.github.io/PGS/).

## **Overview**

Library functionality is split over the following classes:

* `PGS_CirclePacking`
  * Circle packings of shapes, subject to varying constraints and patterns of tangencies
* `PGS_Coloring`
  * Minimal colorings of meshes (or mesh-like shapes)
* `PGS_Construction`
  * Construct uncommon/interesting 2D primitives
* `PGS_Contour`
  * Methods that produce various contours from shapes: medial axes, straight skeletons, offset curves, etc.
* `PGS_Conversion`
  * Conversion between *Processing* PShapes and *JTS* Geometries 
* `PGS_Meshing`
  * Mesh generation (excluding triangulation) from shapes
* `PGS_Morphology`
  * Methods that affect the geometry or topology of shapes (buffering, simplification, smoothing, etc.)
* `PGS_Optimisation`
  * Solve geometric optimisation problems, such as finding the maximum inscribed circle, or the closest vertex to a coordinate
* `PGS_PointSet`
  * Generates sets of 2D points having a variety of different distributions and constraints
* `PGS_Processing`
  * Methods that process a shape in some way: compute hulls, partition, slice, etc. 
* `PGS_ShapeBoolean`
  * Boolean set-operations for 2D shapes
* `PGS_ShapePredicates`
  * Various shape metrics (area, circularity, etc.) and predicates (*"do these shapes intersect?"*)
* `PGS_Tiling`
  * Tiling, tessellation and subdivision of the plane using periodic or non-periodic geometric shapes.
* `PGS_Transformation`
  * Various geometric and affine transformations that affect vertex coordinates
* `PGS_Triangulation`
  * Delaunay triangulation (constrained and refined) and earcut triangulation of shapes and point sets
* `PGS_Voronoi`
  * Voronoi Diagrams of shapes and point sets

## **Installation**

<details><summary>Processing IDE — Quick</summary>
<p>

Download the latest *PGS.jar* from [releases](https://github.com/micycle1/PGS/releases) and simply drag-and-drop it onto the [Processing IDE](https://processing.org/environment).
</p>
</details>

<details><summary>Processing IDE — Permanently</summary>
<p>

Download the latest *PGS.jar* from [releases](https://github.com/micycle1/PGS/releases) and save it to `Documents\Processing\libraries\PGS\library`.

Result: `Documents\Processing\libraries\PGS\library\PGS.jar`.

(Note the *.jar* and the folder **must** be called `PGS` — rename the .jar if this is not the case).
</p>
</details>

<details><summary>Maven/Gradle</summary>
<p>

PGS is hosted as an artifact for use in Maven or Gradle projects via [Jitpack](https://jitpack.io/#micycle1/PGS) — follow the instructions there (very easy). 
</p>
</details>

## **Examples**

A number of example Processing sketches are provided in [examples](https://github.com/micycle1/PGS/tree/master/examples).

<p float="middle">
  <a href="examples/partitionSmooth"><img src="resources/examples/partitionSmooth.png" alt="" width="24%"/></a>
  <a href="examples/textOffsetCurves"><img src="resources/examples/textOffsetCurves.png" alt="" width="24%"/></a>
  <a href="examples/ripplingTriangles"><img src="resources/examples/ripplingTriangles.png" alt="" width="24%"/></a>
  <a href="examples/contourMap"><img src="resources/examples/contourMap.png" alt="" width="24%"/></a>
</p>

## **Illustrations**

Much of the functionality (but by no means all) is demonstrated below:

## *2D Boolean Operations*

<table>
  <tr>
    <td align="center" valign="center"><b>Union</td>
    <td align="center" valign="center"><b>Intersection</td>
    <td align="center" valign="center"><b>Subtraction</td>
    <td align="center" valign="center"><b>Symmetric Difference</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/boolean/union.gif"></td>
    <td valign="top"><img src="resources/boolean/intersect.gif"></td>
    <td valign="top"><img src="resources/boolean/subtract.gif"></td>
    <td valign="top" ><img src="resources/boolean/symDifference.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Complement</td>
    <td align="center" valign="center"><b>Mesh Union</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/boolean/complement.png"></td>
    <td valign="top" width="25%"><img src="resources/boolean/meshUnion.gif"></td>
  </tr>
</table>

## *Transformation*

<table>
  <tr>
    <td align="center" valign="center" colspan="2"><b>Rotate Around</td>
    <td align="center" valign="center"><b>Translate To</td>
    <td align="center" valign="center"><b>Touch Scale</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/transform/rotateCenter.gif"></td>
    <td valign="top" width="25%"><img src="resources/transform/rotate.gif"></td>
    <td valign="top" width="25%"><img src="resources/transform/translateTo.gif"></td>
    <td valign="top" width="25%"><img src="resources/transform/touchScale.gif"></td>
  </tr>
    <tr>
    <td align="center" valign="center" colspan="2">Rotate a shape around its centroid or an arbitrary point.</td>
    <td align="center" valign="center">Translate a shape such that its centroid matches some position.</td>
    <td align="center" valign="center">Scale one shape such that it touches another.</td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Resize</td>
    <td align="center" valign="center"><b>Homothetic Transformation</td>
    <td align="center" valign="center"><b>Shear</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/transform/resize.gif"></td>
    <td valign="top" width="25%"><img src="resources/transform/homothetic.gif"></td>
    <td valign="top" width="25%"><img src="resources/transform/shear.gif"></td>
  </tr>
    <tr>
    <td align="center" valign="center"></td>
    <td align="center" valign="center">Projection-transform a shape with respect to a fixed point.</td>
    <td align="center" valign="center"></td>
  </tr>
</table>

## *Geometric Predicates & Metrics*

<table>
  <tr>
    <td align="center" valign="center"><b>Intersects</td>
    <td align="center" valign="center"><b>Contains Shape</td>
    <td align="center" valign="center"><b>Contains Point</td>
    <td align="center" valign="center"><b>Containing Cell</td>
  </tr>
  <tr>
    <td valign="top" align="center" width="25%"> <img src="resources/predicate/intersect.gif"><br>Do shapes intersect with each other?</td>
    <td valign="top" align="center" width="25%"><img src="resources/predicate/contains.gif"><br>Does one shape fully contain another?</td>
    <td valign="top" align="center" width="25%"><img src="resources/predicate/containsPoint.gif"><br>For individual points and point sets.</td>
    <td valign="top" align="center" width="25%"><img src="resources/predicate/containingCell.gif"><br>Which cell contains the query point?</td>
</table>


### Metrics
* Length/perimeter
* Width & Height
* Diameter
* Circularity
* Similarity
* Sphericity
* Elongation
* Density
* Holes
* Maximum interior angle
* Is simple?
* Is convex?
* Equal? (structural and topological equivalence)
* Distance
* Area
* Centroid

## *Contour*

<table>
  <tr>
    <td align="center" valign="center" colspan="2"><b>Isolines</td>
    <td align="center" valign="center" colspan="2"><b>Offset Curves</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/contour/isolines.gif"></td>
    <td valign="top"><img src="resources/contour/isolines2.gif"></td>
    <td valign="top"><img src="resources/contour/miteredInterior.gif"></td>
    <td valign="top" ><img src="resources/contour/miteredExterior.gif"></td>
  </tr>
    </tr>
    <tr>
    <td align="center" valign="center" colspan="2">Isolines from intra-shape euclidean distance, or point sets.</td>
    <td align="center" valign="center" colspan="2">Inner and exterior offset curves; based on <i>miter, bevel</i> or <i>round</i> offset styles.</td>
  </tr>

  <tr>
    <td align="center" valign="center" colspan="2"><b>Voronoi Diagram</td>
    <td align="center" valign="center" colspan="2"><b>Circle-site Voronoi Diagram</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/contour/voronoi.gif"></td>
    <td valign="top" width="25%"><img src="resources/contour/voronoi1.png"></td>
    <td valign="top" width="25%"><img src="resources/contour/voronoiCircles.gif"></td>
    <td valign="top" width="25%"><img src="resources/contour/voronoiCircles.png"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Straight Skeleton</td>
    <td align="center" valign="center" colspan="2"><b>Medial Axis</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/contour/straightSkeleton.png"></td>
    <td valign="top"><img src="resources/contour/medialAxis.gif"></td>
    <td valign="top"><img src="resources/contour/medialAxis.png"></td>
  </tr>
    </tr>
    <tr>
    <td align="center" valign="center"></td>
    <td align="center" valign="center" colspan="2">Medial axis transform with feature pruning via distance, area or <i>axial angle</i>. </td>
  </tr>
</table>

## *Morphology*

<table>
  <tr>
    <td align="center" valign="center"><b>Buffer</td>
    <td align="center" valign="center"><b>Erosion-Dilation</td>
    <td align="center" valign="center" colspan="2"><b>Minkowski Addition</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/morphology/buffer.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/erosionDilation.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/minkSum.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/minkDiff.gif"></td>
  </tr>
    <tr>
    <td align="center" valign="center"></td>
    <td align="center" valign="center">A negative followed by a positive buffer (in a single operation).</td>
    <td align="center" valign="center" colspan="2">Minkowski sum and difference (a.k.a buffer one shape using another shape; the examples add a rotating & growing triangle).</td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Smoothing</td>
    <td align="center" valign="center"><b>Gaussian Smoothing</td>
    <td align="center" valign="center" colspan="2"><b>Rounding</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/morphology/smooth.gif"></td>
    <td valign="top"><img src="resources/morphology/gaussianSmooth.gif"></td>
    <td valign="top"><img src="resources/morphology/round.gif"></td>
    <td valign="top" ><img src="resources/morphology/round2.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center" colspan="2"><b>Radial Warp</td>
    <td align="center" valign="center" colspan="2"><b>Field Warp</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/morphology/radialWarp1.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/radialWarp2.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/fieldWarp.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/fieldWarp2.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Simplification</td>
    <td align="center" valign="center" colspan="2"><b>Chaikin Cutting</td>
    <td align="center" valign="center" colspan="2"><b>Interpolation</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/morphology/simplifyVW.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/chaikinCut1.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/chaikinCut2.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/morph.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Variable Buffer</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/morphology/variableBuffer.gif"></td>
  </tr>
</table>

## *Geometry Processing*

<table>
  <tr>
    <td align="center" valign="center" colspan="2"><b>Points on Perimeter</td>
    <td align="center" valign="center"><b>Point on Perimeter</td>
    <td align="center" valign="center"><b>Perimeter Extraction</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/geometry_processing/pointsOnPerimeter.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/pointsOnPerimeter2.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/pointOnPerimeter.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/perimeterExtract.gif"></td>
  </tr>
    <tr>
    <td align="center" valign="center" colspan="2">Find <i>N</i> points (evenly distributed) along the perimeter of a shape, or points every <i>D</i> distance (with optional perpendicular offset).</td>
    <td align="center" valign="center">Find a point some fraction along the perimeter of a shape (with perpendicular offset).</td>
    <td align="center" valign="center"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Splitting</td>
    <td align="center" valign="center"><b>Convex Partitioning</td>
    <td align="center" valign="center"><b>Equal Partitioning</td>
    <td align="center" valign="center"><b>Trapezoid Partitioning</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/geometry_processing/split.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/decompose2.png"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/partition.png"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/trapezoidPartition.gif"></td>
    
  </tr>
  <tr>
    <td align="center" valign="center">Subdivide (recursively) a shape into quadrants.</td>
    <td align="center" valign="center">Partition a shape into convex polygons.</td>
    <td align="center" valign="center">Partition a shape into N equal area polygons.</td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Slicing</td>
    <td align="center" valign="center"><b>Constrained Random Point Set</td>
    <td align="center" valign="center" colspan="2"><b>Segment Set Intersection</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/geometry_processing/slice.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/randomGridPoints.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/segmentIntersection.png"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/segmentIntersection2.png"></td>
  </tr>
    <tr>
    <td align="center" valign="center">Slice a shape in two along a given line.</td>
    <td align="center" valign="center">Generate constrained random point sets where all points lie within a shape.</td>
    <td align="center" valign="center" colspan="2">Find all points of intersection between a collection of line segments.</td>
  </tr>

  <tr>
    <td align="center" valign="center" colspan="3"><b>Concave Hull</td>
    <td align="center" valign="center"><b>Convex Hull</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/morphology/concaveHull.gif"></td>
    <td valign="top" width="25%"><img src="resources/morphology/concaveHullBFS.png"></td>
    <td valign="top" width="25%"><img src="resources/morphology/concaveHullDFS.png"></td>
    <td valign="top" width="25%"><img src="resources/morphology/convexHull.png"></td>
  </tr>
    <tr>
    <td align="center" valign="center" colspan="3">Concave hull of point sets via breadth-first or depth-first approaches.</td>
    <td align="center" valign="center"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Snap Hull</td>
    <td align="center" valign="center"><b>Shape Intersection</td>
    <td align="center" valign="center" colspan="2"><b>Polygonize Lines</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/morphology/snapHull.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/shapeIntersection.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/polygonizeLines2.png"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/polygonizeLines.png"></td>
  </tr>
    <tr>
      <td align="center" valign="center">A convex hull with some level of shape-feature snapping.</td>
      <td align="center" valign="center">Find all points of intersection between two shapes.</td>
      <td align="center" valign="center" colspan="2">Find the polygonal faces formed by a set of intersecting line segments.</td>
  </tr>
    <tr>
    <td align="center" valign="center"><b>Densification</td>
    <td align="center" valign="center"><b>Tangent Angle</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/geometry_processing/densify.gif"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/tangentAngle.png"></td>
  </tr>
</table>

## *Triangulation*

<table>
  <tr>
    <td align="center" valign="center" colspan="2"><b>Delaunay Triangulation</td>
    <td align="center" valign="center" colspan="2"><b>Earcut Triangulation</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/triangulation/triangulation2.png"></td>
    <td valign="top" width="25%"><img src="resources/triangulation/triangulation1.png"></td>
    <td valign="top" width="25%"><img src="resources/triangulation/earCut.png"></td>
    <td valign="top" width="25%"><img src="resources/triangulation/earCut2.png"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Poisson Delaunay Triangulation</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/contour/poissonTriangulation.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><i>Delaunay triangulation</i> of shapes where <i>steiner points</i> generated by poisson disk sampling are inserted.</td>
  </tr>
</table>

## *Meshing*

<table>
  <tr>
    <td align="center" valign="center" colspan="2"><b>Urquhart Faces</td>
    <td align="center" valign="center"><b>Gabriel Faces</td>
    <td align="center" valign="center"><b>Triangulation Dual</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/triangulation/urquhart1.png"></td>
    <td valign="top" width="25%"><img src="resources/triangulation/urquhart2.png"></td>
    <td valign="top" width="25%"><img src="resources/triangulation/gabriel1.png"></td>
    <td valign="top" width="25%"><img src="resources/triangulation/triangulationDual.png"></td>
  </tr>
  <tr>
    <td align="center" valign="center" colspan="2">Polygon faces of an <i>Urquhart Graph</i> (derived from a triangulation).</td>
    <td align="center" valign="center">Polygon faces of a <i>Gabriel Graph</i> (derived from a triangulation).</td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Centroid Quadrangulation</td>
    <td align="center" valign="center"><b>Edge Collapse Quadrangulation</td>
    <td align="center" valign="center"><b>Split Quadrangulation</td>
    <td align="center" valign="center"><b>Spiral Quadrangulation</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/meshing/centroidQuadrangulation.png"></td>
    <td valign="top" width="25%"><img src="resources/meshing/ecQuadrangulation.png"></td>
    <td valign="top" width="25%"><img src="resources/meshing/quadrangulation.png"></td>
    <td valign="top" width="25%"><img src="resources/meshing/spiralQuadrangulation.gif"></td>
  </tr>
</table>

## *Geometric Optimisation*

<table>
  <tr>
    <td align="center" valign="center"><b>Maximum Inscribed Circle</td>
    <td align="center" valign="center"><b>Minimum Bounding Rectangle</td>
    <td align="center" valign="center" colspan="2"><b>Maximum Inscribed Rectangle</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/pgs/inscribedCircle.gif"></td>
    <td valign="top"><img src="resources/pgs/minimumBoundingRectangle.png"></td>
    <td valign="top"><img src="resources/optimisation/mir.png"></td>
    <td valign="top" ><img src="resources/optimisation/mir2.png"></td>
  </tr>

  <tr>
    <td align="center" valign="center" colspan="2"><b>Minimum Bounding Circle</td>
    <td align="center" valign="center" colspan="2"><b>Minimum Bounding Ellipse</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/pgs/minimumBoundingCircle.png"></td>
    <td valign="top"><img src="resources/optimisation/mbc2.png"></td>
    <td valign="top"><img src="resources/optimisation/mbe1.png"></td>
    <td valign="top" ><img src="resources/optimisation/mbe2.png"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Minimum Bounding Triangle</td>
    <td align="center" valign="center"><b>Envelope</td>
    <td align="center" valign="center" colspan="2"><b>Problem of Apollonius</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/optimisation/mbt.png"></td>
    <td valign="top" width="25%"><img src="resources/geometry_processing/envelope.png"></td>
    <td valign="top"><img src="resources/optimisation/apollonius1.gif"></td>
    <td valign="top"><img src="resources/optimisation/apollonius2.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Largest Empty Circle</td>
    <td align="center" valign="center"><b>Closest Vertex</td>
    <td align="center" valign="center"><b>Closest Point Pair</td>
    <td align="center" valign="center"><b>Farthest Point Pair</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/optimisation/lec.png"></td>
    <td valign="top" width="25%"><img src="resources/pgs/closestVertex.gif"></td>
    <td valign="top" width="25%"><img src="resources/optimisation/closestPair.png"></td>
    <td valign="top" width="25%"><img src="resources/optimisation/farthestPair.png"></td>
  </tr>
  
</table>

## *Circle Packing*

<table>
  <tr>
    <td align="center" valign="center" colspan="2"><b>Front Chain</td>
    <td align="center" valign="center" colspan="2"><b>Trinscribed</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/circle_packing/frontChain1.png"></td>
    <td valign="top"><img src="resources/circle_packing/frontChain2.png"></td>
    <td valign="top"><img src="resources/circle_packing/inscribed1.png"></td>
    <td valign="top"><img src="resources/circle_packing/inscribed2.png"></td>
  </tr>

  <tr>
    <td align="center" valign="center" colspan="2"><b>Maximum Inscribed</td>
    <td align="center" valign="center" colspan="2"><b>Stochastic</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/circle_packing/maxInscribed2.png"></td>
    <td valign="top"><img src="resources/circle_packing/maxInscribed1.png"></td>
    <td valign="top"><img src="resources/circle_packing/pack1.png"></td>
    <td valign="top"><img src="resources/circle_packing/pack2.png"></td>
  </tr>

  <tr>
    <td align="center" valign="center" colspan="2"><b>Square Lattice</td>
    <td align="center" valign="center" colspan="2"><b>Hex Lattice</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/circle_packing/squareLattice1.png"></td>
    <td valign="top"><img src="resources/circle_packing/squareLattice2.png"></td>
    <td valign="top"><img src="resources/circle_packing/hexLattice1.png"></td>
    <td valign="top"><img src="resources/circle_packing/hexLattice2.png"></td>
  </tr>
</table>

## *Coloring*

<table>
  <tr>
    <td valign="top" width="33%"><img src="resources/coloring/leaf.png"></td>
    <td valign="top" width="33%"><img src="resources/coloring/subdivision.png"></td>
    <td valign="top" width="33%"><img src="resources/coloring/dual.png"></td>
  </tr>
  <tr>
    <td valign="top" width="33%"><img src="resources/coloring/voro.png"></td>
    <td valign="top" width="33%"><img src="resources/coloring/swirl.png"></td>
    <td valign="top" width="33%"><img src="resources/coloring/dino.png"></td>
  </tr>
</table>

## *Construction*

<table>
  <tr>
    <td align="center" valign="center"><b>Supercircle</td>
    <td align="center" valign="center"><b>Supershape</td>
    <td align="center" valign="center" colspan="2"><b>Star</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/pgs/superCircle.gif"></td>
    <td valign="top"><img src="resources/pgs/supershape.gif"></td>
    <td valign="top"><img src="resources/pgs/star.gif"></td>
    <td valign="top" ><img src="resources/pgs/star2.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Random Convex Polygon</td>
    <td align="center" valign="center"><b>Heart</td>
    <td align="center" valign="center"><b>Ring</td>
    <td align="center" valign="center"><b>Linear Spiral</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/pgs/randomPolygon.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/heart.png"></td>
    <td valign="top" width="25%"><img src="resources/pgs/ring.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/spiral.gif"></td>
  </tr>
  <tr>
    <td align="center" valign="center"><b>Fermat Spiral</td>
    <td align="center" valign="center"><b>Rectangular Spiral</td>
    <td align="center" valign="center"><b>Koch Snowflake</td>
    <td align="center" valign="center"><b>Hilbert Curve</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/pgs/fermatSpiral.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/rectSpiral.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/kochSnowflake.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/hilbertCurve.gif"></td>
  </tr>
  <tr>
    <td align="center" valign="center"><b>Sierpinski Carpet</td>
    <td align="center" valign="center"><b>Sierpinski Curve</td>
    <td align="center" valign="center"><b>Sierpinski Tri-Curve A</td>
    <td align="center" valign="center"><b>Sierpinski Tri-Curve B</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/pgs/sierpinskiCarpet.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/sierpinskiCurve.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/sierpinskiThreeSteps.gif"></td>
    <td valign="top" width="25%"><img src="resources/pgs/sierpinskiFourSteps.gif"></td>
  </tr>
</table>

## *Point Sets*

<table>
  <tr>
    <td align="center" valign="center"><b>Random</td>
    <td align="center" valign="center"><b>Gaussian</td>
    <td align="center" valign="center"><b>Square Grid</td>
    <td align="center" valign="center"><b>Hex Grid</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/point_set/random.png"></td>
    <td valign="top" width="25%"><img src="resources/point_set/gaussian.png"></td>
    <td valign="top" width="25%"><img src="resources/point_set/squareGrid.png"></td>
    <td valign="top" width="25%"><img src="resources/point_set/hexGrid.png"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Phyllotaxis</td>
    <td align="center" valign="center"><b>Poisson</td>
    <td align="center" valign="center"><b>Hexagon</td>
    <td align="center" valign="center"><b>Ring</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/point_set/phyllotaxis.gif"></td>
    <td valign="top" width="25%"><img src="resources/point_set/poisson.png"></td>
    <td valign="top" width="25%"><img src="resources/point_set/hexaPoints.gif"></td>
    <td valign="top" width="25%"><img src="resources/point_set/ring.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>Halton LDS</td>
    <td align="center" valign="center"><b>Hammersley LDS</td>
    <td align="center" valign="center"><b>Plastic LDS</td>
    <td align="center" valign="center"><b>Jittered Plastic LDS</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/point_set/haltonLDS.gif"></td>
    <td valign="top" width="25%"><img src="resources/point_set/hammersleyLDS.gif"></td>
    <td valign="top" width="25%"><img src="resources/point_set/plasticLDS.gif"></td>
    <td valign="top" width="25%"><img src="resources/point_set/plasticJitteredLDS.gif"></td>
  </tr>

  <tr>
    <td align="center" valign="center"><b>N-Rooks LDS</td>
    <td align="center" valign="center"><b>Hilbert Sort</td>
    <td align="center" valign="center"><b>Distance Prune</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/point_set/nRooksLDS.png"></td>
    <td valign="top" width="25%"><img src="resources/point_set/hilbertSort.gif"></td>
    <td valign="top" width="25%"><img src="resources/point_set/removeWithinDistance.gif"></td>
  </tr>
</table>

## *Tiling & Subdivision*

<table>
  <tr>
    <td align="center" valign="center"><b>Random Quad Subdivision</td>
    <td align="center" valign="center"><b>Random Rect Subdivision</td>
    <td align="center" valign="center"><b>Random Triangle Subdivision</td>
    <td align="center" valign="center"><b>Islamic Tiling</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/tiling/randomSubdivision.png"></td>
    <td valign="top" width="25%"><img src="resources/tiling/rectSubdivision.png"></td>
    <td valign="top" width="25%"><img src="resources/tiling/triangleSubdivision.png"></td>
    <td valign="top" width="25%"><img src="resources/tiling/islamic.png"></td>
  </tr>
  <tr>
    <td align="center" valign="center" colspan="2"><b>Doyle Spiral</td>
    <td align="center" valign="center"><b>Hexagon Tiling</td>
    <td align="center" valign="center"><b>Penrose Tiling</td>
  </tr>
  <tr>
    <td valign="top" width="25%"><img src="resources/tiling/doyeSpiral1.png"></td>
    <td valign="top" width="25%"><img src="resources/tiling/doyeSpiral2.png"></td>
    <td valign="top" width="25%"><img src="resources/tiling/hex.png"></td>
    <td valign="top" width="25%"><img src="resources/tiling/penrose.png"></td>
  </tr>
</table>