[![](https://jitpack.io/v/micycle1/PGS.svg)](https://jitpack.io/#micycle1/PGS) [![Lines of Code](https://sonarcloud.io/api/project_badges/measure?project=micycle1_PTS&metric=ncloc)](https://sonarcloud.io/dashboard?id=micycle1_PTS)

# Processing Geometry Suite

*Processing Geometry Suite* is a software project that provides easy access to 2D geometric algorithms in the form of a [Processing](https://processing.org/) library.

The focus of the library is on visualisation rather than providing underlying data structures. To this end all methods in the library are static and most of them take in and return [`PShapes`](https://processing.org/reference/PShape.html) or [`PVectors`](https://processing.org/reference/PVector.html).

Docs are hosted via *GitHub Pages* [here](https://micycle1.github.io/PGS/).

## **Overview**

Library functionality is split over the following classes:

* `PGS_Construction`
  * Construct uncommon 2D primitives
* `PGS_Contour`
  * Methods that produce various contours from shapes: medial axes, straight skeletons, offset curves, etc.
* `PGS_Conversion`
  * Conversion between PShapes and JTS Geometries 
* `PGS_Morphology`
  * Methods that affect the geometry or topology of shapes (buffering, simplification, smoothing, etc.)
* `PGS_Optimsation`
  * Solve geometric optimisation problems, such as finding the maximum inscribed circle, or the closest vertex to a coordinate
* `PGS_Processing`
  * Methods that process a shape in some way: compute hulls, partition, slice, etc. 
* `PGS_ShapeBoolean`
  * Boolean set-operations for 2D shapes
* `PGS_ShapePredicates`
  * Various shape metrics (area, circularity, etc.) and predicates (*"do these shapes intersect?"*)
* `PGS_Transformation`
  * Various geometric and affine transformations that affect vertex coordinates
* `PGS_Triangulation`
  * Delaunay triangulation (constrained and refined) and earcut triangulation of shapes and point sets
* `PGS_Voronoi`
  * Voronoi Diagrams of shapes and point sets

## **Installation**

<details><summary>Processing IDE — Quick</summary>
<p>

Download the latest *PGS.jar* from [releases](https://github.com/micycle1/PGS/releases) and simply drag-and-drop it onto the [Processing IDE](https://processing.org/reference/environment/).
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

A handful of example Processing sketches are provided in [examples](https://github.com/micycle1/PGS/tree/master/examples).

<p float="middle">
  <img src="resources/examples/drawOffsetCurves.png" alt="" width="24%"/>
  <img src="resources/examples/partitionSmooth.png" alt="" width="24%"/>
  <img src="resources/examples/triangulate.png" alt="" width="24%"/>
  <img src="resources/examples/minkShearLetters.png" alt="" width="24%"/>
</p>


## **Illustrations**

Much of the functionality (but by no means all) is demonstrated below:

## *2D Boolean Operations*

<table>
  <tr>
    <td align="center" valign="center"><h3>Union</td>
    <td align="center" valign="center"><h3>Intersection</td>
    <td align="center" valign="center"><h3>Subtraction</td>
    <td align="center" valign="center"><h3>Symmetric Difference</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/boolean/union.gif"></td>
    <td valign="top"><img src="resources/boolean/intersect.gif"></td>
    <td valign="top"><img src="resources/boolean/subtract.gif"></td>
    <td valign="top" ><img src="resources/boolean/symDifference.gif"></td>
  </tr>
</table>

### Complement
<img src="resources/boolean/complement.png" alt="" width="25%"/>

## *Transformation*

*Only methods beyond those offered in Processing are illustrated below.*

### Rotate Around
Rotate a shape around its centroid or an arbitrary point.

<p float="middle">
  <img src="resources/transform/rotateCenter.gif" alt="" width="25%"/>
  <img src="resources/transform/rotate.gif" alt="" width="25%"/>
</p>

### Translate To
Translate a shape such that its centroid matches some position.

<img src="resources/transform/translateTo.gif" alt="" width="25%"/>

### Touch Scale
Scale one shape such that it touches another.

<img src="resources/transform/touchScale.gif" alt="" width="25%"/>

### Homothetic Transformation
Projection-transform a shape with respect to a fixed point.

<img src="resources/transform/homothetic.gif" alt="" width="25%"/>

## *Geometric Predicates & Metrics*

<table table-layout="fixed">
  <tr>
    <td align="center" valign="center"><h3>Intersects</td>
    <td align="center" valign="center"><h3>Contains Shape</td>
    <td align="center" valign="center"><h3>Contains Point</td>
  </tr>
  <tr>
    <td valign="top" align="center" width="33%"> <img src="resources/predicate/intersect.gif"><br>Do shapes intersect with each other?</td>
    <td valign="top" align="center" width="33%"><img src="resources/predicate/contains.gif"><br>Does one shape fully contain another?</td>
    <td valign="top" align="center" width="33%"><img src="resources/predicate/containsPoint.gif"><br>For individual points and point sets.</td>
</table>


### Metrics
* Length
* Circularity
* Similarity
* Holes
* Is simple?
* Is convex?
* Distance
* Area
* Centroid

## *Contour*

### Medial Axis
Medial axis transform with feature pruning via distance, area or *axial angle*. 

<p float="middle">
  <img src="resources/contour/medialAxis.gif" alt="" width="25%"/>
  <img src="resources/contour/medialAxis.png" alt="" width="25%"/>
</p>

### Straight Skeleton
<p float="middle">
  <img src="resources/contour/straightSkeleton.png" alt="" width="25%"/>
</p>

### Isolines (topographic contour lines)
Isolines from intra-shape euclidean distance, or point sets.

<p float="middle">
  <img src="resources/contour/isolines.gif" alt="" width="25%"/>
  <img src="resources/contour/isolines2.gif" alt="" width="25%"/>
</p>


### Offset Curves
Inner and exterior offset curves; based on *miter*, *bevel* or *round* offset styles. 

<p float="middle">
  <img src="resources/contour/miteredInterior.gif" alt="" width="25%"/>
  <img src="resources/contour/miteredExterior.gif" alt="" width="25%"/>
</p>

### Voronoi Diagram
Point site voronoi diagrams from shape vertices.

<p float="middle">
  <img src="resources/contour/voronoi.gif" alt="" width="25%"/>
  <img src="resources/contour/voronoi1.png" alt="" width="25%"/>
</p>

### Circle Site Voronoi Diagram
Circle-site voronoi diagrams via point-site approximation.
<p float="middle">
  <img src="resources/contour/voronoiCircles.gif" alt="" width="25%"/>
  <img src="resources/contour/voronoiCircles.png" alt="" width="25%"/>
</p>

### Delaunay Triangulation
Constrained & refined *Delaunay triangulation* of shapes and point sets.

<p float="middle">
  <img src="resources/contour/triangulation1.png" alt="" width="25%"/>
  <img src="resources/contour/triangulation2.png" alt="" width="25%"/>
</p>

### Poisson Delaunay Triangulation
*Delaunay triangulation* of shapes where *steiner points* generated by poisson disk sampling are inserted.

<img src="resources/contour/poissonTriangulation.gif" alt="" width="25%"/>

### Earcut Triangulation
<p float="middle">
  <img src="resources/contour/earCut.png" alt="" width="25%"/>
  <img src="resources/contour/earCut2.png" alt="" width="25%"/>
</p>

## *Morphology*

### Buffer
<img src="resources/morphology/buffer.gif" alt="" width="25%"/>


### Erosion-Dilation
A negative followed by a positive buffer (in a single operation).

<img src="resources/morphology/erosionDilation.gif" alt="" width="25%"/>

### Minkowski Addition
Minkowski sum and difference (a.k.a buffer one shape using another shape; pictured: buffering using a rotating & growing triangle).
<p float="middle">
  <img src="resources/morphology/minkSum.gif" alt="" width="25%"/>
  <img src="resources/morphology/minkDiff.gif" alt="" width="25%"/>
</p>

### Simplification
<img src="resources/morphology/simplifyVW.gif" alt="" width="25%"/>

### Smoothing
<img src="resources/morphology/smooth.gif" alt="" width="25%"/>

### Gaussian Smoothing
<img src="resources/morphology/gaussianSmooth.gif" alt="" width="25%"/>

### Rounding

<p float="middle">
  <img src="resources/morphology/round.gif" alt="" width="25%"/>
  <img src="resources/morphology/round2.gif" alt="" width="25%"/>
</p>

## *Geometry Processing*

### Point on Perimeter
Find a point some fraction along the perimeter of a shape (with perpendicular offset).

<img src="resources/geometry_processing/pointOnPerimeter.gif" alt="" width="25%"/>

### Points on Perimeter
Find *N* points (evenly distributed) along the perimeter of a shape, or points every *D* distance (with optional perpendicular offset).

<p float="middle">
  <img src="resources/geometry_processing/pointsOnPerimeter.gif" alt="" width="25%"/>
  <img src="resources/geometry_processing/pointsOnPerimeter2.gif" alt="" width="25%"/>
</p>

### Partitioning
Partition a shape into simple (convex) polygons.

<p float="middle">
  <img src="resources/geometry_processing/decompose1.png" alt="" width="25%"/>
  <img src="resources/geometry_processing/decompose2.png" alt="" width="25%"/>
</p>

### Splitting
Subdivide (recursively) a shape into quadrants

<img src="resources/geometry_processing/split.gif" alt="" width="25%"/>

### Slicing
Slice a shape in two along a given line

<img src="resources/geometry_processing/slice.gif" alt="" width="25%"/>

### Densification
<img src="resources/geometry_processing/densify.gif" alt="" width="25%"/>

### Constrained Random Point Set
Generate constrained random point sets where all points lie within a shape. Points can be distributed entirely randomly or according to grid with configurable tightness.

<p float="middle">
  <img src="resources/geometry_processing/randomGridPoints.gif" alt="" width="25%"/>
    <img src="resources/geometry_processing/randomGridPoints2.gif" alt="" width="25%"/>
</p>

### Envelope
<img src="resources/geometry_processing/envelope.png" alt="" width="25%"/>

### Concave Hull
Concave hull of point sets.
<p float="middle">
  <img src="resources/morphology/concaveHull.gif" alt="" width="25%"/>
  <img src="resources/morphology/concaveHull2.gif" alt="" width="25%"/>
</p>

### Convex Hull
<img src="resources/morphology/convexHull.png" alt="" width="25%"/>

### Snap Hull
A convex hull with some level of shape-feature snapping.

<img src="resources/morphology/snapHull.gif" alt="" width="25%"/>

### Shape Intersection
Find all points of intersection between two shapes.

<img src="resources/geometry_processing/shapeIntersection.gif" alt="" width="25%"/>

### Segment Set Intersection

Find all points of intersection between a collection of line segments.

<p float="middle">
  <img src="resources/geometry_processing/segmentIntersection.png" alt="" width="25%"/>
  <img src="resources/geometry_processing/segmentIntersection2.png" alt="" width="25%"/>
</p>

## *Geometric Optimisation*

### Closest Point
<img src="resources/pgs/closestVertex.gif" alt="" width="25%"/>

### Maximum Inscribed Circle
<img src="resources/pgs/inscribedCircle.gif" alt="" width="25%"/>

### Maximum Inscribed Rectangle
Maximum inscribed axis-aligned rectangle of convex shapes.

<p float="middle">
  <img src="resources/optimisation/mir1.png" alt="" width="25%"/>
  <img src="resources/optimisation/mir2.png" alt="" width="25%"/>
</p>

### Minimum Bounding Circle

<p float="middle">
  <img src="resources/pgs/minimumBoundingCircle.png" alt="" width="25%"/>
  <img src="resources/optimisation/mbc2.png" alt="" width="25%"/>
</p>

### Minimum Bounding Ellipse

<p float="middle">
  <img src="resources/optimisation/mbe1.png" alt="" width="25%"/>
  <img src="resources/optimisation/mbe2.png" alt="" width="25%"/>
</p>

### Minimum Bounding Rectangle
<img src="resources/pgs/minimumBoundingRectangle.png" alt="" width="25%"/>

### Problem of Apollonius

<p float="middle">
  <img src="resources/optimisation/apollonius1.gif" alt="" width="25%"/>
  <img src="resources/optimisation/apollonius2.gif" alt="" width="25%"/>
</p>

## *Construction*

### Supercircle
Generate *supercircles*, with a configurable constant.

<img src="resources/pgs/superCircle.gif" alt="" width="25%"/>

### Supershape

<img src="resources/pgs/supershape.gif" alt="" width="25%"/>


### Star

<p float="middle">
  <img src="resources/pgs/star.gif" alt="" width="25%"/>
  <img src="resources/pgs/star2.gif" alt="" width="25%"/>
</p>

### Random Polygon
Generate a random convex n-gon

<img src="resources/pgs/randomPolygon.gif" alt="" width="25%"/>