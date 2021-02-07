# Processing Topology Suite

<h3 align="center"> 🚧 Under Construction 🚧 </h3>

---

A library for manipulating PShape geometry

PTS wraps JTS, enabling its methods to be applied to Processing's `PShape` objects. Beyond that, PTS provides other geometry __ such as splines.

[Contents from https://doc.cgal.org/latest/Manual/packages.html]

## Morphology

### PShapes

- Boolean operations: union, difference, intersection, etc. (/OVERLAY OPERATIONS)
- Shape Boundaries
- Smoothing
- Simplification
- Area, centroid, etc.

### SPATIAL RELATIONSHIPS
###  BUFFERS

### POLYGONIZATION

### Geometry methods
- Spatial Predicates, relate()
- Overlay ops, buffer(), convexHull()
- Metrics

### Geometry Processing
- Line Merging
- Noding & Polygonization
- Simplification
- Linear Referencing



### Fields/ Point Sets

- Voronoi
- Poisson-Disc

## Libraries

- [OS_Minkowski_Sum_Diff
](https://github.com/OrdnanceSurvey/OS_Minkowski_Sum_Diff)
- [JTS](https://github.com/locationtech/jts)

Shortcomings with JTS native triangulation (`DelaunayTriangulationBuilder`):

- Doesn't respect concave shapes/holes (which arises from computing triangulation of the vertices only, not edges) (effectively triangulates the convex hull)
- No refinement: comparison [here](http://www.cs.cmu.edu/~quake/triangle.quality.html)
  - Long, thin triangles (bad angles)
  - Large difference (non-uniform) in triangle areas
  - No way to sub-divide without 
  - Many triangles may share one boundary vertex