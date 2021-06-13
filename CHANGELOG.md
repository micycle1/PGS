# Changelog

All notable changes to PGS will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Dates are *YYYY-MM-DD*.

## [1.1.0] – 2021-06-xx

### Added

- **`PGS_CirclePacking`** — a class for circle packings of shapes, subject to varying constraints and patterns of tangencies
- `closestPointPair()` to `PGS_Optimisation`. Efficiently computes the **closest** pair of points in a set of points.
- `farthestPointPair()` to `PGS_Optimisation`. Efficiently computes the **farthest** pair of points in a set of points.
- `chaikinCut()` to `PGS_Morphology`. Smoothes shapes via iterated corner cuts.
- `createHeart()` to `PGS_Construction`. Generates heart-shaped PShapes.
- `urquhartFaces()` to `PGS_Triangulation`. Tesselates a triangulation into polygons corresponding to the faces of an _Urquhart graph_.
- `gabrielFaces()` to `PGS_Triangulation`. Tesselates a triangulation into polygons corresponding to the faces of an _Gabriel graph_.
- A new `earCutTriangulation()` method signature that takes in a PShape argument (previously it accepted a list of points only)
- Expand `PGS_Conversion` to support conversion between:
  - `PATH` PShape<->JTS `LineString`
  - `POINTS` PShape<->JTS `MultiPoint`
  - `LINES` PShape<->JTS `MultiLineString`

### Changed 

- Split `PGS_Processing.concaveHull()` into `concaveHullDFS()` and `concaveHullBFS()` (the method previously used the BFS approach only).
- Compute (rather than ignore) circle sites of radius 0 (these are effectively points) in `PGS_Voronoi.voronoiCirclesDiagram()`
- Replaced the algorithm used by `PGS_Processing.generateRandomPoints()` with a triangulation-based approach. The new approach is ~10x faster!
- Renamed `delaunayTriangulationTin()` to `delaunayTriangulationMesh()`.

### Fixed
- Error when `concaveHull2()` was called with alpha > 1.
- Concave hull methods no longer mutate the input point set.
- PShapes marked as closed and having less than 3 vertices could cause an error during conversion ([#22](https://github.com/micycle1/PGS/issues/22)).
- `PGS_Conversion.toPVector()` now handles [primitive](https://processing.org/examples/shapeprimitives.html) PShapes
- Constrained delaunay triangulations now respect shape holes

### Removed
- `PGS_Processing.concaveHull()` (see *Changed*)

## [1.0.0] – 2021-05-06