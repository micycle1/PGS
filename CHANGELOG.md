# Changelog

All notable changes to PGS will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Dates are *YYYY-MM-DD*.

## [1.1.0] – 2021-00-00

### Added


- `closestPointPair()` to `PGS_Optimisation`. The method efficiently computes the **closest** pair of points in a set of points.
- `farthestPointPair()` to `PGS_Optimisation`. The method efficiently computes the **farthest** pair of points in a set of points.
- `chaikinCut()` to `PGS_Morphology`. The method smooths shapes via iterated corner cuts.
- Expand `PGS_Conversion` to support conversion between:
  - `PATH` PShape<->JTS `LineString`
  - `POINTS` PShape<->JTS `MultiPoint`
  - `LINES` PShape<->JTS `MultiLineString`

### Changed 

- Split `PGS_Processing.concaveHull()` into `concaveHullDFS()` and `concaveHullBFS()` (the method previously used the BFS approach only).

### Fixed
- Error when `concaveHull2()` was called with alpha > 1.
- Concave hull methods no longer mutate the input point set.
- PShapes marked as closed and having less than 3 vertices could cause an error during conversion ([#22](https://github.com/micycle1/PGS/issues/22)).

### Removed
- `PGS_Processing.concaveHull()` (see *Changed* above)

## [1.0.0] – 2021-05-06