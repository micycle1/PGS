import processing.javafx.*;
import micycle.pgs.*;
import micycle.pgs.PGS_Coloring.ColoringAlgorithm;
import java.util.List;
import java.util.Collection;
import micycle.pgs.color.Palette;

void setup() {
  size(1000, 1000, FX2D);
  smooth();
}

void draw() {
  quilt(4, 2);
}

void quilt(int grid, long seed) { // grid is NxN
  
  // create blob
  float buffer = 50;
  var points = PGS_PointSet.poisson(buffer, buffer, width-buffer, height-buffer, 35, seed);
  var s = PGS_Hull.concaveHullBFS(points, 0);
  s = PGS_Morphology.simplifyHobby(s, 1); // make bendy/"organic"

  // create grid
  float cX = width/2;
  float cY = height/2;
  var segz = PGS_SegmentSet.parallelSegments(cX, cY, width * 1.5, width / grid, 0, grid - 1);
  segz.addAll(PGS_SegmentSet.parallelSegments(cX, cY, height * 1.5, height / grid, PI / 2, grid - 1));
  var segs = PGS_SegmentSet.toPShape(segz);
  
  // warp the grid
  float time = 0;
  segs = PGS_Morphology.fieldWarp(segs, 50, 1, time, true, seed);

  // create quilt patches
  var plane = PGS_Construction.createRect(0, 0, width, height, 0);
  var patches = PGS_ShapeBoolean.unionLines(segs, plane);

  Palette p = Palette.getPalette(38);
  final var hobby = s;

  PGS_Coloring.colorMesh(patches, ColoringAlgorithm.RLF).forEach((patch, col) -> {
    int patchBG = col == 0 ? p.get(0) : p.get(2);
    int curveColor = col == 0 ? p.get(1) : p.get(3);

    patch.setFill(patchBG);
    patch = PGS_Conversion.setAllStrokeToFillColor(patch, 0.5);
    PGS_Conversion.disableAllStroke(patch);
    shape(patch);
    try {
      var hobbySlice = PGS_ShapeBoolean.intersect(patch, hobby);
      hobbySlice = PGS_Conversion.setAllFillColor(hobbySlice, curveColor);
      PGS_Conversion.disableAllStroke(hobbySlice);
      shape(hobbySlice);
    }
    catch (Exception e) {
      // just in case
    }
  }
  );

  PGS_Conversion.setAllStrokeColor(segs, 0, 2);
  shape(segs);
}
