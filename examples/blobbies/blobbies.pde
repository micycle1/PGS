import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;
import java.util.Collection;
import micycle.pgs.color.Palette;
import micycle.pgs.PGS_Contour.OffsetStyle;

void setup() {
  size(1000, 1000, FX2D);
  smooth();
}

void draw() {
  background(0);
  blobbiesArt(3 + mouseX/100d);
}

void blobbiesArt(double offset) {
  var b1 = blobbie(250, 250, 66);
  var b2 = blobbie(450, 750, 10);
  var b3 = blobbie(700, 225, 1110);
  var b4 = blobbie(1000, 1000, 11111);

  var a = PGS_ShapeBoolean.union(b1, b2, b3, b4); // union rather than flatten, in case of overlap

  var outerCurves = PGS_Contour.offsetCurvesOutward(a, OffsetStyle.ROUND, offset, (int) (400 / offset));
  strokeWeight((float) offset);
  var palette = Palette.FREAK;
  var colOffset = 0;
  PGS_Processing.applyWithIndex(outerCurves, (i, c) -> {
    stroke(palette.get(i + colOffset));
    var points = PGS_Processing.pointsOnExterior(c, offset, 0);
    shape(points);
  }
  );

  var innerCurves = PGS_Contour.offsetCurvesInward(a, OffsetStyle.ROUND, offset);
  PGS_Processing.applyWithIndex(innerCurves, (i, c) -> {
    stroke(Palette.SHEETS.get(i + colOffset));
    var points = PGS_Processing.pointsOnExterior(c, offset, 0);
    shape(points);
  }
  );
}

PShape blobbie(double cx, double cy, long seed) {
  final double GOLDEN_RATIO = (1 + Math.sqrt(5)) / 2;
  var s = PGS_Construction.createStar(cx, cy, 5, 125, 400, 0);
  s = PGS_Morphology.simplify(s, 10);
  s = PGS_Morphology.fieldWarp(s, 100, 1, 0, false, seed);
  s = PGS_Morphology.simplifyHobby(s, 0.9);
  s = PGS_Morphology.fieldWarp(s, 50, 2, 0, false, seed);
  s = PGS_Transformation.scale(s, 0.6);
  s = PGS_Transformation.rotateAroundCenter(s, seed * GOLDEN_RATIO);
  return s;
}

public void shape(Collection<PVector> points) {
  points.forEach(p -> circle(p));
}

public void circle(PVector v) {
  if (v.z == 0) {
    point(v);
  } else {
    ellipse(v.x, v.y, v.z * 2, v.z * 2);
  }
}

public void point(PVector v) {
  super.point(v.x, v.y);
}
