import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;
import java.util.Collection;
import micycle.pgs.color.Palette;

void setup() {
  size(1000, 1000, FX2D);
  smooth();
}

void draw() {
  var palette = Palette.FREAK;

  int i = 3;

  background(palette.get(i));

  long seed = 110;
  PShape s = PGS_Construction.createSponge(width, height, 20, 20, 80, 5, seed);
  s = PGS_Transformation.scale(s, 1.15);
  // sponge can create an arrangement with nested holes; if so, pick the first child
  if (s.getChildCount() > 0) {
    s = PGS_Conversion.reorderChildren(s, (a, b) -> Double.compare(PGS_ShapePredicates.area(b), PGS_ShapePredicates.area(a)));
    s = s.getChild(0);
  }

  // circle packing of the structure
  var c = PGS_CirclePacking.frontChainPack(s, 3, 15, 0);
  fill(palette.get(i + 1));
  strokeWeight(2);
  stroke(palette.get(i + 2));
  shape(c);

  boolean dense = mouseX > width/2;
  // create a dense packing of structure (SLOW)
  if (dense) {
    var gaps = PGS_ShapeBoolean.subtract(s, PGS_Conversion.toCircles(c));
    gaps = PGS_Morphology.simplify(gaps, 2);
    var c3 = PGS_CirclePacking.maximumInscribedPack(gaps, 3d, 1);
    shape(c3);
  }

  var holes = PGS_Processing.extractHoles(s);
  holes = PGS_Morphology.buffer(holes, -50);

  // circle packing of the structure's holes
  var c2 = PGS_CirclePacking.frontChainPack(holes, 6, 20, 1003130);
  fill(palette.get(i + 3));
  stroke(palette.get(i + 4));
  shape(c2);
}

public void shape(Collection<PVector> points) {
  points.forEach(p -> circle(p));
}

public void point(PVector v) {
  super.point(v.x, v.y);
}
public void circle(PVector v) {
  if (v.z == 0) {
    point(v);
  } else {
    ellipse(v.x, v.y, v.z * 2, v.z * 2);
  }
}
