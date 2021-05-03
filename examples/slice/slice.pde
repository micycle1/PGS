import micycle.pgs.*;
import java.util.List;

PShape polygon;

void setup() {
  size(800, 800);
  smooth();
}

void draw() {
  background(0, 0, 40);

  if ((frameCount-1) % 180 == 0) {
    polygon = PGS_Construction.createRandomPolygon((int) random(3, 7), width, height);
  }

  slice(polygon);
}

void slice(PShape shape) {
  PShape mec = PGS_Optimisation.minimumBoundingCircle(shape);
  PVector pv1 = PGS_Processing.pointOnExterior(mec, frameCount * 0.006f, -1);
  PVector pv2 = PGS_Processing.pointOnExterior(mec, frameCount * 0.003f + 0.5f, -1);

  List<PShape> slice = PGS_Processing.slice(shape, pv1, pv2);
  PShape s = slice.get(0);
  PShape s2 = slice.get(1);
  
  PGS_Conversion.setAllFillColor(s, color(237, 50, 162));
  PGS_Conversion.setAllStrokeColor(s, color(237, 50, 162), 4);
  PGS_Conversion.setAllFillColor(s2, color(255));

  shape(s);
  shape(s2);
  
  strokeWeight(5);
  stroke(0, 0, 40);
  line(pv1.x, pv1.y, pv2.x, pv2.y);
  stroke(color(237, 50, 162));
  strokeWeight(15);
  point(pv1.x, pv1.y);
  point(pv2.x, pv2.y);
}
