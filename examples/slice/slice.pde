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
    polygon = PGS_Construction.createSuperRandomPolygon(width, choice(5,50) , random(0.05, 0.2), choice(0, 3), choice(1, 3), false, true, millis());
    polygon = PGS_Transformation.resizeByMajorAxis(polygon, 666);
    polygon = PGS_Transformation.translateEnvelopeTo(polygon, width/2, height/2);
  }

  slice(polygon);
}

void slice(PShape shape) {
  PShape mec = PGS_Optimisation.minimumBoundingCircle(shape);
  PVector pv1 = PGS_Processing.pointOnExterior(mec, frameCount * 0.006f, -1);
  PVector pv2 = PGS_Processing.pointOnExterior(mec, frameCount * 0.003f + 0.5f, -1);

  PShape slices = PGS_Processing.slice(shape, pv1, pv2);
  PShape s = slices.getChild(0);
  PShape s2 = slices.getChild(1);

  PGS_Conversion.setAllFillColor(s, color(237, 50, 162));
  PGS_Conversion.setAllStrokeColor(s, color(237, 50, 162), 3);
  PGS_Conversion.setAllFillColor(s2, color(255));
  //PGS_Conversion.setAllStrokeColor(s2, color(237, 50, 162), 3);

  shape(s);
  shape(s2);

  strokeWeight(5);
  stroke(255, 255, 0, 100);
  line(pv1.x, pv1.y, pv2.x, pv2.y);
  stroke(color(237, 50, 162));
  strokeWeight(15);
  point(pv1.x, pv1.y);
  point(pv2.x, pv2.y);
}
