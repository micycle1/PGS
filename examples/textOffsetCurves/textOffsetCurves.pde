import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;

PShape pgs;
PFont font;

void setup() {
  size(1000, 1000, FX2D);
  smooth();
  colorMode(HSB, 1, 1, 1, 1);
  font = createFont(PFont.list()[(int) random(PFont.list().length)], 300, true);
}

void draw() {
  background(0, 0, 0.1);
  prepareText();
  PShape offsetCurves = PGS_Contour.offsetCurvesOutward(pgs, PGS_Contour.OffsetStyle.BEVEL, 15 + sin(frameCount*0.02f)*8, 25);

  float hue = 0;
  final float inc = 1f / offsetCurves.getChildCount();
  for (PShape offsetCurve : offsetCurves.getChildren()) {
    offsetCurve.setStroke(color((hue+frameCount*0.01f) % 1, 0.8f, 1, 1));
    offsetCurve.setStrokeWeight(3);
    hue += inc;
  }

  shape(offsetCurves);
  shape(pgs);
}

void prepareText() {
  PShape p = font.getShape('P');
  PShape g = font.getShape('G');
  PShape s = font.getShape('S');
  float pWidth = PGS_ShapePredicates.width(p);
  float gWidth = PGS_ShapePredicates.width(g);
  p = PGS_Morphology.fieldWarp(p, 50, 0.6, frameCount*0.01, false, 81);
  g = PGS_Transformation.translate(g, pWidth + 30, 0);
  g = PGS_Morphology.fieldWarp(g, 60, 0.75, frameCount*0.01, false, 1337);
  s = PGS_Transformation.translate(s, pWidth + gWidth + 60, 0);
  s = PGS_Morphology.fieldWarp(s, 30, 0.5, frameCount*0.01, false, 123081);
  pgs = createShape(GROUP);
  pgs.addChild(p);
  pgs.addChild(g);
  pgs.addChild(s);

  pgs = PGS_Transformation.translateEnvelopeTo(pgs, width/2f, height/2f);
  PGS_Conversion.setAllFillColor(pgs, color(0.5f, 0, 1));
}

void keyPressed() {
  font = createFont(PFont.list()[(int) random(PFont.list().length)], 300, true);
}
