import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;

PShape shape;
boolean inward = true;

void setup() {
  size(800, 800, FX2D);
  smooth();
  textAlign(CENTER, CENTER);

  shape = new PShape();
}

void draw() {
  background(50, 127, 168);

  if (shape.getVertexCount() > 2 || shape.getChildCount() > 0) {
    shape(shape);

    PShape offsetCurves;

    if (inward) {
      offsetCurves = PGS_Contour.offsetCurvesInward(shape, PGS_Contour.OffsetStyle.ROUND, map(mouseX, 0, width, 10, 30));
    } else {
      offsetCurves = PGS_Contour.offsetCurvesOutward(shape, PGS_Contour.OffsetStyle.ROUND, map(mouseX, 0, width, 10, 30), 20);
    }

    colorMode(HSB, 1);
    for (int i = 0; i < offsetCurves.getChildCount(); i++) {
      offsetCurves.getChild(i).setStroke(color((float) mouseY/height, ((float) i / (offsetCurves.getChildCount()-1)) + 0.1, 1));
    }
    colorMode(RGB, 255);

    shape(offsetCurves);
  } else {
    textSize(36);
    text("Use mouse to draw!", width/2, height/2);
    textSize(14);
    text("Press any key to switch between inward/outward style.", width/2, height/2 + 30);
  }
}

void mouseDragged() {
  shape = PGS_ShapeBoolean.union(shape, createShape(ELLIPSE, mouseX, mouseY, 100, 100));
  shape = PGS_Morphology.simplify(shape, 1);
  PGS_Conversion.disableAllStroke(shape);
}

void mouseClicked() {
  mouseDragged();
}

void keyPressed() {
  inward = !inward;
}
