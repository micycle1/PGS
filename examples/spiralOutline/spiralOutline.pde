import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;
import micycle.uniformnoise.*;
import micycle.pgs.PGS_Contour.OffsetStyle;
import micycle.pgs.PGS_Morphology.CapStyle;

PShape polygon;
PShape triangles;
UniformNoise noise;

void setup() {
  size(1000, 1000, FX2D);
  smooth();
  colorMode(HSB, 1, 1, 1);
  rectMode(CORNER);
  noStroke();
  noise = new UniformNoise();
}

void draw() {
  fill(0, 0, 0.05, 64);
  rect(0, 0, width, height);

  PShape spiral = PGS_Construction.createLinearSpiral(width/2, height/2, 0.5+mouseX/200f, 250+mouseY/5f);
  spiral = PGS_Transformation.rotate(spiral, new PVector(width/2, height/2), frameCount/100f);
  spiral = PGS_Processing.extractPerimeter(spiral, 0.005, 1);
  spiral.setFill(false);
  shape(spiral);

  spiral = PGS_Morphology.buffer(spiral, 20, OffsetStyle.ROUND, CapStyle.ROUND);

  int perimeters = 30; // perimeter sections
  for (double i = 0; i < 1; i += 1f/perimeters) {
    PShape perimeter = PGS_Processing.extractPerimeter(spiral, (i +frameCount * 0.001) % 1, (i + 0.5f/perimeters + frameCount * 0.001) % 1);
    perimeter.setStrokeWeight(6);
    perimeter.setStroke(color(noise.uniformNoise(frameCount / 200f, i), 1, 1));
    perimeter.setFill(false);
    perimeter.setStrokeCap(ROUND);
    shape(perimeter);
  }
}
