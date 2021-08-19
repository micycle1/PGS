import processing.javafx.*;
import java.util.List;
import micycle.pgs.*;
import micycle.pgs.utility.PoissonDistribution;
import micycle.uniformnoise.*;

PShape triangles; // triangle lattice arranged in a hexagon pattern
boolean warp = true;
UniformNoise noise;

void setup() {
  size(800, 800, FX2D);
  smooth();
  textAlign(LEFT, TOP);

  noise = new UniformNoise(1337);
  triangles = PGS_Triangulation.delaunayTriangulation(hexaPoints(width, height, 9, 45));
  triangles = PGS_Morphology.buffer(triangles, -2);
}

void draw() {
  colorMode(RGB, 255, 255, 255, 255);
  background(0, 0, 40);

  fill(0, 255, 255);
  text(frameRate, 2, 2); // fps

  colorMode(HSB, 1, 1, 1, 1);

  int z = 0;
  for (PShape triangle : triangles.getChildren()) {
    final PVector c = PGS_ShapePredicates.centroid(triangle);
    final float hue = noise.uniformNoise(c.x / 500f, c.y / 500f, frameCount * 0.005f); // uniform noise
    float scale = noise(c.x, c.y, frameCount * 0.01f) * 2f; // processing noise

    // consider toggling the lines below:
    triangle = PGS_Transformation.scale(triangle, scale);
    triangle = PGS_Transformation.rotateAroundCenter(triangle, hue * TWO_PI * (z++ % 2 == 0 ? 1 : -1));
    triangle = PGS_Transformation.translate(triangle, (hue - 0.5f) * 100, (-hue + 0.5f) * 100);

    if (warp) {
      triangle = PGS_Morphology.fieldWarp(triangle, 11, .11, frameCount * 0.005f, true, 1337);
      triangle = PGS_Morphology.simplify(triangle, 1);
    }

    PGS_Conversion.setAllFillColor(triangle, color(hue, 1, 0.6f, 0.6f));
    PGS_Conversion.setAllStrokeColor(triangle, color(hue, 1, 1), 4);

    shape(triangle);
  }
}

static List<PVector> hexaPoints(int w, int h, int size, float d) {
  List<PVector> outList = new ArrayList<PVector>();
  int i, j;
  float y, x;
  for (i = 0; i <= (size - 1); i++) {
    y = (sqrt(3) * i * d) / 2.0f;
    for (j = 0; j < (2 * size - 1 - i); j++) {
      x = (-(2 * size - i - 2) * d) / 2.0f + j * d;
      outList.add(new PVector(x + w/2f, y + h/2f));
      if (y != 0) {
        outList.add(new PVector(x + w/2f, -y + h/2f));
      }
    }
  }
  return outList;
}

void mousePressed() {
  warp = !warp;
}

void keyPressed() {
  warp = !warp;
}
