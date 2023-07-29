import processing.javafx.*;
import java.util.List;
import micycle.pgs.*;
import micycle.uniformnoise.*;
import java.util.concurrent.atomic.AtomicInteger;

PShape triangles; // triangle lattice arranged in a hexagon pattern
boolean warp = true;
UniformNoise noise;

void setup() {
  size(800, 800, FX2D);
  smooth();
  textAlign(LEFT, TOP);

  noise = new UniformNoise(1337);
  triangles = PGS_Triangulation.delaunayTriangulation(PGS_PointSet.hexagon(width/2, height/2, 9, 45));
  triangles = PGS_Morphology.buffer(triangles, -2);
  frameRate(600);
}

void draw() {
  colorMode(RGB, 255, 255, 255, 255);
  background(0, 0, 40);

  fill(0, 255, 255);
  text(frameRate, 2, 2); // fps

  colorMode(HSB, 1, 1, 1, 1);

  final AtomicInteger z = new AtomicInteger();
  PGS_Conversion.getChildren(triangles).stream().map(triangle -> {
      
    final PVector c = PGS_ShapePredicates.centroid(triangle);
    final float hue = noise.uniformNoise(c.x / 500f, c.y / 500f, millis() * 0.0005f); // uniform noise
    float scale = noise(c.x, c.y, frameCount * 0.01f) * 2f; // processing noise

    // consider toggling the lines below:
    triangle = PGS_Transformation.scale(triangle, scale);
    triangle = PGS_Transformation.rotateAroundCenter(triangle, hue * TWO_PI * (z.getAndIncrement() % 2 == 0 ? 1 : -1));
    triangle = PGS_Transformation.translate(triangle, (hue - 0.5f) * 100, (-hue + 0.5f) * 100);

    if (warp) {
      triangle = PGS_Morphology.fieldWarp(triangle, 11, .11, millis() * 0.0005f, true, 1337);
      triangle = PGS_Morphology.simplify(triangle, 1);
    }

    PGS_Conversion.setAllFillColor(triangle, color(hue, 1, 0.6f, 0.6f));
    PGS_Conversion.setAllStrokeColor(triangle, color(hue, 1, 1), 4);

    return triangle;
  }).sequential().forEach(t -> shape(t));
}

void mousePressed() {
  warp = !warp;
}

void keyPressed() {
  warp = !warp;
}
