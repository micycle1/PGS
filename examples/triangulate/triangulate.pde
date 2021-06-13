import micycle.pgs.*;
import java.util.List;

PShape polygon;

void setup() {
  size(800, 800, FX2D);
  smooth();
  colorMode(HSB, 1, 1, 1);
  polygon = PGS_Construction.createRandomPolygonExact((int) random(3, 7), width, height);
}

void draw() {
  background(0, 0, 1);

  if (frameCount % 180 == 0) {
    polygon = PGS_Construction.createRandomPolygonExact((int) random(3, 7), width, height);
  }

  List<PVector> trianglePoints = PGS_Triangulation.poissonTriangulationPoints(polygon, 40);

  beginShape(TRIANGLES);
  strokeWeight(2);
  stroke(0);
  for (int i = 0; i < trianglePoints.size(); i += 3) {
    PVector origin = PGS_Processing.pointOnExterior(polygon, frameCount*0.004f, 0);
    float dist = map((origin.dist(new PVector(trianglePoints.get(i).x, trianglePoints.get(i).y)) + frameCount) % width, 0, width, 0, 1);

    fill(dist, map(trianglePoints.get(i).x, 0, width, 0.5, 1), 1);
    vertex(trianglePoints.get(i).x, trianglePoints.get(i).y);
    vertex(trianglePoints.get(i + 1).x, trianglePoints.get(i + 1).y);
    vertex(trianglePoints.get(i + 2).x, trianglePoints.get(i + 2).y);
  }
  endShape();
}
