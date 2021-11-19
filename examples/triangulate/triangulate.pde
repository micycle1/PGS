import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;

PShape polygon;
PShape triangles;

void setup() {
  size(800, 800, FX2D);
  smooth();
  colorMode(HSB, 1, 1, 1);
  makeTriangulation();
}

void draw() {
  background(0, 0, 1);

  if (frameCount % 180 == 0) {
    makeTriangulation();
  }

  final PVector origin = PGS_Processing.pointOnExterior(polygon, frameCount*0.004f, 0);
  for (PShape triangle : triangles.getChildren()) {
    PVector centroid = PGS_ShapePredicates.centroid(triangle);
    float dist = map((origin.dist(centroid) + frameCount) % width, 0, width, 0, 1);
    triangle.setFill(color(dist, map(centroid.x, 0, width, 0.5, 1), round(centroid.x) % 2 == 0 ? 0.9 : 0.7));
    triangle.setStrokeWeight(1);
  }
  PGS_Conversion.setAllStrokeToFillColor(triangles);

  shape(polygon);
  shape(triangles);
}

void makeTriangulation() {
  polygon = PGS_Construction.createRandomPolygon((int) random(4, 9), width, height);
  polygon.setStroke(0);
  polygon.setStrokeWeight(10);
  triangles = PGS_Triangulation.delaunayTriangulation(polygon, PGS_PointSet.poisson(0, 0, width, height, random(33, 99)), true, 1, true);
}
