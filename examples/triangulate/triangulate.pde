import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;
import micycle.pgs.color.ColorUtils;

PShape polygon;
PShape polygonSmooth;
PShape triangles;

void setup() {
  size(800, 800, FX2D);
  smooth();
  colorMode(HSB, 1, 1, 1);
  frameRate(60);
  makeTriangulation();
}

void draw() {
  background(0, 0, 1);

  if (frameCount % (60*4) == 0) {
    makeTriangulation();
  }

  PVector origin = PGS_Processing.pointOnExterior(polygonSmooth, frameCount*0.004f, 0);
  for (PShape triangle : triangles.getChildren()) {
    PVector centroid = PGS_ShapePredicates.centroid(triangle);
    float dist = map((origin.dist(centroid) + frameCount) % width, 0, width, 0, 1);
    dist+=(noise(centroid.x, centroid.y)-0.5)/3;
    int fill = color(dist, map(centroid.x, 0, width, 0.5, 0.75), noise(centroid.x/10, centroid.y/10)+0.2);
    
    triangle.setFill(fill);
    triangle.setStrokeWeight(1);
  }

  PGS_Conversion.setAllStrokeToFillColor(triangles, 1);

  shape(triangles);
  shape(polygon);
}

void makeTriangulation() {
  polygon = PGS_Construction.createSuperRandomPolygon(width*0.85, (int)random(3, 50), random(0.01f, 0.25f), (int)random(1, 4), (int)random(0, 3), random(1) > 0.5, random(1) > 0.5, 1337);
  polygon = PGS_Transformation.translateEnvelopeTo(polygon, width/2, height/2);
  polygon = PGS_Morphology.simplify(polygon, 2);
  polygon.setStroke(0);
  polygon.setStrokeWeight(4);
  polygon.setFill(false);

  polygonSmooth = PGS_Morphology.smoothEllipticFourier(polygon, 6);

  triangles = PGS_Triangulation.delaunayTriangulation(polygon, PGS_PointSet.poisson(0, 0, width, height, random(33, 99)), true, 1, true);
}
