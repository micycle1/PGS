import processing.javafx.*;
import micycle.pgs.*;
import micycle.uniformnoise.UniformNoise;
import java.util.List;
import org.tinfour.standard.IncrementalTin;

ArrayList<PShape> shapes;
UniformNoise noise;

void setup() {
  size(800, 800, FX2D);
  smooth();
  colorMode(HSB, 1, 1, 1);
  noise = new UniformNoise();
  run();
}

void draw() {
  background(0, 0, 1);

  for (PShape shape : shapes) {
    shape(shape);
  }

  if (frameCount % 90 == 0) {
    run();
  }
}

void run() {
  shapes = new ArrayList<PShape>();
  shapes.add(prepareFaces(randomPoints(round(random(150, 350)), 5, 0, width/2, 0, height/2)));
  shapes.add(prepareFaces(randomPoints(round(random(150, 350)), 5, width/2, width, 0, height/2)));
  shapes.add(prepareFaces(randomPoints(round(random(150, 350)), 5, 0, width/2, height/2, height)));
  shapes.add(prepareFaces(randomPoints(round(random(150, 350)), 5, width/2, width, height/2, height)));
}

ArrayList<PVector> randomPoints(int n, int buffer, float xMin, float xMax, float yMin, float yMax) {
  ArrayList<PVector> randomPoints = new ArrayList<PVector>();
  for (int i = 0; i < n; i++) {
    randomPoints.add(new PVector(random(xMin + buffer, xMax - buffer), random(yMin+buffer, yMax - buffer)));
  }
  return randomPoints;
}

PShape prepareFaces(ArrayList<PVector> points) {
  PShape hull = PGS_Hull.concaveHullBFS2(points, random(0.15, 0.3));
  IncrementalTin mesh = PGS_Triangulation.delaunayTriangulationMesh(hull, points, true, random(1) > 0.5 ? 0 : random(1) > 0.25 ? 1 : 2, true);

  PShape faces;
  if (random(1) > 0.2) {
    faces = PGS_Meshing.urquhartFaces(mesh, true);
  } else {
    faces = PGS_Meshing.gabrielFaces(mesh, true);
  }


  float hueOffset = random(1);
  boolean outline = random(1) > 0.33;
  for (int i = 0; i < faces.getChildCount(); i++) {
    PShape face = faces.getChild(i);
    PVector centroid = PGS_ShapePredicates.centroid(face);
    int fill = color((noise.uniformNoise(centroid.x*0.015, centroid.y*0.015)+hueOffset)%1, random(0.75, 1), random(0.9, 1));
    face.setStroke(outline ? 0 : fill);
    face.setFill(fill);
    face.setStrokeWeight(2);
  }
  return faces;
}
