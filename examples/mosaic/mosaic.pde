import processing.javafx.*;
import micycle.pgs.*;
import micycle.uniformnoise.UniformNoise;
import java.util.List;
import org.tinfour.common.IIncrementalTin;

List<PShape> shapes;
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
  long seed = System.currentTimeMillis();
  float gutter = 5;

  shapes.add(prepareFaces(PGS_PointSet.random(0+gutter, 0+gutter, width/2-gutter, height/2-gutter, round(random(150, 350)), seed)));
  shapes.add(prepareFaces(PGS_PointSet.random(0+gutter, height/2+gutter, width/2-gutter, height-gutter, round(random(150, 350)), seed+1)));
  shapes.add(prepareFaces(PGS_PointSet.random(width/2+gutter, 0+gutter, width-gutter, height/2-gutter, round(random(150, 350)), seed+2)));
  shapes.add(prepareFaces(PGS_PointSet.random(width/2+gutter, height/2+gutter, width-gutter, height-gutter, round(random(150, 350)), seed+3)));
}

PShape prepareFaces(List<PVector> points) {
  PShape hull = PGS_Hull.concaveHullBFS(points, random(0.25, 0.5));
  IIncrementalTin mesh = PGS_Triangulation.delaunayTriangulationMesh(hull, points, true, random(1) > 0.5 ? 0 : random(1) > 0.25 ? 1 : 2, true);

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
