import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;

PShape polygon;
List<PShape> subPartitions;

void setup() {
  size(800, 800, FX2D);
  smooth();

  List<PVector> randomPoints = PGS_PointSet.poisson(30, 30, width - 30, height - 30, 40);
  polygon = PGS_Hull.concaveHullBFS(randomPoints, 0.1);

  PShape partitions = PGS_Processing.convexPartition(polygon);
  subPartitions = new ArrayList<PShape>();
  for (PShape p : partitions.getChildren()) {
    PShape split = PGS_Processing.split(p);
    subPartitions.add(split.getChild(0));
    subPartitions.add(split.getChild(1));
  }
}

void draw() {
  background(0, 0, 40);
  for (PShape p : subPartitions) {

    PVector centroid = PGS_ShapePredicates.centroid(p);
    if (centroid == null) {
      continue;
    }

    float smooth= map(centroid.x+centroid.y, 0, width+height, 0, 01);
    smooth = triangleWave(1, (smooth + frameCount*0.01) % 2);
    int fill = color(triangleWave(255, frameCount), 255*centroid.x/width, 255*centroid.y/height);
    
    p = PGS_Transformation.scale(p, 1-smooth*0.4);
    p = PGS_Morphology.smooth(p, smooth);

    p.setStrokeWeight(2);
    p.setStroke(255);
    p.setFill(fill);
    shape(p);
  }
}

static float triangleWave(float max, float time) {
  return max - abs(time % (2 * max) - max);
}
