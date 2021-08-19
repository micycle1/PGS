import processing.javafx.*;
import micycle.pgs.*;
import micycle.pgs.utility.PoissonDistribution;
import java.util.List;

PShape polygon;
List<PShape> subPartitions;

void setup() {
  size(800, 800, FX2D);
  smooth();

  List<PVector> randomPoints = new PoissonDistribution().generate(30, 30, width - 30, height - 30, 35, 7);
  polygon = PGS_Processing.concaveHullBFS(randomPoints, 25);

  List<PShape> partitions = PGS_Processing.partition(polygon);
  subPartitions = new ArrayList<PShape>();
  for (PShape p : partitions) {
    subPartitions.addAll(PGS_Processing.split(p));
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
