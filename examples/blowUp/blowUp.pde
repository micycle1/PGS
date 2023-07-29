import processing.javafx.*;
import micycle.pgs.*;
import micycle.pgs.PGS_Coloring.ColoringAlgorithm;

final int[] palette = new int[]{-3196779, -6237555, -16337744, -68513, -77824, -15985089};
PShape mesh;
int millisLast;
int timePeriod = 4000; // ms

void setup() {
  size(1000, 1000, FX2D);
  frameRate(60);
  mesh = initMesh();
  millisLast = millis();
}

void draw() {
  background(0, 0, 40);
//shape(mesh);
  randomSeed(0);
  var step = (float) smoothStep5(triangleWave(0, 1, millis(), timePeriod)); // oscillating alignment in [0...1]
  PGS_Conversion.getChildren(mesh).forEach(c -> {
    PShape transform;
    transform = PGS_Transformation.translate(c, random(-300, 300), random(-300, 300));
    transform = PGS_Transformation.scale(transform, random(0.5, 2f));
    transform = PGS_Transformation.rotateAroundCenter(transform, random(0, TWO_PI));
    transform = PGS_Transformation.align(transform, c, step);
    transform.setFill(setAlpha(PGS_Conversion.getFillColor(transform), map(step, 1, 0, 0, 255)));
    transform.setStroke(color(0, map(step, 1, 0, 0, 255)));
    shape(transform);
  }
  );

  if (millis() - millisLast >= timePeriod && step > 0.99) {
     mesh = initMesh();
    millisLast = millis();
  }
}

PShape initMesh() {
  var segs = PGS_SegmentSet.graphMatchedSegments(PGS_PointSet.poisson(50, 50, width-50, height-50, random(15, 35)));
  segs = PGS_SegmentSet.filterAxisAligned(segs, radians(2)); // filter highly horizontal or vertical segments
  var segsShape = PGS_SegmentSet.toPShape(segs);

  var mesh = PGS_Voronoi.compoundVoronoi(segsShape);

  mesh = PGS_Meshing.simplifyMesh(mesh, 2, false);
  mesh = PGS_Meshing.stochasticMerge(mesh, choice(3, 5), choice(9999));
  mesh = PGS_Processing.removeSmallHoles(mesh, 1e9); // remove occasional artifact
  mesh = PGS_Coloring.colorMesh(mesh, ColoringAlgorithm.RLF, palette);
  mesh = PGS_Conversion.setAllStrokeColor(mesh, color(0), 2);
  mesh = PGS_Meshing.smoothMesh(mesh, choice(20,200), true);
  

  return mesh;
}

static double triangleWave(double min, double max, double time, double timePeriod) {
  double amplitude = (max - min);
  double phase = (time / timePeriod) % 1.0;
  if (phase < 0.5) {
    return amplitude * (2 * phase) + min;
  } else {
    return amplitude * (2 - 2 * phase) + min;
  }
}

static double smoothStep5(double x) {
  return x * x * x * (x * (x * 6. - 15.) + 10.);
}

public static int setAlpha(int c, float alpha) {
  return (c & 16777215) | (int) alpha << 24;
}
