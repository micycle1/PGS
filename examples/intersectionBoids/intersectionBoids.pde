// Boids example adapted from https://processing.org/examples/flocking.html
// This sketch computes a triangulation from points of intersection of a segment set, where segments are formed between pairs of boids
import micycle.pgs.*;
import java.util.List;

final int BOIDS = 100;

Flock flock;
List<PVector> segments;

void setup() {
  size(800, 800, FX2D);
  smooth();
  stroke(0, 200);
  strokeWeight(3);
  textAlign(LEFT, TOP);

  segments = new ArrayList<PVector>();
  flock = new Flock();
  for (int i = 0; i < BOIDS; i++) {
    Boid b = new Boid(random(width), random(height));
    flock.addBoid(b);
    segments.add(b.position);
  }
}

void draw() {
  background(255);

  flock.run();

  List<PVector> intersections = PGS_Processing.lineSegmentIntersections(segments); // compute line intersection points

  PShape t = PGS_Triangulation.delaunayTriangulation(new PShape(), intersections, false, 0, false); // triangulate intersection points
  PGS_Conversion.setAllFillColor(t, color(0, 50, 100, 50));
  PGS_Conversion.setAllStrokeColor(t, color(0, 50, 100), 1);
  shape(t);

  for (int i = 0; i < BOIDS; i+=2) {
    stroke(0, 200, 0, 100);
    line(segments.get(i).x, segments.get(i).y, segments.get(i+1).x, segments.get(i+1).y);
    stroke(0, 100);
    point(segments.get(i).x, segments.get(i).y);
    point(segments.get(i+1).x, segments.get(i+1).y);
  }

  noStroke();
  fill(237, 50, 162, 150);
  for (PVector x : intersections) {
    ellipse(x.x, x.y, 10, 10);
  }

  fill(0);
  text("Intersections: " + intersections.size(), 0, 0);
}
