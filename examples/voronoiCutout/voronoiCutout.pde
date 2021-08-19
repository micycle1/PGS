// This example shows concave hull, constrained voronoi, shape subtract & intersection, slicing, pointsOnExterior, shape intersection points
import processing.javafx.*;
import micycle.pgs.*;
import micycle.pgs.utility.PoissonDistribution;
import java.util.List;

PShape polygon;
List<PShape> subPartitions;
PVector center;

void setup() {
  size(800, 800, FX2D);
  smooth();
  colorMode(HSB, 1);
  stroke(color(324f/360f, 0.79, 0.93));
  strokeWeight(8);
  center = new PVector(width/2, height/2);

  List<PVector> randomPoints = new PoissonDistribution(1337).generate(30, 30, width - 30, height - 30, 15, 7);
  polygon = PGS_Processing.concaveHull2(randomPoints, 0);

  polygon.setFill(false);
  polygon.setStroke(color(1));
  polygon.setStrokeWeight(3);

  List<PShape> partitions = PGS_Processing.partition(polygon);
  subPartitions = new ArrayList<PShape>();
  for (PShape p : partitions) {
    subPartitions.addAll(PGS_Processing.split(p));
  }
}

void draw() {
  background(0, 0, 0.075);

  PShape star = PGS_Construction.createStar(width/2, height/2, 11, mouseX*0.4, mouseX*0.4+100, 1);
  star = PGS_Transformation.rotateAroundCenter(star, frameCount*0.01f);
  star = PGS_Morphology.simplify(star, 1);
  PShape smallStar = PGS_Transformation.scale(star, 0.5);

  PShape outer =  PGS_ShapeBoolean.subtract(polygon, star);
  PShape inner = PGS_ShapeBoolean.subtract(star, polygon);
  inner = PGS_ShapeBoolean.subtract(inner, smallStar);
  PGS_Conversion.setAllFillColor(inner, color(.45, .5, .9));
  shape(inner);

  for (int i = 0; i < outer.getChildCount(); i++) {
    try {
      PShape p = PGS_Voronoi.voronoiDiagram(outer.getChild(i), true);
      PGS_Conversion.setAllStrokeColor(p, color(i / ((float) outer.getChildCount() - 1), 1, 1), 2);
      shape(p);
    }
    catch (Exception e) {
    }
  }

  PShape polygon2 = PGS_ShapeBoolean.subtract(polygon, smallStar);

  PGS_Conversion.disableAllFill(polygon2);
  PGS_Conversion.setAllStrokeColor(polygon2, color(1), 3);
  shape(polygon2);

  PShape innerInner = PGS_ShapeBoolean.intersect(polygon, smallStar);
  List<PVector> ps = PGS_Processing.pointsOnExterior(star, 2, 0);
  innerInner = PGS_Processing.slice(innerInner, ps.get(0), ps.get(1)).get(0);
  PGS_Conversion.disableAllStroke(innerInner);
  shape(innerInner);

  List<PVector> intersections = PGS_Processing.shapeIntersection(polygon, inner);
  for (PVector x : intersections) {
    point(x.x, x.y);
  }
}
