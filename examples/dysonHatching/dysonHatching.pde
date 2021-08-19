import processing.javafx.*;
import micycle.pgs.*;
import micycle.pgs.utility.PoissonDistribution;
import java.util.List;

PShape polygon;
PShape cells;
boolean drawCellBorders = false;

void setup() {
  size(800, 800, FX2D);
  smooth();
  colorMode(HSB, 1, 1, 1, 1);

  List<PVector> randomPoints = new PoissonDistribution().generate(30, 30, width - 30, height - 30, 35, 7);
  polygon = PGS_Processing.concaveHullBFS(randomPoints, 25);
  cells = PGS_Voronoi.voronoiCells(PGS_Processing.generateRandomPoints(polygon, 200));
  PGS_Conversion.disableAllFill(cells);
  PGS_Conversion.setAllStrokeColor(cells, color(0, 0, 1, 1), 1);
}

void draw() {
  background(240/360f, 1, .16);
  text(frameRate, 5, 15);

  for (PShape cell : cells.getChildren()) {
    PVector center = PGS_ShapePredicates.centroid(cell);
    float angle = noise(cell.getVertex(0).x/width, cell.getVertex(2).y/height, frameCount*0.01f)*TWO_PI;

    PShape envelope = PGS_Optimisation.envelope(cell);
    float maxL = envelope.getVertex(0).dist(envelope.getVertex(2));
    float l = maxL/2f; // line length
    float d = 12; // distance between lines
    d-=map(dist(center.x, center.y, mouseX, mouseY), 0, width * 0.6f, 6, 0);

    final int linesN = ceil(l/d);
    float dx = cos(angle + HALF_PI) * d;
    float dy = sin(angle + HALF_PI) * d;

    PShape lines = createShape(); // cell hatching
    lines.setStroke(true);
    lines.setStrokeWeight(1);
    lines.setStrokeCap(ROUND);
    lines.beginShape(PConstants.LINES);
    PVector a, b;
    for (int i = 0; i < linesN; i++) {
      a = PVector.add(center, new PVector(dx * i, dy * i));
      b = PVector.add(a, new PVector(cos(angle) * l, sin(angle) * l));
      a.add(cos(angle) * -l, sin(angle) * -l);
      lines.vertex(a.x, a.y);
      lines.vertex(b.x, b.y);

      a = PVector.add(center, new PVector(dx * -i, dy * -i));
      b = PVector.add(a, new PVector(cos(angle) * l, sin(angle) * l));
      a.add(cos(angle) * -l, sin(angle) * -l);
      lines.vertex(a.x, a.y);
      lines.vertex(b.x, b.y);
    }
    lines.endShape();

    PShape o = PGS_ShapeBoolean.intersect(cell, lines); // crop lines to the cell
    final int col = color((noise(2*cell.getVertex(1).y/width, 2*cell.getVertex(2).x/height) + frameCount*0.002f - map(dist(center.x, center.y, mouseX, mouseY), 0, width * 0.75f, 0, .2f)) % 1, 1, .9f, 1);
    PGS_Conversion.setAllStrokeColor(o, col, noise(cell.getVertex(0).x/width, cell.getVertex(2).y/height, frameCount*0.01f) > 0.5 ? 2 : 4);
    shape(o);
  }

  if (drawCellBorders) {
    shape(cells);
  }
}

void keyPressed() {
  drawCellBorders = !drawCellBorders;
}
