import processing.javafx.*;
import micycle.pgs.*;
import java.util.List;
import java.util.Collection;
import micycle.pgs.color.Palette;
import micycle.pgs.PGS_Coloring.ColoringAlgorithm;

void setup() {
  size(1000, 1000, FX2D);
  smooth();
}

void draw() {
  shape(interlock(50+mouseX/2, 1337));
}

PShape interlock(double r, long seed) {
  var palette = Palette.DAILY;
  int k = 1;
  background(palette.get(k));

  var plane = PGS_Construction.createRect(0, 0, width, height, 0);
  var sites = PGS_Optimisation.circleCoverage(plane, 15, seed);
  var circles = PGS_PointSet.applyRandomWeights(sites, r, r + r/2, seed);
  var shapes = PGS_Conversion.toCircles(circles);
  shapes = PGS_ShapeBoolean.unionLines(shapes, null); // union and polygonise its own linework

  var classes = PGS_Coloring.colorMesh(shapes, ColoringAlgorithm.RLF); // colour classes

  var interlock = classes.entrySet().parallelStream().map(e -> {
    var s = e.getKey(); // shape
    var i = e.getValue().intValue(); // colour class

    s = PGS_Morphology.buffer(s, -5);
    if (PGS_ShapePredicates.area(s) < 200) {
      return null;
    }
    s = PGS_Morphology.simplifyDCE(s, (a, b, c) -> {
      return b > 30;
    }
    );

    s = PGS_Voronoi.innerVoronoi(s, PGS_Processing.generateRandomPoints(s, 33, seed), 10);
    s = PGS_Meshing.areaMerge(s, 4);
    s = PGS_Meshing.smoothMesh(s, 1d, true);
    s = PGS_Morphology.chaikinCut(s, 0.5, 3);

    PGS_Conversion.setAllFillColor(s, palette.get(k+i * 2));
    PGS_Conversion.setAllStrokeColor(s, palette.get(k+i * 2 + 1), 2);
    return s;
  }
  ).toList();

  return PGS_Conversion.flatten(interlock);
}
