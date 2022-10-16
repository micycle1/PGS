import processing.javafx.*;
import micycle.pgs.*;
import micycle.trapmap.TrapMap;
import org.jgrapht.traverse.BreadthFirstIterator;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

TrapMap trapMap;
SimpleGraph<PShape, DefaultEdge> tileGraph;

void setup() {
  size(1000, 1000, FX2D);
  colorMode(RGB, 1, 1, 1, 1);

  initTiling();
}

void draw() {
  background(0, 0, 0.15);

  var root = trapMap.findContainingPolygon(mouseX, mouseY);
  var tileIterator = new BreadthFirstIterator<>(tileGraph, root);

  float inc = abs(sin(millis() / 10000f)) / 200f;
  float z = 0;
  while (tileIterator.hasNext()) {
    PShape tile = tileIterator.next();
    float b = (z += inc )%1; // brightness
    tile.setFill(color(b));
    tile.setStroke(false);
    shape(tile);
  }
}

void initTiling() {
  var tiling = PGS_Tiling.squareTriangleTiling(width*0.8, height*0.8, random(15, 25));
  tiling = PGS_Transformation.translateEnvelopeTo(tiling, width/2, height/2);
  tiling = PGS_Transformation.shear(tiling, 1e-5, 0); // shear slightly to make valid for TrapMap

  trapMap = new TrapMap(PGS_Conversion.getChildren(tiling));
  tileGraph = PGS_Conversion.toDualGraph(tiling);
}

void keyPressed() {
 initTiling(); 
}
