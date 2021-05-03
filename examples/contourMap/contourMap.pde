import micycle.pgs.*;
import java.util.List;
import java.util.Map;

PShape polygon;

List<PVector> heights;

float max = -1, min = 9999;

void setup() {
  size(800, 800, FX2D);
  smooth();
}

void draw() {
  background(0, 0, 40);
  populateHeightMap();
  Map<PShape, Float> isolines = PGS_Contour.isolines(heights, 0.08, min, max);

  for (Map.Entry<PShape, Float> entry : isolines.entrySet()) {
    PShape isoline = entry.getKey();
    float isoHeight = entry.getValue();
    isoline.setStroke(color(
        map(isoHeight, min, max, 50, 255), 
        map(isoline.getVertex(0).x, 0, width, 50, 255), 
        map(isoline.getVertex(0).y, 0, height, 50, 255)
       ));
    shape(isoline);
  }
}

void populateHeightMap() {
  heights = new ArrayList<PVector>();

  final int resolution = 15;
  final float animSpeed = 0.005;

  for (int x = 0; x <= width; x+=resolution) {
    for (int y = 0; y <= height; y+=resolution) {
      float z = noise(x*0.01 + frameCount*animSpeed, y*0.01 + frameCount*animSpeed);
      PVector h = new PVector(x, y, 0);
      
      z+=h.dist(new PVector(mouseX, mouseY))*0.005;
      h.z = z;

      heights.add(h);
      max = max(max, z);
      min = min(min, z);
    }
  }
}
