import micycle.pgs.*;
import micycle.uniformnoise.*;
import processing.javafx.*;

UniformNoise noise;

void setup() {
  size(1000, 1000, FX2D);
  fill(50, 125, 50);
  stroke(0, 255, 0);

  noise = new UniformNoise(1337);
  colorMode(HSB, 1, 1, 1, 1);
}

void draw() {
  background(0, 0, 0);
  var origin = new PVector(-width/1.5f, height);
  //origin = new PVector(mouseX,mouseY);
  var points = PGS_Tiling.doyleSpiral(origin.x, origin.y, 9, 99, width*2);
  var s = PGS_Conversion.toPointsPShape(points);
  s = PGS_Transformation.rotate(s, origin, frameCount/500f);
  var points2 = PGS_Conversion.toPVector(s);

  int i = 0;
  for (var c : points2) {
    var cB = points.get(i++);
    float z = cB.z;
    if (c.x+z < 0 || c.y+z < 0 || c.x-z > width || c.y-z > height || z < 10) {
      continue;
    }

    //circle(c.x, c.y, z*2);
    var spiral = PGS_Construction.createLinearSpiral(c.x, c.y, 7, z);
    spiral = PGS_Morphology.simplify(spiral, 0.5f);
    if (noise.uniformNoise(cB.x, cB.y) > 0.5) {
      spiral = PGS_Transformation.flipVertical(spiral, c.x);
    }
    spiral = PGS_Transformation.rotate(spiral, c, noise.uniformNoise(cB.x/500f, cB.y/500f, frameCount/100f)*TWO_PI*2);
    spiral.setStrokeWeight(2);
    float dist = map(c.dist(new PVector(0, height)), 0, width*sqrt(2), 1, 0.5f);
    spiral.setStroke(color(noise.uniformNoise(cB.x/500f, cB.y/500f, frameCount/100f), dist, 1));
    shape(spiral);
  }
}
