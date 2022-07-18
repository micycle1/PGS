// This example is compatible with Processing 4 only (uses Java 11+ syntax).

import static micycle.pgs.PGS_Conversion.fromPShape;
import static micycle.pgs.PGS_Conversion.toPShape;

import processing.javafx.*;
import micycle.pgs.*;
import micycle.uniformnoise.UniformNoise;
import java.util.List;
import org.locationtech.jts.geom.util.GeometryFixer;
import micycle.pgs.color.Palette;

int particles = 666; // a ceiling value (since points are later pruned by distance)

PFont font;
PShape letter;
List<PVector> targetParticles;
List<PVector> actualParticles;
int[] palette;
UniformNoise noise = new UniformNoise();

void setup() {
  size(1000, 1000, FX2D);
  smooth();
  noStroke();
  textAlign(CENTER, CENTER);

  font = createFont("Kristen ITC", 1100, true);
  actualParticles = PGS_PointSet.gaussian(500, 500, 1, particles); // init
  changeLetter();
}

void draw() {
  fill(palette[0], 50);
  rect(0, 0, width, height);

  if (millis() - time > 2000) { // change letter every 2s
    changeLetter();
  }

  float density = map(noise.uniformNoise(millis()/5000f, 0), 0, 1, 5, 40);
  letterParticles(density);
}

void letterParticles(float density) {
  for (int i = 0; i < targetParticles.size(); i++) { // lerp
    actualParticles.get(i).lerp(targetParticles.get(i), 0.075f);
  }

  var ps = PGS_Conversion.toPointsPShape(actualParticles);
  var warped = PGS_Conversion.toPVector(PGS_Morphology.fieldWarp(ps, 30, 2, millis() / 1000f, false, 1337));

  PGS_PointSet.prunePointsWithinDistance(warped, density).forEach(p -> {
    var character = PGS_Transformation.resizeByMajorAxis(letter, map(noise.uniformNoise(p.x / 50, p.y / 50), 0, 1, 15, 45));
    character = PGS_Transformation.translateTo(character, p.x, p.y);
    character = PGS_Transformation.rotate(character, p, map(noise.uniformNoise((p.x + 1337) / 50, (p.y + 1337) / 50), 0, 1, -PI / 3, PI / 3));
    character.setFill(palette[ceil(noise.uniformNoise(p.x / 350, p.y / 350) * (palette.length - 2) + 0.5f)]);
    shape(character);
  });
}

int time = 0;
void changeLetter() {
  time = millis();
  //noise = new UniformNoise();
  palette = Palette.values()[(int) random(Palette.values().length)].intValue();

  char c = (char) ((int) (random(1000)) % ('z'-'0') + '0'); // random lower case char
  letter = PGS_Transformation.translateTo(font.getShape(c), width / 2, height / 2);
  letter = toPShape(GeometryFixer.fix(fromPShape(letter))); // holes are malformed
  if (letter.getVertexCount() == 0) {
    changeLetter();
    return;
  }
  letter = PGS_Morphology.simplify(letter, 5);
  letter = PGS_Transformation.resizeByMajorAxis(letter, max(width, height)*0.8);
  letter.setStroke(false);

  targetParticles = PGS_Processing.generateRandomPoints(letter, particles);
}
