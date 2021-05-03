// This example demos: medial axis, mink sum, shear, intersection
import micycle.pgs.*;
import java.util.List;

Letter l1, l2;

void setup() {
  size(800, 800, FX2D);
  colorMode(HSB, 1);

  l1 = new Letter('M');
  l2 = new Letter('L');
}

void draw() {
  fill(0.1, 0.2);
  rect(-1, -1, width+1, height+1);

  try {
    l1.update();
    l2.update();

    PShape mink = PGS_Morphology.minkSum(l1.letter, l2.letter);
    mink = PGS_Transformation.translateTo(mink, (l1.pos.x+l2.pos.x) / 2, (l1.pos.y+l2.pos.y) / 2);
    shape(mink);

    shape(l1.letter);
    shape(l2.letter);

    PShape z =  PGS_Contour.medialAxis(mink, 0.3, 0, 0.1);
    shape(z);

    PShape intersect = PGS_ShapeBoolean.intersect(l1.letter, mink);
    PGS_Conversion.setAllFillColor(intersect, color(0, 0.5));
    shape(intersect);
    intersect = PGS_ShapeBoolean.intersect(l2.letter, mink);
    PGS_Conversion.setAllFillColor(intersect, color(0, 0.5));
    shape(intersect);
  }
  catch (Exception e) {
  }

  if (frameCount % 120 == 0) {
    l1.randomise();
    l2.randomise();
  }
}

class Letter {

  PVector pos;
  PShape letter;
  private PShape l;
  float xn = random(4096);
  float yn = random(4096); 

  private float hue; 

  Letter(char c) {
    pos = new PVector(random(width), random(height));
    String randomFont = PFont.list()[round(random(PFont.list().length))];
    PFont font = createFont(randomFont, 96, true);
    l = font.getShape(c);
    hue = random(1);
  }

  void update() {

    xn+=.005;
    yn+=.005;
    pos.x = noise(xn)*width;
    pos.y = noise(yn)*height;

    letter = PGS_Transformation.translateTo(l, 0, 0 );
    letter = PGS_Transformation.shear(letter, map(pos.x, 0, width, -TWO_PI, TWO_PI), map(pos.y, 0, height, -TWO_PI, TWO_PI));
    letter = PGS_Transformation.translateTo(letter, pos.x, pos.y);
    letter = PGS_Morphology.simplify(letter, 1); // as some fonts have very dense vertices 
    letter.setStroke(color(hue, 1, 1));
  }

  void randomise() {
    hue = random(1);
    String randomFont = PFont.list()[(int)random(PFont.list().length)];
    PFont font = createFont(randomFont, 128, true);
    l = font.getShape((char)random(48, 91));
  }
}
