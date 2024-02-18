import processing.javafx.*;
import java.util.List;
import micycle.pgs.*;
import micycle.uniformnoise.*;
import java.util.concurrent.atomic.AtomicInteger;
import net.jafama.FastMath;


PShape triangles; // triangle lattice arranged in a hexagon pattern
boolean warp = true;
UniformNoise noise;

void setup() {
  size(1000, 1000, FX2D);
  smooth();
  textAlign(LEFT, TOP);

  noise = new UniformNoise(1337);
  triangles = PGS_Triangulation.delaunayTriangulation(PGS_PointSet.hexagon(width/2, height/2, 19, 30));
  triangles = PGS_Morphology.buffer(triangles, -2);
  frameRate(600);
}

void draw() {
  colorMode(RGB, 255, 255, 255, 255);
  //background(0, 0, 40);
  fill(0,0,40,33);
  rect(0,0,width,height);

  fill(0, 255, 255);
  //text(frameRate, 2, 2); // fps

  colorMode(HSB, 1, 1, 1, 1);

  final AtomicInteger z = new AtomicInteger();
  PGS_Conversion.getChildren(triangles).stream().map(triangle -> {
      
    final PVector c = PGS_ShapePredicates.centroid(triangle);
    final float hue = noise.uniformNoise(c.x / 500f, c.y / 500f, millis() * 0.0005f); // uniform noise
    float scale = noise(c.x, c.y, frameCount * 0.01f) * 2f; // processing noise

    // consider toggling the lines below:
    triangle = PGS_Transformation.scale(triangle, scale);
    triangle = PGS_Transformation.rotateAroundCenter(triangle, hue * TWO_PI * (z.getAndIncrement() % 2 == 0 ? 1 : -1));
    triangle = PGS_Transformation.translate(triangle, (hue - 0.5f) * 100, (-hue + 0.5f) * 100);

    if (warp) {
      triangle = PGS_Morphology.fieldWarp(PGS_Processing.densify(triangle, 5), 11, .11, millis() * 0.0005f, !true, 1337);
      triangle = PGS_Morphology.simplify(triangle, 0.75);
    }


    PGS_Conversion.setAllFillColor(triangle, setAlpha(sinebow(hue), 85)); // color(hue, 1, 0.6f, 0.6f)
    PGS_Conversion.setAllStrokeColor(triangle, sinebow(hue), 3);

    return triangle;
  }).sequential().forEach(t -> shape(t));
}

void mousePressed() {
  warp = !warp;
}

void keyPressed() {
  warp = !warp;
}

  public static int setAlpha(int c, int alpha) {
    return (c & 16777215) | alpha << 24;
  }

  public static int sinebow(double t) {
    t = 0.5f - t;
    return rgbColor(255 * (sin2(t + 0 / 3f)),
        255 * (sin2(t + 1 / 3d)), 255 * (sin2(t + 2 / 3d)));
  }
  
    private static double sin2(double t) {
    double z = FastMath.sin(PI * t);
    return z * z;
  }
  
    public static int rgbColor(double red, double green, double blue) {
    return -16777216 | FastMath.roundToInt(red) << 16 | FastMath.roundToInt(green) << 8 | FastMath.roundToInt(blue);
  }
