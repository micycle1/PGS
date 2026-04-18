import processing.javafx.*;
import micycle.pgs.*;
import micycle.uniformnoise.UniformNoise;

UniformNoise noise;

final int RESOLUTION = 5;
final float CONTOUR_INTERVAL = 0.09;
final float ANIM_SPEED = 0.0035;

// Match these to the actual field range
final float ISOLINE_MIN = -1;
final float ISOLINE_MAX = 1;

float colorMinSmooth = Float.NaN;
float colorMaxSmooth = Float.NaN;

void setup() {
  size(1000, 1000, FX2D);
  smooth(2);
  noise = new UniformNoise();

  colorMode(HSB, 360, 100, 100, 100);
}

void draw() {
  background(230, 65, 10);

  double[] bounds = {0, 0, width, height};

  PShape isolines = PGS_Contour.isolinesFromFunction(
    bounds,
    RESOLUTION,
    CONTOUR_INTERVAL,
    (x, y) -> heightAt((float) x, (float) y),
    ISOLINE_MIN,
    ISOLINE_MAX
  );

  float frameMin = Float.MAX_VALUE;
  float frameMax = -Float.MAX_VALUE;

  for (int i = 0; i < isolines.getChildCount(); i++) {
    PShape line = isolines.getChild(i);
    if (line.getVertexCount() < 2) continue;

    PVector center = averageVertex(line);
    float isoHeight = quantizedHeight(center.x, center.y);

    frameMin = min(frameMin, isoHeight);
    frameMax = max(frameMax, isoHeight);
  }

  if (Float.isNaN(colorMinSmooth)) {
    colorMinSmooth = frameMin;
    colorMaxSmooth = frameMax;
  } else {
    colorMinSmooth = lerp(colorMinSmooth, frameMin, 0.12);
    colorMaxSmooth = lerp(colorMaxSmooth, frameMax, 0.12);
  }

  if (abs(colorMaxSmooth - colorMinSmooth) < 0.0001) {
    colorMaxSmooth = colorMinSmooth + CONTOUR_INTERVAL;
  }

  for (int i = 0; i < isolines.getChildCount(); i++) {
    PShape line = isolines.getChild(i);
    if (line.getVertexCount() < 2) continue;

    PVector center = averageVertex(line);
    float isoHeight = quantizedHeight(center.x, center.y);

    line.setStroke(contourColor(isoHeight, center, colorMinSmooth, colorMaxSmooth));
    line.setStrokeWeight(3);
    shape(line);
  }
}

float heightAt(float x, float y) {
  float t = frameCount * ANIM_SPEED;

  float z = noise.uniformNoise(
    x * 0.0055 + t,
    y * 0.0055 + t,
    2,
    0.5
  );

  float d = dist(x, y, mouseX, mouseY);
  float envelope = exp(-sq(d / 260.0));
  z += 0.15 * envelope * sin(d * 0.025 - frameCount * 0.04);

  return z;
}

float quantizedHeight(float x, float y) {
  float h = heightAt(x, y);
  h = round(h / CONTOUR_INTERVAL) * CONTOUR_INTERVAL;
  return constrain(h, ISOLINE_MIN, ISOLINE_MAX);
}

PVector averageVertex(PShape s) {
  float sx = 0;
  float sy = 0;
  int n = s.getVertexCount();

  for (int i = 0; i < n; i++) {
    PVector v = s.getVertex(i);
    sx += v.x;
    sy += v.y;
  }

  return new PVector(sx / n, sy / n);
}

int contourColor(float isoHeight, PVector center, float colorMin, float colorMax) {
  float t = norm(isoHeight, colorMin, colorMax);
  t = constrain(t, 0, 1);

  float hue = lerp(0, 55, t);

  // subtle spatial drift
  hue += map(center.x, 0, width, -8, 8);
  hue += map(center.y, 0, height, -6, 6);

  hue = (hue % 360 + 360) % 360;

  float sat = lerp(70, 90, t);
  float bri = lerp(75, 100, t);

  return color(hue, sat, bri, 95);
}
