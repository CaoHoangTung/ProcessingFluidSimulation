final int N = 128;
final int SCALE = 4;
final int iter = 16;

Fluid fluid;

void settings() {
  size(N*SCALE, N*SCALE);
}

void setup() {
  fluid = new Fluid(0.1, 0, 0.0001);
}

void mouseDragged() {
  fluid.addDensity(mouseX/SCALE, mouseY/SCALE, 1000);
  float atmX = mouseX - pmouseX;
  float atmY = mouseY - pmouseY;
  fluid.addVelocity(mouseX/SCALE, mouseY/SCALE, atmX, atmY);
}

void draw() {
  background(0);
  fluid.step();
  fluid.render();
}
