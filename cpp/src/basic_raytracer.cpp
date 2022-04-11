
#include <fstream>
#include <stdlib.h>
#include "shapes.hpp"



void clamp255(Vec3& col) {
  col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
  col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
  col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}

int main() {

  const int H = 500;
  const int W = 500;
  const int brightness = 1;


  const Vec3 white(255, 255, 255);
  const Vec3 black(0, 0, 0);
  const Vec3 red(255, 0, 0);

  Sphere *spheres = (Sphere*) malloc(2 * sizeof(Sphere));
  spheres[0] = Sphere(Vec3(W*0.5, H*0.5, 50), 50);
  spheres[1] = Sphere(Vec3(W*0.8, H*0.8, 50), 30);
  
  const Sphere light(Vec3(W, H, 0), 1);

  std::ofstream out("out2.ppm");
  out << "P3\n" << W << ' ' << H << ' ' << "255\n";

  double t;
  Vec3 pix_col(black);

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      pix_col = black;

      const Ray ray(Vec3(x,y,0), Vec3(0,0,1));

      for (int z = 0; z < 2; z++) {
        if (spheres[z].intersect(ray, t)) {
          const Vec3 pi = ray.o + ray.d*t;
          const Vec3 L = light.c - pi;
          const Vec3 N = spheres[z].getNormal(pi);
          const double dt = dot(L.normalize(), N.normalize());

          pix_col = (red + white*dt) * brightness;
          clamp255(pix_col);
        }
      }

      out << (int)pix_col.x << ' '
          << (int)pix_col.y << ' '
          << (int)pix_col.z << '\n';
    }
  }
  free(spheres);
}
