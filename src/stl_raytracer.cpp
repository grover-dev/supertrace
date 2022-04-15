#include <fstream>
#include <stdlib.h>
#include "shapes.hpp"
#include <string>
#include <iostream>

// void write_to_file()
// {



//   const Vec3 white(255, 255, 255);
//   const Vec3 black(0, 0, 0);
//   const Vec3 red(255, 0, 0);

//   Sphere *spheres = (Sphere*) malloc(2 * sizeof(Sphere));
//   spheres[0] = Sphere(Vec3(W*0.5, H*0.5, 50), 50);
//   spheres[1] = Sphere(Vec3(W*0.8, H*0.8, 50), 30);
  

//   std::ofstream out("out2.ppm");
//   out << "P3\n" << W << ' ' << H << ' ' << "255\n";

//   double t;
//   Vec3 pix_col(black);

//   for (int y = 0; y < H; ++y) {
//     for (int x = 0; x < W; ++x) {
//       pix_col = black;

//       const Ray ray(Vec3(x,y,0), Vec3(0,0,1));

//       for (int z = 0; z < 2; z++) {
//         if (spheres[z].intersect(ray, t)) {
//           const Vec3 pi = ray.o + ray.d*t;
//           const Vec3 L = light.c - pi;
//           const Vec3 N = spheres[z].getNormal(pi);
//           const double dt = dot(L.normalize(), N.normalize());

//           pix_col = (red + white*dt) * brightness;
//           clamp255(pix_col);
//         }
//       }

//       out << (int)pix_col.x << ' '
//           << (int)pix_col.y << ' '
//           << (int)pix_col.z << '\n';
//     }
//   }
//   free(spheres);
// }


void clamp255(Vec3& col)
{
  col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
  col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
  col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}



bool ray_triangle_intersect(struct Ray * ray, struct Triangle * tri, struct Vec3 * intersection_point){
  // error bound for 0
  const float epsilon = 0.000001;

  struct Vec3 c_a_vector = tri->v2 - tri->v0;
  struct Vec3 b_a_vector = tri->v1 - tri->v0;

  struct Vec3 * d_cross_c_a = cross_vec3(ray->d, c_a_vector);

  // first calculating determinant, 
  //  if its ~0 then the ray is parallel to the triangle
  //  if it is <0, then we are hitting the back of the triangle (counting as not intersecting for now)
  //  this will need to be adjusted in the future (especially with refraction) (i.e. use absolute value of det)
  // can therefore ignore it  
  double det = dot_vec3(*d_cross_c_a, b_a_vector); 
  if (det < epsilon){
    return false;
  }
  double inv_det = 1.0 / det;

  struct Vec3 o_a_vector = ray->o - tri->v0;
  struct Vec3 * o_a_cross_b_a = cross_vec3(o_a_vector, b_a_vector);

  // start calculating barycentric coord vectors
  double u = dot_vec3(*d_cross_c_a, o_a_vector) * inv_det;
  // since the vectors are normalized, anything < 0 or > 1 means that the intersection
  // is not in the bounds of the triangle 
  if (u < 0 || u > 1){
    return false;
  }
  double v = dot_vec3(*o_a_cross_b_a, ray->d) * inv_det;
  if (v < 0 || v > 1){
    return false;
  }
  double t = dot_vec3(*o_a_cross_b_a, c_a_vector);
  *intersection_point = ray->o + (ray->d * t);
  return true;
}


#define H 100
#define W 100
#define BRIGHTNESS 0.5

void raytrace(struct STL * stl){
  // creating light source point
  const Sphere light(Vec3(W, H, 0), 1);

  std::ofstream out("out.ppm");
  out << "P3\n" << W << ' ' << H << ' ' << "255\n";

  const Vec3 white(255, 255, 255);
  const Vec3 black(0, 0, 0);
  const Vec3 red(255, 0, 0);

  int length = stl->length;
  Vec3 pix_col(black);

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      pix_col = black;

      Ray ray(Vec3(x,y,0), Vec3(0,0,1));
      for (int i = 0; i < length; i++){
        Vec3 * pi = (Vec3 *)malloc(sizeof(Vec3));
        if(ray_triangle_intersect(&ray, &(stl->triangles[i]), pi)){
            const Vec3 L = light.c - *pi;
            const Vec3 N = stl->triangles->normal;
            const double dt = dot_vec3(L.normalize(), N.normalize());
            pix_col = (red + white*dt) * BRIGHTNESS;
            clamp255(pix_col); 
        }
        free(pi);
      }
      out << (int)pix_col.x << ' '
          << (int)pix_col.y << ' '
          << (int)pix_col.z << '\n';
    }
  }
}


// using Möller-Trumbore algorithm for raytracing w/ triangles 
int main() 
{
  //todo: multiple file loading?
  struct STL * stl;
  std::string filename = "pyramid.stl";
  stl = load_stl(filename);
  std::cout<<"Successfully loaded " <<  filename << "%s\n";
  printf("Number of triangles: %i\n", stl->length);

  raytrace(stl);

  
}
