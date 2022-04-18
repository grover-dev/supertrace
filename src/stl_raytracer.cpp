#include <fstream>
#include <stdlib.h>
#include "shapes.hpp"
#include <string>
#include <iostream>
#include <cmath>

// Update both or find a macro trick
#define FILE_LIST {"sphere.stl"}
#define NUMBER_OF_FILES 1


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
  const float epsilon = 0.0000001;

  struct Vec3 c_a_vector = tri->v2 - tri->v0; //edge2
  struct Vec3 b_a_vector = tri->v1 - tri->v0; //edge1 

  struct Vec3 * d_cross_c_a = cross_vec3(ray->d, c_a_vector); // h = ray cross edge2

  // first calculating determinant, 
  //  if its ~0 then the ray is parallel to the triangle
  //  if it is <0, then we are hitting the back of the triangle (counting as not intersecting for now)
  //  this will need to be adjusted in the future (especially with refraction) (i.e. use absolute value of det)
  // can therefore ignore it  
  double det = dot_vec3(*d_cross_c_a, b_a_vector); // a = (ray cross edge2) => h dot edge1
  if (det < epsilon && det > -epsilon){
    return false;
  }
  double inv_det = 1.0 / det; // f = 1/a

  struct Vec3 o_a_vector = ray->o - tri->v0; // s = ray origin - vertex0

  // start calculating barycentric coord vectors
  double u = dot_vec3(o_a_vector, *d_cross_c_a) * inv_det; // u = s dot h * f
  // since the vectors are normalized, anything < 0 or > 1 means that the intersection
  // is not in the bounds of the triangle 
  if (u < 0.0 || u > 1.0){
    return false;
  }
  struct Vec3 * o_a_cross_b_a = cross_vec3(o_a_vector, b_a_vector); // q = s cross edge1
  double v = dot_vec3(*o_a_cross_b_a, ray->d) * inv_det;
  if (v < 0.0 || (v+u) > 1.0){
    return false;
  }
  double t = dot_vec3(*o_a_cross_b_a, c_a_vector) * inv_det;
  if (t > epsilon){
    *intersection_point = ray->o + (ray->d * t);
    return true;
  }
  return false;
}


#define H 500 // pixel height
#define W 500 // pixel width
#define BRIGHTNESS 0.4
#define SCALING 5.0
#define OFFSET 10.0
#define ZOOM 1

#define STEPS 10

void raytrace(struct STL *stl[], const int number_of_stls, const std::string& filename, float light_angle){
  // creating light source point
  const Sphere light(Vec3(W/2+W*cos(light_angle)/2,H/2+H*sin(light_angle)/2, 0), 1);

  std::ofstream out(filename);
  out << "P3\n" << W << ' ' << H << ' ' << "255\n";

  const Vec3 white(255, 255, 255); // the red will likely need to substituted with surface parameters
  const Vec3 black(0, 0, 0);
  const Vec3 red(0, 255, 0);

  int length;
  Vec3 pix_col(black);
  Vec3 * pi = (Vec3 *)malloc(sizeof(Vec3));
  
  // printf("light_angle: %f\n", light_angle);

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      for (int z = 0; z < number_of_stls; z++) {
        length = stl[z]->length;
        pix_col = black;

        Ray ray(Vec3(x/ZOOM,y/ZOOM,0), Vec3(0,0,1));
        for (int i = 0; i < length; i++){
          if(ray_triangle_intersect(&ray, &(stl[z]->triangles[i]), pi)){
              const Vec3 L = light.c - *pi;
              const Vec3 N = stl[z]->triangles[i].normal;
              const double dt = dot_vec3(L.normalize(), N.normalize());
              pix_col = (red + white*dt) * BRIGHTNESS;
              clamp255(pix_col); 
          }
        }
        // debugging highlighting origin with red square
        if (x <= 10 && y <= 10 ){
          pix_col = Vec3(255,0,0);
        } else if (fabs(x - W*cos(light_angle)) <= 1.1 && fabs(y- H*sin(light_angle) <= 1)){
          pix_col = white;
        }
        // paint y axis green
        if (y == 0){
          pix_col = Vec3(0,255,0);
        // paint x axis blue
        } else if (x == 0){
          pix_col = Vec3(0,0,255);
        }
        out << pix_col.x << ' '
            << pix_col.y << ' '
            << pix_col.z << '\n';
        
      }
    }
  }
  out.close();
  free(pi);
  // for (int z = 0; z<number_of_stls; z++){

  //   // for (int i = 0; i < stl[z]->length; i++){
  //   //   struct Triangle * tri = &(stl[z]->triangles[i]); 
  //   //   free(tri);
  //   // }
  //   free(stl[z]);
  // }
  
}


// using Möller-Trumbore algorithm for raytracing w/ triangles 
int main() 
{
  struct STL *stl[NUMBER_OF_FILES];
  struct STL *objects[NUMBER_OF_FILES];
  const std::string filenames[NUMBER_OF_FILES] = FILE_LIST;
  
  struct Parameters params = Parameters(SCALING, OFFSET, H, W);

  for (int i = 0; i < NUMBER_OF_FILES; i++) {
    stl[i] = load_stl(filenames[i], params);
    std::cout << "Successfully loaded " <<  filenames[i] << "%s\n";
    printf("Number of triangles: %i\n", stl[i]->length);
  }

  std::string output_filename = "output/out.ppm";
  const int start = 0;
  for (int i = start; i < STEPS+start; i++){
    std::string appended_info = std::to_string(i);
    raytrace(stl, NUMBER_OF_FILES, output_filename.insert(10,appended_info), i* M_PI/(float)STEPS);
    output_filename = "output/out.ppm";
  }
}
