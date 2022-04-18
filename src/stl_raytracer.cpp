#include <fstream>
#include <stdlib.h>
#include "shapes.hpp"
#include <string>
#include <iostream>
#include <cmath>
#include <time.h>

// Update both or find a macro trick
#define FILE_LIST {"sphere.stl"}
#define NUMBER_OF_FILES 1
#define DEBUG_MODE true

double CLOCK() {
        struct timespec t;
        clock_gettime(CLOCK_MONOTONIC,  &t);
        return (t.tv_sec * 1000)+(t.tv_nsec*1e-6);
}


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


void clamp_pixels(Vec3& col)
{
  col.x = (col.x > 1023) ? 1023 : (col.x < 0) ? 0 : col.x;
  col.y = (col.y > 1023) ? 1023 : (col.y < 0) ? 0 : col.y;
  col.z = (col.z > 1023) ? 1023 : (col.z < 0) ? 0 : col.z;
}



bool ray_triangle_intersect(struct Ray * ray, struct Triangle * tri, struct Vec3 * intersection_point){
  // error bound for 0
  const float epsilon = 0.0000001;

  struct Vec3 c_a_vector = tri->v2 - tri->v0; //edge2
  struct Vec3 b_a_vector = tri->v1 - tri->v0; //edge1 

  struct Vec3 *d_cross_c_a = cross_vec3(ray->d, c_a_vector); // h = ray cross edge2

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
    free(o_a_cross_b_a);
    free(d_cross_c_a);
    return false;
  }
  double t = dot_vec3(*o_a_cross_b_a, c_a_vector) * inv_det;
  if (t > epsilon){
    *intersection_point = ray->o + (ray->d * t);
    free(o_a_cross_b_a);
    free(d_cross_c_a);
    return true;
  }
  free(o_a_cross_b_a);
  free(d_cross_c_a);
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
  double light_source_x = W/2+W*cos(light_angle)/2;
  double light_source_y = H/2+H*sin(light_angle)/2;
  double light_source_z = 500.0;
  const Sphere light(Vec3(light_source_x,light_source_y,light_source_z ), 1);

  std::ofstream out(filename);
  out << "P3\n" << W << ' ' << H << ' ' << "1023\n";

  const Vec3 white(1023, 1023, 1023); // the red will likely need to substituted with surface parameters
  const Vec3 black(0, 0, 0);
  const Vec3 red(1023, 0, 0);

  Vec3 pix_col(black);
  Vec3 *pi = (Vec3 *)malloc(sizeof(Vec3));
  

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      for (int z = 0; z < number_of_stls; z++) {
        pix_col = black;

        Ray ray(Vec3(x/ZOOM,y/ZOOM,0), Vec3(0,0,1));
        for (int i = 0; i < stl[z]->length; i++){
          if(ray_triangle_intersect(&ray, &(stl[z]->triangles[i]), pi)){
              const Vec3 L = light.c - *pi;
              const Vec3 N = stl[z]->triangles[i].normal;
              const double dt = dot_vec3(L.normalize(), N.normalize());
              pix_col = (red + white*dt) * BRIGHTNESS;
              clamp_pixels(pix_col); 
          }
        }
        // debugging highlighting origin with red square
        if (x <= 10 && y <= 10 ){
          pix_col = Vec3(1023,0,0);
        } else if (fabs(x - light_source_x) <= 1 && fabs(y-light_source_y ) <= 1){
          pix_col = white;
        }
        // paint y axis green
        if (y == 0){
          pix_col = Vec3(0,1023,0);
        // paint x axis blue
        } else if (x == 0){
          pix_col = Vec3(0,0,1023);
        }
        out << (int) pix_col.x << ' '
            << (int) pix_col.y << ' '
            << (int) pix_col.z << '\n';
        
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

  double start_time, finish_time, total_time, current_time, increment_point;
  
  struct Parameters params = Parameters(SCALING, OFFSET, H, W);

  for (int i = 0; i < NUMBER_OF_FILES; i++) {
    stl[i] = load_stl(filenames[i], params);
    std::cout << "Successfully loaded " <<  filenames[i] << "%s\n";
    printf("Number of triangles: %i\n", stl[i]->length);
  }

  std::string output_filename = "output/out.ppm";
  const int start = 0;
  start_time = CLOCK();
  increment_point = start_time;
  for (int i = start; i < STEPS+start; i++){
    std::string appended_info = std::to_string(i);
    raytrace(stl, NUMBER_OF_FILES, output_filename.insert(10,appended_info), i* M_PI/(float)STEPS);
    output_filename = "output/out.ppm";
    if (DEBUG_MODE) {
      current_time = CLOCK();
      printf("Frame %d processed in %f ms\n", i, current_time-increment_point);
      increment_point = current_time;
    }
  }
  if (DEBUG_MODE) {
    finish_time = CLOCK();
    total_time = finish_time-start_time;
    printf("The total time to raytrace was: %f ms\n", total_time);
  }
}
