#include <fstream>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdint.h>
#include "shapes.hpp"
#include <cstring>
#include <omp.h>

 double CLOCK() {
        struct timespec t;
        clock_gettime(CLOCK_MONOTONIC,  &t);
        return (t.tv_sec * 1000)+(t.tv_nsec*1e-6);
}

#define MAX_PIXEL 1023
void clamp_pixels(Vec3& col)
{
  col.x = (col.x > MAX_PIXEL) ? MAX_PIXEL : (col.x < 0) ? 0 : col.x;
  col.y = (col.y > MAX_PIXEL) ? MAX_PIXEL : (col.y < 0) ? 0 : col.y;
  col.z = (col.z > MAX_PIXEL) ? MAX_PIXEL : (col.z < 0) ? 0 : col.z;
}



bool ray_triangle_intersect(struct Ray * ray, struct Triangle * tri, struct Vec3 * intersection_point){
  // error bound for 0
  bool intersect = true;
  const float epsilon = 0.00000001;

  struct Vec3 c_a_vector = tri->v2 - tri->v0; //edge2
  struct Vec3 b_a_vector = tri->v1 - tri->v0; //edge1 

  struct Vec3 d_cross_c_a = cross_vec3(ray->d, c_a_vector); // h = ray cross edge2

  // first calculating determinant, 
  //  if its ~0 then the ray is parallel to the triangle
  //  if it is <0, then we are hitting the back of the triangle (counting as not intersecting for now)
  //  this will need to be adjusted in the future (especially with refraction) (i.e. use absolute value of det)
  // can therefore ignore it  
  double det = dot_vec3(d_cross_c_a, b_a_vector); // a = (ray cross edge2) => h dot edge1
  if (det < epsilon && det > -epsilon){
    // free(d_cross_c_a);
    return false;
  }
  double inv_det = 1.0 / det; // f = 1/a

  struct Vec3 o_a_vector = ray->o - tri->v0; // s = ray origin - vertex0

  // start calculating barycentric coord vectors
  double u = dot_vec3(o_a_vector, d_cross_c_a) * inv_det; // u = s dot h * f
  // since the vectors are normalized, anything < 0 or > 1 means that the intersection
  // is not in the bounds of the triangle 
  if (u < 0.0 || u > 1.0){
    intersect = intersect * 0;
  }
  struct Vec3 o_a_cross_b_a = cross_vec3(o_a_vector, b_a_vector); // q = s cross edge1
  double v = dot_vec3(o_a_cross_b_a, ray->d) * inv_det;
  if (v < 0.0 || (v+u) > 1.0){
    // free(o_a_cross_b_a);
    // free(d_cross_c_a);
    intersect = intersect * 0;
  }
  double t = dot_vec3(o_a_cross_b_a, c_a_vector) * inv_det;
  if (t > epsilon){
    *intersection_point = ray->o + (ray->d * t);
    // free(o_a_cross_b_a);
    // free(d_cross_c_a);
    return intersect;
  }
  // free(o_a_cross_b_a);
  // free(d_cross_c_a);
  return false;
}

// Update both or find a macro trick
#define FILE_LIST {"sphere.stl"}//,"sphere.stl"}
#define NUMBER_OF_FILES 1
struct Vec3 file_offsets[NUMBER_OF_FILES] = {Vec3(0,0,500)};//, Vec3(100,0,0)};
#define DEBUG_MODE true

#define H 500 // pixel height
#define W 500 // pixel width
#define BRIGHTNESS 0.5
#define SCALING 2.5
#define OFFSET 0.0
#define ZOOM 1


// generate a raytraced framed
// requires an array of stls, the number of stls, the output file name, the light angle (angle of the light source, this is ABSOLUTE)
// and the angle of the object (this is INCREMENTING, each frame generation with a given object angle MODIFIES THE STL)
void raytrace(struct STL *stl[], struct Triangle * tri_d, const int number_of_stls, Vec3 *output, float light_angle, float object_angle)
{

  // creating light source point
  double light_source_x = W/2+W*cos(light_angle)/2;
  double light_source_y = H/2+H*sin(light_angle)/2;
  double light_source_z = 1000.0;
  const Sphere light(Vec3(light_source_x,light_source_y,light_source_z ), 1);

  const struct Vec3 white(MAX_PIXEL, MAX_PIXEL, MAX_PIXEL); // the red will likely need to substituted with surface parameters
  const struct Vec3 black(0, 0, 0);
  const struct Vec3 red(MAX_PIXEL, 0, 0);

  struct Vec3 pix_col(black);
  struct Vec3 *pi;
  pi = (Vec3 *) malloc(sizeof(Vec3));
  int i = 0;
  int j = 0;
  int z = 0;
  static struct Ray * ray;
  #pragma omp threadprivate(ray) 
    for(i = 0; i < W; i++){
  #pragma omp parallel for shared(output,i) private(j,z)
      for( j = 0; j < H; j++){
        ray = new Ray(Vec3(i/ZOOM, j/ZOOM, 400), Vec3(0,0,1));
        for ( z = 0; z < number_of_stls; z++) {
          pix_col = black;
          // ray = Ray(Vec3(i/ZOOM,j/ZOOM,400),Vec3(0,0,1));
          for (int ind = 0; ind < stl[z]->length; ind++){
            if(ray_triangle_intersect(ray, &(tri_d[ind]), pi)){
                const Vec3 L = light.c - *pi;
                const Vec3 N = stl[z]->triangles[ind].normal;
                const double dt = abs(dot_vec3(L.normalize(), N.normalize()));
                pix_col = (red + white*dt) * BRIGHTNESS;
                clamp_pixels(pix_col);
            }
          }
          // debugging highlighting origin with red square
          if (i <= 10 && j <= 10 ){
            pix_col = Vec3(MAX_PIXEL,0,0);
          } else if (fabs(i - light_source_x) <= 1 && fabs(j-light_source_y ) <= 1){
            pix_col = white;
          }
          // paint y axis green
          if (j == 0){
            pix_col = Vec3(0,MAX_PIXEL,0);
          // paint x axis blue
          } else if (i == 0){
            pix_col = Vec3(0,0,MAX_PIXEL);
          }
          int pixel_index = i*W +j;
          output[pixel_index] = pix_col;
        }
      }
    }
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
void stl_raytracer_main(int frame_arr [], int frame_arr_length, int total_steps) 
{
  struct STL *stl[NUMBER_OF_FILES];
  const std::string filenames[NUMBER_OF_FILES] = FILE_LIST;

  double start_time, finish_time, total_time, current_time, increment_point;
  
  struct Parameters params = Parameters(SCALING, OFFSET, H, W);

  for (int i = 0; i < NUMBER_OF_FILES; i++) {
    stl[i] = load_stl(filenames[i], params, file_offsets[i]);
    std::cout << "Successfully loaded " <<  filenames[i] << "%s\n";
    printf("Number of triangles: %i\n", stl[i]->length);
  }

  std::string output_filename = "output/out.ppm";
  start_time = CLOCK();
  increment_point = start_time;


  uint32_t stl_size = NUMBER_OF_FILES * (sizeof(struct STL)); //stl[0]->length * sizeof(struct Triangle) + sizeof(uint32_t) + sizeof(struct Vec3));
  uint32_t output_size = H*W*sizeof(struct Vec3);
  Vec3 *output_values;
  output_values = (Vec3 *)malloc(output_size);
  struct STL **stl_d;
  uint32_t tri_size = sizeof(struct Triangle) * stl[0]->length;
  struct Triangle * tri_d = (struct Triangle *)malloc(tri_size);

  struct Vec3 *output_values_d = (Vec3 *)malloc(sizeof(Vec3)*H*W);
  
  int last_frame = 0;
  for (int i = 0; i < frame_arr_length; i++){
    std::string appended_info = std::to_string(i+1);
    // copy values to the gpu kernel
    

    float object_angle;
    if (i > 0){
      object_angle =  (frame_arr[i] - last_frame) * M_PI/(float)total_steps;
    } else {
      object_angle = frame_arr[i ] * M_PI/(float)total_steps; 
    }
    last_frame = frame_arr[i];
    rotate_stl(ROT_Z, stl[0], -object_angle/2);
    // rotate_stl(ROT_X, stl[0], -object_angle*2);
    // rotate_stl(ROT_Y, stl[0], object_angle);


    memcpy(tri_d, stl[0]->triangles, tri_size);

    // stl_d[0]->triangles = tri_d;





    raytrace(stl, tri_d, NUMBER_OF_FILES, output_values_d, frame_arr[i]*M_PI/(float)total_steps, M_PI/(float)total_steps);

    // copy values back out
    memcpy(output_values, output_values_d, output_size);

    // Save values locally
    output_filename.insert(10,appended_info);
    std::ofstream out(output_filename);
    out << "P3\n" << W << ' ' << H << ' ' << MAX_PIXEL<<"\n";

    for (int j = 0; j < H*W; j++) {
      out << (int) output_values[j].x << ' '
          << (int) output_values[j].y << ' '
          << (int) output_values[j].z << '\n';
    }
    out.close();
    
    output_filename = "output/out.ppm";
    if (DEBUG_MODE) {
      current_time = CLOCK();
      printf("Frame %d processed in %f ms\n", frame_arr[i], current_time-increment_point);
      increment_point = current_time;
    }
  }
  free(output_values_d);
  // free(stl);
  if (DEBUG_MODE) {
    finish_time = CLOCK();
    total_time = finish_time-start_time;
    printf("The total time to raytrace was: %f ms\n", total_time);
  }
}

#define STEPS 100
// #define MPI
#ifdef MPI
  #include "mpi.h"
#endif
int main(int argc, char *argv[]){
  #ifdef MPI
    int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    // mpi initialization 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);
    int block_size = STEPS/numprocs;
  #else
    int block_size = STEPS;
    int rank = 0;
  #endif
  int * frame_arr = (int *)malloc(sizeof(int)*block_size);
  
  for (int i = 0; i< block_size; i++){
    frame_arr[i] = block_size * rank + i;
  }

  stl_raytracer_main(frame_arr, block_size, STEPS);
  #ifdef MPI
    MPI_Finalize();
  #endif
}