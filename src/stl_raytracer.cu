#include <fstream>
#include <stdlib.h>
#include "shapes.cuh"
#include <string>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdint.h>


double CLOCK() {
        struct timespec t;
        clock_gettime(CLOCK_MONOTONIC,  &t);
        return (t.tv_sec * 1000)+(t.tv_nsec*1e-6);
}

#define MAX_PIXEL 1023
__device__ void clamp_pixels(Vec3& col)
{
  col.x = (col.x > MAX_PIXEL) ? MAX_PIXEL : (col.x < 0) ? 0 : col.x;
  col.y = (col.y > MAX_PIXEL) ? MAX_PIXEL : (col.y < 0) ? 0 : col.y;
  col.z = (col.z > MAX_PIXEL) ? MAX_PIXEL : (col.z < 0) ? 0 : col.z;
}



__device__ bool ray_triangle_intersect(struct Ray * ray, struct Triangle * tri, struct Vec3 * intersection_point){
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
    free(d_cross_c_a);
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

// Update both or find a macro trick
#define FILE_LIST {"pyramid.stl"}//,"sphere.stl"}
#define NUMBER_OF_FILES 1
struct Vec3 file_offsets[NUMBER_OF_FILES] = {Vec3(0,0,100)};//, Vec3(100,0,0)};
#define DEBUG_MODE true

#define H 500 // pixel height
#define W 500 // pixel width
#define BRIGHTNESS 0.5
#define SCALING 5.0
#define OFFSET 0.0
#define ZOOM 1

#define STEPS 10

// generate a raytraced framed
// requires an array of stls, the number of stls, the output file name, the light angle (angle of the light source, this is ABSOLUTE)
// and the angle of the object (this is INCREMENTING, each frame generation with a given object angle MODIFIES THE STL)
__global__ void raytrace(struct STL *stl[], const int number_of_stls, Vec3 *output, float light_angle, float object_angle)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  if((i >= W) || (j >= H)) return;
  int pixel_index = j*W + i;

  // creating light source point
  double light_source_x = W/2+W*cos(light_angle)/2;
  double light_source_y = H/2+H*sin(light_angle)/2;
  double light_source_z = 500.0;
  const Sphere light(Vec3(light_source_x,light_source_y,light_source_z ), 1);

  const struct Vec3 white(MAX_PIXEL, MAX_PIXEL, MAX_PIXEL); // the red will likely need to substituted with surface parameters
  const struct Vec3 black(0, 0, 0);
  const struct Vec3 red(MAX_PIXEL, 0, 0);

  struct Vec3 pix_col(black);
  struct Vec3 *pi;
  cudaMalloc(&pi, sizeof(Vec3));

  
  struct Vec3 pix_col_tmp = black;

  rotate_stl(ROT_Z, stl[0], object_angle);
  rotate_stl(ROT_X, stl[0], object_angle);
  rotate_stl(ROT_Y, stl[0], -object_angle);
  
  
  for (int z = 0; z < number_of_stls; z++) {
    pix_col = black;
    Ray ray(Vec3(i/ZOOM,j/ZOOM,0), Vec3(0,0,1));
    for (int ind = 0; ind < stl[z]->length; ind++){
      if(ray_triangle_intersect(&ray, &(stl[z]->triangles[ind]), pi)){
          const Vec3 L = light.c - *pi;
          const Vec3 N = stl[z]->triangles[ind].normal;
          const double dt = dot_vec3(L.normalize(), N.normalize());
          pix_col = (red + white*dt) * BRIGHTNESS;
          clamp_pixels(pix_col); 
          if(ind > 0){
            pix_col = pix_col_tmp.max(pix_col);
          }
          pix_col_tmp = pix_col;
      }
    }
    pix_col_tmp = black;
    output[pixel_index] = pix_col;
    cudaFree(pi);
  }
  //   // debugging highlighting origin with red square
  //   if (i <= 10 && j <= 10 ){
  //     pix_col = Vec3(MAX_PIXEL,0,0);
  //   } else if (fabs(i - light_source_x) <= 1 && fabs(j-light_source_y ) <= 1){
  //     pix_col = white;
  //   }
  //   // paint y axis green
  //   if (j == 0){
  //     pix_col = Vec3(0,MAX_PIXEL,0);
  //   // paint x axis blue
  //   } else if (i == 0){
  //     pix_col = Vec3(0,0,MAX_PIXEL);
  //   }
  // }
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
  const std::string filenames[NUMBER_OF_FILES] = FILE_LIST;

  double start_time, finish_time, total_time, current_time, increment_point;
  
  struct Parameters params = Parameters(SCALING, OFFSET, H, W);

  for (int i = 0; i < NUMBER_OF_FILES; i++) {
    stl[i] = load_stl(filenames[i], params, file_offsets[i]);
    std::cout << "Successfully loaded " <<  filenames[i] << "%s\n";
    printf("Number of triangles: %i\n", stl[i]->length);
  }

  std::string output_filename = "output/out.ppm";
  const int start = 0;
  start_time = CLOCK();
  increment_point = start_time;
  for (int i = start; i < STEPS+start; i++){
    std::string appended_info = std::to_string(i+1);

    // copy values to the gpu kernel
    uint32_t stl_size = NUMBER_OF_FILES * sizeof(struct STL);
    uint32_t output_size = H*W*sizeof(struct Vec3);
    Vec3 *output_values;
    cudaMallocHost(&output_values, output_size);

    struct STL **stl_d;
    cudaMalloc(&stl_d, stl_size);
    cudaMemcpy(stl_d, stl, stl_size, cudaMemcpyHostToDevice);
    struct Vec3 *output_values_d;
    cudaMalloc(&output_values_d, output_size);
    cudaMemcpy(output_values_d, output_values, output_size, cudaMemcpyHostToDevice);


    raytrace<<<H, W>>>(stl_d, NUMBER_OF_FILES, output_values_d, i*2*M_PI/(float)STEPS, M_PI/(float)STEPS);
    //raytrace(stl, NUMBER_OF_FILES, i*2*M_PI/(float)STEPS,M_PI/(float)STEPS);

    // copy values back out
    cudaMemcpy(output_values_d, output_values, output_size, cudaMemcpyDeviceToHost);
    cudaFree(output_values_d);
    cudaFree(stl_d);

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
