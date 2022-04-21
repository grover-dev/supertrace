#include <fstream>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdint.h>
#include "shapes.hpp"
#include <cstring>
#include <unistd.h>

 #define MPI
 #ifdef MPI
   #include "mpi.h"
 #endif

 double CLOCK() {
        struct timespec t;
        clock_gettime(CLOCK_MONOTONIC,  &t);
        return (t.tv_sec * 1000)+(t.tv_nsec*1e-6);
}


#define FILE_LIST {"chair.stl"}
struct Vec3 file_offsets[1] = {Vec3(0,0,500)};
#define DEBUG_MODE true

#define H 500 // pixel height
#define W 500 // pixel width
#define BRIGHTNESS 0.5
#define SCALING 2.5
#define OFFSET 0.0
#define ZOOM 1

#define NUM_THREADS 1
static pthread_t threads[NUM_THREADS];
static pthread_mutex_t total_mutex;

float light_angle = 0.0; 

struct Vec3 *output_values;
const uint32_t output_size = H*W*sizeof(struct Vec3);
struct STL *stl[1];
struct Triangle * tris;

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

  struct Vec3 d_cross_c_a = cross_vec3(ray->d, c_a_vector); 

  // first calculating determinant, 
  //  if its ~0 then the ray is parallel to the triangle
  //  if it is <0, then we are hitting the back of the triangle (counting as not intersecting for now)
  //  this will need to be adjusted in the future (especially with refraction) (i.e. use absolute value of det)
  // can therefore ignore it  
  double det = dot_vec3(d_cross_c_a, b_a_vector);
  if (det < epsilon && det > -epsilon){
    return false;
  }
  double inv_det = 1.0 / det; // f = 1/a

  struct Vec3 o_a_vector = ray->o - tri->v0; // s = ray origin - vertex0

  // start calculating barycentric coord vectors
  double u = dot_vec3(o_a_vector, d_cross_c_a) * inv_det;
  // since the vectors are normalized, anything < 0 or > 1 means that the intersection
  // is not in the bounds of the triangle 
  if (u < 0.0 || u > 1.0){
    return false;
  }
  struct Vec3 o_a_cross_b_a = cross_vec3(o_a_vector, b_a_vector); 
  double v = dot_vec3(o_a_cross_b_a, ray->d) * inv_det;
  if (v < 0.0 || (v+u) > 1.0){
    return false;
  }
  double t = dot_vec3(o_a_cross_b_a, c_a_vector) * inv_det;
  if (t > epsilon){
    *intersection_point = ray->o + (ray->d * t);
    return true;
  }
  return false;
}

void print_frame(int *pixels, std::string filename)
{
  std::ofstream out(filename);
    out << "P3\n" << W << ' ' << H << ' ' << MAX_PIXEL<<"\n";

    for (int j = 0; j < H*W; j++) {
      out << (int) pixels[j*3] << ' '
          << (int) pixels[j*3+1] << ' '
          << (int) pixels[j*3+2] << '\n';
    }
    out.close();
}

// generate a raytraced framed
// requires an array of stls, the number of stls, the output file name, the light angle (angle of the light source, this is ABSOLUTE)
// and the angle of the object (this is INCREMENTING, each frame generation with a given object angle MODIFIES THE STL)
void * raytrace(void)
{
  // intptr_t tmp = (intptr_t)args;
  // // int i = (int)tmp;

  // printf("index: %f\n", light_angle);
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
  int j = 0;
  int z = 0;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int *output_int = (int *) malloc(W*H*3*sizeof(int));
  for (int i = 0; i < W; i++){
    for( j = 0; j < H; j++){
      for ( z = 0; z < 1; z++) {
        pix_col = black;
        Ray ray = Ray(Vec3(i/ZOOM,j/ZOOM,400),Vec3(0,0,1));
        for (int ind = 0; ind < stl[z]->length; ind++){
          if(ray_triangle_intersect(&ray, &(tris[ind]), pi)){
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
        int pixel_index = i*W+j;
        output_int[pixel_index*3] = pix_col.x;
        output_int[pixel_index*3 + 1] = pix_col.y;
        output_int[pixel_index*3 + 2] = pix_col.z;
      }
    }
  }
  MPI_Send(&output_int, 1,  MPI_INT,0,0,MPI_COMM_WORLD);
  printf("frame finished\n");
  MPI_Send(&output_int, 3*W*H, MPI_INT, 0, 0, MPI_COMM_WORLD);
  // free(output_int);
  free(pi);
}


// using Möller-Trumbore algorithm for raytracing w/ triangles 
void stl_raytracer_main(int frame_arr [], int frame_arr_length, int total_steps) 
{
  const std::string filenames[1] = FILE_LIST;

  double start_time, finish_time, total_time, current_time, increment_point;
  
  struct Parameters params = Parameters(SCALING, OFFSET, H, W);

  for (int i = 0; i < 1; i++) {
    stl[i] = load_stl(filenames[i], params, file_offsets[i]);
    std::cout << "Successfully loaded " <<  filenames[i] << "%s\n";
    printf("Number of triangles: %i\n", stl[i]->length);
  }

  start_time = CLOCK();
  increment_point = start_time;


  uint32_t stl_size = (sizeof(struct STL)); 
  output_values = (Vec3 *)malloc(output_size);
  uint32_t tri_size = sizeof(struct Triangle) * stl[0]->length;

  tris = (struct Triangle *)malloc(tri_size);
  memcpy(tris, stl[0]->triangles, tri_size);
  
  int last_frame = 0;
  for (int i = 0; i < frame_arr_length; i++){
    //std::string appended_info = std::to_string(i+1);
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


    light_angle = frame_arr[i]*M_PI/(float)total_steps;

    for (int ind = 0; ind < 1; ind++) {
      raytrace();
      //if (pthread_create(&threads[i], NULL, raytrace, (void *)&p)){
      //    printf("failed to create pthread %i\n", i);
      //    exit(-1);
      //}
    }
    // thread joining 
    //for(int i = 0; i < NUM_THREADS; i++) {
    //    pthread_join(threads[i], NULL);
   // }





    // Save values locally
    //output_filename.insert(10,appended_info);
    //std::ofstream out(output_filename);
    //out << "P3\n" << W << ' ' << H << ' ' << MAX_PIXEL<<"\n";

    //for (int j = 0; j < H*W; j++) {
    //  out << (int) output_values[j].x << ' '
    //      << (int) output_values[j].y << ' '
    //      << (int) output_values[j].z << '\n';
    //}
    //out.close();
    
    //output_filename = "output/out.ppm";
    if (DEBUG_MODE) {
      current_time = CLOCK();
      printf("Frame %d processed in %f ms\n", frame_arr[i], current_time-increment_point);
      increment_point = current_time;
    }
  }
  // free(stl);
  if (DEBUG_MODE) {
    finish_time = CLOCK();
    total_time = finish_time-start_time;
    printf("The total time to raytrace was: %f ms\n", total_time);
  }
}


#define STEPS 10
int main(int argc, char *argv[]){
  #ifdef MPI
    int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    // mpi initialization 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);
    int block_size = STEPS/(numprocs-1);
  #else
    int block_size = STEPS;
    int rank = 0;
  #endif
  
  int * frame_arr = (int *)malloc(sizeof(int)*block_size);
  if (rank > 0) {
    for (int i = 0; i< block_size; i++){
      frame_arr[i] = block_size * (rank-1) + i;
    }
    stl_raytracer_main(frame_arr, block_size, STEPS);
  } else {
    std::string filename = "output/out.ppm";
    int *output_int = (int *) malloc(3*W*H*sizeof(int));
    for (int i =0; i < STEPS; i++) {
      // The expression for rank is wrong
      MPI_Recv(&output_int, 1, MPI_INT, 1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%i, numprocs: %i\n", i, numprocs);

      MPI_Recv(&output_int, 3*W*H, MPI_INT, i % (numprocs-1) + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("Receivee from %d from node %d\n", i, i % (numprocs-1) + 1);
      std::string appended_info = std::to_string(i+1);
      filename.insert(10,appended_info);
      print_frame(output_int, filename);
      filename = "output/out.ppm";
    }
  }

  #ifdef MPI
    MPI_Finalize();
  #endif
}