#ifndef SHAPES_CUH
#define SHAPES_CUH

#include <cmath>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>



struct Vec3 {
  double x,y,z;
  __device__ __host__ Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
  __device__ __host__ Vec3(const Vec3& v): x(v.x), y(v.y), z(v.z) {}
  __device__ __host__ Vec3 operator + (const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
  __device__ __host__ Vec3 operator - (const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
  __device__ __host__ Vec3 operator * (float d) const { return Vec3(x*d, y*d, z*d); }
  __device__ __host__ Vec3 operator / (float d) const { return Vec3(x/d, y/d, z/d); }
  __device__ __host__ Vec3 max(const Vec3& v)
  { 
    return Vec3((x>v.x) ? x : v.x,(y>v.y) ? y : v.y,(z>v.z) ? z : v.z);
  }
  __device__ __host__ Vec3 normalize() const {
    double mg = sqrt(x*x + y*y + z*z);
    return Vec3(x/mg,y/mg,z/mg);
  }
  __device__ __host__ void print() const {
    printf("x: %f, y: %f, z: %f\n",x,y,z);
  }
};


__device__ __host__ double dot_vec3(const Vec3& a, const Vec3& b);

__device__ __host__ struct Vec3 cross_vec3(const Vec3& a, const Vec3& b);


struct Triangle {
  Vec3 v0, v1, v2;
  Vec3 normal;
  uint32_t * attributes;
  uint32_t attribute_length;
};

struct STL {
  struct Triangle * triangles;
  uint32_t length;
  Vec3 center;
};

enum ROT_MATRIX_TYPE {ROT_X, ROT_Y, ROT_Z};
 __host__ struct Vec3 rotate_vec3(enum ROT_MATRIX_TYPE matrix, const Vec3& v, double theta_rad);
 __host__ void rotate_triangle(enum ROT_MATRIX_TYPE matrix, Triangle& tri, double theta_rad);
 __host__ void rotate_stl(enum ROT_MATRIX_TYPE matrix, struct STL * stl, double theta_rad);

struct STL* load_stl(const std::string& filename, struct Parameters params, struct Vec3 file_offsets);

void dump_stl(struct STL* polygon, std::string filename);

struct Parameters {
  double s; // scaling
  double c_o; // constant offset
  unsigned int h; // height
  unsigned int w; // width
  Parameters(const double s, const double c_o, const unsigned int h, const unsigned int w) :
             s(s), c_o(c_o), h(h), w(w) {}
};


struct Ray {
  Vec3 o,d;
  __device__ __host__ Ray(const Vec3& o, const Vec3& d) : o(o), d(d) {}
};


struct Sphere {
  Vec3 c;
  double r;
  __device__ __host__ Sphere(const Vec3& c, float r) : c(c), r(r) {}
  __device__ __host__ Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
  __device__ __host__ bool intersect(const Ray& ray, double &t) const {
    const Vec3 o = ray.o;
    const Vec3 d = ray.d;
    const Vec3 oc = o - c;
    const double b = 2 * dot_vec3(oc, d);
    const double c = dot_vec3(oc, oc) - r*r;
    double disc = b*b - 4 * c;
    if (disc < 1e-4) return false;
    disc = sqrt(disc);
    const double t0 = -b - disc;
    const double t1 = -b + disc;
    t = (t0 < t1) ? t0 : t1;
    return true;
  }
};
#endif
