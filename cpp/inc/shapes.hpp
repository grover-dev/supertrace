#ifndef SHAPES_HPP
#define SHAPES_HPP

#include <cmath>
#include <stdint.h>
#include <string>


struct Vec3 {
  float x,y,z;
  Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
  Vec3 operator + (const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
  Vec3 operator - (const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
  Vec3 operator * (float d) const { return Vec3(x*d, y*d, z*d); }
  Vec3 operator / (float d) const { return Vec3(x/d, y/d, z/d); }
  Vec3 normalize() const {
    double mg = sqrt(x*x + y*y + z*z);
    return Vec3(x/mg,y/mg,z/mg);
  }
};


double dot(const Vec3& a, const Vec3& b);

struct STL* load_stl(const std::string& filename);

struct Ray {
  Vec3 o,d;
  Ray(const Vec3& o, const Vec3& d) : o(o), d(d) {}
};


struct Sphere {
  Vec3 c;
  double r;
  Sphere(const Vec3& c, float r) : c(c), r(r) {}
  Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
  bool intersect(const Ray& ray, double &t) const {
    const Vec3 o = ray.o;
    const Vec3 d = ray.d;
    const Vec3 oc = o - c;
    const double b = 2 * dot(oc, d);
    const double c = dot(oc, oc) - r*r;
    double disc = b*b - 4 * c;
    if (disc < 1e-4) return false;
    disc = sqrt(disc);
    const double t0 = -b - disc;
    const double t1 = -b + disc;
    t = (t0 < t1) ? t0 : t1;
    return true;
  }
};

struct Triangle {
  Vec3 v0, v1, v2;
  Vec3 normal;
  uint32_t * attributes;
  uint32_t attribute_length;
};

struct STL {
  struct Triangle * triangles;
  uint32_t length;
};






#endif