#include "shapes.hpp"

#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

uint8_t bit_lookup[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};

double dot(const Vec3& a, const Vec3& b)
{
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

inline bool file_exists(const std::string& filename)
{
  return(access(filename.c_str(), F_OK) != -1);
}

void map_to_vector(uint32_t * buffer, Vec3 &vector)
{
  vector.x = *(float *) buffer;
  vector.y = *(float *) (buffer + 1);
  vector.z = *(float *) (buffer + 2);
}

// stl structure:80 bytes of header, 4 bytes of length, then sets of 50 bytes (at least) 
// that describe the verticies and the normal vector
// there are also 2 bytes in each triangle that have the length of the attribute (usually set to 0),
// however, we may actually use these (might be in ply?) so this assumes that there may be additional bytes of 
// attributes

// BINARY stl loader (not meant for ascii)
struct STL* load_stl(const std::string& filename)
{
  if (file_exists(filename)){
    //printf("Found %s\n", filename.c_str());
    std::ifstream ifs(filename);
    char c;

    uint32_t index = 0;
    uint32_t length = 0;
    uint32_t temp_4byte_reader = 0;
    struct STL * stl_struct;
    uint32_t tricount = 0;

    bool at_start = true;

    uint32_t buffer[3] = {0};
    
    struct Triangle *tri;
    uint32_t triangle_index = 0;

    uint32_t start_index = 0;
    if(ifs.is_open()){
      while(ifs.good()) {
        ifs.get(c);
        if (index >= 80 && index < 84){
          temp_4byte_reader = 0x000000ff & c;
          length |= (temp_4byte_reader << (8 * (index - 80)));

          if (index == 83) {
            stl_struct = (struct STL *) malloc((sizeof(struct STL)));
            stl_struct->triangles = (struct Triangle *) malloc(sizeof(struct Triangle) * length);
            stl_struct->length = length;
          }
        } else if (index >= 84) {
          if (at_start){
              triangle_index = 0;
              start_index = index;
              at_start = false;
          }
          
          // read value in and store it in appropriate byte 
          // stls are little endian, so bit shift here converts to big endian
          // dividing index by 3 also puts the 12 bytes into 3 uint32s
          temp_4byte_reader = 0x000000ff & c;
          //printf("Inedex:    %u, Bytestring %x:      ",index, temp_4byte_reader);
          buffer[bit_lookup[(index-start_index) % 12]] |= temp_4byte_reader << (8 * (((index-start_index) % 12) % 4));
          if ((index-start_index) % 12 == 11 && triangle_index < 4) {
            if (triangle_index == 0) {
              map_to_vector(buffer, stl_struct->triangles[tricount].normal);
            } else if (triangle_index == 1) {
              map_to_vector(buffer, stl_struct->triangles[tricount].v0);
            } else if (triangle_index == 2) {
              map_to_vector(buffer, stl_struct->triangles[tricount].v1);
            } else if (triangle_index == 3) {
              map_to_vector(buffer, stl_struct->triangles[tricount].v2);
            }

            buffer[0] = 0;
            buffer[1] = 0;
            buffer[2] = 0;
            temp_4byte_reader = 0;
            triangle_index++;
          } else if (triangle_index == 4 || index-start_index > 48) {// the first condition should be sufficient
            //todo: handle the length bytes, then check how long attributes are

            // Once all data for a single triangle is processed, move on to the next
            if (index-start_index == 49) {
              // printf("Normal: %f, %f, %f\n", stl_struct->triangles[tricount].normal.x, stl_struct->triangles[tricount].normal.y, stl_struct->triangles[tricount].normal.z);
              // printf("Vertex1: %f, %f, %f\n", stl_struct->triangles[tricount].v0.x, stl_struct->triangles[tricount].v0.y, stl_struct->triangles[tricount].v0.z);
              // printf("Vertex2: %f, %f, %f\n", stl_struct->triangles[tricount].v1.x, stl_struct->triangles[tricount].v1.y, stl_struct->triangles[tricount].v1.z);
              // printf("Vertex3: %f, %f, %f\n", stl_struct->triangles[tricount].v2.x, stl_struct->triangles[tricount].v2.y, stl_struct->triangles[tricount].v2.z);
              tricount++;
              //printf("Processed triangle: %d\n", tricount);
              at_start = true;
            }
          }
        }
        index++;
      }
    }
    ifs.close();

    return stl_struct;
  }

  return nullptr; // check for this "exception"
}