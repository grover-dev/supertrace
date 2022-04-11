#include "shapes.hpp"

#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>


double dot(const Vec3& a, const Vec3& b) {
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

inline bool file_exists(const std::string& filename){
  return(access(filename.c_str(), F_OK) != -1);
}

void map_to_vector(uint32_t * buffer, Vec3 * vector){

}

// stl structure:80 bytes of header, 4 bytes of length, then sets of 50 bytes (at least) 
// that describe the verticies and the normal vector
// there are also 2 bytes in each triangle that have the length of the attribute (usually set to 0),
// however, we may actually use these (might be in ply?) so this assumes that there may be additional bytes of 
// attributes

// BINARY stl loader (not meant for ascii)
void load_stl(const std::string& filename){
  if (file_exists(filename)){
    std::ifstream ifs(filename);
    char c;

    uint32_t index = 0;
    
    uint32_t length = 0;
    struct STL * stl_struct;

    bool at_start = true;

    uint32_t * buffer = (uint32_t *)malloc(3*sizeof(uint32_t));
    for (int i = 0; i < 3; i++){
      buffer[i] = 0;
    }
    struct Triangle * tri;
    uint32_t triangle_index = 0;

    uint32_t start_index = 0;
    if(ifs.is_open()){
      while(ifs.good()) {
        ifs.get(c);
        if (index == 84){
          stl_struct = (struct STL *) malloc(sizeof(struct Triangle)*length+sizeof(uint32_t));
        }
        if (index >= 80 && index < 84){
          length |= c << (index - 80);
        } else if (index >= 84){
          if (at_start){
              start_index = index;
              tri = (struct Triangle *)malloc(sizeof(struct Triangle));
              at_start = false;
          } else {
            // read value in and store it in appropriate byte 
            // stls are little endian, so bit shift here converts to big endian
            // dividing index by 3 also puts the 12 bytes into 3 uint32s
            buffer[((index-start_index)%12)/3] |= c << ((index-start_index)%4);
            if ((index-start_index) % 12 == 11 && triangle_index < 3){
              struct Vec3 * output_vector = (struct Vec3 *)malloc(sizeof(struct Vec3));
              map_to_vector(buffer, output_vector);
              if (triangle_index == 0){
                tri->v0 = *output_vector;
              } else if (triangle_index == 1){
                tri->v1 = *output_vector;
              } else if (triangle_index == 2){
                tri->v2 = *output_vector;
              }
              triangle_index++;
            } else {
              //todo: handle the length bytes, then check how long attributes are 

            }
          }
        }


        index++;
      }
    }

    ifs.close();
  }

}