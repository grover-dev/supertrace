#include "shapes.hpp"
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

const uint8_t index_mapping[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};

double dot_vec3(const Vec3& a, const Vec3& b)
{
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

struct Vec3 * cross_vec3(const Vec3& a, const Vec3& b)
{
  struct Vec3 * cross_product = (struct Vec3 *)malloc(sizeof(struct Vec3));
  //todo: check this
  cross_product->x = (a.y * b.z) - (b.y * a.z); 
  cross_product->y = (a.z * b.x) - (b.z * a.x);
  cross_product->z = (a.x * b.y) - (b.x * a.y);
  return cross_product;
} 


inline bool file_exists(const std::string& filename)
{
  return(access(filename.c_str(), F_OK) != -1);
}

inline void big_to_little_endian_unint_32(uint32_t big_endian_value, uint8_t* little_endian_buffer)
{
  little_endian_buffer[0] = (big_endian_value) & 0xff;
  little_endian_buffer[1] = (big_endian_value >> 8) & 0xff;
  little_endian_buffer[2] = (big_endian_value >> 16) & 0xff;
  little_endian_buffer[3] = (big_endian_value >> 24) & 0xff;
}

inline void big_to_little_endian_float(float big_endian_value, uint8_t* little_endian_buffer)
{
  uint32_t* byte_reader;
  byte_reader = (uint32_t*) &big_endian_value;
  big_to_little_endian_unint_32(*byte_reader, little_endian_buffer);
}

inline void write_vec3_to_file(struct Vec3* vec3, std::ofstream &output_file)
{
  uint8_t vec3_byte_buffer[4] = {0};

  big_to_little_endian_float(vec3->x, vec3_byte_buffer);
  output_file.write((char *) vec3_byte_buffer, 4);
  big_to_little_endian_float(vec3->y, vec3_byte_buffer);
  output_file.write((char *) vec3_byte_buffer, 4);
  big_to_little_endian_float(vec3->z, vec3_byte_buffer);
  output_file.write((char *) vec3_byte_buffer, 4);
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
    std::ifstream ifs(filename);
    char c;

    float max_value = 0.0;
    float min_value = 99999.99;

    uint32_t index = 0;
    uint32_t length = 0;
    uint32_t temp_4byte_reader = 0;
    uint32_t buffer[3] = {0};
    struct STL * stl_struct;
    uint32_t tricount = 0;

    bool at_start = true;
    uint32_t triangle_index = 0;
    uint32_t start_index = 0;

    Vec3 zeroing_offset = Vec3(0,0,0);


    if(ifs.is_open()){
      while(ifs.good()) {
        ifs.get(c);
        if (index >= 80 && index < 84){
          temp_4byte_reader = 0xff & c;
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
          temp_4byte_reader = 0xff & c;
          buffer[index_mapping[(index-start_index) % 12]] |= temp_4byte_reader << (8 * (((index-start_index) % 12) % 4));
          if ((index-start_index) % 12 == 11 && triangle_index < 4) {
            
            if (*(float *)buffer < zeroing_offset.x){
              zeroing_offset.x = *(float *)buffer;
            }
            if (*(float *)(buffer +1) < zeroing_offset.y){
              zeroing_offset.y = *(float *)(buffer +1);
            }
            if (*(float *)(buffer + 2) < zeroing_offset.z){
              zeroing_offset.z = *(float *)(buffer +2);
            }

            if (*(float *)buffer > max_value){
              max_value = *(float *)buffer;
            }
            if (*(float *)(buffer +1) > max_value){
              max_value = *(float *)(buffer +1);
            }
            if (*(float *)(buffer + 2) > max_value){
              max_value = *(float *)(buffer +2);
            }

            if (*(float *)buffer < min_value){
              min_value = *(float *)buffer;
            }
            if (*(float *)(buffer +1) < min_value){
              min_value = *(float *)(buffer +1);
            }
            if (*(float *)(buffer + 2) < min_value){
              min_value = *(float *)(buffer +2);
            }

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
            triangle_index++;
          } else if (triangle_index == 4) {
            //todo: handle the length bytes, then check how long attributes are

            // Once all data for a single triangle is processed, move on to the next
            if (index-start_index == 49) {
              tricount++;
              at_start = true;
            }
          }
        }
        index++;
      }
    }
    ifs.close();

    if (zeroing_offset.x != 0 || zeroing_offset.y != 0 || zeroing_offset.z != 0){
      Vec3 offset = Vec3(250,250,250);
      for (int i = 0; i < stl_struct->length; i++){
        // stl_struct->triangles[i].normal = stl_struct->triangles[i].normal - zeroing_offset;
        // stl_struct->triangles[i].v0 = stl_struct->triangles[i].v0 - zeroing_offset;
        // stl_struct->triangles[i].v1 = stl_struct->triangles[i].v1 - zeroing_offset;
        // stl_struct->triangles[i].v2 = stl_struct->triangles[i].v2 - zeroing_offset;
        stl_struct->triangles[i].normal = stl_struct->triangles[i].normal + offset;
        stl_struct->triangles[i].v0 = stl_struct->triangles[i].v0 + offset;
        stl_struct->triangles[i].v1 = stl_struct->triangles[i].v1 + offset;
        stl_struct->triangles[i].v2 = stl_struct->triangles[i].v2 + offset;
      }
    }

    printf("max value: %f and min value: %f\n", max_value, min_value);
    return stl_struct;
  }

  return nullptr; // check for this "exception"
}

// BINARY stl dumper (not meant for ascii)
// Create a binary stl file from a C struct
void dump_stl(struct STL* polygon, std::string filename)
{
  std::ofstream output_file (filename);
  if (output_file.is_open()){
    char c;
    int i, j;
    uint8_t byte_buffer[4] = {0};


    // Write empty header: 80 Bytes
    c = '\0';
    for (i = 0; i < 80; i++) {
      output_file.write(&c, 1);
    }

    // Write length
    big_to_little_endian_unint_32(polygon->length, byte_buffer);
    output_file.write((char *) byte_buffer, 4);
    

    // Write each triangle
    for (i = 0; i < polygon->length; i++) {
      // Write normal
      write_vec3_to_file(&polygon->triangles[i].normal, output_file);
      
      // Write each vertex
      write_vec3_to_file(&polygon->triangles[i].v0, output_file);
      write_vec3_to_file(&polygon->triangles[i].v1, output_file);
      write_vec3_to_file(&polygon->triangles[i].v2, output_file);

      // Write attributes as empty
      // todo: add attribute support
      for (j = 0; j < 2; j++) {
        output_file.write(&c, 1);
      }
    }
  }
}