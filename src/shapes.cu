#include "shapes.hpp"
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

const uint8_t index_mapping[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};

__device__ double dot_vec3(const Vec3& a, const Vec3& b)
{
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

__device__ struct Vec3 * cross_vec3(const Vec3& a, const Vec3& b)
{
  struct Vec3 * cross_product = (struct Vec3 *)malloc(sizeof(struct Vec3));
  cross_product->x = (a.y * b.z) - (b.y * a.z); 
  cross_product->y = (a.z * b.x) - (b.z * a.x);
  cross_product->z = (a.x * b.y) - (b.x * a.y);
  return cross_product;
}

__device__ struct Vec3 rotate_vec3(enum ROT_MATRIX_TYPE matrix, const Vec3& v, double theta_rad){
  Vec3 rotation_matrix [3] = {Vec3(0,0,0), Vec3(0,0,0), Vec3(0,0,0)};
  // printf("cos: %f, sin: %f\n", cos(theta_rad), sin(theta_rad));
  if (matrix == ROT_X){
    rotation_matrix[0] = Vec3(1, 0, 0);
    rotation_matrix[1] = Vec3(0,cos(theta_rad), -sin(theta_rad));
     rotation_matrix[2] = Vec3(0,sin(theta_rad),cos(theta_rad));
  } else if (matrix == ROT_Y){
    rotation_matrix[0] = Vec3(cos(theta_rad), 0, sin(theta_rad));
    rotation_matrix[1] = Vec3(0,1, 0);
    rotation_matrix[2] = Vec3(-sin(theta_rad),0,cos(theta_rad));
  } else if (matrix == ROT_Z){
    rotation_matrix[0] = Vec3(cos(theta_rad), -sin(theta_rad), 0);
    rotation_matrix[1] = Vec3(sin(theta_rad), cos(theta_rad), 0);
    rotation_matrix[2] = Vec3(0,0,1);
  }

  // rotation_matrix[0].print();
  // rotation_matrix[1].print();
  // rotation_matrix[2].print();

  return Vec3(dot_vec3(rotation_matrix[0], v),
              dot_vec3(rotation_matrix[1], v),
              dot_vec3(rotation_matrix[2], v));
}

__device__ void rotate_triangle(enum ROT_MATRIX_TYPE matrix, Triangle& tri, double theta_rad){
  tri.normal = rotate_vec3(matrix, tri.normal, theta_rad);
  tri.v0 = rotate_vec3(matrix, tri.v0, theta_rad);
  tri.v1 = rotate_vec3(matrix, tri.v1, theta_rad);
  tri.v2 = rotate_vec3(matrix, tri.v2, theta_rad);
}

__device__ void shift_triangle(Triangle &tri, struct Vec3 shift){
  tri.v0 = tri.v0 + shift;
  tri.v1 = tri.v1 + shift;
  tri.v2 = tri.v2 + shift;
}

__device__ void rotate_stl(enum ROT_MATRIX_TYPE matrix, struct STL * stl, double theta_rad){
  for (int i = 0; i < stl->length; i++){
    shift_triangle(stl->triangles[i], stl->center * -1);
    rotate_triangle(matrix, stl->triangles[i] , theta_rad);
    shift_triangle(stl->triangles[i], stl->center);
  }
}


double magnitude_vec3(const Vec3*vec){
  return (sqrt(vec->x* vec->x + vec->y * vec->y + vec->z * vec->z));
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

// move data from a buffer to a 3 element vector
void map_to_vector(uint32_t * buffer, Vec3 &vector)
{
  vector.x = *(float *) buffer;
  vector.y = *(float *) (buffer + 1);
  vector.z = *(float *) (buffer + 2);
}

struct Vec3 get_triangle_center(struct Triangle tri){
  Vec3 center = Vec3((tri.v0.x + tri.v1.x + tri.v2.x),
                    (tri.v0.y + tri.v1.y + tri.v2.y),
                    (tri.v0.z + tri.v1.z + tri.v2.z)); 
  return (center / 3);
}

struct Vec3 get_model_center(STL * stl){
  double total_area = 0.0;
  struct Vec3 total_weighted_area = Vec3(0,0,0);
  for (int i = 0; i < stl->length; i++){
    struct Vec3 edge_0 = stl->triangles[i].v1 - stl->triangles[i].v0;
    struct Vec3 edge_1 = stl->triangles[i].v2 - stl->triangles[i].v0;
    struct Vec3 *cross_prod = cross_vec3(edge_0, edge_1);
    double area = 0.5 * magnitude_vec3(cross_prod);
    total_weighted_area = total_weighted_area + (get_triangle_center(stl->triangles[i]) * area); 
    total_area+=area;
    free(cross_prod);
  }
  return total_weighted_area/total_area;
}


// stl structure:80 bytes of header, 4 bytes of length, then sets of 50 bytes (at least) 
// that describe the verticies and the normal vector
// there are also 2 bytes in each triangle that have the length of the attribute (usually set to 0),
// however, we may actually use these (might be in ply?) so this assumes that there may be additional bytes of 
// attributes

// BINARY stl loader (not meant for ascii)
struct STL* load_stl(const std::string& filename, struct Parameters params, struct Vec3 file_offsets)
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

    struct Vec3 center = get_model_center(stl_struct);
    center.print();
    Vec3 model_centering_dif =  Vec3(params.h/2, params.w/2, 0.0)- center;
    for (int i = 0; i < stl_struct->length; i++){
      stl_struct->triangles[i].v0 = stl_struct->triangles[i].v0 + model_centering_dif;
      stl_struct->triangles[i].v1 = stl_struct->triangles[i].v1 + model_centering_dif;
      stl_struct->triangles[i].v2 = stl_struct->triangles[i].v2 + model_centering_dif;
    }
    for (int i = 0; i < stl_struct->length; i++){
      stl_struct->triangles[i].v0 = stl_struct->triangles[i].v0 + file_offsets;
      stl_struct->triangles[i].v1 = stl_struct->triangles[i].v1 + file_offsets;
      stl_struct->triangles[i].v2 = stl_struct->triangles[i].v2 + file_offsets;
    }

    if (params.c_o){

      Vec3 offset = Vec3(params.c_o,params.c_o,params.c_o);
      for (int i = 0; i < stl_struct->length; i++){
        stl_struct->triangles[i].v0 = stl_struct->triangles[i].v0 + offset;
        stl_struct->triangles[i].v1 = stl_struct->triangles[i].v1 + offset;
        stl_struct->triangles[i].v2 = stl_struct->triangles[i].v2 + offset;
      }
    }
    center = get_model_center(stl_struct);
    stl_struct->center = center;

    if(params.s != 1.0){
      center.print();
      for (int i = 0; i < stl_struct->length; i++){
        stl_struct->triangles[i].v0 = center + ((stl_struct->triangles[i].v0 - center) * params.s) ;
        stl_struct->triangles[i].v1 = center + ((stl_struct->triangles[i].v1 - center) * params.s) ;
        stl_struct->triangles[i].v2 = center + ((stl_struct->triangles[i].v2 - center) * params.s) ;
      }
    }
    
    // if (params.h)


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
    output_file.close();
  }
}