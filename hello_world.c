#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define STEP 1.0

#define BOUNDS 256



// structure holding ray parameters
struct ray_t{
    float x;
    float y;
    float z;

    float vec_x;
    float vec_y;
    float vec_z;

    float red;
    float green;
    float blue;
};

struct image_t{
    unsigned char * red;
    unsigned char * green;
    unsigned char * blue;
};


//__global__ void cuda_hello(){
//    printf("hello world from gpu!\n");
//}

float * vector_norm(float vec_x, float vec_y, float vec_z){
    float mag = sqrt(pow(vec_x, 2) + pow(vec_y,2)+ pow(vec_z, 2));    
    float * vec_n = malloc(3* sizeof(float));    
    vec_n[0] = vec_x / mag;
    vec_n[1] = vec_y / mag;
    vec_n[2] = vec_z / mag;
    return vec_n;
}


// initialize new ray
struct ray_t * new_ray(float x, float y, float z){
    struct ray_t * ray = malloc(sizeof(struct ray_t)); 
    
    //x,y,z coordinates initialzied with non normalized vector 
    ray->x = x;
    ray->y = y;
    ray->z = z;
    float * vec_n;
    vec_n = vector_norm(x, y, z); 
    ray->vec_x = vec_n[0]; 
    ray->vec_y = vec_n[1]; 
    ray->vec_z = vec_n[2]; 
    free(vec_n);

    ray->red = 255;
    ray->green = 255;
    ray->blue = 255;
    
    return ray;
};

// writing ppm image file for debugging
// based on: http://rosettacode.org/wiki/Bitmap/Write_a_PPM_file#C
void image_write(struct image_t * image){
    FILE *fp =  fopen("first.ppm", "wb");
    (void) fprintf(fp, "P6\n%d %d\n255\n", BOUNDS, BOUNDS); 
    for (int x = 0; x < BOUNDS; x++){
        for (int y = 0; y < BOUNDS; y++){
            static unsigned char color[3];
            color[0] = image->red[x+y*BOUNDS];
            color[1] = image->green[x+y*BOUNDS];
            color[3] = image->blue[x+y*BOUNDS];
            (void) fwrite(color, 1, 3, fp);
        }
    } 
    (void) fclose(fp);
}



#define D 10.0
#define PI 3.14159265
#define THETA (PI/2)

int main(){
//    cuda_hello<<<1,1>>>();
    struct ray_t ** rays = malloc(BOUNDS*BOUNDS*sizeof(struct ray_t *));
    // creating BOUNDS x BOUNDS grid of rays
    for(int x = 0; x < BOUNDS; x++){
        for (int y = 0; y < BOUNDS; y++){
                        
        }
    }
     

    return 0;
}
