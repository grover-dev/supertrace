#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define BOUNDS 11

// structure holding ray parameters
struct ray_t{
    float x;
    float y;
    float z;

    float dx;
    float dy;
    float dz;

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

float * vector_norm(float x, float y, float z){
    float mag = sqrt(pow(x, 2) + pow(y,2)+ pow(z, 2));    
    float * vec_n = malloc(3* sizeof(float));    
    vec_n[0] = x / mag;
    vec_n[1] = y / mag;
    vec_n[2] = z / mag;
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
    ray->dx = vec_n[0]; 
    ray->dy = vec_n[1]; 
    ray->dz = vec_n[2]; 
    free(vec_n);

    // starting off with white light
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

void trace_step(struct ray_t ** rays){
    struct ray_t * ray = rays[0];
    printf("ray 0,0: (%f, %f, %f), delta: (%f, %f, %f)\n", ray->x,ray->y,ray->z, ray->dx, ray->dy, ray->dz);
    for (int i = 0; i < BOUNDS * BOUNDS; i++){
        ray = rays[i];
        ray->x += ray->dx;
        ray->y += ray->dy;
        ray->z += ray->dz;
        if(i == 0){
            printf("ray 0,0: (%f, %f, %f)\n", ray->x,ray->y,ray->z);
        }
    }    
}



#define STEPS 10


float D = 10.0;
const float PI = 3.14159265;
float THETA = PI/2;

int main(){
//    cuda_hello<<<1,1>>>();
    struct ray_t ** rays = malloc(BOUNDS*BOUNDS*sizeof(struct ray_t *));
    // calculating bounding coords
    float start_x = - D * tan(THETA/2);
    float start_y = D * tan(THETA/2);
    float start_z = D;
    //float delta_theta = THETA / (float)BOUNDS;
    float delta = D * tan(THETA/2) * 2.0 / ((float)BOUNDS-1.0);
    //printf("theta: %f, D: %f, x size: %f\n", THETA, D, size_x); 
    
    float x = start_x;
    float y = start_y;
    float z = start_z;
    // creating BOUNDS x BOUNDS grid of rays 
    for(int j = 0; j < BOUNDS; j++){
        x = start_x;
        for (int i = 0; i < BOUNDS; i++){
            rays[i+j*BOUNDS] = new_ray(x, y, z); 
            //printf("(%.2f,%.2f)", rays[i+j*BOUNDS]->x, rays[i+j*BOUNDS]->y); 
            x += delta;
        }
        y -= delta;
        //printf("\n");

    }
    
    for (int i = 0; i < STEPS; i ++){
        trace_step(rays);        
    }
    

    // cleaning up
    for (int i = 0; i < BOUNDS*BOUNDS;i++){
        free(rays[i]);   
    }
    free(rays);
    return 0;
}
