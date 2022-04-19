#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include "shapes.hpp"
#include <string.h>
#include <stdio.h>


// Both in pixel count
#define IMAGE_HEIGHT 500
#define IMAGE_WIDTH  500
#define ZOOM 1.0

// Constant for force integration
// #define h 1
#define SOLAR_MASS (1.98 * pow(10,30))
// #define c 299792458
#define G (6.6743*pow(10,-11))  

// system scaling gives the meters/pixel
float SYSTEM_SCALING_FACTOR = 1.0; 

struct body {
    int mass; // mass 
    double radius; // radius
    struct Vec3 center;
    body(int mass, double radius, struct Vec3& center) : mass(mass), radius(radius), center(center) {}
};

struct background {
    unsigned int height;
    unsigned int width;
    struct Vec3 * pix_arr;
    unsigned int pix_arr_length;
    background(const unsigned int height, const unsigned int width,
                struct Vec3 * pix_arr, unsigned int pix_arr_length) :
        height(height), width(width), pix_arr(pix_arr), pix_arr_length(pix_arr_length){}
};

float schwarzchild_radius(body b){
    return (2 * G * b.mass *SOLAR_MASS / pow(c,2));
}

float sqr_norm(struct Vec3& v){
    return dot_vec3(v,v);
} 


// Search script for "RK4" to see what this is supposed to do
struct Vec3 specific_force_field(struct Vec3& r,float h2)
{
    // r = 
    return (r * -1.5f * h2 / pow(sqr_norm(r),2.5));
}


// Approximate the solution to first order ODE
// x0 and y0 are initial condition
// x is function input
// step_size is the step size
void runge_kutta_4_vec3(Vec3& r0, Vec3& v0, body b, float step_size, int steps, float h2)
{
    // Count number of iterations using step size or
    // step height h
    struct Vec3 k1 = Vec3(0,0,0);
    struct Vec3 k2 = Vec3(0,0,0);
    Vec3 k3 = Vec3(0,0,0);
    Vec3 k4 = Vec3(0,0,0);
    // Iterate for number of iterations
    for (int i=1; i<=steps; i++)
    {
        // Apply Runge Kutta Formulas to find
        // next value of y
        k1 = specific_force_field(r0, h2);
        struct Vec3 r0_1 = r0 + (k1*step_size*0.5);
        k2 = specific_force_field(r0_1, h2);
        r0_1 = r0 + (k2*step_size*0.5);
        k3 = specific_force_field(r0_1, h2);
        r0_1 = r0 + (k3*step_size);
        k4 = specific_force_field(r0_1,h2);
        // Update next value of y
        Vec3 inc = ((k1 + k2 * 2+ k3 * 2 + k4)*step_size/6.0);
        r0 = r0 + v0 * step_size; 

        v0 = v0 + inc;
        // Update next value of x
    }
}




void raytrace(struct body * body, struct background * background_field, std::string& filename){
    const Vec3 white(1023, 1023, 1023); // the red will likely need to substituted with surface parameters
    const Vec3 black(0, 0, 0);
    const Vec3 red(1023, 0, 0); 
    
    std::ofstream out(filename);

    Vec3 cam_center = Vec3(IMAGE_HEIGHT/2,IMAGE_WIDTH/2,500);
    const Sphere camera (cam_center, 1);

    out << "P3\n" << IMAGE_WIDTH << ' ' << IMAGE_HEIGHT << ' ' << "1023\n";
    struct Vec3 pix_col(black);


    Vec3 norm_cam = cam_center / sqr_norm(cam_center);
    Ray ray(Vec3(0,200,0), Vec3(,0,0));
    int steps = 1000;
    Vec3 * ray_path = (struct Vec3 *)malloc(steps*sizeof(struct Vec3));

    Vec3 r0 = (ray.o - body->center)* SYSTEM_SCALING_FACTOR;
    Vec3 r0_1 = Vec3(0,0,0);
    Vec3 v0 = ray.d;

    struct Vec3 * cross = cross_vec3(r0,v0);
    Vec3 tmp = *cross;// * ( body->mass * SOLAR_MASS);
    float h2 = sqr_norm(tmp);// SOLAR_MASS;
    printf("h2: %f\n", h2);
    free(cross);
    for(int i = 0; i < steps; i++){
        runge_kutta_4_vec3(r0, v0,*body, 0.00002, 1, h2);
        r0_1 = (r0 / SYSTEM_SCALING_FACTOR) + body->center;
        memcpy((void*)&ray_path[i],(void *)&r0_1,sizeof(struct Vec3));
        // if (r0_1.x)
        // r0_1.print();
    }
    
    for (int y = 0; y < IMAGE_HEIGHT; ++y) {
        for (int x = 0; x < IMAGE_WIDTH; ++x) {
            pix_col = black;
            
            // debug event horizon print
            if ((pow(body->center.x - x, 2) + pow(body->center.y -y,2)) <= pow(body->radius,2)){
                pix_col = white;
            // debug photon sphere print
            } else if (fabs((pow(body->center.x - x, 2) + pow(body->center.y -y,2)) - pow(1.5*body->radius,2)) <= 100){
                pix_col = Vec3(0,1023,0);
            }
            if (x == IMAGE_WIDTH/2 || y ==IMAGE_HEIGHT/2 ){
                pix_col = Vec3(0,255,0);
            }
            for (int i = 0; i < steps; i++){
                if(fabs(ray_path[i].x - x) <=1 && fabs(ray_path[i].y -y ) <= 1){
                    pix_col = red;
                    break;
                }
            }
            out << (int) pix_col.x << ' '
                << (int) pix_col.y << ' '
                << (int) pix_col.z << '\n';
        }
    }
    out.close();
}

int main(void){
    struct Vec3 * tmp_background_pixels = (struct Vec3 *)malloc(sizeof(Vec3)*IMAGE_HEIGHT*IMAGE_WIDTH);
    struct background tmp_background = background(IMAGE_HEIGHT, IMAGE_WIDTH, 
                                        tmp_background_pixels, IMAGE_HEIGHT*IMAGE_WIDTH); 
    std::string output_filename = "output/gravity_out.ppm";
    for (int i = 0; i < IMAGE_WIDTH; i++){
        for (int j = 0; j < IMAGE_HEIGHT; j++){
            tmp_background.pix_arr[i*IMAGE_WIDTH+j] = Vec3(0,0,0);
        }
    }
    struct Vec3 center = Vec3(IMAGE_WIDTH/2, IMAGE_HEIGHT/2, 0);
    struct body blackhole = body(1, 0.0, center);

    float radius = schwarzchild_radius(blackhole);

    if (radius > IMAGE_HEIGHT/10 || radius > IMAGE_WIDTH/10){
        SYSTEM_SCALING_FACTOR = (IMAGE_HEIGHT > IMAGE_WIDTH) ? 10.0*radius/IMAGE_WIDTH : 10.0*radius/IMAGE_WIDTH;
    } else{
        SYSTEM_SCALING_FACTOR = (IMAGE_HEIGHT > IMAGE_WIDTH) ? IMAGE_WIDTH/(10.0*radius): IMAGE_WIDTH/(10.0*radius);

    }
    SYSTEM_SCALING_FACTOR = SYSTEM_SCALING_FACTOR / ZOOM; 
    blackhole.radius = radius/SYSTEM_SCALING_FACTOR;
    
    
    raytrace(&blackhole, &tmp_background, output_filename);
}