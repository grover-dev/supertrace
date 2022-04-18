#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "shapes.hpp"


// Both in pixel count
#define IMAGE_HEIGHT 500
#define IMAGE_WIDTH  500



// Search script for "RK4" to see what this is supposed to do
float specific_force_field(float y, float h2)
{
    float vector = 0.0f; // This should be the position vel vector for the point
    return (-1.5f * h2 * y / pow(vector, 2.5));
}


// Approximate the solution to first order ODE
// x0 and y0 are initial condition
// x is function input
// step_size is the step size
float runge_kutta_4(float x0, float y0, float x, float step_size)
{
    // Count number of iterations using step size or
    // step height h
    int n = (int)((x - x0) / step_size);
 
    float k1, k2, k3, k4, k5;
 
    // Iterate for number of iterations
    float y = y0;
    for (int i=1; i<=n; i++)
    {
        // Apply Runge Kutta Formulas to find
        // next value of y
        k1 = h*specific_force_field(x0, y);
        k2 = h*specific_force_field(x0 + 0.5*step_size, y + 0.5*k1);
        k3 = h*specific_force_field(x0 + 0.5*step_size, y + 0.5*k2);
        k4 = h*specific_force_field(x0 + step_size, y + k3);
 
        // Update next value of y
        y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
 
        // Update next value of x
        x0 = x0 + h;
    }
 
    return y;
}