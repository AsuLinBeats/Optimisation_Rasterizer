#pragma once

#include "Vec4_Optimised.h"
//#include "vec4.h"
// #include "colour.h"
#include "colour_optimised.h"

// keep light straightforward - struct for storing information
struct Light {
    vec4 omega_i; // light direction
    colour L; // light colour
    colour ambient; // ambient light component 
};

