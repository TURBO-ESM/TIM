#pragma once

#include <stdint.h>
#include <stddef.h>

struct RealArray_C {
    double* data;   // Pointer to multidimensional array
    int* shape;     // An array of dimension extents
    int* lb;        // Lower bounds
    int* ub;        // Upper bounds
    int dim;        // The number of dimensions
};
struct Box_C {
    int* idxS;
    int* idxE;
};
