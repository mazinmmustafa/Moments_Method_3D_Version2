#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
//

// Definitions
class shape{
    private:

    public:
        shape(){}
        ~shape(){}

};

// Functions
void call_gmsh(const real_t tol);
void create_vertical_wire_dipole(const real_t length, const real_t port_length);
void create_sphere(const real_t radius);

#endif