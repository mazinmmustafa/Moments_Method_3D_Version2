#ifndef __TESTBENCH_HPP__
#define __TESTBENCH_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
#include "shape.hpp"
#include "projection.hpp"
#include "engine.hpp"

// Definitions

// Functions
void test_utilities();
void test_gmsh();
void test_shape();

void test_engine_1d_1d();
void test_engine_1d_vertical_dipole();

#endif