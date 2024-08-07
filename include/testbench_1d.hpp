#ifndef __TESTBENCH_1D_HPP__
#define __TESTBENCH_1D_HPP__

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
void test_engine_1d_vertical_dipole_input_adminttance();
void test_engine_1d_loop_input_impedance();
void test_engine_1d_vertical_dipole_mutual_impedance();
void test_engine_1d_transmission_line_S_parameters();
void test_engine_1d_RCS_vertical_wire();

#endif