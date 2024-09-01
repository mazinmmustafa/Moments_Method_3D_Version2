#ifndef __NU_INTEGRAND_1D_HPP__
#define __NU_INTEGRAND_1D_HPP__

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
//
#include "shared_definitions.hpp"
#include "R_mn.hpp"
#include "integrand_projection_1d.hpp"
#include "integrand_projection_2d.hpp"
#include "integrand_projection_3d.hpp"

// Definitions

// Functions
complex_t nu_3d_3d(const basis_3d_t b_m, const basis_3d_t b_n, const real_t lambda, const complex_t eps_b);

#endif
