#ifndef __INTEGRAND_PROJECTION_1D_HPP__
#define __INTEGRAND_PROJECTION_1D_HPP__

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
#include "R_mn.hpp"

// Definitions

// Functions

// 1d
void integrand_L1_1d_1d(const real_t alpha, basis_1d_t b_m, basis_1d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp, const real_t a, const real_t lambda);
void integrand_L2_1d_1d(const real_t alpha, basis_1d_t b_m,basis_1d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp, 
    const real_t a);
void integrand_L3_1d_1d(const real_t alpha, basis_1d_t b_m,basis_1d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp, 
    const real_t a, const real_t lambda);

#endif
