#ifndef __INTEGRAND_PROJECTION_3D_HPP__
#define __INTEGRAND_PROJECTION_3D_HPP__

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

// 3d
void integrand_L1_3d_3d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_3d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp);
void integrand_L2_3d_3d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_3d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp);
void integrand_L3_3d_3d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_3d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp);
#endif
