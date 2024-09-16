#ifndef __INTEGRAND_PROJECTION_2D_HPP__
#define __INTEGRAND_PROJECTION_2D_HPP__

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

// 2d
void integrand_L1_2d_2d(const real_t alpha, const real_t beta, basis_2d_t b_m, basis_2d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp);
void integrand_L2_2d_2d(const real_t alpha, const real_t beta, basis_2d_t b_m, basis_2d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp);
void integrand_L3_2d_2d(const real_t alpha, const real_t beta, basis_2d_t b_m, basis_2d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp);
//
void integrand_L1_3d_2d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_2d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp);
void integrand_L2_3d_2d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_2d_t b_n, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p);

#endif
