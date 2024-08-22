#ifndef __R_MN_HPP__
#define __R_MN_HPP__

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

// Definitions

// Functions
// 1d
void R_mn_1d_1d(const real_t alpha, 
    const real_t alpha_, 
    const basis_1d_t b_m, const basis_1d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a);
void R_mn_2d_2d(const real_t alpha, const real_t beta, const real_t alpha_, 
    const real_t beta_, const basis_2d_t b_m, const basis_2d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp);
    
#endif
