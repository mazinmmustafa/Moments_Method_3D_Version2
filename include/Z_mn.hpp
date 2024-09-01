#ifndef __Z_MN_HPP__
#define __Z_MN_HPP__

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
#include "psi_integrand.hpp"
#include "phi_integrand.hpp"
#include "delta_integrand.hpp"
#include "nu_integrand.hpp"

// Definitions

// Functions
complex_t Z_mn_1d_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, const real_t a, 
    const quadl_domain_t quadl, int_t &flag);
complex_t Z_mn_2d_2d(const basis_2d_t b_m, const basis_2d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, 
    const quadl_domain_t quadl, int_t &flag);
complex_t Z_mn_3d_3d(const basis_3d_t b_m, const basis_3d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, const complex_t eps_b, 
    const quadl_domain_t quadl, int_t &flag);

#endif
