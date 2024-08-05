#ifndef __SAHRED_DEFINITIONS_HPP__
#define __SAHRED_DEFINITIONS_HPP__

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

// 1d
struct integrand_1d_1d_args{
    quadl_domain_t quadl;
    basis_1d_t b_m;
    basis_1d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

struct integrand_1d_2d_args{
    quadl_domain_t quadl;
    basis_1d_t b_m;
    basis_2d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

struct integrand_1d_3d_args{
    quadl_domain_t quadl;
    basis_1d_t b_m;
    basis_3d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

// 2d
struct integrand_2d_1d_args{
    quadl_domain_t quadl;
    basis_2d_t b_m;
    basis_1d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

struct integrand_2d_2d_args{
    quadl_domain_t quadl;
    basis_2d_t b_m;
    basis_2d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

struct integrand_2d_3d_args{
    quadl_domain_t quadl;
    basis_2d_t b_m;
    basis_3d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

// 3d
struct integrand_3d_1d_args{
    quadl_domain_t quadl;
    basis_3d_t b_m;
    basis_1d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

struct integrand_3d_2d_args{
    quadl_domain_t quadl;
    basis_3d_t b_m;
    basis_2d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

struct integrand_3d_3d_args{
    quadl_domain_t quadl;
    basis_3d_t b_m;
    basis_3d_t b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

// Functions

#endif
