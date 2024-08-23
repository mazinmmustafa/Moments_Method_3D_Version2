#ifndef __SCATTERED_FIELD_HPP__
#define __SCATTERED_FIELD_HPP__

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
struct scattered_field_args_1d_t{
    basis_1d_t b_m;
    complex_t k=0.0; 
    complex_t eta=0.0;
    real_t a=0.0;
    vector_t<real_t> r;
    vector_t<real_t> unit_vector;
};

struct scattered_field_args_2d_t{
    basis_2d_t b_m;
    complex_t k=0.0; 
    complex_t eta=0.0;
    vector_t<real_t> r;
    vector_t<real_t> unit_vector;
};

// Functions
complex_t compute_E_1d(const basis_1d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, 
    const real_t a, quadl_domain_t quadl);
complex_t compute_H_1d(const basis_1d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, 
    const real_t a, quadl_domain_t quadl);
complex_t compute_E_2d(const basis_2d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl);
complex_t compute_H_2d(const basis_2d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl);
#endif
