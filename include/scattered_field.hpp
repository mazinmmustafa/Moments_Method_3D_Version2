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
struct scattered_field_args_t{
    basis_1d_t b_m;
    complex_t E_TM=0.0; 
    complex_t E_TE=0.0;
    real_t theta_i=0.0;
    real_t phi_i=0.0; 
    real_t k=0.0; 
    real_t eta=0.0;
    vector_t<real_t> r;
    //
    vector_t<real_t> unit_vector;
};

// Functions

#endif
