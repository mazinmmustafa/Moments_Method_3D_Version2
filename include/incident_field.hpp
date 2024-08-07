#ifndef __INCIDNET_FIELD_HPP__
#define __INCIDNET_FIELD_HPP__

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
struct incident_field_t{
    vector_t<complex_t> E, H;
};

struct incident_field_args_t{
    basis_1d_t b_m;
    complex_t E_TM=0.0; 
    complex_t E_TE=0.0;
    real_t theta_i=0.0;
    real_t phi_i=0.0; 
    real_t k=0.0; 
    real_t eta=0.0;
    vector_t<real_t> r, u, n;
};

// Functions
incident_field_t compute_incident_field(const complex_t E_TM, const complex_t E_TE, const real_t theta_i, const real_t phi_i, 
    const real_t k, const real_t eta, const vector_t<real_t> r);
complex_t compute_incident_E_integrand_1d(const complex_t alpha, void *args_);

#endif
