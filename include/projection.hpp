#ifndef __PROJECTION_HPP__
#define __PROJECTION_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
#include "shape.hpp"
//

// Definitions
struct projection_1d_para{
    real_t l_m, l_p, P_m, P_p;
    real_t P_0;
    vector_t<real_t> P_0_unit, l_unit, p_0;
    projection_1d_para(){}
};

struct projection_2d_para{
    projection_1d_para para_1d[3];
    real_t R_0[3], R_m[3], R_p[3];
    real_t d;
    vector_t<real_t> u[3], n, p_0;
    projection_2d_para(){}
};

// Functions
projection_1d_para prjection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2, 
    const vector_t<real_t> p);
projection_2d_para prjection_2d(const vector_t<real_t> v1, const vector_t<real_t> v2, 
    const vector_t<real_t> v3, const vector_t<real_t> p);

#endif
