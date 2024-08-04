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
    vector_t<real_t> P_0_unit, u, l_unit, p_0;
    projection_1d_para(){}
};

// Functions
projection_1d_para prjection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2, 
    const vector_t<real_t> p);

#endif
