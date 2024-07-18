#ifndef __PROJECTION_HPP__
#define __PROJECTION_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
//

// Definitions
struct projection_1d_t{
    vector_t<real_t> p0, P0_u;
    real_t l_m, l_p;
    real_t P0, P_m, P_p;
};

struct projection_2d_t{
    projection_1d_t para_1d[3];
    vector_t<real_t> u[3], p0[3];
    real_t R0[3], R_m[3], R_p[3], d;
    vector_t<real_t> n;
};

struct projection_3d_t{
    projection_2d_t para_2d[4];
    vector_t<real_t> n[4];
};

// Functions
projection_1d_t get_projection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> p, const real_t lambda);
projection_2d_t get_projection_2d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, const real_t lambda);
projection_3d_t get_projection_3d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> v4, const vector_t<real_t> p, const real_t lambda);

void get_projection_2d_edge(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, projection_1d_t &para_1d, vector_t<real_t> &p0, 
    real_t &d, real_t &R0, real_t &R_m, real_t &R_p, vector_t<real_t> &u, vector_t<real_t> &n, const real_t lambda);
void get_projection_3d_triangle(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, vector_t<real_t> &n,
    projection_2d_t &para_2d, const real_t lambda);

#endif