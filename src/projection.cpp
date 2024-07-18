//
#include "projection.hpp"

const real_t eps_tol=1.0E-4;

projection_1d_t get_projection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> p, const real_t lambda){
    projection_1d_t para_1d;
    vector_t<real_t> v21=v2-v1;
    real_t alpha=v21*(p-v1)/(v21*v21);
    vector_t<real_t> p0=v1+alpha*v21;
    vector_t<real_t> l=unit(v2-v1);
    para_1d.p0 = p0;
    para_1d.P0 = mag(p0-p);
    para_1d.P0_u = para_1d.P0<(eps_tol*lambda) ? vector_t<real_t>(0.0, 0.0, 0.0) : unit(p0-p);
    para_1d.l_m = (v1-p0)*l;
    para_1d.l_p = (v2-p0)*l;
    para_1d.P_m = mag(p-v1);
    para_1d.P_p = mag(p-v2);
    return para_1d;
}

void get_projection_2d_edge(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, projection_1d_t &para_1d, vector_t<real_t> &p0_, 
    real_t &d, real_t &R0, real_t &R_m, real_t &R_p, vector_t<real_t> &u, vector_t<real_t> &n, 
    const real_t lambda){
    vector_t<real_t> v21=v2-v1;
    vector_t<real_t> v31=v3-v1;
    vector_t<real_t> v32=v3-v2;
    real_t a1=v21*v21; real_t a2=v21*v31; real_t a3=v21*(v1-p);
    real_t b1=v31*v31; real_t b2=v31*v21; real_t b3=v31*(v1-p);
    real_t alpha=(a2*b3-b1*a3)/(a1*b1-a2*b2);
    real_t beta=(b2*a3-a1*b3)/(a1*b1-a2*b2);
    vector_t<real_t> p0=v1+alpha*v21+beta*v31;
    p0_ = p0;
    para_1d = get_projection_1d(v1, v2, p0, lambda);
    d = mag(p-p0);
    R0 = mag(p-para_1d.p0);
    R_m = mag(p-v1);
    R_p = mag(p-v2);
    n = unit(v21^v32);
    u = unit(v21^(v21^v32));
}

projection_2d_t get_projection_2d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, const real_t lambda){
    projection_2d_t para_2d;
    get_projection_2d_edge(v2, v3, v1, p, para_2d.para_1d[0], para_2d.p0[0],
        para_2d.d, para_2d.R0[0], para_2d.R_m[0], 
        para_2d.R_p[0], para_2d.u[0], para_2d.n, lambda);
    get_projection_2d_edge(v3, v1, v2, p, para_2d.para_1d[1], para_2d.p0[1],
        para_2d.d, para_2d.R0[1], para_2d.R_m[1], 
        para_2d.R_p[1], para_2d.u[1], para_2d.n, lambda);
    get_projection_2d_edge(v1, v2, v3, p, para_2d.para_1d[2], para_2d.p0[2],
        para_2d.d, para_2d.R0[2], para_2d.R_m[2], 
        para_2d.R_p[2], para_2d.u[2], para_2d.n, lambda);
    return para_2d;
}

void get_projection_3d_triangle(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, vector_t<real_t> &n,
    projection_2d_t &para_2d, const real_t lambda){
    para_2d = get_projection_2d(v1, v2, v3, p, lambda);
    n = para_2d.n;
}

projection_3d_t get_projection_3d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> v4, const vector_t<real_t> p, 
    const real_t lambda){
    projection_3d_t para_3d;
    get_projection_3d_triangle(v2, v3, v4, p, para_3d.n[0], para_3d.para_2d[0], lambda);
    get_projection_3d_triangle(v3, v1, v4, p, para_3d.n[1], para_3d.para_2d[1], lambda);
    get_projection_3d_triangle(v1, v2, v4, p, para_3d.n[2], para_3d.para_2d[2], lambda);
    get_projection_3d_triangle(v3, v2, v1, p, para_3d.n[3], para_3d.para_2d[3], lambda);
    return para_3d;
}