//
#include "projection.hpp"

projection_1d_para prjection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2, 
    const vector_t<real_t> p){
    vector_t<real_t> v21=v2-v1;
    const real_t alpha=v21*(p-v1)/(v21*v21);
    vector_t<real_t> p_0=v1+alpha*v21;
    projection_1d_para para;
    para.p_0 = p_0;
    vector_t<real_t> P_0=para.p_0-p;
    para.P_0 = mag(P_0);
    para.P_0_unit = unit(P_0);
    para.l_unit = unit(v2-v1);
    para.l_m = (v1-para.p_0)*para.l_unit;
    para.l_p = (v2-para.p_0)*para.l_unit;
    para.P_m = mag(v1-p);
    para.P_p = mag(v2-p);
    return para;
}

projection_2d_para prjection_2d(const vector_t<real_t> v1, const vector_t<real_t> v2, 
    const vector_t<real_t> v3, const vector_t<real_t> p){
    vector_t<real_t> v21=v2-v1;
    vector_t<real_t> v31=v3-v1;
    real_t a1=v21*v21;
    real_t a2=v21*v31;
    real_t a3=v21*(v1-p);
    real_t b1=v31*v31;
    real_t b2=v31*v21;
    real_t b3=v31*(v1-p);
    const real_t alpha=(a2*b3-b1*a3)/(a1*b1-a2*b2);
    const real_t beta =(b2*a3-a1*b3)/(a1*b1-a2*b2);
    vector_t<real_t> p_0=v1+alpha*v21+beta*v31;
    projection_2d_para para;
    para.para_1d[0] = prjection_1d(v2, v3, p_0);
    para.para_1d[1] = prjection_1d(v3, v1, p_0);
    para.para_1d[2] = prjection_1d(v1, v2, p_0);
    const vector_t<real_t> n=unit(v21^v31);
    para.d = (p-p_0)*n;
    para.R_0[0] = mag(p-para.para_1d[0].p_0);
    para.R_0[1] = mag(p-para.para_1d[1].p_0);
    para.R_0[2] = mag(p-para.para_1d[2].p_0);
    para.R_m[0] = mag(p-v2);
    para.R_p[0] = mag(p-v3);
    para.R_m[1] = mag(p-v3);
    para.R_p[1] = mag(p-v1);
    para.R_m[2] = mag(p-v1);
    para.R_p[2] = mag(p-v2);
    para.u[0] = unit((v3-v2)^n);
    para.u[1] = unit((v1-v3)^n);
    para.u[2] = unit((v2-v1)^n);
    para.n = n;
    para.p_0 = p_0;
    return para;
}