//
#include "integrand_projection_1d.hpp"

const real_t projection_tol=1.0E-4;

real_t sign_function(const real_t l_m, const real_t l_p){
    if ((l_m*l_p)>0.0){
        return 0.0;
    }else{
        return 1.0;
    }
}

// 1d
void integrand_L1_1d_1d(const real_t alpha, basis_1d_t b_m, basis_1d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp, const real_t a, const real_t lambda){
    assert(lambda>0.0);
    projection_1d_para para;
    vector_t<real_t> p;
    real_t ans;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    // I_mm = para.P_0<projection_tol*lambda ? sign_function(para.l_m, para.l_p) : ans;
    I_mm = ans;
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.e[0], b_n.r_p , p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    // I_mp = para.P_0<projection_tol*lambda ? sign_function(para.l_m, para.l_p) : ans;
    I_mp = ans;
    // pm
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    // I_pm = para.P_0<projection_tol*lambda ? sign_function(para.l_m, para.l_p) : ans;
    I_pm = ans;
    // pp
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.e[0], b_n.r_p , p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    // I_pp = para.P_0<projection_tol*lambda ? sign_function(para.l_m, para.l_p) : ans;
    I_pp = ans;
}

void integrand_L2_1d_1d(const real_t alpha, basis_1d_t b_m,basis_1d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp, 
    const real_t a){
    projection_1d_para para;
    vector_t<real_t> p;
    vector_t<real_t> ans;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_mm = ans;
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_mp = ans;
    // pm
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_pm = ans;
    // pp
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_pp = ans;
}

void integrand_L3_1d_1d(const real_t alpha, basis_1d_t b_m,basis_1d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp, 
    const real_t lambda){
    projection_1d_para para;
    vector_t<real_t> p;
    vector_t<real_t> ans;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_mm = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_mp = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // pm
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_pm = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // pp
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_pp = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
}
