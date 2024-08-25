//
#include "integrand_projection_2d.hpp"

const real_t projection_tol=1.0E-4;

// 2d
void integrand_L1_2d_2d(const real_t alpha, const real_t beta, basis_2d_t b_m, basis_2d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp){
    projection_2d_para para;
    vector_t<real_t> p;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1];
    para = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], p);
    I_mm = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_mm+=A*(B-abs(para.d)*(C-D));
    }
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1];
    para = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], p);
    I_mp = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_mp+=A*(B-abs(para.d)*(C-D));
    }
    // pm
    p = b_m.r_p+alpha*b_m.L_p[1]+beta*b_m.L_p[0];
    para = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], p);
    I_pm = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_pm+=A*(B-abs(para.d)*(C-D));
    }
    // pp
    p = b_m.r_p+alpha*b_m.L_p[1]+beta*b_m.L_p[0];
    para = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], p);
    I_pp = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_pp+=A*(B-abs(para.d)*(C-D));
    }
}

void integrand_L2_2d_2d(const real_t alpha, const real_t beta, basis_2d_t b_m, basis_2d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp){
    projection_2d_para para;
    vector_t<real_t> p;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1];
    para = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], p);
    I_mm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        I_mm = I_mm+0.5*(A+B)*para.u[i];
    }
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1];
    para = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], p);
    I_mp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        I_mp = I_mp+0.5*(A+B)*para.u[i];
    }
    // pm
    p = b_m.r_p+alpha*b_m.L_p[1]+beta*b_m.L_p[0];
    para = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], p);
    I_pm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        I_pm = I_pm+0.5*(A+B)*para.u[i];
    }
    // pp
    p = b_m.r_p+alpha*b_m.L_p[1]+beta*b_m.L_p[0];
    para = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], p);
    I_pp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        I_pp = I_pp+0.5*(A+B)*para.u[i];
    }
}

void integrand_L3_2d_2d(const real_t alpha, const real_t beta, basis_2d_t b_m, basis_2d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp){
    projection_2d_para para;
    vector_t<real_t> p;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1];
    para = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], p);
    I_mm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_mm = I_mm+A*para.u[i]+sign(para.d)*(B-C)*para.n;
    }
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1];
    para = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], p);
    I_mp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_mp = I_mp+A*para.u[i]+sign(para.d)*(B-C)*para.n;
    }
    // pm
    p = b_m.r_p+alpha*b_m.L_p[1]+beta*b_m.L_p[0];
    para = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], p);
    I_pm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_pm = I_pm+A*para.u[i]+sign(para.d)*(B-C)*para.n;
    }
    // pp
    p = b_m.r_p+alpha*b_m.L_p[1]+beta*b_m.L_p[0];
    para = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], p);
    I_pp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_pp = I_pp+A*para.u[i]+sign(para.d)*(B-C)*para.n;
    }
}
