//
#include "integrand_projection_3d.hpp"

const real_t projection_tol=1.0E-4;

// 3d
void integrand_L1_3d_3d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_3d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp){
    projection_3d_para para;
    vector_t<real_t> p;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2];
    para = prjection_3d(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], p);
    I_mm = 0.0;
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_mm = I_mm+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2];
    para = prjection_3d(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], p);
    I_mp = 0.0;
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_mp = I_mp+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
    // pm
    p = b_m.r_p+alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0];
    para = prjection_3d(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], p);
    I_pm = 0.0;
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_pm = I_pm+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
    // pp
    p = b_m.r_p+alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0];
    para = prjection_3d(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], p);
    I_pp = 0.0;
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_pp = I_pp+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
}

void integrand_L2_3d_3d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_3d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp){
    projection_3d_para para;
    vector_t<real_t> p;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2];
    para = prjection_3d(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], p);
    I_mm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=0.5*para.para_2d[j].para_1d[i].P_0*(pow(para.para_2d[j].R_0[i], 2.0)+2.0*pow(para.para_2d[j].d, 2.0));
            real_t C=
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t D=0.5*para.para_2d[j].para_1d[i].P_0*(
                    para.para_2d[j].R_p[i]*para.para_2d[j].para_1d[i].l_p-
                    para.para_2d[j].R_m[i]*para.para_2d[j].para_1d[i].l_m
            );
            real_t E=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t F=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_mm = I_mm+(1.0/3.0)*A*(B*C+D-pow(abs(para.para_2d[j].d), 3.0)*(E-F))*para.para_2d[j].n;
        }
    }
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2];
    para = prjection_3d(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], p);
    I_mp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=0.5*para.para_2d[j].para_1d[i].P_0*(pow(para.para_2d[j].R_0[i], 2.0)+2.0*pow(para.para_2d[j].d, 2.0));
            real_t C=
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t D=0.5*para.para_2d[j].para_1d[i].P_0*(
                    para.para_2d[j].R_p[i]*para.para_2d[j].para_1d[i].l_p-
                    para.para_2d[j].R_m[i]*para.para_2d[j].para_1d[i].l_m
            );
            real_t E=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t F=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_mp = I_mp+(1.0/3.0)*A*(B*C+D-pow(abs(para.para_2d[j].d), 3.0)*(E-F))*para.para_2d[j].n;
        }
    }
    // pm
    p = b_m.r_p+alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0];
    para = prjection_3d(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], p);
    I_pm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=0.5*para.para_2d[j].para_1d[i].P_0*(pow(para.para_2d[j].R_0[i], 2.0)+2.0*pow(para.para_2d[j].d, 2.0));
            real_t C=
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t D=0.5*para.para_2d[j].para_1d[i].P_0*(
                    para.para_2d[j].R_p[i]*para.para_2d[j].para_1d[i].l_p-
                    para.para_2d[j].R_m[i]*para.para_2d[j].para_1d[i].l_m
            );
            real_t E=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t F=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_pm = I_pm+(1.0/3.0)*A*(B*C+D-pow(abs(para.para_2d[j].d), 3.0)*(E-F))*para.para_2d[j].n;
        }
    }
    // pp
    p = b_m.r_p+alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0];
    para = prjection_3d(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], p);
    I_pp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=0.5*para.para_2d[j].para_1d[i].P_0*(pow(para.para_2d[j].R_0[i], 2.0)+2.0*pow(para.para_2d[j].d, 2.0));
            real_t C=
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t D=0.5*para.para_2d[j].para_1d[i].P_0*(
                    para.para_2d[j].R_p[i]*para.para_2d[j].para_1d[i].l_p-
                    para.para_2d[j].R_m[i]*para.para_2d[j].para_1d[i].l_m
            );
            real_t E=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t F=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_pp = I_pp+(1.0/3.0)*A*(B*C+D-pow(abs(para.para_2d[j].d), 3.0)*(E-F))*para.para_2d[j].n;
        }
    }
}

void integrand_L3_3d_3d(const real_t alpha, const real_t beta, const real_t gamma, basis_3d_t b_m, basis_3d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp){
    projection_3d_para para;
    vector_t<real_t> p;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2];
    para = prjection_3d(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], p);
    I_mm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_mm = I_mm+A*(B-abs(para.para_2d[j].d)*(C-D))*para.para_2d[j].n;
        }
    }
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2];
    para = prjection_3d(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], p);
    I_mp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_mp = I_mp+A*(B-abs(para.para_2d[j].d)*(C-D))*para.para_2d[j].n;
        }
    }
    // pm
    p = b_m.r_p+alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0];
    para = prjection_3d(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], p);
    I_pm = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_pm = I_pm+A*(B-abs(para.para_2d[j].d)*(C-D))*para.para_2d[j].n;
        }
    }
    // pp
    p = b_m.r_p+alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0];
    para = prjection_3d(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], p);
    I_pp = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            I_pp = I_pp+A*(B-abs(para.para_2d[j].d)*(C-D))*para.para_2d[j].n;
        }
    }
}
