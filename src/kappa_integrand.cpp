//
#include "kappa_integrand.hpp"

// 3d 3d
complex_t kappa_3d_3d_singular_integrand_outer(const complex_t alpha, const complex_t beta, const complex_t gamma, void *args_){
    integrand_3d_3d_args *args=(integrand_3d_3d_args*)args_;
    basis_3d_t b_m=args->b_m;
    basis_3d_t b_n=args->b_n;
    const complex_t k=args->k;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_m, R_p;
    vector_t<real_t> rho_m=b_n.r_m+b_n.L_m[0]/3.0+b_n.L_m[1]/3.0+b_n.L_m[2]/3.0;
    R_m = mag(b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]+real(gamma)*b_m.L_m[2]-rho_m);
    R_p = mag(b_m.r_p+real(alpha)*b_m.L_p[2]+real(beta)*b_m.L_p[1]+real(gamma)*b_m.L_p[0]-rho_m);
    real_t factor=b_n.nA*(b_n.L_m[0]/3.0+b_n.L_m[1]/3.0+b_n.L_m[2]/3.0);
    I_m = +1.0*(-j*k*exp(-j*k*R_m/2.0))*sinc(k*R_m/2.0)*factor;
    I_p = -1.0*(-j*k*exp(-j*k*R_p/2.0))*sinc(k*R_p/2.0)*factor;
    return 2.0*b_m.A*b_n.A*(I_m/b_n.V_m+I_p/b_n.V_m)/(4.0*pi);
}


void integrand_L1_3d(basis_3d_t b_m, basis_3d_t b_n, real_t &I_m, real_t &I_p){
    projection_3d_para para;
    vector_t<real_t> p;
    // m
    p = b_n.r_m+(b_n.L_m[0]+b_n.L_m[1]+b_n.L_m[2])/3.0;
    para = prjection_3d(b_m.r_m, b_m.e[0], b_m.e[1], b_m.e[2], p);
    I_m = 0.0;
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
            I_m = I_m+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
    // p
    p = b_n.r_p+(b_n.L_p[2]+b_n.L_p[1]+b_n.L_p[0])/3.0;
    para = prjection_3d(b_m.r_p, b_m.e[2], b_m.e[1], b_m.e[0], p);
    I_p = 0.0;
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
            I_p = I_p+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
}

complex_t kappai_3d_3d_integrand_1(void *args_){
    integrand_3d_3d_args *args=(integrand_3d_3d_args*)args_;
    basis_3d_t b_m=args->b_m;
    basis_3d_t b_n=args->b_n;
    real_t I_m, I_p;
    integrand_L1_3d(b_m, b_n, I_m, I_p);
    vector_t<real_t> rho_m=b_n.r_m+b_n.L_m[0]/3.0+b_n.L_m[1]/3.0+b_n.L_m[2]/3.0;
    real_t factor=b_n.nA*(b_n.L_m[0]/3.0+b_n.L_m[1]/3.0+b_n.L_m[2]/3.0);
    I_m = +1.0*I_m*factor/(3.0*b_n.V_m);
    I_p = -1.0*I_p*factor/(3.0*b_n.V_m);
    return b_m.A*b_n.A*(I_m/b_m.V_m+I_p/b_m.V_p)/(4.0*pi);
}

complex_t kappa_3d_3d(const basis_3d_t b_m, const basis_3d_t b_n, const complex_t k, const real_t lambda, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_3d_3d_args args={quadl, b_m, b_n};
    args.k = k;
    args.lambda = lambda;
    complex_t I1=0.0, I2=0.0;
    tetrahedron_domain_t tetrahedron={vector_t<real_t>(0.0, 0.0, 0.0), 
        vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0), 
        vector_t<real_t>(0.0, 0.0, 1.0)};
    I1 = args.quadl.integral_3d(kappa_3d_3d_singular_integrand_outer, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    I2 = kappai_3d_3d_integrand_1(&args);
    complex_t kappa_m=(b_n.eps_m-1.0)/b_n.eps_m;
    complex_t kappa_p=(b_n.eps_p-1.0)/b_n.eps_p;
    return (kappa_p-kappa_m)*(I1+I2);
}
