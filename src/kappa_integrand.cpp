//
#include "kappa_integrand.hpp"

// 3d 3d
complex_t kappa_3d_3d_singular_integrand_inner(const complex_t alpha_, const complex_t beta_, void *args_){
    integrand_3d_3d_args *args=(integrand_3d_3d_args*)args_;
    basis_3d_t b_m=args->b_m;
    basis_3d_t b_n=args->b_n;
    const complex_t k=args->k;
    complex_t alpha=args->alpha;
    complex_t beta=args->beta;
    complex_t gamma=args->gamma;
    complex_t I_mm, I_mp, I_pm, I_pp;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mm, R_mp, R_pm, R_pp;
    complex_t gamma_=1.0-alpha_-beta_;
    R_mn_3d_3d(real(alpha), real(beta), real(gamma), real(alpha_), real(beta_), real(gamma_), 
        b_m, b_n, R_mm, R_mp, R_pm, R_pp);
    vector_t<real_t> rho_n_m=real(alpha_)*b_n.L_m[0]+real(beta_)*b_n.L_m[1]+real(gamma_)*b_n.L_m[2];
    vector_t<real_t> rho_n_p=real(alpha_)*b_n.L_p[2]+real(beta_)*b_n.L_p[1]+real(gamma_)*b_n.L_p[0];
    I_mm = +1.0*b_n.nA*rho_n_m*(-j*k*exp(-j*k*R_mm/2.0))*sinc(k*R_mm/2.0);
    I_mp = -1.0*b_n.nA*rho_n_p*(-j*k*exp(-j*k*R_mp/2.0))*sinc(k*R_mp/2.0);
    I_pm = -1.0*b_n.nA*rho_n_m*(-j*k*exp(-j*k*R_pm/2.0))*sinc(k*R_pm/2.0);
    I_pp = +1.0*b_n.nA*rho_n_p*(-j*k*exp(-j*k*R_pp/2.0))*sinc(k*R_pp/2.0);
    return 4.0*b_m.A*b_n.A*b_n.A*(I_mm/b_n.V_m+I_mp/b_n.V_p+I_pm/b_n.V_m+I_pp/b_n.V_p)/(4.0*pi);
}

complex_t kappa_3d_3d_singular_integrand_outer(const complex_t alpha, const complex_t beta, const complex_t gamma, void *args_){
    integrand_3d_3d_args *args=(integrand_3d_3d_args*)args_;
    args->alpha = alpha;
    args->beta = beta;
    args->gamma = gamma;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    complex_t ans=args->quadl.integral_2d(kappa_3d_3d_singular_integrand_inner, args, triangle, flag);
    return ans;
}

complex_t kappa_3d_3d_integrand_1(const complex_t alpha, const complex_t beta, const complex_t gamma, void *args_){
    integrand_3d_3d_args *args=(integrand_3d_3d_args*)args_;
    basis_3d_t b_m=args->b_m;
    basis_3d_t b_n=args->b_n;
    real_t I_mm, I_mp, I_pm, I_pp;
    basis_2d_t b_n_=basis_2d_t(b_n.e[0], b_n.e[1], b_n.e[2], b_n.e[0], b_n.pg_m, b_n.pg_p);
    integrand_L1_3d_2d(real(alpha), real(beta), real(gamma), b_m, b_n_, I_mm, I_mp, I_pm, I_pp);
    complex_t ans=0.0;  
    ans+= -1.0*(b_n.nA*(b_n.e[0]-b_m.r_m))*I_mm/(b_n.V_m);
    ans+= +1.0*(b_n.nA*(b_n.e[0]-b_m.r_m))*I_mp/(b_n.V_p);
    ans+= +1.0*(b_n.nA*(b_n.e[0]-b_m.r_p))*I_pm/(b_n.V_m);
    ans+= -1.0*(b_n.nA*(b_n.e[0]-b_m.r_p))*I_pp/(b_n.V_p);
    return 2.0*b_m.A*b_n.A*ans/(4.0*pi);
}

complex_t kappa_3d_3d_integrand_2(const complex_t alpha, const complex_t beta, const complex_t gamma, void *args_){
    integrand_3d_3d_args *args=(integrand_3d_3d_args*)args_;
    basis_3d_t b_m=args->b_m;
    basis_3d_t b_n=args->b_n;
    vector_t<real_t> I_mm, I_mp, I_pm, I_pp;
    basis_2d_t b_n_=basis_2d_t(b_n.e[0], b_n.e[1], b_n.e[2], b_n.e[0], b_n.pg_m, b_n.pg_p);
    integrand_L2_3d_2d(real(alpha), real(beta), real(gamma), b_m, b_n_, I_mm, I_mp, I_pm, I_pp);
    complex_t ans=0.0;
    ans+= +1.0*b_n.nA*I_mm/(b_n.V_m);
    ans+= -1.0*b_n.nA*I_mp/(b_n.V_p);
    ans+= -1.0*b_n.nA*I_pm/(b_n.V_m);
    ans+= +1.0*b_n.nA*I_pp/(b_n.V_p);
    return 2.0*b_m.A*b_n.A*ans/(4.0*pi);
}

complex_t kappa_3d_3d(const basis_3d_t b_m, const basis_3d_t b_n, const complex_t k, const real_t lambda, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_3d_3d_args args={quadl, b_m, b_n};
    args.k = k;
    args.lambda = lambda;
    complex_t I1=0.0, I2=0.0, I3=0.0;
    tetrahedron_domain_t tetrahedron={vector_t<real_t>(0.0, 0.0, 0.0), 
        vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0), 
        vector_t<real_t>(0.0, 0.0, 1.0)};
    I1 = args.quadl.integral_3d(kappa_3d_3d_singular_integrand_outer, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    I2 = args.quadl.integral_3d(kappa_3d_3d_integrand_1, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    I3 = args.quadl.integral_3d(kappa_3d_3d_integrand_2, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    return (b_n.eps_p-b_n.eps_m)*(I1+I2+I3);
}
