//
#include "psi_integrand.hpp"

// 1d 1d
complex_t psi_1d_1d_singular_integrand_inner(const complex_t alpha_, void *args_){
    integrand_1d_1d_args *args=(integrand_1d_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    const complex_t k=args->k;
    const real_t a=args->a;
    complex_t alpha=args->alpha;
    complex_t I_mm, I_mp, I_pm, I_pp;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mm, R_mp, R_pm, R_pp;
    R_mn_1d_1d(real(alpha), real(alpha_), b_m, b_n, R_mm, R_mp, R_pm, R_pp, a);
    I_mm = +1.0*(b_m.L_m[0]*b_n.L_m[0])*alpha*alpha_*(-j*k*exp(-j*k*R_mm/2.0))*sinc(k*R_mm/2.0);
    I_mp = -1.0*(b_m.L_m[0]*b_n.L_p[0])*alpha*alpha_*(-j*k*exp(-j*k*R_mp/2.0))*sinc(k*R_mp/2.0);
    I_pm = -1.0*(b_m.L_p[0]*b_n.L_m[0])*alpha*alpha_*(-j*k*exp(-j*k*R_pm/2.0))*sinc(k*R_pm/2.0);
    I_pp = +1.0*(b_m.L_p[0]*b_n.L_p[0])*alpha*alpha_*(-j*k*exp(-j*k*R_pp/2.0))*sinc(k*R_pp/2.0);
    return (I_mm+I_mp+I_pm+I_pp)/(4.0*pi);
}

complex_t psi_1d_1d_singular_integrand_outer(const complex_t alpha, void *args_){
    integrand_1d_1d_args *args=(integrand_1d_1d_args*)args_;
    args->alpha = alpha;
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
    int flag;
    complex_t ans=args->quadl.integral_1d(psi_1d_1d_singular_integrand_inner, args, line, flag);
    return ans;
}

complex_t psi_1d_1d_integrand_1(const complex_t alpha, void *args_){
    integrand_1d_1d_args *args=(integrand_1d_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    real_t a=args->a;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*real(alpha)*b_m.L_m[0];
    rho_p = -1.0*real(alpha)*b_m.L_p[0];
    integrand_L1_1d_1d(real(alpha), b_m, b_n, I_mm, I_mp, I_pm, I_pp, a);
    projection_1d_para para_m_m = prjection_1d(b_n.r_m, b_n.e[0], b_m.r_m+real(alpha)*b_m.L_m[0]);
    projection_1d_para para_m_p = prjection_1d(b_n.e[0], b_n.r_p, b_m.r_m+real(alpha)*b_m.L_m[0]);
    projection_1d_para para_p_m = prjection_1d(b_n.r_m, b_n.e[0], b_m.r_p+real(alpha)*b_m.L_p[0]);
    projection_1d_para para_p_p = prjection_1d(b_n.e[0], b_n.r_p, b_m.r_p+real(alpha)*b_m.L_p[0]);
    complex_t ans=0.0;
    ans+= -1.0*(rho_m*(para_m_m.l_m*para_m_m.l_unit))*I_mm/mag(b_n.L_m[0]);
    ans+= +1.0*(rho_m*(para_m_p.l_p*para_m_p.l_unit))*I_mp/mag(b_n.L_p[0]);
    ans+= -1.0*(rho_p*(para_p_m.l_m*para_p_m.l_unit))*I_pm/mag(b_n.L_m[0]);
    ans+= +1.0*(rho_p*(para_p_p.l_p*para_p_p.l_unit))*I_pp/mag(b_n.L_p[0]);
    return ans/(4.0*pi);
}

complex_t psi_1d_1d_integrand_2(const complex_t alpha, void *args_){
    integrand_1d_1d_args *args=(integrand_1d_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    real_t a=args->a;
    vector_t<real_t> I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*real(alpha)*b_m.L_m[0];
    rho_p = -1.0*real(alpha)*b_m.L_p[0];
    integrand_L2_1d_1d(real(alpha), b_m, b_n, I_mm, I_mp, I_pm, I_pp, a);
    complex_t ans=0.0;
    ans+= +1.0*rho_m*I_mm/mag(b_n.L_m[0]);
    ans+= -1.0*rho_m*I_mp/mag(b_n.L_p[0]);
    ans+= +1.0*rho_p*I_pm/mag(b_n.L_m[0]);
    ans+= -1.0*rho_p*I_pp/mag(b_n.L_p[0]);
    return ans/(4.0*pi);
}

complex_t psi_1d_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, const real_t lambda, 
    const real_t a, const quadl_domain_t quadl, int_t &flag){
    integrand_1d_1d_args args={quadl, b_m, b_n};
    args.a = a;
    args.k = k;
    args.lambda = lambda;
    complex_t I1=0.0, I2=0.0, I3=0.0;
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    I1 = args.quadl.integral_1d(psi_1d_1d_singular_integrand_outer, &args, line, flag);
    assert_error(!flag, "no convergence");
    I2 = args.quadl.integral_1d(psi_1d_1d_integrand_1, &args, line, flag);
    assert_error(!flag, "no convergence");
    I3 = args.quadl.integral_1d(psi_1d_1d_integrand_2, &args, line, flag);
    assert_error(!flag, "no convergence");
    return I1+I2+I3;
}

// 1d 2d
// 1d 3d

// 2d 1d
// 2d 2d
// 2d 3d

// 3d 1d
// 3d 2d
// 3d 3d