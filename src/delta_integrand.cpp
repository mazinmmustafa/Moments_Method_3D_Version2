//
#include "delta_integrand.hpp"

// 1d 1d
// 1d 2d
// 1d 3d

// 2d 1d
// 2d 2d
complex_t delta_2d_2d_singular_integrand_inner(const complex_t alpha_, const complex_t beta_, void *args_){
    integrand_2d_2d_args *args=(integrand_2d_2d_args*)args_;
    basis_2d_t b_m=args->b_m;
    basis_2d_t b_n=args->b_n;
    const complex_t k=args->k;
    complex_t alpha=args->alpha;
    complex_t beta=args->beta;
    complex_t I_mm, I_mp, I_pm, I_pp;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mm, R_mp, R_pm, R_pp;
    R_mn_2d_2d(real(alpha), real(beta), real(alpha_), real(beta_), b_m, b_n, R_mm, R_mp, R_pm, R_pp);
    vector_t<real_t> rho_m_m=real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1];
    vector_t<real_t> rho_m_p=real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0];
    vector_t<real_t> rho_n_m=real(alpha_)*b_n.L_m[0]+real(beta_)*b_n.L_m[1];
    vector_t<real_t> rho_n_p=real(alpha_)*b_n.L_p[1]+real(beta_)*b_n.L_p[0];
    I_mm = +1.0*rho_m_m*(rho_n_m^unit((b_n.r_m+rho_n_m)-(b_m.r_m+rho_m_m)))*-0.5*k*k*exp(-j*k*R_mm/2.0)*(sinc(k*R_mm/2.0)+j*sinc_dx(k*R_mm/2.0));
    I_mp = -1.0*rho_m_m*(rho_n_p^unit((b_n.r_p+rho_n_p)-(b_m.r_m+rho_m_m)))*-0.5*k*k*exp(-j*k*R_mp/2.0)*(sinc(k*R_mp/2.0)+j*sinc_dx(k*R_mp/2.0));
    I_pm = -1.0*rho_m_p*(rho_n_m^unit((b_n.r_m+rho_n_m)-(b_m.r_p+rho_m_p)))*-0.5*k*k*exp(-j*k*R_pm/2.0)*(sinc(k*R_pm/2.0)+j*sinc_dx(k*R_pm/2.0));
    I_pp = +1.0*rho_m_p*(rho_n_p^unit((b_n.r_p+rho_n_p)-(b_m.r_p+rho_m_p)))*-0.5*k*k*exp(-j*k*R_pp/2.0)*(sinc(k*R_pp/2.0)+j*sinc_dx(k*R_pp/2.0));
    return b_m.L*b_n.L*(I_mm+I_mp+I_pm+I_pp)/(4.0*pi);
}

complex_t delta_2d_2d_singular_integrand_outer(const complex_t alpha, const complex_t beta, void *args_){
    integrand_2d_2d_args *args=(integrand_2d_2d_args*)args_;
    args->alpha = alpha;
    args->beta = beta;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    complex_t ans=args->quadl.integral_2d(delta_2d_2d_singular_integrand_inner, args, triangle, flag);
    return ans;
}

complex_t delta_2d_2d_integrand_1(const complex_t alpha, const complex_t beta, void *args_){
    integrand_2d_2d_args *args=(integrand_2d_2d_args*)args_;
    basis_2d_t b_m=args->b_m;
    basis_2d_t b_n=args->b_n;
    vector_t<real_t> I_mm, I_mp, I_pm, I_pp;
    integrand_L3_2d_2d(real(alpha), real(beta), b_m, b_n, I_mm, I_mp, I_pm, I_pp);
    projection_2d_para para_m_m = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]);
    projection_2d_para para_m_p = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]);
    projection_2d_para para_p_m = prjection_2d(b_n.r_m, b_n.e[0], b_n.e[1], b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0]);
    projection_2d_para para_p_p = prjection_2d(b_n.r_p, b_n.e[1], b_n.e[0], b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0]);
    complex_t ans=0.0;  
    ans+= -1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1])*((b_n.r_m-para_m_m.p_0)^I_mm)/(2.0*b_n.A_m[0]);
    ans+= +1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1])*((b_n.r_p-para_m_p.p_0)^I_mp)/(2.0*b_n.A_p[0]);
    ans+= +1.0*(alpha*b_m.L_p[1]+beta*b_m.L_p[0])*((b_n.r_m-para_p_m.p_0)^I_pm)/(2.0*b_n.A_m[0]);
    ans+= -1.0*(alpha*b_m.L_p[1]+beta*b_m.L_p[0])*((b_n.r_p-para_p_p.p_0)^I_pp)/(2.0*b_n.A_p[0]);
    return b_m.L*b_n.L*ans/(4.0*pi);
}

complex_t delta_2d_2d(const basis_2d_t b_m, const basis_2d_t b_n, const complex_t k, const real_t lambda, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_2d_2d_args args={quadl, b_m, b_n};
    args.k = k;
    args.lambda = lambda;
    complex_t I1=0.0, I2=0.0;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    I1 = args.quadl.integral_2d(delta_2d_2d_singular_integrand_outer, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    I2 = args.quadl.integral_2d(delta_2d_2d_integrand_1, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    return I1+I2;
}

// 2d 3d

// 3d 1d
// 3d 2d
// 3d 3d