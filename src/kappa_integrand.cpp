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
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_m, R_p;
    vector_t<real_t> dist=b_n.e[0]+real(alpha_)*b_n.e[1]+real(beta_)*b_n.e[2];
    R_m = mag(b_m.r_m+real(alpha)*b_m.e[0]+real(beta)*b_m.e[1]+real(gamma)*b_m.e[2]-dist);
    R_p = mag(b_m.r_p+real(alpha)*b_m.e[2]+real(beta)*b_m.e[1]+real(gamma)*b_m.e[0]-dist);
    real_t factor=b_n.nA*(real(alpha_)*b_n.e[0]+real(beta_)*b_n.e[1]+real(1.0-alpha_-beta_)*b_n.e[2]);
    I_m = +1.0*(-j*k*exp(-j*k*R_m/2.0))*sinc(k*R_m/2.0)*factor;
    I_p = -1.0*(-j*k*exp(-j*k*R_p/2.0))*sinc(k*R_p/2.0)*factor;
    return 4.0*b_m.A*b_n.A*b_n.A*(I_m/b_n.V_m+I_p/b_n.V_p)/(4.0*pi);
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
    vector_t<real_t> I_m, I_p;
    basis_2d_t b_n_=basis_2d_t(b_n.e[0], b_n.e[1], b_n.e[2], b_n.e[0], b_n.pg_m, b_n.pg_p);
    integrand_L2_3d_2d(real(alpha), real(beta), real(gamma), b_m, b_n_, I_m, I_p);
    complex_t ans=0.0;  
    ans+= +1.0*b_n.nA*I_m/b_n.V_m;
    ans+= -1.0*b_n.nA*I_p/b_n.V_p;
    return 2.0*b_m.A*b_n.A*ans/(4.0*pi);
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
    I2 = args.quadl.integral_3d(kappa_3d_3d_integrand_1, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    complex_t kappa_m=(b_n.eps_m-1.0)/b_n.eps_m;
    complex_t kappa_p=(b_n.eps_p-1.0)/b_n.eps_p;
    return (kappa_p-kappa_m)*(I1+I2);
}
