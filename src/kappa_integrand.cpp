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
    vector_t<real_t> dist=b_n.e[0]+(b_n.e[1]-b_n.e[0])/3.0+(b_n.e[2]-b_n.e[0])/3.0;
    R_m = mag(b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]+real(gamma)*b_m.L_m[2]-dist);
    R_p = mag(b_m.r_p+real(alpha)*b_m.L_p[2]+real(beta)*b_m.L_p[1]+real(gamma)*b_m.L_p[0]-dist);
    real_t factor=b_n.nA*(b_n.r_m+b_n.L_m[0]/3.0+b_n.L_m[1]/3.0+b_n.L_m[2]/3.0);
    I_m = +1.0*(exp(-j*k*R_m)/R_m)*factor;
    I_p = -1.0*(exp(-j*k*R_p)/R_p)*factor;
    return 36.0*b_m.A*b_n.A*(I_m/b_n.V_m+I_p/b_n.V_p)/(4.0*pi);
}

complex_t kappa_3d_3d(const basis_3d_t b_m, const basis_3d_t b_n, const complex_t k, const real_t lambda, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_3d_3d_args args={quadl, b_m, b_n};
    args.k = k;
    args.lambda = lambda;
    complex_t I1=0.0;
    tetrahedron_domain_t tetrahedron={vector_t<real_t>(0.0, 0.0, 0.0), 
        vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0), 
        vector_t<real_t>(0.0, 0.0, 1.0)};
    I1 = args.quadl.integral_3d(kappa_3d_3d_singular_integrand_outer, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    complex_t kappa_m=(b_n.eps_m-1.0)/b_n.eps_m;
    complex_t kappa_p=(b_n.eps_p-1.0)/b_n.eps_p;
    return (kappa_p-kappa_m)*(I1);
}
