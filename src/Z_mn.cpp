//
#include "Z_mn.hpp"

// 1d 1d
complex_t Z_mn_1d_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, const real_t a, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_1d_1d_args args={quadl, b_m, b_n};
    args.a = a;
    args.k = k;
    args.eta = eta;
    args.lambda = lambda;
    const complex_t j=complex_t(0.0, 1.0);
    complex_t psi, phi;
    psi = psi_1d_1d(b_m, b_n, k, lambda, a, quadl, flag);
    phi = phi_1d_1d(b_m, b_n, k, lambda, a, quadl, flag);
    return +j*k*eta*psi-j*(eta/k)*phi;
}