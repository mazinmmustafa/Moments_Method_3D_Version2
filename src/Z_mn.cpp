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
    if (flag){print("psi\n");}
    phi = phi_1d_1d(b_m, b_n, k, lambda, a, quadl, flag);
    if (flag){print("phi\n");}
    return +j*k*eta*psi-j*(eta/k)*phi;
}

// 2d 2d
complex_t Z_mn_2d_2d(const basis_2d_t b_m, const basis_2d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_2d_2d_args args={quadl, b_m, b_n};
    args.k = k;
    args.eta = eta;
    args.lambda = lambda;
    const complex_t j=complex_t(0.0, 1.0);
    complex_t psi, phi;
    psi = psi_2d_2d(b_m, b_n, k, lambda, quadl, flag);
    if (flag){print("psi\n");}
    phi = phi_2d_2d(b_m, b_n, k, lambda, quadl, flag);
    if (flag){print("phi\n");}
    return +j*k*eta*psi-j*(eta/k)*phi;
}

// 3d 3d
complex_t Z_mn_3d_3d(const basis_3d_t b_m, const basis_3d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, const complex_t eps_b, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_3d_3d_args args={quadl, b_m, b_n};
    args.k = k;
    args.eta = eta;
    args.lambda = lambda;
    const complex_t j=complex_t(0.0, 1.0);
    complex_t psi, phi, nu, kappa;
    psi = psi_3d_3d(b_m, b_n, k, lambda, quadl, flag);
    if (flag){print("psi\n");}
    phi = phi_3d_3d(b_m, b_n, k, lambda, quadl, flag);
    if (flag){print("phi\n");}
    kappa = kappa_3d_3d(b_m, b_n, k, lambda, quadl, flag);
    if (flag){print("kappa\n");}
    nu = nu_3d_3d(b_m, b_n, lambda, eps_b);
    return +j*k*eta*psi-j*(eta/k)*phi-j*(eta/k)*nu-j*(eta/k)*kappa/eps_b;
}