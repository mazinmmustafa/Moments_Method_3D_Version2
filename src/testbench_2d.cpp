//
#include "testbench_2d.hpp"


void test_engine_2d_2d(){

    const complex_t eta=sqrt(mu_0/eps_0);

    stopwatch_t T;
    quadl_domain_t quadl;   
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set_1d(k_max, tol);
    complex_t ans, psi, phi;
    int_t flag;
    vector_t<real_t> r_m_m, r_m_p, r_n_m, r_n_p, e_m, e_n;
    basis_1d_t b_m, b_n;

    real_t lambda=c_0/1.0E8;
    // real_t lambda=1.0;

    const complex_t k=2.0*pi/lambda;

    const real_t a=1.0E-3;

    r_m_m = vector_t<real_t>(+0.5, +0.5, +0.0);
    e_m   = vector_t<real_t>(+0.0, +0.0, +0.0);
    r_m_p = vector_t<real_t>(+0.0, -0.5, +0.0);

    r_n_m = vector_t<real_t>(-0.5, -0.5, +0.0);
    e_n   = vector_t<real_t>(+0.0, -0.5, +0.0);
    r_n_p = vector_t<real_t>(-0.0, +0.0, +0.0);

    // r_m_m = vector_t<real_t>(+0.20, -1.6, 0.0);
    // e_m   = vector_t<real_t>(+0.40, -1.2, 0.0);
    // r_m_p = vector_t<real_t>(+0.60, -0.0, 0.0);

    // r_n_m = vector_t<real_t>(+0.20, -1.6, 0.0);
    // e_n   = vector_t<real_t>(+0.40, -1.2, 0.0);
    // r_n_p = vector_t<real_t>(+0.60, -0.0, 0.0);

    b_m = basis_1d_t(r_m_m, e_m, r_m_p, 1, 1);
    b_n = basis_1d_t(r_n_m, e_n, r_n_p, 1, 1);

    T.set();
    print(psi_1d_1d(b_m, b_n, k, lambda, a, quadl, flag)); print(flag);
    print(phi_1d_1d(b_m, b_n, k, lambda, a, quadl, flag)); print(flag);
    print(Z_mn_1d_1d(b_m, b_n, k, eta, lambda, a, quadl, flag)); print(flag);
    T.unset();
    
}