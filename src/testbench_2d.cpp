//
#include "testbench_2d.hpp"


void test_engine_2d_2d(){

    const complex_t eta=sqrt(mu_0/eps_0);

    stopwatch_t T;
    quadl_domain_t quadl;   
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set_2d(k_max, tol);
    complex_t ans, psi, phi;
    int_t flag;
    vector_t<real_t> r_m_m, r_m_p, r_n_m, r_n_p, e_m_1, e_m_2, e_n_1, e_n_2;
    basis_2d_t b_m, b_n;

    real_t lambda=1.0;

    const complex_t k=2.0*pi/lambda;


    r_m_m = vector_t<real_t>(+0.3, -0.4, +0.3);
    e_m_1 = vector_t<real_t>(+0.1, +0.0, -0.7);
    e_m_2 = vector_t<real_t>(-0.3, +0.0, +0.0);
    r_m_p = vector_t<real_t>(+0.0, +0.2, +0.0);

    r_n_m = vector_t<real_t>(+0.3, -0.4, +0.3);
    e_n_1 = vector_t<real_t>(+0.1, +0.0, -0.7);
    e_n_2 = vector_t<real_t>(-0.3, +0.0, +0.0);
    r_n_p = vector_t<real_t>(+0.0, +0.2, +0.0);

    b_m = basis_2d_t(r_m_m, e_m_1, e_m_2, r_m_p, 1, 1);
    b_n = basis_2d_t(r_n_m, e_n_1, e_n_2, r_n_p, 1, 1);

    T.set();
    print(psi_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    print(phi_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    print(Z_mn_2d_2d(b_m, b_n, k, eta, lambda, quadl, flag)); print(flag);
    T.unset();
    
}