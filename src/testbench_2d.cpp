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

    // r_n_m = vector_t<real_t>(-0.3, +0.0, +0.0);
    // e_n_1 = vector_t<real_t>(+0.1, +0.0, -0.7);
    // e_n_2 = vector_t<real_t>(+0.0, +0.2, +0.0);
    // r_n_p = vector_t<real_t>(+0.6, +0.2, +0.8);

    b_m = basis_2d_t(r_m_m, e_m_1, e_m_2, r_m_p, 1, 1);
    b_n = basis_2d_t(r_n_m, e_n_1, e_n_2, r_n_p, 1, 1);

    T.set();
    // print(phi_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    // print(psi_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    // print(delta_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    print(Z_mn_2d_2d(b_m, b_n, k, eta, lambda, quadl, flag)); print(flag);
    T.unset();
    
}

void test_engine_2d_sphere_RCS(){

    // problem defintions
    const real_t GHz=1.0E+9;
    const real_t freq=0.75*GHz;
    const real_t clmax=1.0/11.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t radius=0.5;

    const size_t Ns=1001;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;
    range_t theta_s, phi_s;
    theta_s.set(deg2rad(-180.0), deg2rad(+180.0), Ns);
    phi_s.set(deg2rad(-180.0), deg2rad(+180.0), Ns);
    theta_s.linspace();
    phi_s.linspace();

    engine_t engine;
    create_sphere(radius);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, 0, 0);

    engine.compute_Z_mn();
    // engine.load_Z_mn("data/Z_mn.bin");
    file_t file;
    sigma_t sigma;

    // theta

    E_TM = +1.0;
    E_TE = +0.0;
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_solutions();
    engine.export_currents("data/currents.pos", 0);
    exit(0);

    file.open("data/RCS_1.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        sigma = engine.compute_RCS(theta_s(i), phi_i);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(phi_s(i)), sigma.theta, sigma.phi);
    }
    file.close();

    // phi

    E_TM = +0.0;
    E_TE = +1.0;
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_solutions();
    
    file.open("data/RCS_2.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        sigma = engine.compute_RCS(theta_s(i), phi_i);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), sigma.theta, sigma.phi);
    }
    file.close();

    //
    theta_s.unset();
    phi_s.unset();

    engine.unset();
    
}