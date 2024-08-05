//
#include "testbench.hpp"

void test_utilities(){

    stopwatch_t T;
    T.set();

    const real_t R_max=0.1;
    const size_t Ns=101;
    range_t range;
    range.set(-R_max, +R_max, Ns);
    range.linspace();

    file_t file;
    file.open("data/test.dat", 'w');

    for (size_t i=0; i<Ns; i++){
        file.write("%21.14E %21.14E %21.14E\n", 
            range(i), sinc(range(i)), sinc_dx(range(i)));
    }

    file.close();

    range.unset();
    T.unset();

    file.open("data/test.dat", 'r');

    real_t x[3]={0.0, 0.0, 0.0};
    for (size_t i=0; i<Ns; i++){
        file.read("%lf %lf %lf\n", &x[0], &x[1], &x[2]);
        print("%21.14E %21.14E %21.14E\n", x[0], x[1], x[2]);
    }

    file.close();

}

void test_gmsh(){

    // create_vertical_wire_dipole(0.47, 0.1);
    create_sphere(0.5);
    call_gmsh(0.2);

}

void test_shape(){

    const real_t GHz=1.0E+9;
    // const real_t mm=1.0E-3;

    const real_t freq=2.45*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/21.0;

    shape_t shape(freq, 1.0, 1.0);
    // create_vertical_wire_dipole(0.47*lambda, 0.1*lambda);
    // create_sphere(100);
    create_patch_antenna();
    shape.get_basis_functions(clmax, 1.0E-3);

    shape.clear();

}


void test_engine_1d_1d(){

    const complex_t eta=sqrt(mu_0/eps_0);
    const complex_t k=2.0*pi;
    const real_t lambda=1.0;
    
    const real_t a=1.0E-4;
    const real_t L=1.0/11.0;

    stopwatch_t T;
    quadl_domain_t quadl;   
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set_1d(k_max, tol);
    complex_t ans, psi, phi;
    int_t flag;
    vector_t<real_t> r_m_m, r_m_p, r_n_m, r_n_p, e_m, e_n;
    basis_1d_t b_m, b_n;

    r_m_m = vector_t<real_t>(+0.0*L, +0.0*L, +0.0*L);
    e_m = vector_t<real_t>(+0.2*L, +0.0*L, +0.0*L);
    r_m_p = vector_t<real_t>(+0.6*L, +0.0*L, +0.0*L);

    r_n_m = vector_t<real_t>(+0.0*L, +0.0*L, +0.0*L);
    e_n = vector_t<real_t>(+0.2*L, +0.0*L, +0.0*L);
    r_n_p = vector_t<real_t>(+0.6*L, +0.0*L, +0.0*L);

    b_m = basis_1d_t(r_m_m, e_m, r_m_p, 1, 1);
    b_n = basis_1d_t(r_n_m, e_n, r_n_p, 1, 1);
    
    T.set();
    ans = Z_mn_1d_1d(b_m, b_n, k, eta, lambda, a, quadl, flag);
    T.unset();
    print(ans);
    print(flag);
    
}

void test_engine_1d_vertical_dipole(){

    // problem defintions
    const real_t freq=c_0;
    const real_t clmax=1.0/31.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t a=5.0E-3;
    const real_t L=0.47;
    const size_t N_ports=1;
    const complex_t V1=1.0;
    const complex_t Z1=50.0;
    const int_t pg1=1; 
    const vector_t<real_t> p1=vector_t<real_t>(+0.0, +0.0, +1.0);
    const real_t port_length=1.0*clmax;

    engine_t engine;
    create_vertical_wire_dipole(L, port_length, clmax);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, a, N_ports);
    engine.assign_port(0, V1, Z1, pg1, p1, port_length, 0.0);

    engine.compute_Z_mn();
    engine.compute_V_m_ports();
    engine.compute_I_n();
    engine.export_solutions();

    print(engine.compute_Z_in(0));
    
    engine.unset();
    
}

void test_engine_1d_vertical_dipole_input_adminttance(){

    // problem defintions
    const real_t freq=c_0;
    const real_t clmax=1.0/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t a=1.0E-3;
    const size_t N_ports=1;
    const complex_t V1=1.0;
    const complex_t Z1=50.0;
    const int_t pg1=1; 
    const vector_t<real_t> p1=vector_t<real_t>(+0.0, +0.0, +1.0);

    const size_t Ns=801;
    const real_t L_min=0.01;
    const real_t L_max=4.0;
    const real_t port_length=L_min/2.0;
    range_t L;
    L.set(L_min, L_max, Ns);
    L.linspace();

    engine_t engine;
    file_t file;
    file.open("data/Y_in.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        print("\nstep %zu:\n", i);
        create_vertical_wire_dipole(L(i), port_length, clmax);
        engine.set(freq, mu_b, eps_b, clmax, 1.0, a, N_ports);
        engine.assign_port(0, V1, Z1, pg1, p1, port_length, 0.0);

        engine.compute_Z_mn();
        engine.compute_V_m_ports();
        engine.compute_I_n();
        complex_t Y_in=1.0/engine.compute_Z_in(0);
        file.write("%21.14E %21.14E %21.14E\n", L(i), real(Y_in), imag(Y_in));
        engine.unset();
    }
    file.close();

    L.unset();
    
}