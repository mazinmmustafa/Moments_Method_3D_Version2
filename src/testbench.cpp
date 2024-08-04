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


void test_engine_1d(){

    const complex_t eta=sqrt(mu_0/eps_0);
    const complex_t k=2.0*pi;
    const real_t lambda=1.0;
    
    const real_t a=1.0E-4;
    const real_t L=1.0/31.0;

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

    r_n_m = vector_t<real_t>(+0.6*L-0.0*L, +0.0*L, +0.0*L);
    e_n = vector_t<real_t>(+0.2*L-0.0*L, +0.0*L, +0.0*L);
    r_n_p = vector_t<real_t>(+0.0*L-0.0*L, +0.3*L, +0.4*L);

    b_m = basis_1d_t(r_m_m, e_m, r_m_p, 1, 1);
    b_n = basis_1d_t(r_n_m, e_n, r_n_p, 1, 1);
    
    T.set();
    ans = Z_mn_1d(b_m, b_n, k, eta, lambda, a, quadl, flag);
    T.unset();
    print(ans);
    print(flag);
    
    
}

void test_engine_1d_vertical_dipole(){

    // units
    const real_t GHz=1.0E+9;
    const real_t mm=1.0E-3;
    // problem defintions
    const real_t freq=2.45*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t a=1.0E-4*mm;
    const real_t L=0.47*mm;

    engine_1d_t engine;
    engine.set(freq, mu_b, eps_b, clmax, mm, a);
    create_vertical_wire_dipole(L, 0.1*mm);
    
    

    engine.unset();
    
}
