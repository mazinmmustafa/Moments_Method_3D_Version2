//
#include "testbench_1d.hpp"

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

    engine.export_currents("data/currents.txt");

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

    const size_t Ns=401;
    const real_t L_min=0.05;
    const real_t L_max=4.0;
    const real_t port_length=clmax;
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

void test_engine_1d_loop_input_impedance(){

    // problem defintions
    const real_t freq=c_0;
    const real_t clmax=1.0/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t a=1.0E-4;
    const size_t N_ports=1;
    const complex_t V1=1.0;
    const complex_t Z1=50.0;
    const int_t pg1=1; 
    const vector_t<real_t> p1=vector_t<real_t>(+0.0, +1.0, +0.0);

    const size_t Ns=401;
    const real_t S_min=0.05;
    const real_t S_max=4.0;
    const real_t port_length=clmax;
    range_t S;
    S.set(S_min, S_max, Ns);
    S.linspace();

    engine_t engine;
    file_t file;
    file.open("data/Z_in.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        print("\nstep %zu:\n", i);
        create_loop(S(i)/(2.0*pi), clmax);
        engine.set(freq, mu_b, eps_b, clmax, 1.0, a, N_ports);
        engine.assign_port(0, V1, Z1, pg1, p1, port_length, 0.0);

        engine.compute_Z_mn();
        engine.compute_V_m_ports();
        engine.compute_I_n();
        complex_t Z_in=engine.compute_Z_in(0);
        file.write("%21.14E %21.14E %21.14E\n", S(i), real(Z_in), imag(Z_in));
        engine.unset();
    }
    file.close();

    S.unset();
    
}

void test_engine_1d_vertical_dipole_mutual_impedance(){

    // problem defintions
    const real_t freq=c_0;
    const real_t clmax=1.0/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t L=0.5;
    const real_t a=1.0E-5;
    const size_t N_ports=2;
    const int_t pg1=1; 
    const int_t pg2=2; 
    const vector_t<real_t> p1=vector_t<real_t>(+0.0, +0.0, +1.0);
    const vector_t<real_t> p2=vector_t<real_t>(+0.0, +0.0, +1.0);
    const complex_t Z_0=50.0;

    const size_t Ns=401;
    const real_t d_min=0.01;
    const real_t d_max=4.0;
    const real_t port_length=clmax;
    range_t d;
    d.set(d_min, d_max, Ns);
    d.linspace();

    engine_t engine;
    file_t file;
    file.open("data/S_matrix.txt", 'w');
    matrix_t<complex_t> S_matrix;
    S_matrix.set(N_ports, N_ports);
    for (size_t i=0; i<Ns; i++){
        print("\nstep %zu:\n", i);
        create_two_vertical_wire_dipole(L, port_length, d(i), clmax);
        engine.set(freq, mu_b, eps_b, clmax, 1.0, a, N_ports);
        engine.assign_port(0, +0.0, +0.0, pg1, p1, port_length, 0.0);
        engine.assign_port(1, +0.0, +0.0, pg2, p2, port_length, 0.0);

        engine.compute_Z_mn();
        engine.compute_S_matrix(S_matrix, Z_0);
        file.write("%21.14E ", d(i));
        for (size_t m=0; m<N_ports; m++){
            for (size_t n=0; n<N_ports; n++){
                file.write("%21.14E %21.14E ", real(S_matrix(m, n)), imag(S_matrix(m, n)));
            }
        }
        file.write("\n");
        engine.unset();
    }
    S_matrix.unset();
    file.close();

    d.unset();
    
}

void test_engine_1d_transmission_line_S_parameters(){

    const real_t MHz=1.0E+6;
    const real_t cm=1.0E-2;
    const real_t mm=1.0E-3;

    // problem defintions
    real_t clmax, lambda;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t L=250.0*cm;
    const real_t S=3.0*cm;
    const real_t a=1.0*mm;
    const size_t N_ports=2;
    const int_t pg1=1; 
    const int_t pg2=2; 
    const vector_t<real_t> p1=vector_t<real_t>(+0.0, +1.0, +0.0);
    const vector_t<real_t> p2=vector_t<real_t>(+0.0, +1.0, +0.0);
    const complex_t Z_0=50.0;

    const size_t Ns=401;
    const real_t freq_min=1.0*MHz;
    const real_t freq_max=200.0*MHz;
    range_t freq;
    freq.set(freq_min, freq_max, Ns);
    freq.linspace();

    engine_t engine;
    file_t file;
    file.open("data/S_matrix.txt", 'w');
    matrix_t<complex_t> S_matrix;
    S_matrix.set(N_ports, N_ports);
    for (size_t i=0; i<Ns; i++){
        print("\nstep %zu:\n", i);
        lambda = c_0/freq(i);
        clmax = lambda/21.0;
        create_transmission_line(L, S, clmax);
        engine.set(freq(i), mu_b, eps_b, clmax, 1.0, a, N_ports);
        engine.assign_port(0, +0.0, Z_0, pg1, p1, 0.0, 0.0);
        engine.assign_port(1, +0.0, Z_0, pg2, p2, 0.0, 0.0);

        engine.compute_Z_mn();
        engine.compute_S_matrix(S_matrix, Z_0);
        file.write("%21.14E ", freq(i));
        for (size_t m=0; m<N_ports; m++){
            for (size_t n=0; n<N_ports; n++){
                file.write("%21.14E %21.14E ", real(S_matrix(m, n)), imag(S_matrix(m, n)));
            }
        }
        file.write("\n");
        engine.unset(); 

    }
    S_matrix.unset();
    file.close();

    freq.unset();
    
}

void test_engine_1d_RCS_vertical_wire(){

    // problem defintions
    const real_t freq=c_0;
    const real_t clmax=1.0/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const complex_t E_TM=1.0;
    const complex_t E_TE=0.0;
    
    const real_t theta_i=deg2rad(30.0);
    const real_t phi_i=deg2rad(0.0);
    const real_t theta_s=deg2rad(30.0);
    const real_t phi_s=deg2rad(0.0);

    const size_t Ns=201;
    const real_t L_min=0.05;
    const real_t L_max=4.0;
    range_t L;
    L.set(L_min, L_max, Ns);
    L.linspace();

    engine_t engine;
    sigma_t sigma;
    file_t file;
    file.open("data/RCS.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        print("\nstep %zu:\n", i);
        const real_t a=L(i)/200.0;
        create_vertical_wire(L(i), clmax);

        engine.set(freq, mu_b, eps_b, clmax, 1.0, a, 0);
        engine.compute_Z_mn();
        engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
        engine.compute_I_n();

        sigma = engine.compute_RCS(theta_s, phi_s);
        file.write("%21.14E %21.14E %21.14E\n", L(i), sigma.theta, sigma.phi);
        engine.unset();
    }
    file.close();
    
    L.unset();
}

void test_engine_1d_far_field_transmission_line(){

    // problem defintions
    const real_t cm=1.0E-2;
    const real_t mm=1.0E-3;
    const real_t GHz=1.0E+9;

    const real_t freq=1.0*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t L=150.0*cm;
    const real_t h=25.0*cm;
    const real_t a=0.8*mm;

    const size_t N_ports=2;
    const int_t pg1=1; 
    const int_t pg2=2; 
    const vector_t<real_t> p1=vector_t<real_t>(+0.0, +1.0, +0.0);
    const vector_t<real_t> p2=vector_t<real_t>(+0.0, +1.0, +0.0);
    const complex_t Z_0=50.0;

    const real_t theta=deg2rad(90.0);
    const size_t Ns=401;
    const real_t phi_min=deg2rad(-180.0);
    const real_t phi_max=deg2rad(+180.0);
    range_t phi;
    phi.set(phi_min, phi_max, Ns);
    phi.linspace();

    engine_t engine;
    create_transmission_line(L, 2.0*h, clmax);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, a, N_ports);
    engine.assign_port(0, +2.0, 2.0*Z_0, pg1, p1, 0.0, 0.0);
    engine.assign_port(1, +0.0, 0.0, pg2, p2, 0.0, 0.0);

    engine.compute_Z_mn();
    engine.compute_V_m_ports();
    engine.compute_I_n();
    engine.export_solutions();
    engine.export_currents("data/currents.txt");

    far_field_t far_field;
    file_t file;
    file.open("data/far_field.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        far_field = engine.compute_far_field(theta, phi(i));
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(phi(i)), abs(far_field.theta), abs(far_field.phi));
    }
    engine.unset();
    file.close();
    
    phi.unset();
}


struct integrand_args_1d{
    line_t line;
    vector_t<real_t> r, unit_vector;
    real_t a;
};

#define EPS 1.0E-20

complex_t integrand_1_1d(const complex_t alpha, void *args_){
    assert(args_!=null);
    integrand_args_1d* args=(integrand_args_1d*)args_;
    vector_t<real_t> rho=args->line.v[0]
                         +real(alpha)*(args->line.v[1]-args->line.v[0]);
    real_t R=mag(args->r-rho)+EPS;
    R = sqrt(R*R+args->a*args->a);
    return args->line.length/R;
}

complex_t integrand_2_1d(const complex_t alpha, void *args_){
    assert(args_!=null);
    integrand_args_1d* args=(integrand_args_1d*)args_;
    vector_t<real_t> rho=args->line.v[0]
                         +real(alpha)*(args->line.v[1]-args->line.v[0]);
    real_t R=mag(args->r-rho)+EPS;
    R = sqrt(R*R+args->a*args->a);
    return args->line.length*(rho-args->r)*args->unit_vector/R;
}

complex_t integrand_3_1d(const complex_t alpha, void *args_){
    assert(args_!=null);
    integrand_args_1d* args=(integrand_args_1d*)args_;
    vector_t<real_t> rho=args->line.v[0]
                         +real(alpha)*(args->line.v[1]-args->line.v[0]);
    real_t R=mag(args->r-rho)+EPS;
    R = sqrt(R*R+args->a*args->a);
    return args->line.length*(1.0/pow(R, 2.0))*(args->unit_vector*unit(args->r-rho));
}

void test_engine_1d_debug(){

    //
    vector_t<real_t> v1, v2;
    v1 = vector_t<real_t>(-0.5, -0.6, -0.2);
    v2 = vector_t<real_t>(+0.4, -0.6, +0.0);

    const real_t a=1.0E-3;
    line_t line=line_t(v1, v2, 0);
    vector_t<real_t> p=vector_t<real_t>(+0.2, +0.1, +1.0);

    projection_1d_para para;
    para = prjection_1d(v1, v2, p);

    //
    real_t I1=0.0;
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    I1 = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    print(I1);

    vector_t<real_t> I2 = vector_t<real_t>(0.0, 0.0, 0.0);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    I2 = (para.P_p-para.P_m)*para.l_unit+log((para.P_p+para.l_p)/(para.P_m+para.l_m))*para.P_0*para.P_0_unit;
    print(I2);

    vector_t<real_t> I3 = vector_t<real_t>(0.0, 0.0, 0.0);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    I3 = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
         -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    print(I3);

    vector_t<real_t> x=vector_t<real_t>(1.0, 0.0, 0.0);
    vector_t<real_t> y=vector_t<real_t>(0.0, 1.0, 0.0);
    vector_t<real_t> z=vector_t<real_t>(0.0, 0.0, 1.0);

    int_t flag;
    quadl_domain_t quadl;
    quadl.set_1d(15, 1.0E-10);
    integrand_args_1d args={line, p, p, a};
    line_domain_t line_domain={vector_t<real_t>(0.0, 0.0, 0.0), 
                                vector_t<real_t>(1.0, 0.0, 0.0)};
    
    I1 =real(quadl.integral_1d(integrand_1_1d, &args, line_domain, flag)); assert(!flag);
    print(I1);

    args.unit_vector = x;
    I2.x =real(quadl.integral_1d(integrand_2_1d, &args, line_domain, flag)); assert(!flag);
    args.unit_vector = y;
    I2.y =real(quadl.integral_1d(integrand_2_1d, &args, line_domain, flag)); assert(!flag);
    args.unit_vector = z;
    I2.z =real(quadl.integral_1d(integrand_2_1d, &args, line_domain, flag)); assert(!flag);
    print(I2);

    args.unit_vector = x;
    I3.x =real(quadl.integral_1d(integrand_3_1d, &args, line_domain, flag)); assert(!flag);
    args.unit_vector = y;
    I3.y =real(quadl.integral_1d(integrand_3_1d, &args, line_domain, flag)); assert(!flag);
    args.unit_vector = z;
    I3.z =real(quadl.integral_1d(integrand_3_1d, &args, line_domain, flag)); assert(!flag);
    print(I3);

}


void test_engine_1d_near_field_vertical_dipole(){

    // problem defintions
    const real_t mm=1.0E-3;
    const real_t GHz=1.0E+9;

    const real_t freq=1.0*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t L=150.0*mm;
    const real_t a=1.0*mm;

    const size_t N_ports=1;
    const int_t pg=1; 
    const vector_t<real_t> p=vector_t<real_t>(+0.0, +0.0, +1.0);
    const complex_t Z_0=50.0;

    const size_t Ns=401;
    const real_t z_min=-400.0*mm;
    const real_t z_max=+400.0*mm;
    const real_t x=+200.0*mm;
    const real_t y=+200.0*mm;

    range_t z;
    z.set(z_min, z_max, Ns);
    z.linspace();

    engine_t engine;
    create_vertical_wire_dipole(L, clmax, clmax);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, a, N_ports);
    engine.assign_port(0, 1.0, Z_0, pg, p, 0.0, 0.0);

    engine.compute_Z_mn();
    engine.compute_V_m_ports();
    engine.compute_I_n();
    engine.export_solutions();

    vector_t<real_t> r;
    vector_t<complex_t> E, H;
    file_t file;
    file.open("data/near_field.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        r = vector_t<real_t>(x, y, z(i));
        E = engine.compute_near_field_E(r);
        H = engine.compute_near_field_H(r);
        file.write("%21.14E ", z(i));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(E.x), imag(E.x),
            real(E.y), imag(E.y),
            real(E.z), imag(E.z));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(H.x), imag(H.x),
            real(H.y), imag(H.y),
            real(H.z), imag(H.z));
        file.write("\n");
    }
    engine.unset();
    file.close();
    
    z.unset();
}

void test_engine_2d_transmission_line_near_field_1d(){

    // problem defintions
    const real_t cm=1.0E-2;
    const real_t mm=1.0E-3;
    const real_t GHz=1.0E+9;

    const real_t freq=1.0*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/21.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t L=150.0*cm;
    const real_t h=25.0*cm;
    const real_t a=0.8*mm;

    const size_t N_ports=2;
    const int_t pg1=1; 
    const int_t pg2=2; 
    const vector_t<real_t> p1=vector_t<real_t>(+0.0, +1.0, +0.0);
    const vector_t<real_t> p2=vector_t<real_t>(+0.0, +1.0, +0.0);
    const complex_t Z_0=50.0;

    const size_t Ns_x=101, Ns_y=101;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;

    engine_t engine;
    create_transmission_line(L, 2.0*h, clmax);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, a, N_ports);
    engine.assign_port(0, +2.0, 2.0*Z_0, pg1, p1, 0.0, 0.0);
    engine.assign_port(1, +0.0, 0.0, pg2, p2, 0.0, 0.0);

    engine.compute_Z_mn();
    // engine.load_Z_mn("data/Z_mn.bin");
    file_t file_x, file_y, file_data;
    range_t x, y;
    real_t y_min, y_max, x_min, x_max, z;
    const real_t range=2.0;

    x_min = -range;
    x_max = +range;
    z = +0.0;
    y_min = -range;
    y_max = +range;
    x.set(x_min, x_max, Ns_x);
    y.set(y_min, y_max, Ns_y);
    x.linspace();
    y.linspace();

    engine.compute_V_m_ports();
    engine.compute_I_n();
    vector_t<real_t> r;

    vector_t<complex_t> E, H;
    incident_field_t incident_field;
    file_x.open("data/near_field_2d_x.txt", 'w');
    file_y.open("data/near_field_2d_y.txt", 'w');
    file_data.open("data/near_field_2d_data.txt", 'w');

    size_t counter=0;
    for (size_t i=0; i<Ns_x; i++){
        for (size_t j=0; j<Ns_y; j++){
            progress_bar(counter, Ns_x*Ns_y, "computing near fields...");
            counter++;
            r = vector_t<real_t>(x(i), y(j), z);
            incident_field = compute_incident_field(E_TM, E_TE, theta_i, phi_i, 
                                                    2.0*pi/lambda, sqrt(mu_0/eps_0), r);
            E = engine.compute_near_field_E(r);
            H = engine.compute_near_field_H(r);
            file_x.write("%21.14E ", x(i));
            file_y.write("%21.14E ", y(j));
            vector_t<complex_t> E_field=vector_t<complex_t>(E.x, E.y, E.z);
            file_data.write("%21.14E ", sqrt(pow(real(E.x), 2.0)
                                            +pow(real(E.y), 2.0)
                                            +pow(real(E.z), 2.0)));
        }
        file_x.write("\n");
        file_y.write("\n");
        file_data.write("\n");
    }
    file_x.close();
    file_y.close();
    file_data.close();

    //
    x.unset();
    y.unset();

    engine.unset();
    
}