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
