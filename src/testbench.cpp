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


