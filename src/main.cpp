//
#include "lib_basic.hpp"
#include "testbench_1d.hpp"
#include "testbench_2d.hpp"

int main(){

    // test_utilities();
    // test_gmsh();
    // test_shape();
    // test_engine_1d_1d();
    // test_engine_1d_vertical_dipole();
    test_engine_1d_debug();
    // test_engine_1d_vertical_dipole_input_adminttance();
    // test_engine_1d_loop_input_impedance();
    // test_engine_1d_vertical_dipole_mutual_impedance();
    // test_engine_1d_transmission_line_S_parameters();
    // test_engine_1d_RCS_vertical_wire();
    // test_engine_1d_far_field_transmission_line();
    // test_engine_1d_near_field_vertical_dipole();
    // test_engine_2d_transmission_line_near_field_1d();
    //
    // test_engine_2d_2d();
    // test_engine_2d_sphere_RCS();
    // test_engine_2d_sphere_near_field();
    // test_engine_2d_debug();
    // test_engine_2d_sheet_near_field();
    // test_engine_2d_box_near_field();
    // test_engine_2d_sphere_near_field_2d();

    // const real_t y=-1.0;
    // vector_t<real_t> r_m=vector_t<real_t>(-1.0, 0.0, +0.0);
    // vector_t<real_t> e_1=vector_t<real_t>(+0.0, y, +1.0);
    // vector_t<real_t> e_2=vector_t<real_t>(+0.0, y, -1.0);
    // vector_t<real_t> r_p=vector_t<real_t>(+1.0, 0.0, +0.0);
    // basis_2d_t b_n=basis_2d_t(r_m, e_1, e_2, r_p, 0, 0);
    
    // vector_t<real_t> l_n_m=0.5*(b_n.L_m[0]+b_n.L_m[1]);
    // l_n_m = unit(l_n_m);
    // real_t sign=l_n_m*b_n.n_p[0]<0.0 ? -1.0 : +1.0;
    // real_t Omega=2.0*(pi+sign*acos(b_n.n_m[0]*b_n.n_p[0]));
    // print(Omega/(4.0*pi));

    return 0;
}

