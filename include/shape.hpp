#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
//

// Definitions

struct edge_t{
    vector_t<real_t> v[2];
    int_t physical_group=-1;
    real_t length=0.0;
    edge_t(){}
    edge_t(const vector_t<real_t> v1, const vector_t<real_t> v2, const int_t physical_group){
        this->v[0] = v1;
        this->v[1] = v2;
        get_length();
        this->physical_group = physical_group;
    }   
    void get_length(){
        this->length = mag(v[0]-v[1]);
        assert_error(this->length>0.0, "invalid edge");
    }
};

struct triangle_t{
    vector_t<real_t> v[3];
    int_t physical_group=-1;
    real_t area=0.0;
    vector_t<real_t> n;
    triangle_t(){}
    triangle_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const int_t physical_group){
        this->v[0] = v1;
        this->v[1] = v2;
        this->v[2] = v3;
        get_area();
        this->physical_group = physical_group;
    }
    void get_area(){
        vector_t<real_t> vector=(v[1]-v[0])^(v[2]-v[0]);
        this->area = mag(vector)/2.0;
        assert_error(this->area>0.0, "invalid triangle");
        this->n = unit(vector);
    }
};

struct tetrahedron_t{
    vector_t<real_t> v[4];
    int_t physical_group=-1;
    real_t volume=0.0;
    complex_t eps=complex_t(+1.0, -0.0);
    complex_t mu=complex_t(+1.0, -0.0);
    tetrahedron_t(){}
    tetrahedron_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const vector_t<real_t> v4, const int_t physical_group){
        this->v[0] = v1;
        this->v[1] = v2;
        this->v[2] = v3;
        this->v[3] = v4;
        get_volume();
        this->physical_group = physical_group;
    }
    void get_volume(){
        this->volume = ((v[1]-v[0])^(v[2]-v[0]))*(v[3]-v[0])/6.0;
        assert_error(this->volume>0.0, "invalid tetrahedron");
    }
};

struct basis_1d_t{
    vector_t<real_t> r_m, e, r_p;
    vector_t<real_t> L_m, L_p;
    real_t l_m, l_p;
    int_t pg_m, pg_p;
    basis_1d_t(){}
    basis_1d_t(const vector_t<real_t> r_m, const vector_t<real_t> e_1, const vector_t<real_t> r_p, 
        const int_t pg_m, const int_t pg_p){
        this->r_m = r_m;
        this->e = e_1;
        this->r_p = r_p;
        this->pg_m = pg_m;
        this->pg_p = pg_p;
        basis_1d_t::get_parameters();
    }
    void get_parameters(){
        this->L_m = +1.0*(this->e-this->r_m);
        this->L_p = +1.0*(this->e-this->r_p);
        this->l_m = mag(this->L_m);
        this->l_p = mag(this->L_p);
        assert_error(this->l_m>0.0&&this->l_p>0.0, "invalid 1d basis");
    }
};

struct basis_2d_t{
    basis_1d_t b_1, b_2;
    real_t e;
    real_t A_m, A_p;
    vector_t<real_t> n_m, n_p;
    int_t pg_m, pg_p;
    basis_2d_t(){}
    basis_2d_t(const vector_t<real_t> r_m, const vector_t<real_t> e_1, const vector_t<real_t> e_2, const vector_t<real_t> r_p,
        const int_t pg_m, const int_t pg_p){
        this->b_1 = basis_1d_t(r_m, e_1, r_p, pg_m, pg_p);
        this->b_2 = basis_1d_t(r_m, e_2, r_p, pg_m, pg_p);
        this->pg_m = pg_m;
        this->pg_p = pg_p;
        basis_2d_t::get_parameters();
    }
    void get_parameters(){
        this->e = mag(this->b_1.e-this->b_2.e);
        vector_t<real_t> vector_m=(this->b_1.e-this->b_1.r_m)^(this->b_2.e-this->b_2.r_m);
        vector_t<real_t> vector_p=(this->b_2.e-this->b_2.r_p)^(this->b_1.e-this->b_1.r_p);
        this->A_m = mag(vector_m)/2.0;
        this->A_p = mag(vector_p)/2.0;
        assert_error(this->A_m>0.0&&this->A_p>0.0, "invalid 2d basis");
        this->n_m = unit(vector_m);
        this->n_p = unit(vector_p);
    }
};

struct basis_3d_t{
    basis_2d_t b_1, b_2, b_3;
    real_t a;
    real_t V_m, V_p;
    vector_t<real_t> n;
    int_t pg_m, pg_p;
    basis_3d_t(){}
    basis_3d_t(const vector_t<real_t> r_m, const vector_t<real_t> e_1, const vector_t<real_t> e_2, 
        const vector_t<real_t> e_3, const vector_t<real_t> r_p, const int_t pg_m, const int_t pg_p){
        this->b_1 = basis_2d_t(r_m, e_2, e_1, r_p, pg_m, pg_p);
        this->b_2 = basis_2d_t(r_m, e_3, e_2, r_p, pg_m, pg_p);
        this->b_3 = basis_2d_t(r_m, e_1, e_3, r_p, pg_m, pg_p);
        this->pg_m = pg_m;
        this->pg_p = pg_p;
        basis_3d_t::get_parameters();
    }
    void get_parameters(){
        vector_t<real_t> vector=(this->b_1.b_1.e-this->b_1.b_2.e)^(this->b_2.b_1.e-this->b_1.b_2.e);
        this->a = mag(vector)/2.0;
        this->n = unit(vector);
        this->V_m = ((this->b_3.b_2.e-this->b_3.b_1.e)^(this->b_1.b_1.e-this->b_1.b_2.e))*(-1.0*this->b_1.b_2.L_m)/6.0;
        this->V_p = ((this->b_1.b_1.e-this->b_1.b_2.e)^(this->b_3.b_2.e-this->b_3.b_1.e))*(-1.0*this->b_1.b_2.L_p)/6.0;
        assert_error(this->V_m>0.0&&this->V_p>0.0, "invalid 3d basis");
    }
};

class shape_t{
    private:
        int_t is_physical_specified=false;
        int_t is_basis_allocated=false;
        real_t frequency=0.0, lambda=0.0;
        complex_t mu_b=1.0, eps_b=1.0;
        size_t N_0d=0, N_1d=0, N_2d=0, N_3d=0;
        size_t N_basis_1d=0, N_basis_2d=0, N_basis_3d=0;
        basis_1d_t *basis_1d_list=null;
        basis_2d_t *basis_2d_list=null;
        basis_3d_t *basis_3d_list=null;
        const real_t mesh_tol=1.0E-8;
        void load_mesh(const real_t metric_unit);
        void load_basis_functions();
    public:
        shape_t(const real_t frequency, const complex_t mu_b, const complex_t eps_b){
            assert_error(abs(mu_b)>0.0 && abs(eps_b)>0.0, "invalid background medium parameters");
            assert_error(frequency>0.0, "invalid frequency");
            this->frequency = frequency;
            this->lambda = frequency/c_0;
        }
        ~shape_t(){}
        void get_basis_functions(const real_t clmax, const real_t metric_unit);
        void clear();
};

// Functions
void call_gmsh(const real_t tol);
void create_vertical_wire_dipole(const real_t length, const real_t port_length);
void create_sphere(const real_t radius);
void create_patch_antenna();

#endif
