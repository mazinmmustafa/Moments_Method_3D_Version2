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

struct line_t{
    vector_t<real_t> v[2];
    int_t physical_group=-1;
    real_t length=0.0;
    line_t(){}
    line_t(const vector_t<real_t> v1, const vector_t<real_t> v2, const int_t physical_group){
        this->v[0] = v1;
        this->v[1] = v2;
        get_length();
        this->physical_group = physical_group;
    }   
    void get_length(){
        this->length = mag(v[0]-v[1]);
        assert_error(this->length>0.0, "invalid line");
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
    vector_t<real_t> r_m, e[1], r_p;
    vector_t<real_t> L_m[1], L_p[1];
    int_t pg_m, pg_p;
    basis_1d_t(){}
    basis_1d_t(const vector_t<real_t> r_m, const vector_t<real_t> e_1, const vector_t<real_t> r_p, 
        const int_t pg_m, const int_t pg_p){
        this->r_m = r_m;
        this->e[0] = e_1;
        this->r_p = r_p;
        this->pg_m = pg_m;
        this->pg_p = pg_p;
        basis_1d_t::get_parameters();
    }
    void get_parameters(){
        for (size_t i=0; i<1; i++){
            this->L_m[i] = this->e[i]-this->r_m;
            this->L_p[i] = this->e[i]-this->r_p;
            assert_error(mag(this->L_m[i])>0.0&&mag(this->L_p[i])>0.0, "invalid 1d basis");
        }
    }
};

struct basis_2d_t{
    vector_t<real_t> r_m, e[2], r_p;
    vector_t<real_t> L_m[2], L_p[2];
    real_t L;
    real_t A_m[1], A_p[1];
    vector_t<real_t> n_m[1], n_p[1];
    int_t pg_m, pg_p;
    basis_2d_t(){}
    basis_2d_t(const vector_t<real_t> r_m, const vector_t<real_t> e_1, const vector_t<real_t> e_2, const vector_t<real_t> r_p,
        const int_t pg_m, const int_t pg_p){
        this->r_m = r_m;
        this->e[0] = e_1;
        this->e[1] = e_2;
        this->r_p = r_p;
        this->pg_m = pg_m;
        this->pg_p = pg_p;
        basis_2d_t::get_parameters();
    }
    void get_parameters(){
        this->L = mag(this->e[0]-this->e[1]);
        for (size_t i=0; i<2; i++){
            this->L_m[i] = this->e[i]-this->r_m;
            this->L_p[i] = this->e[i]-this->r_p;
            assert_error(mag(this->L_m[i])>0.0&&mag(this->L_p[i])>0.0, "invalid 1d basis");
        }
        vector_t<real_t> vector_m=-1.0*this->L_m[1]^this->L_m[0];
        vector_t<real_t> vector_p=+1.0*this->L_p[1]^this->L_p[0];
        this->A_m[0] = mag(vector_m)/2.0;
        this->A_p[0] = mag(vector_p)/2.0;
        assert_error(this->A_m[0]>0.0&&this->A_p[0]>0.0, "invalid 2d basis");
        this->n_m[0] = unit(vector_m);
        this->n_p[0] = unit(vector_p);
    }
};

struct basis_3d_t{
    vector_t<real_t> r_m, e[3], r_p;
    vector_t<real_t> L_m[3], L_p[3];
    vector_t<real_t> n_m[3], n_p[3];
    real_t A_m[3], A_p[3];
    real_t A;
    vector_t<real_t> nA_m, nA_p;
    real_t V_m, V_p;
    int_t pg_m, pg_p;
    basis_3d_t(){}
    basis_3d_t(const vector_t<real_t> r_m, const vector_t<real_t> e_1, const vector_t<real_t> e_2, 
        const vector_t<real_t> e_3, const vector_t<real_t> r_p, const int_t pg_m, const int_t pg_p){
        this->r_m = r_m;
        this->e[0] = e_1;
        this->e[1] = e_2;
        this->e[2] = e_3;
        this->r_p = r_p;
        this->pg_m = pg_m;
        this->pg_p = pg_p;
        basis_3d_t::get_parameters();
    }
    void get_parameters(){
        vector_t<real_t> vector=(this->e[2]-this->e[0])^(this->e[1]-this->e[0]);
        this->A = mag(vector)/2.0;
        this->nA_m = -1.0*unit(vector);
        this->nA_p = +1.0*unit(vector);
        for (size_t i=0; i<3; i++){
            this->L_m[i] = this->e[i]-this->r_m;
            this->L_p[i] = this->e[i]-this->r_p;
            assert_error(mag(this->L_m[i])>0.0&&mag(this->L_p[i])>0.0, "invalid 1d basis");
        }
        for (size_t i=0; i<3; i++){
            vector_t<real_t> vector_m=this->L_m[(i+1)%3]^this->L_m[(i+2)%3];
            vector_t<real_t> vector_p=this->L_p[(i+1)%3]^this->L_p[(i+2)%3];
            this->n_m[i] = -1.0*unit(vector_m);
            this->n_p[i] = +1.0*unit(vector_p);
            this->A_m[i] = mag(vector_m)/2.0;
            this->A_p[i] = mag(vector_p)/2.0;
            assert_error(this->A_m[0]>0.0&&this->A_p[0]>0.0, "invalid 2d basis");
        }
        this->V_m = -1.0*((this->L_m[2]-this->L_m[0])^(this->L_m[1]-this->L_m[0]))*this->L_m[0]/6.0;
        this->V_p = +1.0*((this->L_p[2]-this->L_p[0])^(this->L_p[1]-this->L_p[0]))*this->L_p[0]/6.0;
        assert_error(this->V_m>0.0&&this->V_p>0.0, "invalid 3d basis");
    }
};

class shape_t{
    private:
        int_t is_physical_specified=false;
        int_t is_basis_allocated=false;
        real_t freq=0.0, lambda=0.0;
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
        shape_t(){}
        ~shape_t(){}
        shape_t(const real_t freq, const complex_t mu_b, const complex_t eps_b){
            shape_t::set(freq, mu_b, eps_b);
        }
        void set(const real_t freq, const complex_t mu_b, const complex_t eps_b){
            assert_error(abs(mu_b)>0.0 && abs(eps_b)>0.0, "invalid background medium parameters");
            assert_error(freq>0.0, "invalid freq");
            this->freq = freq;
            this->mu_b = mu_b;
            this->eps_b = eps_b;
            this->lambda = c_0/(real(sqrt(mu_b*eps_b))*freq);
        }
        void get_basis_functions(const real_t clmax, const real_t metric_unit);
        void clear();
        void get_info(size_t &N_basis_1d, size_t &N_basis_2d, size_t &N_basis_3d){
            N_basis_1d = this->N_basis_1d;
            N_basis_2d = this->N_basis_2d;
            N_basis_3d = this->N_basis_3d;
        }
        basis_1d_t get_basis_1d(const size_t i){
            assert_error(this->is_basis_allocated, "no basis functions available");
            assert_error(i<this->N_basis_1d, "index is out of range");
            return this->basis_1d_list[i];
        }
        basis_2d_t get_basis_2d(const size_t i){
            assert_error(this->is_basis_allocated, "no basis functions available");
            assert_error(i<this->N_basis_2d, "index is out of range");
            return this->basis_2d_list[i];
        }
        basis_3d_t get_basis_3d(const size_t i){
            assert_error(this->is_basis_allocated, "no basis functions available");
            assert_error(i<this->N_basis_3d, "index is out of range");
            return this->basis_3d_list[i];
        }
};

// Functions
void call_gmsh(const real_t tol);
void create_vertical_wire_dipole(const real_t length, const real_t port_length, const real_t clmax);
void create_loop(const real_t radius, const real_t clmax);
void create_two_vertical_wire_dipole(const real_t length, const real_t port_length, const real_t d, 
    const real_t clmax);
void create_vertical_wire(const real_t length, const real_t clmax);
void create_transmission_line(const real_t L, const real_t S, const real_t clmax);
void create_sphere(const real_t radius);
void create_patch_antenna();

#endif
