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

// struct point_t{
//     vector_t<real_t> r;
//     point_t(){}
//     point_t(const vector_t<real_t> r){
//         this->r = r;
//     }
// };

// struct edge_t{
//     point_t point[2];
//     vector_t<real_t> L;
//     real_t l=0.0;
//     int_t physical_group=-1;
//     edge_t(){}
//     edge_t(const point_t p_1, const point_t p_2){
//         this->point[0] = p_1;
//         this->point[1] = p_2;
//         //
//         this->L = p_2.r-p_1.r;
//         this->l = mag(this->L);
//     }
// };

// struct triangle_t{
//     edge_t edge[3];
//     real_t area=0.0;
//     vector_t<real_t> n;
//     int_t physical_group=-1;
//     triangle_t(){}
//     triangle_t(const point_t p_1, const point_t p_2, const point_t p_3){
//         this->edge[0] = edge_t(p_1, p_2);
//         this->edge[1] = edge_t(p_2, p_3);
//         this->edge[2] = edge_t(p_3, p_1);
//         //
//         this->n = -1.0*unit(edge[1].L^edge[0].L);
//         this->area = 0.5*(edge[1].L^edge[0].L)*this->n;
//         assert(this->area>0.0);
//     }
// };

// struct tetrahedron_t{
//     triangle_t triangle[4];
//     real_t volume=0.0;
//     int_t physical_group=-1;
//     tetrahedron_t(){}
//     tetrahedron_t(const point_t p_1, const point_t p_2, const point_t p_3, const point_t p_4){
//         this->triangle[0] = triangle_t(p_1, p_2, p_4);
//         this->triangle[1] = triangle_t(p_1, p_3, p_2);
//         this->triangle[2] = triangle_t(p_1, p_4, p_3);
//         this->triangle[3] = triangle_t(p_2, p_3, p_4);
//         //
//         this->volume = (triangle[0].edge[0].L^triangle[1].edge[0].L)*triangle[3].edge[0].L/6.0;
//         assert(this->volume>0.0);
//     }
// };

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
    vector_t<real_t> r_m, e_1, r_p;
    vector_t<real_t> L_m, L_p;
    real_t l_m, l_p;
    int_t physical_group=-1;
    basis_1d_t(const vector_t<real_t> r_m, const vector_t<real_t> e_1, const vector_t<real_t> r_p){
        this->r_m = r_m;
        this->e_1 = e_1;
        this->r_p = r_p;
        basis_1d_t::get_parameters();
    }
    void get_parameters(){
        this->L_m = +1.0*(this->e_1-this->r_m);
        this->L_p = +1.0*(this->e_1-this->r_p);
        this->l_m = mag(this->L_m);
        this->l_p = mag(this->L_p);
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
        // basis_2d_t *basis_2d_list=null;
        // basis_3d_t *basis_3d_list=null;
        const real_t mesh_tol=1.0E-8;
        void load_mesh(const real_t metric_unit);
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

#endif
