//
#include "nu_integrand.hpp"

// 1d 1d


// 1d 2d
// 1d 3d

// 2d 1d
// 2d 2d

// 2d 3d

// 3d 1d
// 3d 2d
// 3d 3d

complex_t nu_3d_3d(const basis_3d_t b_m, const basis_3d_t b_n, const real_t lambda, const complex_t eps_b){
    const real_t tol=lambda*1.0E-6;
    vector_t<real_t> L_m1, L_m2, L_m3;
    vector_t<real_t> L_n1, L_n2, L_n3;
    complex_t factor=0.0;
    complex_t chi=0.0;
    complex_t ans=0.0;
    // 
    size_t counter=0;
    tetrahedron_t tetrahedron_m, tetrahedron_n;
    // mm
    tetrahedron_m = tetrahedron_t(b_m.r_m, b_m.e[0], b_m.e[1], b_m.e[2], b_m.pg_m);
    tetrahedron_n = tetrahedron_t(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], b_n.pg_m);
    counter = 0;
    for (size_t i=0; i<4; i++){
        for (size_t j=0; j<4; j++){
            if (is_equal(tetrahedron_m.v[i], tetrahedron_n.v[j], tol)){
                counter++;
            }
        }
    }
    if (counter==4){
        chi = (b_n.eps_m/eps_b)-1.0;
        factor = +1.0*(2.0*b_m.A*b_n.A)/(3.0*b_n.V_m);
        factor*=1.0/chi;
        for (size_t i=0; i<3; i++){
            for (size_t j=0; j<3; j++){
                if (is_equal(b_m.L_m[i]-b_n.L_m[j], -1.0*(b_m.r_m-b_n.r_m), tol)){
                    ans+=2.0*factor*(b_m.L_m[i]*b_n.L_m[j])/120.0;
                }else{
                    ans+=1.0*factor*(b_m.L_m[i]*b_n.L_m[j])/120.0;
                }
            }
        }
    }
    // mp
    tetrahedron_m = tetrahedron_t(b_m.r_m, b_m.e[0], b_m.e[1], b_m.e[2], b_m.pg_m);
    tetrahedron_n = tetrahedron_t(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], b_n.pg_p);
    counter = 0;
    for (size_t i=0; i<4; i++){
        for (size_t j=0; j<4; j++){
            if (is_equal(tetrahedron_m.v[i], tetrahedron_n.v[j], tol)){
                counter++;
            }
        }
    }
    if (counter==4){
        chi = (b_n.eps_p/eps_b)-1.0;
        factor = -1.0*(2.0*b_m.A*b_n.A)/(3.0*b_n.V_p);
        factor*=1.0/chi;
        for (size_t i=0; i<3; i++){
            for (size_t j=0; j<3; j++){
                if (is_equal(b_m.L_m[i]-b_n.L_p[j], -1.0*(b_m.r_m-b_n.r_p), tol)){
                    ans+=2.0*factor*(b_m.L_m[i]*b_n.L_p[j])/120.0;
                }else{
                    ans+=1.0*factor*(b_m.L_m[i]*b_n.L_p[j])/120.0;
                }
            }
        }
    }
    // pm
    tetrahedron_m = tetrahedron_t(b_m.r_p, b_m.e[2], b_m.e[1], b_m.e[0], b_m.pg_p);
    tetrahedron_n = tetrahedron_t(b_n.r_m, b_n.e[0], b_n.e[1], b_n.e[2], b_n.pg_m);
    counter = 0;
    for (size_t i=0; i<4; i++){
        for (size_t j=0; j<4; j++){
            if (is_equal(tetrahedron_m.v[i], tetrahedron_n.v[j], tol)){
                counter++;
            }
        }
    }
    if (counter==4){
        chi = (b_n.eps_m/eps_b)-1.0;
        factor = -1.0*(2.0*b_m.A*b_n.A)/(3.0*b_n.V_m);
        factor*=1.0/chi;
        for (size_t i=0; i<3; i++){
            for (size_t j=0; j<3; j++){
                if (is_equal(b_m.L_p[i]-b_n.L_m[j], -1.0*(b_m.r_p-b_n.r_m), tol)){
                    ans+=2.0*factor*(b_m.L_p[i]*b_n.L_m[j])/120.0;
                }else{
                    ans+=1.0*factor*(b_m.L_p[i]*b_n.L_m[j])/120.0;
                }
            }
        }
    }
    // pp
    tetrahedron_m = tetrahedron_t(b_m.r_p, b_m.e[2], b_m.e[1], b_m.e[0], b_m.pg_p);
    tetrahedron_n = tetrahedron_t(b_n.r_p, b_n.e[2], b_n.e[1], b_n.e[0], b_n.pg_p);
    counter = 0;
    for (size_t i=0; i<4; i++){
        for (size_t j=0; j<4; j++){
            if (is_equal(tetrahedron_m.v[i], tetrahedron_n.v[j], tol)){
                counter++;
            }
        }
    }
    if (counter==4){
        chi = (b_n.eps_p/eps_b)-1.0;
        factor = +1.0*(2.0*b_m.A*b_n.A)/(3.0*b_n.V_p);
        factor*=1.0/chi;
        for (size_t i=0; i<3; i++){
            for (size_t j=0; j<3; j++){
                if (is_equal(b_m.L_p[i]-b_n.L_p[j], -1.0*(b_m.r_p-b_n.r_p), tol)){
                    ans+=2.0*factor*(b_m.L_p[i]*b_n.L_p[j])/120.0;
                }else{
                    ans+=1.0*factor*(b_m.L_p[i]*b_n.L_p[j])/120.0;
                }
            }
        }
    }
    //
    return ans;
}