#ifndef __ENGINE_1D_HPP__
#define __ENGINE_1D_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
#include "shape.hpp"
#include "projection.hpp"
//

// Definitions
class engine_1d_t{
    private:
        int_t is_engine_set=false;
        int_t is_Z_mn_allocated=false;
        int_t is_V_m_allocated=false;
        int_t is_I_n_allocated=false;
        int_t is_Z_mn_calculated=false;
        int_t is_V_m_calculated=false;
        int_t is_I_n_calculated=false;
        quadl_domain_t quadl;
        shape_t shape;
        matrix_t<complex_t> Z_mn, V_m, I_n;
        real_t freq=0.0, lambda=0.0;
        complex_t mu_b=1.0, eps_b=1.0;
        complex_t k_b=0.0, eta_b=0.0;
        real_t a=0.0;
    public:         
        engine_1d_t(){}
        ~engine_1d_t(){}
        void set(const real_t freq, const complex_t mu_b, const complex_t eps_b, 
            const real_t clmax, const real_t unit_metric, const real_t a);
        void compute_Z_mn();
        void unset();
};

// Functions
complex_t phi_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, const real_t lambda, 
    const real_t a, const quadl_domain_t quadl, int_t &flag);
complex_t psi_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, const real_t lambda, 
    const real_t a, const quadl_domain_t quadl, int_t &flag);
complex_t Z_mn_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, const real_t a, 
    const quadl_domain_t quadl, int_t &flag);

#endif
