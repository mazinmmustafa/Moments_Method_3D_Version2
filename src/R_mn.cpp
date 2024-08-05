//
#include "R_mn.hpp"

// 1d
void R_mn_1d_1d(const real_t alpha, const real_t alpha_, 
    const basis_1d_t b_m, const basis_1d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

void R_mn_1d_2d(const real_t alpha, const real_t alpha_, const real_t beta_, 
    const basis_1d_t b_m, const basis_2d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m)+beta_*(b_n.e[1]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p)+beta_*(b_n.e[1]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

void R_mn_1d_3d(const real_t alpha, const real_t alpha_, const real_t beta_, const real_t gamma_, 
    const basis_1d_t b_m, const basis_3d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha* (b_m.e[0]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha* (b_m.e[0]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m)+beta_*(b_n.e[1]-b_n.r_m)+gamma_*(b_n.e[2]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p)+beta_*(b_n.e[1]-b_n.r_p)+gamma_*(b_n.e[2]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

// 2d
void R_mn_2d_1d(const real_t alpha, const real_t beta, const real_t alpha_, 
    const basis_2d_t b_m, const basis_1d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m)+beta *(b_m.e[1]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p)+beta *(b_m.e[1]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

void R_mn_2d_2d(const real_t alpha, const real_t beta, const real_t alpha_, const real_t beta_, 
    const basis_2d_t b_m, const basis_2d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m)+beta *(b_m.e[1]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p)+beta *(b_m.e[1]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m)+beta_*(b_n.e[1]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p)+beta_*(b_n.e[1]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

void R_mn_2d_3d(const real_t alpha, const real_t beta, const real_t alpha_, const real_t beta_, const real_t gamma_, 
    const basis_2d_t b_m, const basis_3d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m)+beta *(b_m.e[1]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p)+beta *(b_m.e[1]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m)+beta_*(b_n.e[1]-b_n.r_m)+gamma_*(b_n.e[2]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p)+beta_*(b_n.e[1]-b_n.r_p)+gamma_*(b_n.e[2]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

// 3d
void R_mn_3d_1d(const real_t alpha, const real_t beta, const real_t gamma, const real_t alpha_, 
    const basis_3d_t b_m, const basis_1d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m)+beta *(b_m.e[1]-b_m.r_m)+gamma *(b_m.e[2]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p)+beta *(b_m.e[1]-b_m.r_p)+gamma *(b_m.e[2]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

void R_mn_3d_2d(const real_t alpha, const real_t beta, const real_t gamma, const real_t alpha_, const real_t beta_, 
    const basis_3d_t b_m, const basis_2d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m)+beta *(b_m.e[1]-b_m.r_m)+gamma *(b_m.e[2]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p)+beta *(b_m.e[1]-b_m.r_p)+gamma *(b_m.e[2]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m)+beta_*(b_n.e[1]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p)+beta_*(b_n.e[1]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}

void R_mn_3d_3d(const real_t alpha, const real_t beta, const real_t gamma, const real_t alpha_, const real_t beta_, const real_t gamma_, 
    const basis_3d_t b_m, const basis_3d_t b_n,  
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, const real_t a){  
    vector_t<real_t> rho_m_m = b_m.r_m+alpha *(b_m.e[0]-b_m.r_m)+beta *(b_m.e[1]-b_m.r_m)+gamma *(b_m.e[2]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha *(b_m.e[0]-b_m.r_p)+beta *(b_m.e[1]-b_m.r_p)+gamma *(b_m.e[2]-b_m.r_p);
    vector_t<real_t> rho_n_m = b_n.r_m+alpha_*(b_n.e[0]-b_n.r_m)+beta_*(b_n.e[1]-b_n.r_m)+gamma_*(b_n.e[2]-b_n.r_m);
    vector_t<real_t> rho_n_p = b_n.r_p+alpha_*(b_n.e[0]-b_n.r_p)+beta_*(b_n.e[1]-b_n.r_p)+gamma_*(b_n.e[2]-b_n.r_p);
    R_mm = mag(rho_m_m-rho_n_m); 
    R_mp = mag(rho_m_m-rho_n_p);
    R_pm = mag(rho_m_p-rho_n_m);
    R_pp = mag(rho_m_p-rho_n_p);
    R_mm = sqrt(R_mm*R_mm+a*a);
    R_mp = sqrt(R_mp*R_mp+a*a);
    R_pm = sqrt(R_pm*R_pm+a*a);
    R_pp = sqrt(R_pp*R_pp+a*a);
}