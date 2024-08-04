//
#include "engine_1d.hpp"

void R_mn_1d(const real_t alpha, const real_t alpha_, 
    real_t &R_mm, real_t &R_mp, real_t &R_pm, real_t &R_pp, 
    basis_1d_t b_m, basis_1d_t b_n, const real_t a){   
    vector_t<real_t> rho_m_m = b_m.r_m+alpha*(b_m.e[0]-b_m.r_m);
    vector_t<real_t> rho_m_p = b_m.r_p+alpha*(b_m.e[0]-b_m.r_p);
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

struct integrand_1d_args{
    quadl_domain_t quadl;
    basis_1d_t b_m, b_n;
    complex_t k=0.0, eta=0.0;
    real_t lambda=0.0;
    real_t a=0.0;
    complex_t alpha=0.0;
};

const real_t projection_tol=1.0E-4;

void integrand_L1_internal(const real_t alpha, basis_1d_t b_m, basis_1d_t b_n, 
    real_t &I_mm, real_t &I_mp, real_t &I_pm, real_t &I_pp, const real_t a){
    projection_1d_para para;
    vector_t<real_t> p;
    real_t ans;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_mm = ans;
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_mp = ans;
    // pm
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_pm = ans;
    // pp
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_pp = ans;
}

void integrand_L2_internal(const real_t alpha, basis_1d_t b_m,basis_1d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp, 
    const real_t a){
    projection_1d_para para;
    vector_t<real_t> p;
    vector_t<real_t> ans;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_mm = ans;
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_mp = ans;
    // pm
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_pm = ans;
    // pp
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_pp = ans;
}

void integrand_L3_internal(const real_t alpha, basis_1d_t b_m,basis_1d_t b_n, 
    vector_t<real_t> &I_mm, vector_t<real_t> &I_mp, vector_t<real_t> &I_pm, vector_t<real_t> &I_pp, 
    const real_t lambda){
    projection_1d_para para;
    vector_t<real_t> p;
    vector_t<real_t> ans;
    // mm
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_mm = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // mp
    p = b_m.r_m+alpha*b_m.L_m[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_mp = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // pm
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.r_m, b_n.e[0], p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_pm = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // pp
    p = b_m.r_p+alpha*b_m.L_p[0];
    para = prjection_1d(b_n.e[0], b_n.r_p, p);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
        -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    I_pp = para.P_0 < projection_tol*lambda ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
}

complex_t psi_1d_singular_integrand_inner(const complex_t alpha_, void *args_){
    integrand_1d_args *args=(integrand_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    const complex_t k=args->k;
    const real_t a=args->a;
    complex_t alpha=args->alpha;
    complex_t I_mm, I_mp, I_pm, I_pp;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mm, R_mp, R_pm, R_pp;
    R_mn_1d(real(alpha), real(alpha_), R_mm, R_mp, R_pm, R_pp, b_m, b_n, a);
    I_mm = +1.0*(b_m.L_m[0]*b_n.L_m[0])*alpha*alpha_*(-j*k*exp(-j*k*R_mm/2.0))*sinc(k*R_mm/2.0);
    I_mp = -1.0*(b_m.L_m[0]*b_n.L_p[0])*alpha*alpha_*(-j*k*exp(-j*k*R_mp/2.0))*sinc(k*R_mp/2.0);
    I_pm = -1.0*(b_m.L_p[0]*b_n.L_m[0])*alpha*alpha_*(-j*k*exp(-j*k*R_pm/2.0))*sinc(k*R_pm/2.0);
    I_pp = +1.0*(b_m.L_p[0]*b_n.L_p[0])*alpha*alpha_*(-j*k*exp(-j*k*R_pp/2.0))*sinc(k*R_pp/2.0);
    return (I_mm+I_mp+I_pm+I_pp)/(4.0*pi);
}

complex_t phi_1d_singular_integrand_inner(const complex_t alpha_, void *args_){
    integrand_1d_args *args=(integrand_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    const complex_t k=args->k;
    const real_t a=args->a;
    complex_t alpha=args->alpha;
    complex_t I_mm, I_mp, I_pm, I_pp;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mm, R_mp, R_pm, R_pp;
    R_mn_1d(real(alpha), real(alpha_), R_mm, R_mp, R_pm, R_pp, b_m, b_n, a);
    I_mm = (-j*k*exp(-j*k*R_mm/2.0))*sinc(k*R_mm/2.0);
    I_mp = (-j*k*exp(-j*k*R_mp/2.0))*sinc(k*R_mp/2.0);
    I_pm = (-j*k*exp(-j*k*R_pm/2.0))*sinc(k*R_pm/2.0);
    I_pp = (-j*k*exp(-j*k*R_pp/2.0))*sinc(k*R_pp/2.0);
    return (I_mm-I_mp-I_pm+I_pp)/(4.0*pi);
}

complex_t psi_1d_singular_integrand_outer(const complex_t alpha, void *args_){
    integrand_1d_args *args=(integrand_1d_args*)args_;
    args->alpha = alpha;
    edge_domain_t edge={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
    int flag;
    complex_t ans=args->quadl.integral_1d(psi_1d_singular_integrand_inner, args, edge, flag);
    return ans;
}

complex_t phi_1d_singular_integrand_outer(const complex_t alpha, void *args_){
    integrand_1d_args *args=(integrand_1d_args*)args_;
    args->alpha = alpha;
    edge_domain_t edge={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
    int flag;
    complex_t ans=args->quadl.integral_1d(phi_1d_singular_integrand_inner, args, edge, flag);
    return ans;
}

complex_t psi_1d_integrand_1(const complex_t alpha, void *args_){
    integrand_1d_args *args=(integrand_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    real_t a=args->a;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> R_m_m, R_m_p, R_p_m, R_p_p;
    R_m_m = b_m.r_m+real(alpha)*b_m.L_m[0]-b_n.r_m;
    R_m_p = b_m.r_m+real(alpha)*b_m.L_m[0]-b_n.r_p;
    R_p_m = b_m.r_p+real(alpha)*b_m.L_p[0]-b_n.r_m;
    R_p_p = b_m.r_p+real(alpha)*b_m.L_p[0]-b_n.r_p;
    integrand_L1_internal(real(alpha), b_m, b_n, I_mm, I_mp, I_pm, I_pp, a);
    complex_t ans=0.0;
    ans+= +1.0*(real(alpha)*b_m.L_m[0]*R_m_m)*I_mm/mag(b_n.L_m[0]);
    ans+= -1.0*(real(alpha)*b_m.L_m[0]*R_m_p)*I_mp/mag(b_n.L_p[0]);
    ans+= -1.0*(real(alpha)*b_m.L_p[0]*R_p_m)*I_pm/mag(b_n.L_m[0]);
    ans+= +1.0*(real(alpha)*b_m.L_p[0]*R_p_p)*I_pp/mag(b_n.L_p[0]);
    return ans/(4.0*pi);
}

complex_t psi_1d_integrand_2(const complex_t alpha, void *args_){
    integrand_1d_args *args=(integrand_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    real_t a=args->a;
    vector_t<real_t> I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*real(alpha)*b_m.L_m[0];
    rho_p = -1.0*real(alpha)*b_m.L_p[0];
    integrand_L2_internal(real(alpha), b_m, b_n, I_mm, I_mp, I_pm, I_pp, a);
    complex_t ans=0.0;
    ans+= +1.0*(real(alpha)*b_m.L_m[0])*I_mm/mag(b_n.L_m[0]);
    ans+= -1.0*(real(alpha)*b_m.L_m[0])*I_mp/mag(b_n.L_p[0]);
    ans+= -1.0*(real(alpha)*b_m.L_p[0])*I_pm/mag(b_n.L_m[0]);
    ans+= +1.0*(real(alpha)*b_m.L_p[0])*I_pp/mag(b_n.L_p[0]);
    return ans/(4.0*pi);
}

complex_t phi_1d_integrand_1(const complex_t alpha, void *args_){
    integrand_1d_args *args=(integrand_1d_args*)args_;
    basis_1d_t b_m=args->b_m;
    basis_1d_t b_n=args->b_n;
    real_t a=args->a;
    real_t I_mm, I_mp, I_pm, I_pp;
    integrand_L1_internal(real(alpha), b_m, b_n, I_mm, I_mp, I_pm, I_pp, a);
    complex_t ans=I_mm/mag(b_n.L_m[0])-I_mp/mag(b_n.L_p[0])-I_pm/mag(b_n.L_m[0])+I_pp/mag(b_n.L_p[0]);
    return ans/(4.0*pi);
}

complex_t psi_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, const real_t lambda, 
    const real_t a, const quadl_domain_t quadl, int_t &flag){
    integrand_1d_args args={quadl, b_m, b_n};
    args.a = a;
    args.k = k;
    args.lambda = lambda;
    complex_t I1=0.0, I2=0.0, I3=0.0;
    edge_domain_t edge={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    I1 = args.quadl.integral_1d(psi_1d_singular_integrand_outer, &args, edge, flag);
    I2 = args.quadl.integral_1d(psi_1d_integrand_1, &args, edge, flag);
    I3 = args.quadl.integral_1d(psi_1d_integrand_2, &args, edge, flag);
    return I1+I2+I3;
}

complex_t phi_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, const real_t lambda, 
    const real_t a, const quadl_domain_t quadl, int_t &flag){
    integrand_1d_args args={quadl, b_m, b_n};
    args.a = a;
    args.k = k;
    args.lambda = lambda;
    complex_t I1=0.0, I2=0.0;
    edge_domain_t edge={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    I1 = args.quadl.integral_1d(phi_1d_singular_integrand_outer, &args, edge, flag);
    I2 = args.quadl.integral_1d(phi_1d_integrand_1, &args, edge, flag);
    return I1+I2;
}

complex_t Z_mn_1d(const basis_1d_t b_m, const basis_1d_t b_n, const complex_t k, 
    const complex_t eta, const real_t lambda, const real_t a, 
    const quadl_domain_t quadl, int_t &flag){
    integrand_1d_args args={quadl, b_m, b_n};
    args.a = a;
    args.k = k;
    args.eta = eta;
    args.lambda = lambda;
    const complex_t j=complex_t(0.0, 1.0);
    complex_t psi, phi;
    psi = psi_1d(b_m, b_n, k, lambda, a, quadl, flag);
    phi = phi_1d(b_m, b_n, k, lambda, a, quadl, flag);
    return +j*k*eta*psi-j*(eta/k)*phi;
}

//

void engine_1d_t::set(const real_t freq, const complex_t mu_b, const complex_t eps_b, 
    const real_t clmax, const real_t unit_metric, const real_t a){
    if (this->is_engine_set){engine_1d_t::unset();}
    // check inputs errors
    this->a = a;
    this->mu_b = mu_b;
    this->eps_b = eps_b;
    this->lambda = (c_0/real(sqrt(mu_b*eps_b)))/freq;
    this->k_b = 2.0*pi*freq*sqrt(mu_0*eps_0)*sqrt(mu_b*eps_b);
    this->eta_b = sqrt(mu_0/eps_0)*sqrt(mu_b/eps_b);
    this->shape.set(freq, mu_b, eps_b);
    shape.get_basis_functions(clmax, unit_metric);

}

void engine_1d_t::unset(){
    this->Z_mn.unset();
    this->V_m.unset();
    this->I_n.unset();
    this->is_Z_mn_allocated = false;
    this->is_V_m_allocated = false;
    this->is_I_n_allocated = false;
    this->is_Z_mn_calculated = false;
    this->is_V_m_calculated = false;
    this->is_I_n_calculated = false;
    this->is_engine_set = false;
    shape.clear();
}