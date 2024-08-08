//
#include "scattered_field.hpp"

// 1d

void integrand_L2_1d(basis_1d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p,const real_t a){
    projection_1d_para para;
    vector_t<real_t> ans;
    // m
    para = prjection_1d(b_m.r_m, b_m.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_m = ans;
    // p
    para = prjection_1d(b_m.e[0], b_m.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_p = ans;
}

void integrand_L3_1d(basis_1d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p, const real_t a){
    projection_1d_para para;
    vector_t<real_t> ans;
    // m
    para = prjection_1d(b_m.r_m, b_m.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
         -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    if (isnan(mag(ans))){ans = vector_t<real_t>(0.0, 0.0, 0.0);}
    I_m = ans;
    // p
    para = prjection_1d(b_m.e[0], b_m.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
         -(atan2(para.l_p, para.P_0)-atan2(para.l_m, para.P_0))*para.P_0_unit/para.P_0;
    if (isnan(mag(ans))){ans = vector_t<real_t>(0.0, 0.0, 0.0);}
    I_p = ans;
}

complex_t E_1d_singular_integrand_1(const complex_t alpha, void *args_){
    scattered_field_args_t *args=(scattered_field_args_t*)args_;
    basis_1d_t b_m=args->b_m;
    const complex_t k=args->k;
    const real_t a=args->a;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    R_m = sqrt(R_m*R_m+a*a);
    R_p = sqrt(R_p*R_p+a*a);
    I_m = -j*k*exp(-j*k*R_m/2.0)*sinc(k*R_m/2.0)*(+1.0*alpha*b_m.L_m[0]*args->unit_vector);
    I_p = -j*k*exp(-j*k*R_p/2.0)*sinc(k*R_p/2.0)*(-1.0*alpha*b_m.L_p[0]*args->unit_vector);
    return (I_m+I_p)/(4.0*pi);
}

complex_t E_1d_singular_integrand_2(const complex_t alpha, void *args_){
    scattered_field_args_t *args=(scattered_field_args_t*)args_;
    basis_1d_t b_m=args->b_m;
    const complex_t k=args->k;
    const real_t a=args->a;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    R_m = sqrt(R_m*R_m+a*a);
    R_p = sqrt(R_p*R_p+a*a);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*(unit(R_m_vector-r)*args->unit_vector);
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*(unit(R_p_vector-r)*args->unit_vector);
    return (I_m-I_p)/(4.0*pi);
}

complex_t E_1d_integral_1(void *args_){
    scattered_field_args_t *args=(scattered_field_args_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L2_1d(b_m, r, I_m, I_p, a);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m/mag(b_m.L_m[0]);
    ans+= -1.0*args->unit_vector*I_p/mag(b_m.L_p[0]);
    return ans/(4.0*pi);
}

complex_t E_1d_integral_2(void *args_){
    scattered_field_args_t *args=(scattered_field_args_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_1d(b_m, r, I_m, I_p, a);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m/mag(b_m.L_m[0]);
    ans+= -1.0*args->unit_vector*I_p/mag(b_m.L_p[0]);
    return ans/(4.0*pi);
}

complex_t compute_E_1d(const basis_1d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, 
    const real_t a, quadl_domain_t quadl){
    scattered_field_args_t args;
    const complex_t j=complex_t(0.0, 1.0);
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.a = a;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2, I3, I4;
    int_t flag;
    edge_domain_t edge={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    //
    I1 = quadl.integral_1d(E_1d_singular_integrand_1, &args, edge, flag);
    assert_error(!flag, "no convergence");
    I2 = E_1d_integral_1(&args);
    I3 = quadl.integral_1d(E_1d_singular_integrand_2, &args, edge, flag);
    assert_error(!flag, "no convergence");
    I4 = E_1d_integral_2(&args);
    return -j*k*eta*((1.0*I1+1.0*I2)-(0.0*I3+0.0*I4)/(k*k));
}

// 1d 2d
// 1d 3d

// 2d 1d
// 2d 2d
// 2d 3d

// 3d 1d
// 3d 2d
// 3d 3d