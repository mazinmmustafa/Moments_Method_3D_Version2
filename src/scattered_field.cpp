//
#include "scattered_field.hpp"

// 1d

void integrand_L1_1d(basis_1d_t b_m, const vector_t<real_t> p, 
    real_t &I_m, real_t &I_p, const real_t a){
    projection_1d_para para;
    real_t ans;
    // m
    para = prjection_1d(b_m.r_m, b_m.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_m = ans;
    // p
    para = prjection_1d(b_m.e[0], b_m.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_p = ans;
}

void integrand_L2_1d(basis_1d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p, const real_t a){
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
         -(para.l_p/para.P_p-para.l_m/para.P_m)*para.P_0_unit/para.P_0;
    I_m = ans;
    // p
    para = prjection_1d(b_m.e[0], b_m.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
         -(para.l_p/para.P_p-para.l_m/para.P_m)*para.P_0_unit/para.P_0;
    I_p = ans;
}

complex_t E_1d_singular_integrand_1(const complex_t alpha, void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
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
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
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
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    real_t I_m, I_p;
    integrand_L1_1d(b_m, r, I_m, I_p, a);
    projection_1d_para para_m = prjection_1d(b_m.r_m, b_m.e[0], r);
    projection_1d_para para_p = prjection_1d(b_m.e[0], b_m.r_p, r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*(para_m.l_m*para_m.l_unit)*I_m/mag(b_m.L_m[0]);
    ans+= +1.0*args->unit_vector*(para_p.l_p*para_p.l_unit)*I_p/mag(b_m.L_p[0]);
    return ans/(4.0*pi);
}

complex_t E_1d_integral_2(void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
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

complex_t E_1d_integral_3(void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
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
    scattered_field_args_1d_t args;
    const complex_t j=complex_t(0.0, 1.0);
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.a = a;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2, I3, I4, I5;
    int_t flag;
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    //
    I1 = quadl.integral_1d(E_1d_singular_integrand_1, &args, line, flag);
    assert_error(!flag, "no convergence");
    I2 = E_1d_integral_1(&args);
    I3 = E_1d_integral_2(&args);
    I4 = quadl.integral_1d(E_1d_singular_integrand_2, &args, line, flag);
    assert_error(!flag, "no convergence");
    I5 = E_1d_integral_3(&args);
    return -j*k*eta*(I1+I2+I3)+j*(eta/k)*(I4+I5);
}

complex_t H_1d_singular_integrand_1(const complex_t alpha, void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
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
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*
        (alpha*(+1.0*b_m.L_m[0]^unit(R_m_vector-r))*args->unit_vector);
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*
        (alpha*(-1.0*b_m.L_p[0]^unit(R_p_vector-r))*args->unit_vector);
    return (I_m+I_p)/(4.0*pi);
}

complex_t H_1d_integral_1(void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_1d(b_m, r, I_m, I_p, a);
    projection_1d_para para_m=prjection_1d(b_m.r_m, b_m.e[0], r);
    projection_1d_para para_p=prjection_1d(b_m.e[0], b_m.r_p, r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*((para_m.l_m*para_m.l_unit)^I_m)/mag(b_m.L_m[0]);
    ans+= +1.0*args->unit_vector*((para_p.l_p*para_p.l_unit)^I_p)/mag(b_m.L_p[0]);
    return ans/(4.0*pi);
}

complex_t compute_H_1d(const basis_1d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, 
    const real_t a, quadl_domain_t quadl){
    scattered_field_args_1d_t args;
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.a = a;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2;
    int_t flag;
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    //
    I1 = quadl.integral_1d(H_1d_singular_integrand_1, &args, line, flag);
    assert_error(!flag, "no convergence");
    I2 = H_1d_integral_1(&args);
    return I1+I2;
}

// 2d 

void integrand_L1_2d(basis_2d_t b_m, const vector_t<real_t> p, 
    real_t &I_m, real_t &I_p){
    projection_2d_para para;
    // m
    para = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], p);
    I_m = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_m+=A*(B-abs(para.d)*(C-D));
    }
    // p
    para = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], p);
    I_p = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_p+=A*(B-abs(para.d)*(C-D));
    }
}

void integrand_L2_2d(basis_2d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p){
    projection_2d_para para;
    // m
    para = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], p);
    I_m = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        I_m = I_m+0.5*(A+B)*para.u[i];
    }
    // p
    para = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], p);
    I_p = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        I_p = I_p+0.5*(A+B)*para.u[i];
    }
}

void integrand_L3_2d(basis_2d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p){
    projection_2d_para para;
    // m
    para = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], p);
    I_m = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_m = I_m+A*para.u[i]+sign(para.d)*(B-C)*(para.u[i]*para.para_1d[i].P_0_unit)*para.n;
    }
    // p
    para = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], p);
    I_p = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I_p = I_p+A*para.u[i]+sign(para.d)*(B-C)*(para.u[i]*para.para_1d[i].P_0_unit)*para.n;
    }
}

complex_t E_2d_singular_integrand_1(const complex_t alpha, const complex_t beta, void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    I_m = -j*k*exp(-j*k*R_m/2.0)*sinc(k*R_m/2.0)*b_m.L*(+1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1])*args->unit_vector);
    I_p = -j*k*exp(-j*k*R_p/2.0)*sinc(k*R_p/2.0)*b_m.L*(-1.0*(alpha*b_m.L_p[1]+beta*b_m.L_p[0])*args->unit_vector);
    return (I_m+I_p)/(4.0*pi);
}

complex_t E_2d_singular_integrand_2(const complex_t alpha, const complex_t beta, void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*(unit(R_m_vector-r)*args->unit_vector)*2.0*b_m.L;
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*(unit(R_p_vector-r)*args->unit_vector)*2.0*b_m.L;
    return (I_m-I_p)/(4.0*pi);
}

complex_t E_2d_integral_1(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    real_t I_m, I_p;
    integrand_L1_2d(b_m, r, I_m, I_p);
    projection_2d_para para_m = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], r);
    projection_2d_para para_p = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*(b_m.r_m-para_m.p_0)*I_m*b_m.L/(2.0*b_m.A_m[0]);
    ans+= +1.0*args->unit_vector*(b_m.r_p-para_p.p_0)*I_p*b_m.L/(2.0*b_m.A_p[0]);
    return ans/(4.0*pi);
}

complex_t E_2d_integral_2(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L2_2d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m*b_m.L/(2.0*b_m.A_m[0]);
    ans+= -1.0*args->unit_vector*I_p*b_m.L/(2.0*b_m.A_p[0]);
    return ans/(4.0*pi);
}

complex_t E_2d_integral_3(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_2d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m*b_m.L/b_m.A_m[0];
    ans+= -1.0*args->unit_vector*I_p*b_m.L/b_m.A_p[0];
    return ans/(4.0*pi);
}

complex_t compute_E_2d(const basis_2d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl){
    scattered_field_args_2d_t args;
    const complex_t j=complex_t(0.0, 1.0);
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2, I3, I4, I5;
    int_t flag;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0), vector_t<real_t>(0.0, +1.0, 0.0)};
    //
    I1 = quadl.integral_2d(E_2d_singular_integrand_1, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    I2 = E_2d_integral_1(&args);
    I3 = E_2d_integral_2(&args);
    I4 = quadl.integral_2d(E_2d_singular_integrand_2, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    I5 = E_2d_integral_3(&args);
    return -j*k*eta*(I1+I2+I3)+j*(eta/k)*(I4+I5);
}

complex_t H_2d_singular_integrand_1(const complex_t alpha, const complex_t beta, void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*
        ((+1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1])^unit(R_m_vector-r))*args->unit_vector)*b_m.L;
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*
        ((-1.0*(alpha*b_m.L_p[1]+beta*b_m.L_p[0])^unit(R_p_vector-r))*args->unit_vector)*b_m.L;
    return (I_m+I_p)/(4.0*pi);
}

complex_t H_2d_integral_1(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_2d(b_m, r, I_m, I_p);
    projection_2d_para para_m=prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], r);
    projection_2d_para para_p=prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*((b_m.r_m-para_m.p_0)^I_m)*b_m.L/(2.0*b_m.A_m[0]);
    ans+= +1.0*args->unit_vector*((b_m.r_p-para_p.p_0)^I_p)*b_m.L/(2.0*b_m.A_p[0]);
    return ans/(4.0*pi);
}

complex_t compute_H_2d(const basis_2d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl){
    scattered_field_args_2d_t args;
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2;
    int_t flag;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0), vector_t<real_t>(0.0, +1.0, 0.0)};
    //
    I1 = quadl.integral_2d(H_2d_singular_integrand_1, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    I2 = H_2d_integral_1(&args);
    return I1+I2;
}