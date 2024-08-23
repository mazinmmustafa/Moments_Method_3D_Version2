//
#include "incident_field.hpp"

incident_field_t compute_incident_field(const complex_t E_TM, const complex_t E_TE, const real_t theta_i, const real_t phi_i, 
    const real_t k, const real_t eta, const vector_t<real_t> r){
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> theta_i_u=vector_t<real_t>(cos(theta_i)*cos(phi_i), cos(theta_i)*sin(phi_i), -sin(theta_i));
    vector_t<real_t> phi_i_u=vector_t<real_t>(-sin(phi_i), cos(phi_i), 0.0);
    vector_t<real_t> k_i=k*vector_t<real_t>(sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos(theta_i));
    vector_t<complex_t> E_i=(E_TM*theta_i_u+E_TE*phi_i_u)*exp(+j*(k_i*r));
    vector_t<complex_t> H_i=(E_TM*(k_i^theta_i_u)+E_TE*(k_i^phi_i_u))*exp(+j*(k_i*r))/eta;
    incident_field_t incident_field;
    incident_field.E = E_i;
    incident_field.H = H_i;
    return incident_field;
}

// 1d incident
complex_t compute_incident_E_integrand_1d(const complex_t alpha, void *args_){
    assert(args_!=null);
    incident_field_args_1d_t *args=(incident_field_args_1d_t*)args_;
    vector_t<real_t> r;
    complex_t I_m=0.0, I_p=0.0;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*real(alpha)*args->b_m.L_m[0];
    rho_p = +1.0*real(alpha)*args->b_m.L_p[0];
    incident_field_t incident_field;
    //
    r = args->b_m.r_m+rho_m;
    incident_field = compute_incident_field(args->E_TM, args->E_TE, args->theta_i, args->phi_i, 
        args->k, args->eta, r);
    I_m = +1.0*rho_m*incident_field.E;
    //
    r = args->b_m.r_p+rho_p;
    incident_field = compute_incident_field(args->E_TM, args->E_TE, args->theta_i, args->phi_i, 
        args->k, args->eta, r);
    I_p = -1.0*rho_p*incident_field.E;
    return I_m+I_p;
}

complex_t compute_scattered_far_field_E_theta_integrand_1d(const complex_t alpha, void *args_){
    assert(args_!=null);
    incident_field_args_1d_t *args=(incident_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    const complex_t j=complex_t(0.0, 1.0);
    const real_t theta_i=args->theta_i;
    const real_t phi_i=args->phi_i;
    const real_t k=args->k;
    const real_t eta=args->eta;
    vector_t<real_t> theta_i_u=vector_t<real_t>(cos(theta_i)*cos(phi_i), cos(theta_i)*sin(phi_i), -sin(theta_i));
    vector_t<real_t> phi_i_u=vector_t<real_t>(-sin(phi_i), cos(phi_i), 0.0);
    vector_t<real_t> k_i=k*vector_t<real_t>(sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos(theta_i));
    vector_t<real_t> r;
    complex_t I_m=0.0, I_p=0.0;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*real(alpha)*b_m.L_m[0];
    rho_p = +1.0*real(alpha)*b_m.L_p[0];
    incident_field_t incident_field;
    r = b_m.r_m+rho_m;
    I_m = (k*eta/(4.0*pi))*(+1.0*theta_i_u*rho_m)*exp(+j*(k_i*r));
    r = b_m.r_p+rho_p;
    I_p = (k*eta/(4.0*pi))*(-1.0*theta_i_u*rho_p)*exp(+j*(k_i*r));
    return I_m+I_p;
}

complex_t compute_scattered_far_field_E_phi_integrand_1d(const complex_t alpha, void *args_){
    assert(args_!=null);
    incident_field_args_1d_t *args=(incident_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    const complex_t j=complex_t(0.0, 1.0);
    const real_t theta_i=args->theta_i;
    const real_t phi_i=args->phi_i;
    const real_t k=args->k;
    const real_t eta=args->eta;
    vector_t<real_t> theta_i_u=vector_t<real_t>(cos(theta_i)*cos(phi_i), cos(theta_i)*sin(phi_i), -sin(theta_i));
    vector_t<real_t> phi_i_u=vector_t<real_t>(-sin(phi_i), cos(phi_i), 0.0);
    vector_t<real_t> k_i=k*vector_t<real_t>(sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos(theta_i));
    vector_t<real_t> r;
    complex_t I_m=0.0, I_p=0.0;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*real(alpha)*b_m.L_m[0];
    rho_p = +1.0*real(alpha)*b_m.L_p[0];
    incident_field_t incident_field;
    r = b_m.r_m+rho_m;
    I_m = (k*eta/(4.0*pi))*(+1.0*phi_i_u*rho_m)*exp(+j*(k_i*r));
    r = b_m.r_p+rho_p;
    I_p = (k*eta/(4.0*pi))*(-1.0*phi_i_u*rho_p)*exp(+j*(k_i*r));
    return I_m+I_p;
}

// 2d incident
complex_t compute_incident_E_integrand_2d(const complex_t alpha, const complex_t beta, void *args_){
    assert(args_!=null);
    incident_field_args_2d_t *args=(incident_field_args_2d_t*)args_;
    vector_t<real_t> r;
    complex_t I_m=0.0, I_p=0.0;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*(real(alpha)*args->b_m.L_m[0]+real(beta)*args->b_m.L_m[1]);
    rho_p = +1.0*(real(alpha)*args->b_m.L_p[1]+real(beta)*args->b_m.L_p[0]);
    incident_field_t incident_field;
    //
    r = args->b_m.r_m+rho_m;
    incident_field = compute_incident_field(args->E_TM, args->E_TE, args->theta_i, args->phi_i, 
        args->k, args->eta, r);
    I_m = +1.0*rho_m*incident_field.E;
    //
    r = args->b_m.r_p+rho_p;
    incident_field = compute_incident_field(args->E_TM, args->E_TE, args->theta_i, args->phi_i, 
        args->k, args->eta, r);
    I_p = -1.0*rho_p*incident_field.E;
    return args->b_m.L*(I_m+I_p);
}

complex_t compute_scattered_far_field_E_theta_integrand_2d(const complex_t alpha, const complex_t beta, void *args_){
    assert(args_!=null);
    incident_field_args_2d_t *args=(incident_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t j=complex_t(0.0, 1.0);
    const real_t theta_i=args->theta_i;
    const real_t phi_i=args->phi_i;
    const real_t k=args->k;
    const real_t eta=args->eta;
    vector_t<real_t> theta_i_u=vector_t<real_t>(cos(theta_i)*cos(phi_i), cos(theta_i)*sin(phi_i), -sin(theta_i));
    vector_t<real_t> phi_i_u=vector_t<real_t>(-sin(phi_i), cos(phi_i), 0.0);
    vector_t<real_t> k_i=k*vector_t<real_t>(sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos(theta_i));
    vector_t<real_t> r;
    complex_t I_m=0.0, I_p=0.0;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*(real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]);
    rho_p = +1.0*(real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0]);
    incident_field_t incident_field;
    r = b_m.r_m+rho_m;
    I_m = (k*eta/(4.0*pi))*(+1.0*theta_i_u*rho_m)*exp(+j*(k_i*r));
    r = b_m.r_p+rho_p;
    I_p = (k*eta/(4.0*pi))*(-1.0*theta_i_u*rho_p)*exp(+j*(k_i*r));
    return args->b_m.L*(I_m+I_p);
}

complex_t compute_scattered_far_field_E_phi_integrand_2d(const complex_t alpha, const complex_t beta, void *args_){
    assert(args_!=null);
    incident_field_args_2d_t *args=(incident_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t j=complex_t(0.0, 1.0);
    const real_t theta_i=args->theta_i;
    const real_t phi_i=args->phi_i;
    const real_t k=args->k;
    const real_t eta=args->eta;
    vector_t<real_t> theta_i_u=vector_t<real_t>(cos(theta_i)*cos(phi_i), cos(theta_i)*sin(phi_i), -sin(theta_i));
    vector_t<real_t> phi_i_u=vector_t<real_t>(-sin(phi_i), cos(phi_i), 0.0);
    vector_t<real_t> k_i=k*vector_t<real_t>(sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos(theta_i));
    vector_t<real_t> r;
    complex_t I_m=0.0, I_p=0.0;
    vector_t<real_t> rho_m, rho_p;
    rho_m = +1.0*(real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]);
    rho_p = +1.0*(real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0]);
    incident_field_t incident_field;
    r = b_m.r_m+rho_m;
    I_m = (k*eta/(4.0*pi))*(+1.0*phi_i_u*rho_m)*exp(+j*(k_i*r));
    r = b_m.r_p+rho_p;
    I_p = (k*eta/(4.0*pi))*(-1.0*phi_i_u*rho_p)*exp(+j*(k_i*r));
    return args->b_m.L*(I_m+I_p);
}


