//
#include "engine.hpp"

// #define only_1d
#define only_2d
// #define only_3d

void engine_t::set(const real_t freq, const complex_t mu_b, const complex_t eps_b, 
    const real_t clmax, const real_t unit_metric, const real_t a, const size_t N_ports){
    if (this->is_engine_set){engine_t::unset();}
    // check inputs errors
    assert_error(freq>0.0, "invalid frequency");
    this->freq = freq;
    this->a = a*unit_metric;
    this->mu_b = mu_b;
    this->eps_b = eps_b;
    this->lambda = (c_0/real(sqrt(mu_b*eps_b)))/freq;
    this->k_b = 2.0*pi*freq*sqrt(mu_0*eps_0)*sqrt(mu_b*eps_b);
    this->eta_b = sqrt(mu_0/eps_0)*sqrt(mu_b/eps_b);
    this->shape.set(freq, mu_b, eps_b);
    this->shape.get_basis_functions(clmax, unit_metric);
    this->shape.get_info(this->N_basis_1d, this->N_basis_2d, this->N_basis_3d);
    //
    this->N = 0;
    #ifdef only_1d
    this->N += this->N_basis_1d;
    #endif
    #ifdef only_2d
    this->N += this->N_basis_2d;
    #endif
    #ifdef only_3d
    this->N += this->N_basis_3d;
    #endif
    assert_error(this->N<max_system_size, "too large linear system");
    this->Z_mn.set(this->N, this->N);
    this->V_m.set(this->N, 1);
    this->I_n.set(this->N, 1);
    this->is_Z_mn_allocated = true;
    this->is_V_m_allocated = true;
    this->is_I_n_allocated = true;
    this->quadl.set_1d(this->k_max_1d, this->tol_1d);
    this->quadl.set_2d(this->k_max_2d, this->tol_2d);
    this->quadl.set_3d(this->k_max_3d, this->tol_3d);
    //
    this->N_ports = N_ports;
    this->port_list = (port_t*)calloc(this->N_ports, sizeof(port_t));
    assert(this->port_list!=null);
    this->is_port_list_allocated = true;
    //
    this->is_engine_set = true;
}

void engine_t::unset(){
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
    if (this->is_port_list_allocated){
        free(this->port_list);
    }
    shape.clear();
}

void engine_t::assign_port(const size_t index, const complex_t V, const complex_t Z, const int_t pg,
    const vector_t<real_t> p, const real_t L, const real_t W){
    assert_error(index<this->N_ports, "port index out of range");
    assert_error(this->is_port_list_allocated, "engine is not set yet");
    port_t port;
    port.index = index;
    port.V = V;
    port.Z = Z;
    port.pg = pg;
    port.p = p;
    // assert_error(L>0.0&&W>=0.0, "invalid port description");
    port.L = L;
    port.W = W;
    this->port_list[index] = port;
}

void engine_t::compute_Z_mn(){
    int_t flag;
    basis_2d_t basis_m, basis_n;
    char *msg=(char*)calloc(this->max_line_length, sizeof(char));
    size_t count=0;
    complex_t Z=0.0;
    // Zmn_LL
    #ifdef only_1d
    {   
        size_t i=this->N_basis_3d+this->N_basis_2d;
        size_t j=this->N_basis_3d+this->N_basis_2d;
        basis_1d_t b_m;
        basis_1d_t b_n;
        for (size_t m=0; m<this->N_basis_1d; m++){
            b_m = this->shape.get_basis_1d(m);
            for (size_t n=m; n<this->N_basis_1d; n++){
                b_n = this->shape.get_basis_1d(n);
                sprintf(msg, "LL: Z_mn (%zu, %zu)", m, n);
                progress_bar(count, this->N_basis_1d*(this->N_basis_1d+1)/2, msg);
                Z = Z_mn_1d_1d(b_m, b_n, k_b, eta_b, lambda, a, quadl, flag);
                if (m==n){
                    for (size_t k=0; k<this->N_ports; k++){
                        if ((b_m.pg_m==b_m.pg_p)&&(b_m.pg_m==this->port_list[k].pg)){
                            Z+=this->port_list[k].Z;
                        }
                    }
                }
                assert_error(flag==false, "no convergence");
                assert_error(!isinf(abs(Z)), "inf value for Z_mn");
                assert_error(!isnan(abs(Z)), "nan value for Z_mn");
                this->Z_mn(i+m, j+n) = Z;
                count++;
            }
            for (size_t n=m+1; n<N; n++){
                this->Z_mn(j+n, i+m) = this->Z_mn(i+m, j+n);
            }
        }
    }
    #endif
    // Zmn_SS
    #ifdef only_2d
    {   
        size_t i=this->N_basis_3d;
        size_t j=this->N_basis_3d;
        basis_2d_t b_m;
        basis_2d_t b_n;
        for (size_t m=0; m<this->N_basis_2d; m++){
            b_m = this->shape.get_basis_2d(m);
            for (size_t n=m; n<this->N_basis_2d; n++){
                b_n = this->shape.get_basis_2d(n);
                sprintf(msg, "SS: Z_mn (%zu, %zu)", m, n);
                progress_bar(count, this->N_basis_2d*(this->N_basis_2d+1)/2, msg);
                Z = Z_mn_2d_2d(b_m, b_n, k_b, eta_b, lambda, quadl, flag);
                // if (m==n){
                //     for (size_t k=0; k<this->N_ports; k++){
                //         if ((b_m.pg_m==b_m.pg_p)&&(b_m.pg_m==this->port_list[k].pg)){
                //             Z+=this->port_list[k].Z;
                //         }
                //     }
                // }
                assert_error(flag==false, "no convergence");
                assert_error(!isinf(abs(Z)), "inf value for Z_mn");
                assert_error(!isnan(abs(Z)), "nan value for Z_mn");
                this->Z_mn(i+m, j+n) = Z;
                count++;
            }
            for (size_t n=m+1; n<N; n++){
                this->Z_mn(j+n, i+m) = this->Z_mn(i+m, j+n);
            }
        }
    }
    #endif
    //
    free(msg);
    engine_t::save_Z_mn("data/Z_mn.bin");
    this->is_Z_mn_calculated = true;
}

void engine_t::compute_I_n(){
    assert_error(this->is_engine_set, "egnine is not set yet");
    print("computing I_n...");
    engine_t::load_Z_mn("data/Z_mn.bin");
    this->Z_mn.lup();
    this->Z_mn.solve(this->V_m, this->I_n);
    this->is_I_n_calculated = true;
    // print(", done!\n");
}

complex_t engine_t::compute_Z_in(const size_t port_index){
    assert_error(this->is_engine_set, "egnine is not set yet");
    complex_t Z_in=0.0;
    complex_t V=0.0;
    complex_t I=0.0;
    // 1d
    {   
        basis_1d_t b_m;
        size_t i=this->N_basis_3d+this->N_basis_2d;
        for (size_t m=0; m<this->N_basis_1d; m++){
            b_m = this->shape.get_basis_1d(m);
            if ((b_m.pg_m==b_m.pg_p)&&(b_m.pg_m==this->port_list[port_index].pg)){
                V = this->V_m(i+m, 0);
                I = this->I_n(i+m, 0);
                Z_in = V/I-this->port_list[port_index].Z;
                return Z_in;
            }
        }
    }
    return Z_in;
}

complex_t engine_t::compute_S_mutual(const size_t port_index){
    assert_error(this->is_engine_set, "egnine is not set yet");
    complex_t V=0.0;
    // 1d
    {   
        basis_1d_t b_m;
        size_t i=this->N_basis_3d+this->N_basis_2d;
        for (size_t m=0; m<this->N_basis_1d; m++){
            b_m = this->shape.get_basis_1d(m);
            if ((b_m.pg_m==b_m.pg_p)&&(b_m.pg_m==this->port_list[port_index].pg)){
                V = -2.0*this->I_n(i+m, 0)*this->port_list[port_index].Z;
                V = V*(this->port_list[port_index].p*((unit(b_m.L_m[0])-unit(b_m.L_p[0]))/2.0));
                return V;
            }
        }
    }
    return V;
}

void engine_t::compute_S_matrix(matrix_t<complex_t> &S_matrix, const complex_t Z_0){
    assert_error(this->is_engine_set, "egnine is not set yet");
    // 1d
    {
        for (size_t n=0; n<this->N_ports; n++){
            for (size_t k=0; k<this->N_ports; k++){
                if (k==n){
                    engine_t::assign_port(k, +1.0, Z_0, 
                    this->port_list[k].pg, this->port_list[k].p, 
                    this->port_list[k].L, this->port_list[k].W);
                }else{
                    engine_t::assign_port(k, +0.0, Z_0, 
                    this->port_list[k].pg, this->port_list[k].p, 
                    this->port_list[k].L, this->port_list[k].W);
                }
            }
            engine_t::compute_V_m_ports();
            engine_t::compute_I_n();
            for (size_t m=0; m<this->N_ports; m++){
                if (m==n){
                    complex_t Z=engine_t::compute_Z_in(m);
                    S_matrix(m, n) = (Z-Z_0)/(Z+Z_0);
                }else{
                    S_matrix(m, n) = engine_t::compute_S_mutual(m);
                }
            }
        }  
    }
}

void engine_t::compute_V_m_ports(){
    assert_error(this->is_engine_set, "egnine is not set yet");
    print("computing V_m...");
    // Vm_L
    #ifdef only_1d
    {
        basis_1d_t b_m;
        size_t i=this->N_basis_3d+this->N_basis_2d;
        for (size_t m=0; m<this->N_basis_1d; m++){
            b_m = this->shape.get_basis_1d(m);
            this->V_m(i+m, 0) = 0.0;
            for (size_t k=0; k<this->N_ports; k++){
                if ((b_m.pg_m==b_m.pg_p)&&(b_m.pg_m==this->port_list[k].pg)){
                    const complex_t V=this->port_list[k].V;
                    const vector_t<real_t> p=this->port_list[k].p;
                    this->V_m(i+m, 0) = (unit(b_m.L_m[0])-unit(b_m.L_p[0]))*p*V/2.0;
                }
            }
        }
    }
    #endif
    // Vm_S
    #ifdef only_2d
    {
        size_t i=this->N_basis_3d;
        for (size_t m=0; m<this->N_basis_2d; m++){
            this->V_m(i+m, 0) = 0.0;
        }
    }
    #endif
    print(", done!\n");
}

void engine_t::compute_V_m_incident(const complex_t E_TM, const complex_t E_TE, const real_t theta_i, const real_t phi_i){
    assert_error(this->is_engine_set, "egnine is not set yet");
    print("computing V_m...");
    // Vm_L
    #ifdef only_1d
    {
        basis_1d_t b_m;
        size_t i=this->N_basis_3d+this->N_basis_2d;
        for (size_t m=0; m<this->N_basis_1d; m++){
            b_m = this->shape.get_basis_1d(m);
            line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
            incident_field_args_1d_t args;
            args.b_m = b_m;
            args.E_TM = E_TM;
            args.E_TE = E_TE;
            args.theta_i = theta_i;
            args.phi_i = phi_i;
            args.k = real(this->k_b);
            args.eta = real(this->eta_b);
            int_t flag;
            this->V_m(i+m, 0) =  this->quadl.integral_1d(compute_incident_E_integrand_1d, &args, line, flag);
        }
    }
    #endif
    // Vm_S
    #ifdef only_2d
    {
        basis_2d_t b_m;
        size_t i=this->N_basis_3d;
        for (size_t m=0; m<this->N_basis_2d; m++){
            b_m = this->shape.get_basis_2d(m);
            triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
            incident_field_args_2d_t args;
            args.b_m = b_m;
            args.E_TM = E_TM;
            args.E_TE = E_TE;
            args.theta_i = theta_i;
            args.phi_i = phi_i;
            args.k = real(this->k_b);
            args.eta = real(this->eta_b);
            int_t flag;
            this->V_m(i+m, 0) = this->quadl.integral_2d(compute_incident_E_integrand_2d, &args, triangle, flag);
        }
    }
    #endif
    // Vm_V
    #ifdef only_3d
    {
        size_t i=0;
        for (size_t m=0; m<this->N_basis_3d; m++){
            this->V_m(i+m, 0) = 0.0;
        }
    }
    #endif
    print(", done!\n");
}

void engine_t::save_Z_mn(const char *filename){
    binary_file_t file;
    assert(filename!=null);
    file.open(filename, 'w');
    file.write(&this->N);
    // print("saving Z_mn solutions...");
    for (size_t m=0; m<this->N; m++){
        for (size_t n=0; n<this->N; n++){
            file.write(&this->Z_mn(m, n));
        }
    }
    file.close();
    // print(", done\n");
}

void engine_t::load_Z_mn(const char *filename){
    binary_file_t file;
    assert(filename!=null);
    file.open(filename, 'r');
    file.read(&this->N);
    this->Z_mn.unset();
    this->Z_mn.set(this->N, this->N);
    // print("loading Z_mn solutions...");
    for (size_t m=0; m<this->N; m++){
        for (size_t n=0; n<this->N; n++){
            file.read(&this->Z_mn(m, n));
        }
    }
    file.close();
    // print(", done\n");
    this->is_Z_mn_calculated = true;
}

void engine_t::export_solutions(){
    file_t file;
    file.open("data/Z_mn.txt", 'w');
    for (size_t m=0; m<this->N; m++){
        for (size_t n=0; n<this->N; n++){
            file.write("%zu, %zu: %21.14E, %21.14E\n", m, n, real(this->Z_mn(m, n)), imag(this->Z_mn(m, n)));
        }
        file.write("\n");
    }
    file.close();
    file.open("data/V_m.txt", 'w');
    for (size_t m=0; m<this->N; m++){
        file.write("%zu: %21.14E, %21.14E\n", m, real(this->V_m(m, 0)), imag(this->V_m(m, 0)));
    }
    file.close();
    file.open("data/I_n.txt", 'w');
    for (size_t n=0; n<this->N; n++){
        file.write("%zu: %21.14E, %21.14E\n", n, real(this->I_n(n, 0)), imag(this->I_n(n, 0)));
    }
    file.close();
}

//

sigma_t engine_t::compute_RCS(const real_t theta_i, const real_t phi_i){
    sigma_t sigma;
    complex_t sum_theta=0.0, sum_phi=0.0;
    int_t flag;
    // 1d
    #ifdef only_1d
    incident_field_args_1d_t args;
    args.theta_i = theta_i;
    args.phi_i = phi_i;
    args.k = real(this->k_b);
    args.eta = real(this->eta_b);
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
    for (size_t m=0; m<this->N_basis_1d; m++){
        basis_1d_t b_m=this->shape.get_basis_1d(m);
        args.b_m = b_m;
        sum_theta+=this->quadl.integral_1d(compute_scattered_far_field_E_theta_integrand_1d, &args, line, flag)*this->I_n(m, 0);
        sum_phi+=this->quadl.integral_1d(compute_scattered_far_field_E_phi_integrand_1d, &args, line, flag)*this->I_n(m, 0);
    }
    #endif
    // 2d
    #ifdef only_2d
    incident_field_args_2d_t args;
    args.theta_i = theta_i;
    args.phi_i = phi_i;
    args.k = real(this->k_b);
    args.eta = real(this->eta_b);
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    for (size_t m=0; m<this->N_basis_2d; m++){
        basis_2d_t b_m=this->shape.get_basis_2d(m);
        args.b_m = b_m;
        sum_theta+=this->quadl.integral_2d(compute_scattered_far_field_E_theta_integrand_2d, &args, triangle, flag)*this->I_n(m, 0);
        sum_phi+=this->quadl.integral_2d(compute_scattered_far_field_E_phi_integrand_2d, &args, triangle, flag)*this->I_n(m, 0);
    }
    #endif
    sigma.theta = 4.0*pi*abs(sum_theta)*abs(sum_theta);
    sigma.phi = 4.0*pi*abs(sum_phi)*abs(sum_phi);
    return sigma;
}

far_field_t engine_t::compute_far_field(const real_t theta_i, const real_t phi_i){
    far_field_t far_field;
    incident_field_args_1d_t args;
    args.theta_i = theta_i;
    args.phi_i = phi_i;
    args.k = real(this->k_b);
    args.eta = real(this->eta_b);
    complex_t sum_theta=0.0, sum_phi=0.0;
    // 1d
    int_t flag;
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
    for (size_t m=0; m<this->N_basis_1d; m++){
        basis_1d_t b_m=this->shape.get_basis_1d(m);
        args.b_m = b_m;
        sum_theta+=this->quadl.integral_1d(compute_scattered_far_field_E_theta_integrand_1d, &args, line, flag)*this->I_n(m, 0);
        sum_phi+=this->quadl.integral_1d(compute_scattered_far_field_E_phi_integrand_1d, &args, line, flag)*this->I_n(m, 0);
    }
    far_field.theta = sum_theta;
    far_field.phi = sum_phi;
    return far_field;
}

near_field_t engine_t::compute_near_field_E(const vector_t<real_t> r){
    near_field_t near_field;
    vector_t<real_t> x=vector_t<real_t>(1.0, 0.0, 0.0);
    vector_t<real_t> y=vector_t<real_t>(0.0, 1.0, 0.0);
    vector_t<real_t> z=vector_t<real_t>(0.0, 0.0, 1.0);
    near_field.x = 0.0;
    near_field.y = 0.0;
    near_field.z = 0.0;
    // 1d
    #ifdef only_1d
    for (size_t m=0; m<this->N_basis_1d; m++){
        basis_1d_t b_m=this->shape.get_basis_1d(m);
        near_field.x+= compute_E_1d(b_m, r, x, k_b, eta_b, a, quadl)*this->I_n(m, 0);
        near_field.y+= compute_E_1d(b_m, r, y, k_b, eta_b, a, quadl)*this->I_n(m, 0);
        near_field.z+= compute_E_1d(b_m, r, z, k_b, eta_b, a, quadl)*this->I_n(m, 0);
    }
    #endif
    // 2d
    #ifdef only_2d
    for (size_t m=0; m<this->N_basis_2d; m++){
        basis_2d_t b_m=this->shape.get_basis_2d(m);
        near_field.x+= compute_E_2d(b_m, r, x, k_b, eta_b, quadl)*this->I_n(m, 0);
        near_field.y+= compute_E_2d(b_m, r, y, k_b, eta_b, quadl)*this->I_n(m, 0);
        near_field.z+= compute_E_2d(b_m, r, z, k_b, eta_b, quadl)*this->I_n(m, 0);
    }
    #endif
    return near_field;
}

near_field_t engine_t::compute_near_field_H(const vector_t<real_t> r){
    near_field_t near_field;
    vector_t<real_t> x=vector_t<real_t>(1.0, 0.0, 0.0);
    vector_t<real_t> y=vector_t<real_t>(0.0, 1.0, 0.0);
    vector_t<real_t> z=vector_t<real_t>(0.0, 0.0, 1.0);
    near_field.x = 0.0;
    near_field.y = 0.0;
    near_field.z = 0.0;
    // 1d
    #ifdef only_1d
    for (size_t m=0; m<this->N_basis_1d; m++){
        basis_1d_t b_m=this->shape.get_basis_1d(m);
        near_field.x+= compute_H_1d(b_m, r, x, k_b, eta_b, a, quadl)*this->I_n(m, 0);
        near_field.y+= compute_H_1d(b_m, r, y, k_b, eta_b, a, quadl)*this->I_n(m, 0);
        near_field.z+= compute_H_1d(b_m, r, z, k_b, eta_b, a, quadl)*this->I_n(m, 0);
    }
    #endif
    // 2d
    #ifdef only_2d
    for (size_t m=0; m<this->N_basis_2d; m++){
        basis_2d_t b_m=this->shape.get_basis_2d(m);
        near_field.x+= compute_H_2d(b_m, r, x, k_b, eta_b, quadl)*this->I_n(m, 0);
        near_field.y+= compute_H_2d(b_m, r, y, k_b, eta_b, quadl)*this->I_n(m, 0);
        near_field.z+= compute_H_2d(b_m, r, z, k_b, eta_b, quadl)*this->I_n(m, 0);
    }
    #endif
    return near_field;
}

//

void engine_t::export_currents(const char *filename){
    file_t file;
    binary_file_t binary_file;
    size_t N_0d, N_1d, N_2d, N_3d;
    int_t is_physical_specified;
    const real_t metric_unit=this->shape.get_metric_unit();
    const real_t mesh_tol=this->shape.get_mesh_tol();
    file.open("mesh/mesh/info.txt", 'r');
    file.read("%zu %zu %zu %zu\n", &N_0d, &N_1d, &N_2d, &N_3d);
    file.read("%d\n", &is_physical_specified);
    file.close();
    real_t x, y, z;
    vector_t<real_t> v1, v2, v3, v4, v5;
    int_t pg;
    // 1d bases
    #ifdef only_1d
    line_t *line_list=(line_t*)calloc(N_1d, sizeof(line_t));
    assert(line_list!=null);
    file.open("mesh/mesh/elements_1d.txt", 'r');
    for (size_t i=0; i<N_1d; i++){
        pg = -1;
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(mesh_tol)));
        y = round_m(y, -round(log10(mesh_tol)));
        z = round_m(z, -round(log10(mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(mesh_tol)));
        y = round_m(y, -round(log10(mesh_tol)));
        z = round_m(z, -round(log10(mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        file.read("%d", &pg);
        line_list[i] = line_t(v1, v2, pg);
    }
    file.close();
    //
    file.open(filename, 'w');
    for (size_t i=0; i<N_1d; i++){
        progress_bar(i, N_1d, "exporting currents");
        line_t line=line_list[i];
        vector_t<complex_t> J1=vector_t<complex_t>(0.0, 0.0, 0.0);
        vector_t<complex_t> J2=vector_t<complex_t>(0.0, 0.0, 0.0);
        for (size_t j=0; j<this->N_basis_1d; j++){
            basis_1d_t b=this->shape.get_basis_1d(j);
            // Scenario 1
            if (is_equal(line.v[0], b.e[0], mesh_tol*this->lambda)){
                J1 = J1 + this->I_n(j, 0)*unit(b.L_m[0]);
            }
            // Scenario 2
            if (is_equal(line.v[1], b.e[0], mesh_tol*this->lambda)){
                J2 = J2 + this->I_n(j, 0)*unit(b.L_m[0]);
            }
        }
        real_t I1=mag(J1);
        real_t I2=mag(J2);
        file.write("%21.14E %21.14E %21.14E %21.14E\n", line.v[0].x, line.v[0].y, line.v[0].z, I1);
        file.write("%21.14E %21.14E %21.14E %21.14E\n", line.v[1].x, line.v[1].y, line.v[1].z, I2);
    }
    file.close();
    free(line_list);
    #endif
    // 2d bases
    #ifdef only_2d
    triangle_t *triangle_list=(triangle_t*)calloc(N_2d, sizeof(triangle_t));
    assert(triangle_list!=null);
    file.open("mesh/mesh/elements_2d.txt", 'r');
    for (size_t i=0; i<N_2d; i++){
        pg = -1;
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(mesh_tol)));
        y = round_m(y, -round(log10(mesh_tol)));
        z = round_m(z, -round(log10(mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(mesh_tol)));
        y = round_m(y, -round(log10(mesh_tol)));
        z = round_m(z, -round(log10(mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(mesh_tol)));
        y = round_m(y, -round(log10(mesh_tol)));
        z = round_m(z, -round(log10(mesh_tol)));
        v3 = vector_t<real_t>(x, y, z);
        file.read("%d", &pg);
        triangle_list[i] = triangle_t(v1, v2, v3, pg);
    }
    file.close();
    //    
    file.open(filename, 'w');
    file.write("View \"Current\" {\n");
    for (size_t i=0; i<N_2d; i++){
        progress_bar(i, N_2d, "exporting currents");
        triangle_t triangle=triangle_list[i];
        vector_t<complex_t> J_c=vector_t<complex_t>(0.0, 0.0, 0.0);
        for (size_t j=0; j<this->N_basis_2d; j++){
            basis_2d_t b=this->shape.get_basis_2d(j);
            // Scenario 1
            if (is_equal(triangle.v[0], b.r_m, mesh_tol*this->lambda) &&
                is_equal(triangle.v[1], b.e[0], mesh_tol*this->lambda) &&
                is_equal(triangle.v[2], b.e[1], mesh_tol*this->lambda)){
                J_c = J_c + this->I_n(j, 0)*(1.0/3.0)*(b.L_m[0]+b.L_m[1])*b.L/(2.0*b.A_m[0]);
            }
            // Scenario 2
            if (is_equal(triangle.v[0], b.r_p, mesh_tol*this->lambda) &&
                is_equal(triangle.v[1], b.e[1], mesh_tol*this->lambda) &&
                is_equal(triangle.v[2], b.e[0], mesh_tol*this->lambda)){
                J_c = J_c - this->I_n(j, 0)*(1.0/3.0)*(b.L_p[0]+b.L_p[1])*b.L/(2.0*b.A_p[0]);
            }
            // Scenario 3
            if (is_equal(triangle.v[0], b.e[1], mesh_tol*this->lambda) &&
                is_equal(triangle.v[1], b.r_m, mesh_tol*this->lambda) &&
                is_equal(triangle.v[2], b.e[0], mesh_tol*this->lambda)){
                J_c = J_c + this->I_n(j, 0)*(1.0/3.0)*(b.L_m[0]+b.L_m[1])*b.L/(2.0*b.A_m[0]);
            }
            // Scenario 4
            if (is_equal(triangle.v[0], b.e[0], mesh_tol*this->lambda) &&
                is_equal(triangle.v[1], b.r_p, mesh_tol*this->lambda) &&
                is_equal(triangle.v[2], b.e[1], mesh_tol*this->lambda)){
                J_c = J_c - this->I_n(j, 0)*(1.0/3.0)*(b.L_p[0]+b.L_p[1])*b.L/(2.0*b.A_p[0]);
            }
            // Scenario 5
            if (is_equal(triangle.v[0], b.e[0], mesh_tol*this->lambda) &&
                is_equal(triangle.v[1], b.e[1], mesh_tol*this->lambda) &&
                is_equal(triangle.v[2], b.r_m, mesh_tol*this->lambda)){
                J_c = J_c + this->I_n(j, 0)*(1.0/3.0)*(b.L_m[0]+b.L_m[1])*b.L/(2.0*b.A_m[0]);
            }
            // Scenario 6
            if (is_equal(triangle.v[0], b.e[1], mesh_tol*this->lambda) &&
                is_equal(triangle.v[1], b.e[0], mesh_tol*this->lambda) &&
                is_equal(triangle.v[2], b.r_p, mesh_tol*this->lambda)){
                J_c = J_c - this->I_n(j, 0)*(1.0/3.0)*(b.L_p[0]+b.L_p[1])*b.L/(2.0*b.A_p[0]);
            }
        }
        real_t I_c=mag(J_c);
        // 
        file.write("ST(%21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E){%21.14E, %21.14E, %21.14E};\n", 
            triangle.v[0].x, triangle.v[0].y, triangle.v[0].z, 
            triangle.v[1].x, triangle.v[1].y, triangle.v[1].z, 
            triangle.v[2].x, triangle.v[2].y, triangle.v[2].z, 
            I_c, I_c, I_c);
    }
    file.write("};\n");
    file.close();
    free(triangle_list);
    #endif
}