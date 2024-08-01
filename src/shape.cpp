//
#include "shape.hpp"

void shape_t::clear(){
    this->N_0d = 0;
    this->N_1d = 0;
    this->N_2d = 0;
    this->N_3d = 0;
    this->N_basis_1d = 0;
    this->N_basis_2d = 0;
    this->N_basis_3d = 0;
    this->frequency = 0.0;
    this->lambda = 0.0;
    this->mu_b = 0.0;
    this->eps_b = 0.0;
    if (this->is_basis_allocated){
        free(this->basis_1d_list);
        // free(this->basis_2d_list);
        // free(this->basis_3d_list);
        this->is_basis_allocated = false;
    }
}


void shape_t::get_basis_functions(const real_t clmax, const real_t metric_unit){
    assert_error(clmax>0.0, "invalid maximum element size");
    assert_error(metric_unit>0.0, "invalid mertic unit");
    call_gmsh(clmax/metric_unit);
    shape_t::load_mesh(metric_unit);
}

size_t mod_1d(const size_t a){
    return a==1 ? 0 : 1;
}

size_t mod_2d(const size_t a){
    size_t ans=3-(size_t)(a%3);
    return ans==3 ? 0 : ans;
}

size_t mod_3d(const size_t a){
    size_t ans=6-(size_t)(a%6);
    return ans==6 ? 0 : ans;
}

void shape_t::load_mesh(const real_t metric_unit){
    file_t file;
    file.open("mesh/mesh/info.txt", 'r');
    file.read("%zu %zu %zu %zu\n", &this->N_0d, &this->N_1d, &this->N_2d, &this->N_3d);
    file.read("%d\n", &this->is_physical_specified);
    file.close();
    this->basis_1d_list = (basis_1d_t*)calloc(this->N_1d, sizeof(basis_1d_list));
    assert(this->basis_1d_list!=null);
    this->is_basis_allocated = true;
    real_t x, y, z;
    vector_t<real_t> v1, v2, v3, v4, v5;
    int_t pg;
    // 1d bases
    edge_t *edge_list=(edge_t*)calloc(N_1d, sizeof(edge_t));
    assert(edge_list!=null);
    file.open("mesh/mesh/elements_1d.txt", 'r');
    for (size_t i=0; i<this->N_1d; i++){
        pg = -1;
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        if (is_physical_specified){
            file.read("%d", &pg);
        }
        edge_list[i] = edge_t(v1, v2, pg);
    }
    file.close();
    file.open("mesh/basis/basis_1d.txt", 'w');
    for (size_t i=0; i<this->N_1d; i++){
        progress_bar(i, this->N_1d, "creating 1d basis functions...");
        edge_t edge_s=edge_list[i];
        for (size_t j=(i+1); j<this->N_1d; j++){
            edge_t edge_d=edge_list[j];
            size_t new_edge[1];
            size_t index_edge_s=0, index_edge_d=0;
            size_t counter=0;
            for (size_t ii=0; ii<2; ii++){
                for (size_t jj=0; jj<2; jj++){
                    if (is_equal(edge_s.v[ii], edge_d.v[jj], this->mesh_tol*this->lambda)){
                        new_edge[counter] = ii;
                        index_edge_s+=ii;
                        index_edge_d+=jj;
                        counter++;
                    }
                }
            }
            assert_error(counter<2, "invlid mesh");
            if (counter==1){
                index_edge_s = mod_1d(index_edge_s);
                index_edge_d = mod_1d(index_edge_d);
                v1 = edge_s.v[index_edge_s];
                v2 = edge_s.v[new_edge[0]];
                v3 = edge_d.v[index_edge_d];
                file.write("%21.14E %21.14E %21.14E ", v1.x, v1.y, v1.z);
                file.write("%21.14E %21.14E %21.14E ", v2.x, v2.y, v2.z);
                file.write("%21.14E %21.14E %21.14E ", v3.x, v3.y, v3.z);
                if (edge_s.physical_group==edge_d.physical_group){
                    file.write("%d\n", edge_s.physical_group);
                }else{
                    file.write("%d\n", -1);
                }
                this->N_basis_1d++;
            }
        }
    }
    file.close();
    free(edge_list);
    // 2d bases
    triangle_t *triangle_list=(triangle_t*)calloc(N_2d, sizeof(triangle_t));
    assert(triangle_list!=null);
    file.open("mesh/mesh/elements_2d.txt", 'r');
    for (size_t i=0; i<this->N_2d; i++){
        pg = -1;
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v3 = vector_t<real_t>(x, y, z);
        if (is_physical_specified){
            file.read("%d", &pg);
        }
        triangle_list[i] = triangle_t(v1, v2, v3, pg);
    }
    file.close();
    file.open("mesh/basis/basis_2d.txt", 'w');
    for (size_t i=0; i<this->N_2d; i++){
        progress_bar(i, this->N_2d, "creating 2d basis functions...");
        triangle_t triangle_s=triangle_list[i];
        for (size_t j=(i+1); j<this->N_2d; j++){
            triangle_t triangle_d=triangle_list[j];
            size_t new_triangle[2];
            size_t index_triangle_s=0, index_triangle_d=0;
            size_t counter=0;
            for (size_t ii=0; ii<3; ii++){
                for (size_t jj=0; jj<3; jj++){
                    if (is_equal(triangle_s.v[ii], triangle_d.v[jj], this->mesh_tol*this->lambda)){
                        new_triangle[counter] = ii;
                        index_triangle_s+=ii;
                        index_triangle_d+=jj;
                        counter++;
                    }
                }
            }
            assert_error(counter<3, "invlid mesh");
            if (counter==2){
                index_triangle_s = mod_2d(index_triangle_s);
                index_triangle_d = mod_2d(index_triangle_d);
                v1 = triangle_s.v[index_triangle_s];
                v2 = triangle_s.v[new_triangle[0]];
                v3 = triangle_s.v[new_triangle[1]];
                v4 = triangle_d.v[index_triangle_d];
                vector_t<real_t> n;
                n = unit((v2-v1)^(v3-v1));
                if (triangle_s.n*n<0.0){
                    vector_t<real_t> temp=v2;
                    v2 = v3;
                    v3 = temp;
                }
                file.write("%21.14E %21.14E %21.14E ", v1.x, v1.y, v1.z);
                file.write("%21.14E %21.14E %21.14E ", v2.x, v2.y, v2.z);
                file.write("%21.14E %21.14E %21.14E ", v3.x, v3.y, v3.z);
                file.write("%21.14E %21.14E %21.14E ", v4.x, v4.y, v4.z);
                n = unit((v2-v1)^(v3-v1));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                n = unit((v3-v4)^(v2-v4));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                if (triangle_s.physical_group==triangle_d.physical_group){
                    file.write("%d\n", triangle_s.physical_group);
                }else{
                    file.write("%d\n", -1);
                }
                this->N_basis_2d++;
            }
        }
    }
    file.close();
    free(triangle_list);
    // 3d bases
    tetrahedron_t *tetrahedron_list=(tetrahedron_t*)calloc(N_3d, sizeof(tetrahedron_t));
    assert(tetrahedron_list!=null);
    file.open("mesh/mesh/elements_3d.txt", 'r');
    for (size_t i=0; i<this->N_3d; i++){
        pg = -1;
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v3 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->lambda*this->mesh_tol)));
        y = round_m(y, -round(log10(this->lambda*this->mesh_tol)));
        z = round_m(z, -round(log10(this->lambda*this->mesh_tol)));
        v4 = vector_t<real_t>(x, y, z);
        if (is_physical_specified){
            file.read("%d", &pg);
        }
        tetrahedron_list[i] = tetrahedron_t(v1, v2, v3, v4, pg);
    }
    file.close();
    file.open("mesh/basis/basis_3d.txt", 'w');
    for (size_t i=0; i<this->N_3d; i++){
        progress_bar(i, this->N_3d, "creating 3d basis functions...");
        tetrahedron_t tetrahedron_s=tetrahedron_list[i];
        for (size_t j=(i+1); j<this->N_3d; j++){
            tetrahedron_t tetrahedron_d=tetrahedron_list[j];
            size_t new_tetrahedron[3];
            size_t index_tetrahedron_s=0, index_tetrahedron_d=0;
            size_t counter=0;
            for (size_t ii=0; ii<4; ii++){
                for (size_t jj=0; jj<4; jj++){
                    if (is_equal(tetrahedron_s.v[ii], tetrahedron_d.v[jj], this->mesh_tol*this->lambda)){
                        new_tetrahedron[counter] = ii;
                        index_tetrahedron_s+=ii;
                        index_tetrahedron_d+=jj;
                        counter++;
                    }
                }
            }
            assert_error(counter<4, "invlid mesh");
            if (counter==3){
                index_tetrahedron_s = mod_3d(index_tetrahedron_s);
                index_tetrahedron_d = mod_3d(index_tetrahedron_d);
                v1 = tetrahedron_s.v[index_tetrahedron_s];
                v2 = tetrahedron_s.v[new_tetrahedron[0]];
                v3 = tetrahedron_s.v[new_tetrahedron[1]];
                v4 = tetrahedron_s.v[new_tetrahedron[2]];
                v5 = tetrahedron_d.v[index_tetrahedron_d];
                vector_t<real_t> n, n_ref;
                n = unit((v3-v1)^(v2-v1));
                n_ref = ((v1+v2+v3)/3.0-v4);
                if (n_ref*n<0.0){
                    vector_t<real_t> temp=v2;
                    v2 = v3;
                    v3 = temp;
                }
                file.write("%21.14E %21.14E %21.14E ", v1.x, v1.y, v1.z);
                file.write("%21.14E %21.14E %21.14E ", v2.x, v2.y, v2.z);
                file.write("%21.14E %21.14E %21.14E ", v3.x, v3.y, v3.z);
                file.write("%21.14E %21.14E %21.14E ", v4.x, v4.y, v4.z);
                file.write("%21.14E %21.14E %21.14E ", v5.x, v5.y, v5.z);
                n = unit((v3-v1)^(v2-v1));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                n = unit((v2-v1)^(v4-v1));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                n = unit((v4-v1)^(v3-v1));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                n = unit((v2-v5)^(v3-v5));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                n = unit((v4-v5)^(v2-v5));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                n = unit((v3-v5)^(v4-v5));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                n = unit((v4-v3)^(v2-v3));
                file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                if (tetrahedron_s.physical_group==tetrahedron_d.physical_group){
                    file.write("%d\n", tetrahedron_s.physical_group);
                }else{
                    file.write("%d\n", -1);
                }
                this->N_basis_3d++;
            }
        }
    }
    file.close();
    free(tetrahedron_list);
    //
    print("total number of 1d basis fuctions: %d\n", this->N_basis_1d);
    print("total number of 2d basis fuctions: %d\n", this->N_basis_2d);
    print("total number of 3d basis fuctions: %d\n", this->N_basis_3d);
}

//
void call_gmsh(const real_t tol){
    int_t max_length=200;
    char *cmd=(char*)calloc(max_length, sizeof(char));
    print("calling gmsh...");
    sprintf(cmd, "gmsh mesh/shape.geo -3 -clmax %0.4f -format vtk -save_all -o mesh/shape.vtk > mesh/shape_log.txt", tol);
    assert_error(!system(cmd), "unable to mesh geometry");
    print(", done!\n");
    #ifdef __windows__
    sprintf(cmd, "python mesh/read_vtk.py");
    #endif
    #ifdef __linux__
    sprintf(cmd, "python3 mesh/read_vtk.py");
    #endif
    assert_error(!system(cmd), "unable to generate mesh");
    free(cmd);
}

void create_vertical_wire_dipole(const real_t length, const real_t port_length){
    assert_error(length>0, "invalid legnth");
    assert_error(port_length<=length, "invalid port legnth");
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("Point(1) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, -length/2);
    file.write("Point(2) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, -port_length/2);
    file.write("Point(3) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, +port_length/2);
    file.write("Point(4) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, +length/2);
    file.write("Line(1) = {1, 2};\n");
    file.write("Line(2) = {2, 3};\n");
    file.write("Line(3) = {3, 4};\n");
    file.write("Physical Curve(\"Port\", 1) = {2};\n");
    file.write("Physical Curve(\"Wire\", 2) = {1, 3};\n");
    file.close();
}

void create_sphere(const real_t radius){
    assert_error(radius>0, "invalid radius");
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("SetFactory(\"OpenCASCADE\");\n");
    file.write("Sphere(1) = {%21.14E, %21.14E, %21.14E, %21.14E, -Pi/2, Pi/2, 2*Pi};\n", 0.0, 0.0, 0.0, radius);
    file.write("Physical Surface(\"Surface\", 1) = {1};\n");
    file.write("Physical Volume(\"Volume\", 1) = {1};\n");
    file.close();
}