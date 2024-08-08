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
    this->freq = 0.0;
    this->lambda = 0.0;
    this->mu_b = 0.0;
    this->eps_b = 0.0;
    if (this->is_basis_allocated){
        free(this->basis_1d_list);
        free(this->basis_2d_list);
        free(this->basis_3d_list);
        this->is_basis_allocated = false;
    }
}


void shape_t::get_basis_functions(const real_t clmax, const real_t metric_unit){
    assert_error(clmax>0.0, "invalid maximum element size");
    assert_error(metric_unit>0.0, "invalid mertic unit");
    call_gmsh(clmax);
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
    binary_file_t binary_file;
    file.open("mesh/mesh/info.txt", 'r');
    file.read("%zu %zu %zu %zu\n", &this->N_0d, &this->N_1d, &this->N_2d, &this->N_3d);
    file.read("%d\n", &this->is_physical_specified);
    file.close();
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
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        file.read("%d", &pg);
        edge_list[i] = edge_t(v1, v2, pg);
    }
    file.close();
    file.open("mesh/basis/basis_1d.txt", 'w');
    binary_file.open("mesh/basis/basis_1d.bin", 'w');
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
                if ((this->is_physical_specified&&edge_s.physical_group>0&&edge_d.physical_group>0) || !this->is_physical_specified){
                    index_edge_s = mod_1d(index_edge_s);
                    index_edge_d = mod_1d(index_edge_d);
                    v1 = edge_s.v[index_edge_s];
                    v2 = edge_s.v[new_edge[0]];
                    v3 = edge_d.v[index_edge_d];
                    file.write("%21.14E %21.14E %21.14E ", v1.x, v1.y, v1.z);
                    file.write("%21.14E %21.14E %21.14E ", v2.x, v2.y, v2.z);
                    file.write("%21.14E %21.14E %21.14E ", v3.x, v3.y, v3.z);
                    binary_file.write(&v1);
                    binary_file.write(&v2);
                    binary_file.write(&v3);
                    file.write("%d %d\n", edge_s.physical_group, edge_d.physical_group);
                    binary_file.write(&edge_s.physical_group);
                    binary_file.write(&edge_d.physical_group);
                    this->N_basis_1d++;
                }
            }
        }
    }
    binary_file.close();
    file.close();
    free(edge_list);
    // 2d bases
    triangle_t *triangle_list=(triangle_t*)calloc(N_2d, sizeof(triangle_t));
    assert(triangle_list!=null);
    file.open("mesh/mesh/elements_2d.txt", 'r');
    for (size_t i=0; i<this->N_2d; i++){
        pg = -1;
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v3 = vector_t<real_t>(x, y, z);
        file.read("%d", &pg);
        triangle_list[i] = triangle_t(v1, v2, v3, pg);
    }
    file.close();
    file.open("mesh/basis/basis_2d.txt", 'w');
    binary_file.open("mesh/basis/basis_2d.bin", 'w');
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
                if ((this->is_physical_specified&&triangle_s.physical_group>0&&triangle_d.physical_group>0) || !this->is_physical_specified){
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
                    binary_file.write(&v1);
                    binary_file.write(&v2);
                    binary_file.write(&v3);
                    binary_file.write(&v4);
                    n = unit((v2-v1)^(v3-v1));
                    file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                    n = unit((v3-v4)^(v2-v4));
                    file.write("%21.14E %21.14E %21.14E ", n.x, n.y, n.z);
                    file.write("%d %d\n", triangle_s.physical_group, triangle_d.physical_group);
                    binary_file.write(&triangle_s.physical_group);
                    binary_file.write(&triangle_d.physical_group);
                    this->N_basis_2d++;
                }
            }
        }
    }
    binary_file.close();
    file.close();
    free(triangle_list);
    // 3d bases
    tetrahedron_t *tetrahedron_list=(tetrahedron_t*)calloc(N_3d, sizeof(tetrahedron_t));
    assert(tetrahedron_list!=null);
    file.open("mesh/mesh/elements_3d.txt", 'r');
    for (size_t i=0; i<this->N_3d; i++){
        pg = -1;
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v1 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v2 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v3 = vector_t<real_t>(x, y, z);
        file.read("%lf %lf %lf", &x, &y, &z); x*=metric_unit; y*=metric_unit; z*=metric_unit;
        x = round_m(x, -round(log10(this->mesh_tol)));
        y = round_m(y, -round(log10(this->mesh_tol)));
        z = round_m(z, -round(log10(this->mesh_tol)));
        v4 = vector_t<real_t>(x, y, z);
        file.read("%d", &pg);
        tetrahedron_list[i] = tetrahedron_t(v1, v2, v3, v4, pg);
    }
    file.close();
    file.open("mesh/basis/basis_3d.txt", 'w');
    binary_file.open("mesh/basis/basis_3d.bin", 'w');
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
                if ((this->is_physical_specified&&tetrahedron_s.physical_group>0&&tetrahedron_d.physical_group>0) || !this->is_physical_specified){
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
                    binary_file.write(&v1);
                    binary_file.write(&v2);
                    binary_file.write(&v3);
                    binary_file.write(&v4);
                    binary_file.write(&v5);
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
                    file.write("%d %d\n", tetrahedron_s.physical_group, tetrahedron_d.physical_group);
                    binary_file.write(&tetrahedron_s.physical_group);
                    binary_file.write(&tetrahedron_d.physical_group);
                    this->N_basis_3d++;
                }
            }
        }
    }
    binary_file.close();
    file.close();
    free(tetrahedron_list);
    //
    print("total number of 1d basis fuctions: %d\n", this->N_basis_1d);
    print("total number of 2d basis fuctions: %d\n", this->N_basis_2d);
    print("total number of 3d basis fuctions: %d\n", this->N_basis_3d);
    shape_t::load_basis_functions();
}

void shape_t::load_basis_functions(){
    this->basis_1d_list = (basis_1d_t*)calloc(this->N_basis_1d, sizeof(basis_1d_t)); assert(this->basis_1d_list!=null); 
    this->basis_2d_list = (basis_2d_t*)calloc(this->N_basis_2d, sizeof(basis_2d_t)); assert(this->basis_2d_list!=null);
    this->basis_3d_list = (basis_3d_t*)calloc(this->N_basis_3d, sizeof(basis_3d_t)); assert(this->basis_3d_list!=null);
    this->is_basis_allocated = true;
    vector_t<real_t> v1, v2, v3, v4, v5;
    int_t pg_m, pg_p;
    binary_file_t file;
    file.open("mesh/basis/basis_1d.bin", 'r');
    for (size_t i=0; i<this->N_basis_1d; i++){  
        file.read(&v1);
        file.read(&v2);
        file.read(&v3);
        file.read(&pg_m);
        file.read(&pg_p);
        this->basis_1d_list[i] = basis_1d_t(v1, v2, v3, pg_m, pg_p);
    }
    file.close();
    file.open("mesh/basis/basis_2d.bin", 'r');
    for (size_t i=0; i<this->N_basis_2d; i++){  
        file.read(&v1);
        file.read(&v2);
        file.read(&v3);
        file.read(&v4);
        file.read(&pg_m);
        file.read(&pg_p);
        this->basis_2d_list[i] = basis_2d_t(v1, v2, v3, v4, pg_m, pg_p);
    }
    file.close();
    file.open("mesh/basis/basis_3d.bin", 'r');
    for (size_t i=0; i<this->N_basis_3d; i++){  
        file.read(&v1);
        file.read(&v2);
        file.read(&v3);
        file.read(&v4);
        file.read(&v5);
        file.read(&pg_m);
        file.read(&pg_p);
        this->basis_3d_list[i] = basis_3d_t(v1, v2, v3, v4, v5, pg_m, pg_p);
    }
    file.close();
}

//
void call_gmsh(const real_t tol){
    assert_error(tol>0.0, "invalid mesh tolerance");
    int_t max_length=200;
    char *cmd=(char*)calloc(max_length, sizeof(char));
    print("calling gmsh...");
    sprintf(cmd, "gmsh mesh/shape.geo -3 -clmax %0.4f -format vtk -save_all -o mesh/shape.vtk > mesh/shape_log.txt", tol);
    // sprintf(cmd, "gmsh mesh/shape.geo -3 -format vtk -save_all -o mesh/shape.vtk > mesh/shape_log");
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

void create_vertical_wire_dipole(const real_t length, const real_t port_length, const real_t clmax){
    assert_error(length>0, "invalid length");
    assert_error(port_length<=clmax, "invalid port length");
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("Point(1) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, -length/2.0);
    file.write("Point(2) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, -port_length/2.0);
    file.write("Point(3) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, 0.0);
    file.write("Point(4) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, +port_length/2.0);
    file.write("Point(5) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, +length/2.0);
    file.write("MeshSize {1, 2, 3, 4, 5} = %21.14E;\n", clmax);
    file.write("Line(1) = {1, 2};\n");
    file.write("Line(2) = {2, 3};\n");
    file.write("Line(3) = {3, 4};\n");
    file.write("Line(4) = {4, 5};\n");
    file.write("Physical Curve(\"Port\", 1) = {2, 3};\n");
    file.write("Physical Curve(\"Wire\", 2) = {1, 4};\n");
    file.close();
}

void create_two_vertical_wire_dipole(const real_t length, const real_t port_length, const real_t d, 
    const real_t clmax){
    assert_error(length>0, "invalid length");
    assert_error(port_length<=length, "invalid port length");
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("Point(1) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -d/2.0, 0.0, -length/2.0);
    file.write("Point(2) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -d/2.0, 0.0, -port_length/2.0);
    file.write("Point(3) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -d/2.0, 0.0, 0.0);
    file.write("Point(4) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -d/2.0, 0.0, +port_length/2.0);
    file.write("Point(5) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -d/2.0, 0.0, +length/2.0);
    file.write("MeshSize {1, 2, 3, 4, 5} = %21.14E;\n", clmax);
    file.write("Line(1) = {1, 2};\n");
    file.write("Line(2) = {2, 3};\n");
    file.write("Line(3) = {3, 4};\n");
    file.write("Line(4) = {4, 5};\n");
    file.write("Point(6) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +d/2.0, 0.0, -length/2.0);
    file.write("Point(7) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +d/2.0, 0.0, -port_length/2.0);
    file.write("Point(8) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +d/2.0, 0.0, 0.0);
    file.write("Point(9) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +d/2.0, 0.0, +port_length/2.0);
    file.write("Point(10) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +d/2.0, 0.0, +length/2.0);
    file.write("MeshSize {6, 7, 8, 9, 10} = %21.14E;\n", clmax);
    file.write("Line(5) = {6, 7};\n");
    file.write("Line(6) = {7, 8};\n");
    file.write("Line(7) = {8, 9};\n");
    file.write("Line(8) = {9, 10};\n");
    file.write("Physical Curve(\"Port1\", 1) = {2, 3};\n");
    file.write("Physical Curve(\"Port2\", 2) = {6, 7};\n");
    file.write("Physical Curve(\"Wire1\", 3) = {1, 4};\n");
    file.write("Physical Curve(\"Wire2\", 4) = {5, 8};\n");
    file.close();
}

void create_loop(const real_t radius, const real_t clmax){
    assert_error(radius>0, "invalid radius");
    const real_t port_length= radius/10.0 < clmax ? radius/10.0 : clmax;
    const real_t theta=asin(port_length/(2.0*radius));
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("SetFactory(\"OpenCASCADE\");\n");
    file.write("Circle(1) = {%21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E};\n", 0.0, 0.0, 0.0, radius,
        0.0+theta, pi/2.0);
    file.write("Circle(2) = {%21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E};\n", 0.0, 0.0, 0.0, radius,
        pi/2.0, pi);
    file.write("Circle(3) = {%21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E};\n", 0.0, 0.0, 0.0, radius,
        pi, (3.0/2.0)*pi);
    file.write("Circle(4) = {%21.14E, %21.14E, %21.14E, %21.14E, %21.14E, %21.14E};\n", 0.0, 0.0, 0.0, radius,
        (3.0/2.0)*pi, 2.0*pi-theta);
    file.write("Point(9) = {%21.14E, %21.14E, %21.14E, 1.0};\n", radius*cos(theta), -port_length/2.0, 0.0);
    file.write("Point(10) = {%21.14E, %21.14E, %21.14E, 1.0};\n", radius*cos(theta), 0.0, 0.0);
    file.write("Point(11) = {%21.14E, %21.14E, %21.14E, 1.0};\n", radius*cos(theta), +port_length/2.0, 0.0);
    file.write("MeshSize {1, 3, 5, 7} = %21.14E;\n", clmax);
    file.write("Line(5) = {9, 10};\n");
    file.write("Line(6) = {10, 11};\n");
    file.write("Physical Curve(\"Port\", 1) = {5, 6};\n");
    file.write("Physical Curve(\"Wire\", 2) = {1, 2, 3, 4};\n");
    file.close();
}

void create_vertical_wire(const real_t length, const real_t clmax){
    assert_error(length>0, "invalid length");
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("Point(1) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, -length/2.0);
    file.write("Point(2) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, 0.0);
    file.write("Point(3) = {%21.14E, %21.14E, %21.14E, 1.0};\n", 0.0, 0.0, +length/2.0);
    file.write("MeshSize {1, 2, 3} = %21.14E;\n", clmax);
    file.write("Line(1) = {1, 2};\n");
    file.write("Line(2) = {2, 3};\n");
    file.close();
}

void create_transmission_line(const real_t L, const real_t S, const real_t clmax){
    assert_error(L>0, "invalid length");
    assert_error(S>0, "invalid length");
    const real_t port_length= S/2 < clmax ? S/2 : clmax;
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("Point(1) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -L/2.0, -S/2.0, 0.0);
    file.write("Point(2) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -L/2.0, -port_length/2.0, 0.0);
    file.write("Point(3) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -L/2.0, 0.0, 0.0);
    file.write("Point(4) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -L/2.0, +port_length/2.0, 0.0);
    file.write("Point(5) = {%21.14E, %21.14E, %21.14E, 1.0};\n", -L/2.0, +S/2, 0.0);
    file.write("Point(6) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +L/2.0, -S/2, 0.0);
    file.write("Point(7) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +L/2.0, -port_length/2.0, 0.0);
    file.write("Point(8) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +L/2.0, 0.0, 0.0);
    file.write("Point(9) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +L/2.0, +port_length/2.0, 0.0);
    file.write("Point(10) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +L/2.0, +S/2, 0.0);
    file.write("Point(11) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +0.0, -S/2, 0.0);
    file.write("Point(12) = {%21.14E, %21.14E, %21.14E, 1.0};\n", +0.0, +S/2, 0.0);
    //
    // file.write("MeshSize {11, 12} = %21.14E;\n", clmax<S ? clmax : S);
    file.write("MeshSize {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = %21.14E;\n", clmax);

    file.write("Line(1) = {1, 2};\n");
    file.write("Line(2) = {2, 3};\n");
    file.write("Line(3) = {3, 4};\n");
    file.write("Line(4) = {4, 5};\n");

    file.write("Line(5) = {7, 6};\n");
    file.write("Line(6) = {7, 8};\n");
    file.write("Line(7) = {8, 9};\n");
    file.write("Line(8) = {9 , 10};\n");

    file.write("Line(9) = {1, 11};\n");
    file.write("Line(10) = {11, 6};\n");
    file.write("Line(11) = {5, 12};\n");
    file.write("Line(12) = {12, 10};\n");

    file.write("Physical Curve(\"Port1\", 1) = {2, 3};\n");
    file.write("Physical Curve(\"Port2\", 2) = {6, 7};\n");
    file.write("Physical Curve(\"TL\", 3) = {1, 4, 5, 8, 9, 10, 11, 12};\n");
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

void create_patch_antenna(){
    file_t file;
    file.open("mesh/shape.geo", 'w');
    file.write("Merge \"cad/patch_antenna.brep\";\n");
    file.write("MeshSize {27, 28, 29, 30, 31, 32, 21, 22, 23, 24, 25, 26} = 4.0; // Patch;\n");
    file.write("MeshSize {6, 8, 5, 7, 2, 4, 1, 3} = 6.0; // Substrate Corners\n");
    file.write("MeshSize {15, 16, 17, 18, 19, 20, 9, 10, 11, 12, 13, 14} = 6.0; // GND Patch\n");
    file.write("MeshSize {20, 32, 21, 9} = 4.0; // Port\n");
    file.write("Physical Surface(\"Patch\", 1) = {9};\n");
    file.write("Physical Surface(\"GND\", 2) = {8, 4};\n");
    file.write("Physical Surface(\"Port\", 3) = {10};\n");
    file.write("Physical Volume(\"Substrate\", 4) = {1};\n");
    file.close();
}