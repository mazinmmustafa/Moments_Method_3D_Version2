//
#include "shape.hpp"




//
void call_gmsh(const real_t tol){
    int_t max_length=200;
    char *cmd=(char*)calloc(max_length, sizeof(char));
    sprintf(cmd, "gmsh mesh/shape.geo -3 -clmax %0.4f -format msh -save_all -o mesh/shape.msh", tol);
    system(cmd);
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