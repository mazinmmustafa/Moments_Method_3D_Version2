SetFactory("OpenCASCADE");
Sphere(1) = { 0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,  5.00000000000000E-01, -Pi/2, Pi/2, 2*Pi};
Physical Surface("Surface", 1) = {1};
Physical Volume("Volume", 1) = {-1};
