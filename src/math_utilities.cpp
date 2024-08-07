//
#include "math_utilities.hpp"

const real_t eps_sinc=1.0E-4;

complex_t sinc(const complex_t z){
    return abs(z)<eps_sinc ? 1.0-z*z/6.0+z*z*z*z/120.0 : sin(z)/z;
}

real_t sinc(const real_t x){
    return abs(x)<eps_sinc ? 1.0-x*x/6.0+x*x*x*x/120.0 : sin(x)/x;
}

complex_t sinc_dx(const complex_t z){
    return abs(z)<eps_sinc ? -z/3.0+z*z*z/30.0-z*z*z*z*z/840.0 : (cos(z)-sinc(z))/z;
}

real_t sinc_dx(const real_t x){
    return abs(x)<eps_sinc ? -x/3.0+x*x*x/30.0-x*x*x*x*x/840.0 : (cos(x)-sinc(x))/x;
}

real_t deg2rad(const real_t theta){
    return theta*pi/180.0;
}

real_t rad2deg(const real_t theta){
    return theta*180.0/pi;
}

real_t round_m(const real_t x, const real_t n){
    real_t y=x*pow(10.0, n);
    return round(y)/pow(10.0, n);
}