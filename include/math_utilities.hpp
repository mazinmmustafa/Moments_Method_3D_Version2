#ifndef __MATH_UTILITIES_HPP__
#define __MATH_UTILITIES_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"

// Definitions

real_t deg2rad(const real_t theta);
real_t rad2deg(const real_t theta);

real_t sinc(const real_t x);
complex_t sinc(const complex_t z);
real_t sinc_dx(const real_t x);
complex_t sinc_dx(const complex_t z);
real_t round_m(const real_t x, const real_t n);

#endif