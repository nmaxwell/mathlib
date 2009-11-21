#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H

#include <complex>
#include <vector>

#include <fftw3.h>

#include <mathlib/math/std_math.h>
#include <mathlib/math/kernels/FDkernels.h>


template<class X, class Y >
Y cheap_diff(X x, functor<X,Y > const & f, int p, X h, int m=24);
	// returns pth derivative of f around x, 
	// using finite differences with a step size of h, and using m points
	
double cheap_diff(double x, double (&f)(double), int p, double h, int m=24);
	// returns pth derivative of f around x, 
	// using finite differences with a step size of h, and using m points
		



template<class T >
void differentiate(
	T * in,
    T * & out,
    int n,
    T h,
    int diff_order,
    int in_stride=1,
    int out_stride=1 );

#endif
