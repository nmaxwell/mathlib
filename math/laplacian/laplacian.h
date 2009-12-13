#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include <mathlib/math/std_math.h>

#include <mathlib/math/transforms/fft.h>


void laplacian_2d_kernel_extension ( double *& ker_ext, int n1, int n2, double * kernel1, int n_l1, int n_r1, double * kernel2, int n_l2, int n_r2 );


void laplacian_2d_fft ( double * in, double *& out, int n1, int n2, double L1, double L2 );
    // assuming n1,n2 are even 

void laplacian_3d_fft ( double * in, double *& out, int n1, int n2, int n3, double L1, double L2, double L3 );
    // assuming n1,n2 are even.
    // not yet verified...















#endif
