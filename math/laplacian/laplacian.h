#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include <mathlib/math/std_math.h>

#include <mathlib/math/transforms/fft.h>

//i*n2+j

void laplacian( double * in, double *& out, int n1, int n2, double L1, double L2 )
{
    // assuming n1,n2 are even.
    
    double s1 = ml_2pi/(L1); s1 *= -s1/(n1*n2);
    double s2 = ml_2pi/(L1); s2 *= -s2/(n1*n2);
    
    int n1_ = n1/2;
    int n2_ = n2/2+1;
    
    complex<double> *Q=0;
    fft_2d( in, Q, n1, n2 );
    
    for ( int k2=0; k2<n2_; k2++ )
    {
        for ( int k1=0; k1<n1_; k1++ )
            Q[k1*n2_+k2] *= s2*(k2 << 1) + s1*(k1 << 1);
        
        for ( int k1=n1_; k1<n1; k1++ )
            Q[k1*n2_+k2] *= s2*(k2 << 1) + s1*((k1-n1) << 1);
    }
    
    ifft_2d( Q, out, n1, n2 );
}

















#endif
