#ifndef LAPLACIAN_CPP
#define LAPLACIAN_CPP


#include "laplacian.h"


void laplacian_2d_kernel_extension ( double *& ker_ext, int n1, int n2, double * kernel1, int n_l1, int n_r1, double * kernel2, int n_l2, int n_r2 )
{
    //  i*n2+j
    
    if (ker_ext == 0) ker_ext = (double*)fftw_malloc(sizeof(double)*(n1*n2) );
    
    for (int k=0; k<n1*n2; k++)
        ker_ext[k] = 0;
    
    for (int k1=-n_l1; k1<0; k1++)
    {
        int k = k1;
        while ( k < 0 ) k += n1;
        ker_ext[ k*n2 ] += kernel1[k1];
    }
    
    for (int k1=0; k1<=n_r1; k1++)
        ker_ext[ (k1%n1)*n2 ] += kernel1[k1];
    
    for (int k2=-n_l2; k2<0; k2++)
    {
        int k = k2;
        while ( k < 0 ) k += n2;
        ker_ext[k] += kernel2[k2];
    }
    
    for (int k2=0; k2<=n_r2; k2++)
        ker_ext[ k2%n2 ] += kernel2[k2];
}








void laplacian_2d_fft ( double * in, double *& out, int n1, int n2, double L1, double L2 )
{
    // assuming n1,n2 are even 
    
    double s1 = -ml_4pi2/((L1*L1)*(n1*n2));
    double s2 = -ml_4pi2/((L2*L2)*(n1*n2));
    
    int n1_ = n1/2;
    int n2_ = n2/2+1;
    
    complex<double> *Q=0;
    
    fft_2d( in, Q, n1, n2 );
    
    for ( int k2=0; k2<n2_; k2++ )
    {
        for ( int k1=0; k1<=n1_; k1++ )
            Q[k1*n2_+k2] *= s2*k2*k2 + s1*k1*k1;
        
        for ( int k1=n1_+1; k1<n1; k1++ )
            Q[k1*n2_+k2] *= s2*k2*k2 + s1*(k1-n1)*(k1-n1);
    }
    
    ifft_2d( Q, out, n1, n2 );
    fftw_free( Q );
}

void laplacian_3d_fft ( double * in, double *& out, int n1, int n2, int n3, double L1, double L2, double L3 )
{
    // assuming n1,n2 are even.
    // not yet verified...
    
    double s1 = -ml_4pi2/((L1*L1)*(n1*n2*n3));
    double s2 = -ml_4pi2/((L2*L2)*(n1*n2*n3));
    double s3 = -ml_4pi2/((L3*L3)*(n1*n2*n3));
    
    int n1_ = n1/2;
    int n2_ = n2/2;
    int n3_ = n3/2+1;
    
    complex<double> *Q=0;
    
    fft_3d( in, Q, n1, n2, n3 );
    
    for ( int k3=0; k3<n3_; k3++ )
    {
        for ( int k2=0; k2<=n2_; k2++ )
        {
            for ( int k1=0; k1<=n1_; k1++ )
                Q[k1*n2_+k2] *= s3*k3*k3 + s2*k2*k2 + s1*k1*k1;
            
            for ( int k1=n1_+1; k1<n1; k1++ )
                Q[k1*n2_+k2] *= s3*k3*k3 + s2*k2*k2 + s1*(k1-n1)*(k1-n1);
        }
        
        for ( int k2=n2_+1; k2<=n2; k2++ )
        {
            for ( int k1=0; k1<=n1_; k1++ )
                Q[k1*n2_+k2] *= s3*k3*k3 + s2*k2*k2 + s1*k1*k1;
            
            for ( int k1=n1_+1; k1<n1; k1++ )
                Q[k1*n2_+k2] *= s3*k3*k3 + s2*(k2-n2)*(k2-n2) + s1*(k1-n1)*(k1-n1);
        }
    }
    
    ifft_3d( Q, out, n1, n2, n3 );
    fftw_free( Q );
}












#endif
