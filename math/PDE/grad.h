#ifndef GRAD_H
#define GRAD_H

#include <mathlib/math/std_math.h>

#include <mathlib/math/transforms/fft.h>


void grad_2d_kernel_extension ( double *& ker_ext, int n1, int n2, double * kernel, int n_l, int n_r, int direction )
{
    //  i*n2+j
    
    if (ker_ext == 0) ker_ext = (double*)fftw_malloc(sizeof(double)*(n1*n2) );
    
    for (int k=0; k<n1*n2; k++)
        ker_ext[k] = 0;
    
    if (direction == 1)
    {
        for (int k1=-n_l; k1<0; k1++)
        {
            int k = k1;
            while ( k < 0 ) k += n1;
            ker_ext[ k*n2 ] += kernel[k1];
        }
        
        for (int k1=0; k1<=n_r; k1++)
            ker_ext[ (k1%n1)*n2 ] += kernel[k1];
    }
    
    if (direction == 2)
    {
        for (int k2=-n_l; k2<0; k2++)
        {
            int k = k2;
            while ( k < 0 ) k += n2;
            ker_ext[k] += kernel[k2];
        }
        
        for (int k2=0; k2<=n_r; k2++)
            ker_ext[ k2%n2 ] += kernel[k2];
    }
}



void grad_3d_kernel_extension ( double *& ker_ext, int n1, int n2, int n3, double * kernel, int n_l, int n_r, int direction )
{
    //  n3*(i*n2+j)+k
    
    if (ker_ext == 0) ker_ext = (double*)fftw_malloc(sizeof(double)*(n1*n2*n3) );
    
    for (int k=0; k<n1*n2*n3; k++)
        ker_ext[k] = 0;
    
    if (direction == 1)
    {
        for (int k1=-n_l; k1<0; k1++)
        {
            int k = k1;
            while ( k < 0 ) k += n1;
            ker_ext[ k*n2*n3 ] += kernel[k1];
        }
        
        for (int k1=0; k1<=n_r; k1++)
            ker_ext[ (k1%n1)*n2*n3 ] += kernel[k1];
    }
    
    if (direction == 2)
    {
        for (int k2=-n_l; k2<0; k2++)
        {
            int k = k2;
            while ( k < 0 ) k += n2;
            ker_ext[k*n3] += kernel[k2];
        }
        
        for (int k2=0; k2<=n_r; k2++)
            ker_ext[ (k2%n2)*n3 ] += kernel[k2];
    }
    
    if (direction == 3)
    {
        for (int k3=-n_l; k3<0; k3++)
        {
            int k = k3;
            while ( k < 0 ) k += n3;
            ker_ext[ k ] += kernel[k3];
        }
        
        for (int k3=0; k3<=n_r; k3++)
            ker_ext[ k3%n3 ] += kernel[k3];
    }
}








#endif
