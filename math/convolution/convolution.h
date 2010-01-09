#ifndef ML_CONVOLUTION_H
#define ML_CONVOLUTION_H

#include <mathlib/math/std_math.h>
#include <mathlib/math/transforms/fft.h>



/*
    periodic_convolution:
        just applies fft convolution theorem.
    
    non_periodic_convolution_slow:
    
    idea:
        out_k += sum_{j=-n_r}^{n_l} ker_j * in_{k-j}
        for 0 <= k < n
    
    does:
        for (int k=0; k<n; k++)
        {
            out[k] = 0;
            //for ( int j = min(-n_r,k-n-1); j <= max(n_l,k); j++ ) 
            //        out[k] += ker[j] * in[k-j];
            for ( int j = -n_r; j <= n_l; j++ ) 
                    if ( k-j >= 0 && k-j < n ) out[k] += ker[j] * in[k-j];
        }
    
    ker[0] us the center point in the kernel, so ker[j] for -n_r <= j <= n_l must be accesible.
    
    non_periodic_convolution:
        exact same functionality as non_periodic_convolution_slow, but 
        manipulates fft convolution theorem, to fix end effects to be non-periodic.
        much faster than non_periodic_convolution_slow for large jobs.
*/

template<class T >
void periodic_convolution(T * ker, T * in, T *& out, int n_l, int n_r, int n, int in_stride=1, int out_stride=1 );

template<class T >
void non_periodic_convolution(T * ker, T * in, T *& out, int n_l, int n_r, int n, int in_stride=1, int out_stride=1 );

template<class T >
void non_periodic_convolution_slow(T * ker, T * in, T *& out, int n_l, int n_r, int n, int in_stride=1, int out_stride=1 );


#endif

