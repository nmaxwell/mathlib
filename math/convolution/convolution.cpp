#ifndef ML_CONVOLUTION_CPP
#define ML_CONVOLUTION_CPP

#include "convolution.h"

template<class T >
void periodic_convolution(T * ker, T * in, T *& out, int n_l, int n_r, int n, int in_stride, int out_stride )
{
    if (1+n_l+n_r <=n)
    {
        complex<T > * c_in=0;
        FFT( in, c_in, n, in_stride, 1 );
        T * extended_ker = new T [n];
        for (int k=0; k<n; k++ )
            extended_ker[k] = 0;        
        for (int k=0; k<=n_l; k++ )
            extended_ker[k] = ker[k];
        for (int k=1; k<=n_r; k++ )
            extended_ker[n-k] = ker[-k];
        complex<T > * c_ker=0;
        
        FFT( extended_ker, c_ker, n, 1, 1 );        
        for (int k=0; k<n/2+1; k++ )
            c_in[k] *= c_ker[k]/((T)n);        
        IFFT( c_in, out, n, 1, out_stride );
                
        delete [] extended_ker;
    }
}

template<class T >
void non_periodic_convolution(T * ker, T * in, T *& out, int n_l, int n_r, int n, int in_stride, int out_stride )
{
    periodic_convolution( ker, in, out, n_l, n_r, n, in_stride, out_stride );
    
    for (int k=0; k<n_l; k++)
    {
        out[k*out_stride] = 0;
        for ( int j = -n_r; j <= n_l; j++ ) 
                if ( k-j >= 0 && k-j < n ) out[k*out_stride] += ker[j] * in[(k-j)*in_stride];
    }
    
    for (int k=n-n_r; k<n; k++)
    {
        out[k*out_stride] = 0;
        for ( int j = -n_r; j <= n_l; j++ ) 
                if ( k-j >= 0 && k-j < n ) out[k*out_stride] += ker[j] * in[(k-j)*in_stride];
    }
}


template<class T >
void non_periodic_convolution_slow(T * ker, T * in, T *& out, int n_l, int n_r, int n, int in_stride, int out_stride )
{
    if (!out) out = new T [n];
    
    for (int k=0; k<n; k++)
    {
        out[k*out_stride] = 0;
        for ( int j = -n_r; j <= n_l; j++ ) 
                if ( k-j >= 0 && k-j < n ) out[k*out_stride] += ker[j] * in[(k-j)*in_stride];
    }
}

#endif
