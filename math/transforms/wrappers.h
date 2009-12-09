#ifndef FFT_WRAPPERS_H
#define FFT_WRAPPERS_H

#include "fft.h"




void fft(
	complex<double > * in,
    complex<double > * & out,
    int n,
    int in_stride=1,
    int out_stride=1 )
{
    static fftw_mz_1d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n );
    
    int plan = mz.plan( n, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_mz_1d.execute returned null in fft.\n";
}

void ifft(
	complex<double > * in,
    complex<double > * & out,
    int n,
    int in_stride=1,
    int out_stride=1 )
{
    static ifftw_mz_1d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n );
    
    int plan = mz.plan( n, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: ifftw_mz_1d.execute returned null in fft.\n";
}

void fft_2d(
	complex<double > * in,
    complex<double > * & out,
    int n1,
    int n2,
    int in_stride=1,
    int out_stride=1 )
{
    static fftw_mz_2d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n1*n2 );
    
    int plan = mz.plan( n1, n2, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_mz_2d.execute returned null in fft_2d.\n";
}

void ifft_2d(
	complex<double > * in,
    complex<double > * & out,
    int n1,
    int n2,
    int in_stride=1,
    int out_stride=1 )
{
    static ifftw_mz_2d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n1*n2 );
    
    int plan = mz.plan( n1, n2, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: ifftw_mz_2d.execute returned null in ifft_2d.\n";
}


void fft_3d(
	complex<double > * in,
    complex<double > * & out,
    int n1,
    int n2,
    int n3 )
{
    static fftw_mz_3d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n1*n2*n3 );
    
    int plan = mz.plan( n1, n2, n3, 1, 1 );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_mz_3d.execute returned null in fft.\n";
}

void ifft_3d(
	complex<double > * in,
    complex<double > * & out,
    int n1,
    int n2,
    int n3 )
{
    static ifftw_mz_3d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n1*n2*n3 );
    
    int plan = mz.plan( n1, n2, n3, 1, 1 );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: ifftw_mz_3d.execute returned null in fft.\n";
}







void fft(
	double * in,
    complex<double > * & out,
    int n,
    int in_stride=1,
    int out_stride=1 )
{
    static fftw_r2c_mz_1d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n/2+1 );
    
    int plan = mz.plan( n, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_r2c_mz_1d.execute returned null in fft.\n";
}

void ifft(
	complex<double > * in,
    double * & out,
    int n,
    int in_stride=1,
    int out_stride=1 )
{
    static fftw_c2r_mz_1d mz;
    
    if (out == 0) out = (double *)fftw_malloc(sizeof(double) * n );
    
    int plan = mz.plan( n, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_r2c_mz_1d.execute returned null in fft.\n";
}

void fft_2d(
	double * in,
    complex<double > * & out,
    int n1,
    int n2,
    int in_stride=1,
    int out_stride=1 )
{
    static fftw_r2c_mz_2d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n1*(n2/2+1) );
    
    int plan = mz.plan( n1, n2, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_r2c_mz_2d.execute returned null in fft_2d.\n";
}

void ifft_2d(
	complex<double > * in,
    double * & out,
    int n1,
    int n2,
    int in_stride=1,
    int out_stride=1 )
{
    static fftw_c2r_mz_2d mz;
    
    if (out == 0) out = (double *)fftw_malloc(sizeof(double) * n1*n2 );
    
    int plan = mz.plan( n1, n2, in_stride, out_stride );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_c2r_mz_2d.execute returned null in ifft_2d.\n";
}


void fft_3d(
	double * in,
    complex<double > * & out,
    int n1,
    int n2,
    int n3 )
{
    static fftw_r2c_mz_3d mz;
    
    if (out == 0) out = (complex<double> *)fftw_malloc(sizeof(fftw_complex) * n1*n2*(n3/2+1) );
    
    int plan = mz.plan( n1, n2, n3, 1, 1 );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_r2c_mz_3d.execute returned null in fft.\n";
}

void ifft_3d(
	complex<double > * in,
    double * & out,
    int n1,
    int n2,
    int n3 )
{
    static fftw_c2r_mz_3d mz;
    
    if (out == 0) out = (double *)fftw_malloc(sizeof(double) * n1*n2*n3 );
    
    int plan = mz.plan( n1, n2, n3, 1, 1 );
    
    int res = mz.execute( plan, in, out );
    
    if ( res )
        std::cerr << "error: fftw_c2r_mz_3d.execute returned null in fft.\n";
}







/*



// FFT: 

// computes real to complex fft, in must be of lengh n, out must be of
//    length "floor(n/2) + 1", if out is left unallocated, pass it as
//    out = 0, and the function will allocate it. 
//    scale the output by 1/n for a unitary transformation.
//    input is preserved 

void FFT(
	double * in,
    complex<double > * & out,
    int n,
    int in_stride=1,
    int out_stride=1 );

void FFT(
	float * in,
    complex<float > * & out,
    int n,
    int in_stride=1,
    int out_stride=1 );



void fft(
	double * in,
    complex<double > * & out,
    int n,
    int in_stride=1,
    int out_stride=1 );

// IFFT: 

// computes complex to real inverse fft, in must be of lengh 
//    "floor(n/2) + 1", out must be of
//    length n, if out is left unallocated, pass it as
//    out = 0, and the function will allocate it. 
//    input is NOT guaranteed to be preserved

void IFFT(
    complex<double > * in,
	double * & out,
    int n,
    int in_stride=1,
    int out_stride=1 );

void IFFT(
    complex<float > * in,
	float * & out,
    int n,
    int in_stride=1,
    int out_stride=1 );




*/




#endif
