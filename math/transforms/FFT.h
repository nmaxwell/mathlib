#ifndef FFT_H
#define FFT_H

#include <fftw3.h>

#include <mathlib/math/std_math.h>

#ifndef N_FFT_THREADS
	#define N_FFT_THREADS  N_NATIVE_THREADS
#endif

#ifndef FFTW_ESTIMATE_MODE
	#define FFTW_ESTIMATE_MODE  FFTW_ESTIMATE
	// #define FFTW_ESTIMATE_MODE FFTW_PATIENT
#endif

//  maps index to omega or k value, for fftw type mappings
double fftw_map(int i, int N, double period) 
{
	if (i <= N/2 )
		return ((double)i/(period));
	else
		return (double)(i-N)/(period);				
}

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





// house keeping stuff
static bool fftw_initialized = false;

void initialize_fftw()
{
	if (!fftw_initialized)
		assert(fftw_init_threads());
	
	fftw_initialized = true;
}

void release_fftw()
{
	if (fftw_initialized)	
		fftw_cleanup_threads();
		
	fftw_initialized = false;
}



#endif




