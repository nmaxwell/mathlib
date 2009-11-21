#ifndef GRID1D_FFT_H
#define GRID1D_FFT_H

#include "grid.h"
#include "grid1D.h"

#include <vector>
#include <complex>
#include <fftw3.h>

template<class St1,class Xt1, class St2,class Xt2 >
void FFT(
	grid1D<double,St1,Xt1 > & in,
	grid1D<complex<double>,St2,Xt2 > & out );

template<class St1,class Xt1, class St2,class Xt2 >
void IFFT(
	grid1D<complex<double>,St1,Xt1 > & in,
	grid1D<double,St2,Xt2 > & out );

#endif
