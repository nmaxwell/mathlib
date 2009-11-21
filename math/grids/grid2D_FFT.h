#ifndef GRID2D_FFT_H
#define GRID2D_FFT_H

#include "grid.h"
#include "grid2D.h"

#include <vector>
#include <complex>
#include <fftw3.h>

template<class St1,class Xt1, class St2,class Xt2 >
void FFT(
	grid2D<double,St1,Xt1 > & in,
	grid2D<complex<double>,St2,Xt2 > & out );

template<class St1,class Xt1, class St2,class Xt2 >
void IFFT(
	grid2D<complex<double>,St1,Xt1 > & in,
	grid2D<double,St2,Xt2 > & out );

#endif
