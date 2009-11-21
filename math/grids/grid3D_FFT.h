#ifndef GRID3D_FFT_H
#define GRID3D_FFT_H

#include "grid.h"
#include "grid3D.h"

#include <vector>
#include <complex>
#include <fftw3.h>

template<class St1,class Xt1, class St2,class Xt2 >
void FFT(
	grid3D<double,St1,Xt1 > & in,
	grid3D<complex<double>,St2,Xt2 > & out );

template<class St1,class Xt1, class St2,class Xt2 >
void IFFT(
	grid3D<complex<double>,St1,Xt1 > & in,
	grid3D<double,St2,Xt2 > & out );

#endif
