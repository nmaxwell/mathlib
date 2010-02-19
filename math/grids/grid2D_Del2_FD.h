#ifndef GRID2D_DEL2_FD_H
#define GRID2D_DEL2_FD_H

#include "grid2D.h"
#include "grid2D_FFT.h"

#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/grids/grid_convolution.cpp>

template<class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid2D<float,St1,Xt1 > & in,
	grid2D<float,St2,Xt2 > & out,
	int order,
	int bdcond,
	int BD_order=24 );

template<class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid2D<double,St1,Xt1 > & in,
	grid2D<double,St2,Xt2 > & out,
	int order,
	int bdcond,
	int BD_order=24 );
	
#endif

