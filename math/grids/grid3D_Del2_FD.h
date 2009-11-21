#ifndef GRID3D_DEL2_FD_H
#define GRID3D_DEL2_FD_H

#include "grid3D_FFT.h"
#include "grid3D.h"

#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/grids/grid_convolution.cpp>

template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid3D<Ytype,St1,Xt1 > & in,
	grid3D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond,
	int BD_order=24 );

#endif
