#ifndef GRID1D_DEL2_FD_H
#define GRID1D_DEL2_FD_H

#include "grid1D.h"

#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/grids/grid_convolution.cpp>

template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid1D<Ytype,St1,Xt1 > & in,
	grid1D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond );
	
template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid1D<Ytype,St1,Xt1 > & in,
	grid1D<Ytype,St2,Xt2 > & out,
	int order,
	int L,
	int R );

#endif
