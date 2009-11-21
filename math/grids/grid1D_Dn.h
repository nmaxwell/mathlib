#ifndef GRID1D_Dn_H
#define GRID1D_Dn_H

#include "grid1D.h"

#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/grids/grid_convolution.cpp>

/*
#define GRID1D_D1_FD
#define GRID1D_D2_FD
#define GRID1D_D3_FD
#define GRID1D_D4_FD
*/

#ifdef GRID1D_Del1_FD

#define INCLUDE_FD_D1_LR

template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del1_FD(
	grid1D<Ytype,St1,Xt1 > & in,
	grid1D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond );
	
#endif









#endif
