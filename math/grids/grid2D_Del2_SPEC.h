#ifndef GRID2D_DEL2_SPEC_H
#define GRID2D_DEL2_SPEC_H

#include "grid2D.h"

#include <mathlib/math/hdaf/hdaf.h>

template<class Xtype, class YtypeScalar >
functor<Xtype,YtypeScalar> 
	hdaf_autofilter_fftw_complex(
		grid2D<complex<double>,double,Xtype > & in,
		int order,
		double eps,
		double fudge );

/*
template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid2D<Ytype,St1,Xt1 > & in,
	grid2D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond );
	*/

#endif
