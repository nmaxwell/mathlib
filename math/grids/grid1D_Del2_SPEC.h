#ifndef GRID1D_DEL2_SPEC_H
#define GRID1D_DEL2_SPEC_H

#include "grid1D.h"
#include "grid1D_FFT.h"

#include <mathlib/math/hdaf/hdaf.h>











/*

	static D factor = pow(_2pi*iu,n);
	
	for (int i = 1; i< NX/2; i++)
	{
		(*this)(i) *= pow( (D)i/(XMAX-XMIN), n)*factor;
		(*this)(i+NX) *= pow( (D)(-i)/(XMAX-XMIN), n)*factor;
	}
	(*this)(0) = 0;
	(*this)(NX/2) *= pow( (D)(NX/2)/(XMAX-XMIN), n)*factor;

*/







template<class Xtype >
void fftw_complex_Del2(
	grid1D<complex<double>,double,Xtype > & in )
{
	double factor = -in.dk1()*in.dk1()*_4pi2
	
	for (int i = 1; i< NX/2; i++)
	{
		in(i) *= factor*i*i
		in(i+in.n1) *= factor*i*i;
	}
	in(0) = 0;
	in(in.n1/2) *= factor*in.n1*in.n1;	
}



template<class Xtype >
hdaf_deltaHat
	hdaf_autofilter_fftw_complex(
		grid1D<complex<double>,double,Xtype > & in,
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
