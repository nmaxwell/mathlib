#ifndef GRID1D_NORMS_CPP
#define GRID1D_NORMS_CPP

#ifndef GRID1D_NORMS_H
	#include "grid1D_norms.h"
#endif

template<classgridPars_ >
YtypeScalar maxNorm(grid1D<Ytype,YtypeScalar,Xtype  > & rhs, int O)
{
	YtypeScalar s = norm(rhs(O));
	for (int i = O; i<rhs.n1-O; i++)
		if ( norm(rhs(i)) > s ) s = norm(rhs(i));
	return s;
}

template<classgridPars_ >
YtypeScalar L2norm(grid1D<Ytype,YtypeScalar,Xtype  > & rhs, int O)
{
	YtypeScalar s = 0.0;
	for (int i = O; i<rhs.n1-O; i++)
		s += norm2(rhs(i));
	return sqrt(s*rhs.dx1());
}


template<classgridPars_ >
YtypeScalar norm(grid1D<Ytype,YtypeScalar,Xtype  > & rhs)
{
	return L2norm(rhs);
}

template<classgridPars_ >
YtypeScalar norm2(grid1D<Ytype,YtypeScalar,Xtype  > & rhs)
{
	YtypeScalar nrm = rhs.L2norm();
	return nrm*nrm;
}

#endif
