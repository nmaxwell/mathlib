#ifndef GRID2D_NORMS_CPP
#define GRID2D_NORMS_CPP

#ifndef GRID2D_NORMS_H
	#include "grid2D_norms.h"
#endif

template<classgridPars_ >
YtypeScalar maxNorm(grid2D<Ytype,YtypeScalar,Xtype  > & rhs, int O)
{
	YtypeScalar s = norm(rhs(O,O));
	for (int i = O; i<rhs.n1-O; i++)
	for (int j = O; j<rhs.n2-O; j++)
		if ( norm(rhs(i,j)) > s ) s = norm(rhs(i,j)); 
	return s;
}

template<class Ytype, class YtypeScalar, class Xtype >
YtypeScalar max_diff(grid2D<Ytype, YtypeScalar, Xtype > & A, grid2D<Ytype, YtypeScalar, Xtype > & B, int O )
{
    YtypeScalar s = 0;
	for (int i = O; i<A.n1-O; i++)
	for (int j = O; j<A.n2-O; j++)
		if ( norm(A(i,j) - B(i,j) ) > s ) s = norm(A(i,j) - B(i,j)); 
	return s;
}

template<classgridPars_ >
YtypeScalar L2norm(grid2D<Ytype,YtypeScalar,Xtype  > & rhs, int O)
{
	YtypeScalar s = 0.0;
	for (int i = O; i<rhs.n1-O; i++)
	for (int j = O; j<rhs.n2-O; j++)
		s += norm2(rhs(i,j));
	return sqrt(s*rhs.dx1()*rhs.dx2());
}


template<classgridPars_ >
YtypeScalar norm(grid2D<Ytype,YtypeScalar,Xtype  > & rhs)
{
	return L2norm(rhs);
}

template<classgridPars_ >
YtypeScalar norm2(grid2D<Ytype,YtypeScalar,Xtype  > & rhs)
{
	YtypeScalar nrm = L2norm(rhs);
	return nrm*nrm;
}

#endif
