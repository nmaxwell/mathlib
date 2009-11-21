#ifndef GRID3D_NORMS_CPP
#define GRID3D_NORMS_CPP

#ifndef GRID3D_NORMS_H
	#include "grid3D_norms.h"
#endif

template<classgridPars_ >
YtypeScalar maxNorm(grid3D<Ytype,YtypeScalar,Xtype  > & rhs, int O)
{
	YtypeScalar s = norm(rhs(O,O,O));
	for (int i = O; i<rhs.n1-O; i++)
	for (int j = O; j<rhs.n2-O; j++)
	for (int k = O; k<rhs.n3-O; k++)
		if ( norm(rhs(i,j,k)) > s ) s = norm(rhs(i,j,k)); 
	return s;
}

template<classgridPars_ >
YtypeScalar L2norm(grid3D<Ytype,YtypeScalar,Xtype  > & rhs, int O)
{
	YtypeScalar s = 0.0;
	for (int i = O; i<rhs.n1-O; i++)
	for (int j = O; j<rhs.n2-O; j++)
	for (int k = O; k<rhs.n3-O; k++)
		s += norm2(rhs(i,j,k));
	return sqrt(s*rhs.dx1()*rhs.dx2()*rhs.dx3());
}


template<classgridPars_ >
YtypeScalar norm(grid3D<Ytype,YtypeScalar,Xtype  > & rhs)
{
	return L2norm(rhs);
}

template<classgridPars_ >
YtypeScalar norm2(grid3D<Ytype,YtypeScalar,Xtype  > & rhs)
{
	YtypeScalar nrm = rhs.L2norm();
	return nrm*nrm;
}

#endif
