#ifndef GRID2D_NORMS_H
#define GRID2D_NORMS_H

#include "mathlib/math/norms.h"

#include "grid2D.h"

template<classgridPars_ >
YtypeScalar maxNorm(grid2D<Ytype,YtypeScalar,Xtype  > & rhs, int O=0);

template<class Ytype, class YtypeScalar, class Xtype >
Ytype max_diff(grid2D<Ytype, YtypeScalar, Xtype > & A, grid2D<Ytype, YtypeScalar, Xtype > & B, int O=0 );

template<classgridPars_ >
YtypeScalar L2norm(grid2D<Ytype,YtypeScalar,Xtype  > & rhs, int O=0);

template<classgridPars_ >
YtypeScalar norm(grid2D<Ytype,YtypeScalar,Xtype  > & rhs);

template<classgridPars_ >
YtypeScalar norm2(grid2D<Ytype,YtypeScalar,Xtype  > & rhs);

#endif
