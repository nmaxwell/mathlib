#ifndef GRID1D_NORMS_H
#define GRID1D_NORMS_H

#include "mathlib/math/norms.h"

#include "grid1D.h"

template<classgridPars_ >
YtypeScalar maxNorm(grid1D<Ytype,YtypeScalar,Xtype  > & rhs, int O=0);

template<classgridPars_ >
YtypeScalar L2norm(grid1D<Ytype,YtypeScalar,Xtype  > & rhs, int O=0);

template<classgridPars_ >
YtypeScalar norm(grid1D<Ytype,YtypeScalar,Xtype  > & rhs);

template<classgridPars_ >
YtypeScalar norm2(grid1D<Ytype,YtypeScalar,Xtype  > & rhs);






#endif
