#ifndef GRID3D_NORMS_H
#define GRID3D_NORMS_H

#include "mathlib/math/norms.h"

#include "grid3D.h"

template<classgridPars_ >
YtypeScalar maxNorm(grid3D<Ytype,YtypeScalar,Xtype  > & rhs, int O=0);

template<classgridPars_ >
YtypeScalar L2norm(grid3D<Ytype,YtypeScalar,Xtype  > & rhs, int O=0);

template<classgridPars_ >
YtypeScalar norm(grid3D<Ytype,YtypeScalar,Xtype  > & rhs);

template<classgridPars_ >
YtypeScalar norm2(grid3D<Ytype,YtypeScalar,Xtype  > & rhs);

#endif
