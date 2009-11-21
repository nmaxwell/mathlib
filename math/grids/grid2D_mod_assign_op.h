#ifndef GRID2D_MOD_ASSIGN_OP_H
#define GRID2D_MOD_ASSIGN_OP_H

#include "grid2D.h"

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator+=( grid2D<Yt1,St1,Xt1 > & lhs, grid2D<Yt2,St2,Xt2 > & rhs );

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator-=( grid2D<Yt1,St1,Xt1 > & lhs, grid2D<Yt2,St2,Xt2 > & rhs );

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator*=( grid2D<Yt1,St1,Xt1 > & lhs, grid2D<Yt2,St2,Xt2 > & rhs );

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator/=( grid2D<Yt1,St1,Xt1 > & lhs, grid2D<Yt2,St2,Xt2 > & rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator+=( grid2D<gridPars_ > & lhs, const Ytype & rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator-=( grid2D<gridPars_ > & lhs, const Ytype & rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator*=( grid2D<gridPars_ > & lhs, const Ytype & rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator/=( grid2D<gridPars_ > & lhs, const Ytype & rhs );

template<classgridPars_ >
void operator+=( grid2D<gridPars_ > & lhs, const YtypeScalar & rhs );

template<classgridPars_ >
void operator-=( grid2D<gridPars_ > & lhs, const YtypeScalar & rhs );

template<classgridPars_ >
void operator*=( grid2D<gridPars_ > & lhs, const YtypeScalar & rhs );

template<classgridPars_ >
void operator/=( grid2D<gridPars_ > & lhs, const YtypeScalar & rhs );

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator+=( grid2D<Yt1,St1,Xt1 > & lhs, const functor2<Xt2, T2 > & f );

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator-=( grid2D<Yt1,St1,Xt1 > & lhs, const functor2<Xt2, T2 > & f );

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator*=( grid2D<Yt1,St1,Xt1 > & lhs, const functor2<Xt2, T2 > & f );

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator/=( grid2D<Yt1,St1,Xt1 > & lhs, const functor2<Xt2, T2 > & f );


#endif
