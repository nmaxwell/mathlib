#ifndef GRID1D_MOD_ASSIGN_OP_CPP
#define GRID1D_MOD_ASSIGN_OP_CPP

#ifndef GRID1D_MOD_ASSIGN_OP_H
	#include "grid1D_mod_assign_op.h"
#endif

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator+=( grid1D<Yt1,St1,Xt1 > & lhs, grid1D<Yt2,St2,Xt2 > & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] += rhs.array[i];  
}

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator-=( grid1D<Yt1,St1,Xt1 > & lhs, grid1D<Yt2,St2,Xt2 > & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] -= rhs.array[i];  
}

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator*=( grid1D<Yt1,St1,Xt1 > & lhs, grid1D<Yt2,St2,Xt2 > & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] *= rhs.array[i];  
}

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator/=( grid1D<Yt1,St1,Xt1 > & lhs, grid1D<Yt2,St2,Xt2 > & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] /= rhs.array[i];  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator+=( grid1D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] += rhs;  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator-=( grid1D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] -= rhs;  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator*=( grid1D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] *= rhs;  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator/=( grid1D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] /= rhs;  
}

template<classgridPars_ >
void operator+=( grid1D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] += rhs;
}

template<classgridPars_ >
void operator-=( grid1D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] -= rhs;  
}

template<classgridPars_ >
void operator*=( grid1D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] *= rhs;
}

template<classgridPars_ >
void operator/=( grid1D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID1D_LINLOOP( lhs )
		lhs.array[i] /= rhs;  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator+=( grid1D<Yt1,St1,Xt1 > & lhs, const functor<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	
	GRID1D_LOOP( lhs )
		lhs.array[i] += f(DX1*i+lhs.a1);  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator-=( grid1D<Yt1,St1,Xt1 > & lhs, const functor<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	
	GRID1D_LOOP( lhs )
		lhs.array[i] -= f(DX1*i+lhs.a1);  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator*=( grid1D<Yt1,St1,Xt1 > & lhs, const functor<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	
	GRID1D_LOOP( lhs )
		lhs.array[i] *= f(DX1*i+lhs.a1);  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator/=( grid1D<Yt1,St1,Xt1 > & lhs, const functor<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	
	GRID1D_LOOP( lhs )
		lhs.array[i] /= f(DX1*i+lhs.a1);  
}



#endif
