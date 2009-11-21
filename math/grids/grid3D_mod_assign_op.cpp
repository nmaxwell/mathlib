#ifndef GRID3D_MOD_ASSIGN_OP_CPP
#define GRID3D_MOD_ASSIGN_OP_CPP

#ifndef GRID3D_MOD_ASSIGN_OP_H
	#include "grid3D_mod_assign_op.h"
#endif




template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator+=( grid3D<Yt1,St1,Xt1 > & lhs, grid3D<Yt2,St2,Xt2 > & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] += rhs.array[i];  
}

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator-=( grid3D<Yt1,St1,Xt1 > & lhs, grid3D<Yt2,St2,Xt2 > & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] -= rhs.array[i];  
}

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator*=( grid3D<Yt1,St1,Xt1 > & lhs, grid3D<Yt2,St2,Xt2 > & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] *= rhs.array[i];  
}

template<class Yt1,class St1,class Xt1, class Yt2,class St2,class Xt2  >
void operator/=( grid3D<Yt1,St1,Xt1 > & lhs, grid3D<Yt2,St2,Xt2 > & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] /= rhs.array[i];  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator+=( grid3D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] += rhs;  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator-=( grid3D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] -= rhs;  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator*=( grid3D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] *= rhs;  
}

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator/=( grid3D<gridPars_ > & lhs, const Ytype & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] /= rhs;  
}

template<classgridPars_ >
void operator+=( grid3D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] += rhs;
}

template<classgridPars_ >
void operator-=( grid3D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] -= rhs;  
}

template<classgridPars_ >
void operator*=( grid3D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] *= rhs;
}

template<classgridPars_ >
void operator/=( grid3D<gridPars_ > & lhs, const YtypeScalar & rhs )
{
	GRID3D_LINLOOP( lhs )
		lhs.array[i] /= rhs;  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator+=( grid3D<Yt1,St1,Xt1 > & lhs, const functor3<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	Xt2 DX2 = lhs.dx2();
	Xt2 DX3 = lhs.dx3();
	
	GRID3D_LOOP( lhs )
		lhs(i,j,k) += f(DX1*i+lhs.a1,DX2*j+lhs.a2,DX3*k+lhs.a3);  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator-=( grid3D<Yt1,St1,Xt1 > & lhs, const functor3<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	Xt2 DX2 = lhs.dx2();
	Xt2 DX3 = lhs.dx3();
	
	GRID3D_LOOP( lhs )
		lhs(i,j,k) -= f(DX1*i+lhs.a1,DX2*j+lhs.a2,DX3*k+lhs.a3);  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator*=( grid3D<Yt1,St1,Xt1 > & lhs, const functor3<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	Xt2 DX2 = lhs.dx2();
	Xt2 DX3 = lhs.dx3();
	
	GRID3D_LOOP( lhs )
		lhs(i,j,k) *= f(DX1*i+lhs.a1,DX2*j+lhs.a2,DX3*k+lhs.a3);  
}

template<class Yt1,class St1,class Xt1,class Xt2, class T2 >
void operator/=( grid3D<Yt1,St1,Xt1 > & lhs, const functor3<Xt2, T2 > & f )
{
	Xt2 DX1 = lhs.dx1();
	Xt2 DX2 = lhs.dx2();
	Xt2 DX3 = lhs.dx3();
	
	GRID3D_LOOP( lhs )
		lhs(i,j,k) /= f(DX1*i+lhs.a1,DX2*j+lhs.a2,DX3*k+lhs.a3);  
}



#endif
