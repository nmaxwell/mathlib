#ifndef GRID_H
#define GRID_H


// Design objectives:
/*

1. A grid is a, discetized, truncated, representation of a function
	from 1,2, 3 or 4 -dimensional real space (represented by class type Xtype)
	to a vector space, represented by the class type, Ytype. The vector 
	space is over the scalar field YtypeScalar.
	
	dimensions will be refered to as (x,y,z,t)
	dimensions will be refered to as (x0,x1,x2,x3)
	
	template parameter dim for dimension
	
2. Xtype, Ytype, YtypeScalar will be passed as template type parameters.
	As little restriction as possible will be put on these types.
	
	Xtype should be a real number type, like float or double, or an abstract
	number class that has the same behaviours as float or double. It is not
	expected for Grid to take complex types as an Xtype.
	
	Ytype will be expected to be a vector type, at a minimum it should meet
	the conditions of being a vector space, explicitly, it will not be expected to have 
	multiplicative inverses, and may not be expected to have a * operation.
	
	YtypeScalar should be a real or complex number type, possibly quarternions.
	

	
3.	For this version of Grid, the discretization will be a regular one,
	with regular fixed spacing between grid points.
	
	for 1-D this will be of the form,
	
		f_i = f(x_i), x_i = i*dx+xa, dx = Lx/nx, Lx = xa-xb. 0<=i<N-1,
		
	so a loop through the grid points will look like,
	
		for (int i = 0;i<N;i++)
			cout << dx*i << endl; // output x_i
			
	so x_i will not quite reach x_b.
	
	Then, mechanisms will be implemented to map to k-space via fourier transforms.
	
	the nyquist frequency will be 1/(2*dx). let g = FFT(f); g is discretized as,
	
	g_j = g(k_j), k_j = dk*j, dk = 1/Lx = 1/(xa-xb).
	
4. Array information.
	Array as 4th template parameter
	5th template parameter, 
	

5. euclidian vector type as template parameter, vecType.


The following are be the specific requirements on template parameter types:
(probably not a complete list...)
	
Xtype - 
	1. type cast from double
	2.	
	
Ytype - 
	1. cout << Ytype;
	2. YtypeScalar norm<Ytype, YtypeScalar > ( Ytype );
	3. operator* (double), /,+.-,+=,-=,*=,/=
	
YtypeScalar - 
	1. sqrt
	3. operator* (double), /,+.-,+=,-=,*=,/= 

vecType - 
	1. vecType(x,y) constructor
	2.
	
arrayType
	1.
		a. for 1-D, Ytype & operator[int i]
		b. for n-D, Ytype & operator(i,j,k,...)
		c. for all, Ytype & operator[int i]; contiguous
		
	2. void resize(n1,n2,n3,...)


*/
//-----------------------------------------------------------






//-----ideas:------------------------------------------------------


/*
grid G1,G2,G3;

FFT(G1,G3);
G3 *= -k^2;
G3 *= autoFilter;
IFFT(G3,G2);
/------------

grid G2 = IFFT*autoFilter*(-k^2)*FT*G1


/---------
#define Del2 IFFT*autoFilter*(-k^2)*FT

grid G2 = Del2*G1

//---------

gridOperator Del2 = IFFT*autoFilter*(-k^2)*FT;

*/


















/*-----------------------

	for double FFT:
		#define	GRID1D_DOUBLE_FFT 
		
	for double DCT:
		#define	GRID1D_DOUBLE_DCT 
		
	
	
		
		
	

	for HDAF del2
		#define GRID1D_DOUBLE_HDAF_DEL2
		#define GRID1D_MP_REAL_HDAF_DEL2

	for FD del2
		#define GRID1D_DOUBLE_FD_DEL2
		#define GRID1D_MP_REAL_FD_DEL2
		
	for 	spectral del2
		#define GRID1D_SPECTRAL_DEL2

*/



#define GRIDLOOP_1D( G ) for (int i = 0; i<G.nx; i++)
#define GRIDLOOP_2D( G ) for (int i = 0; i<G.nx; i++) for (int j = 0; j<G.ny; j++)
#define GRIDLOOP_2D( G ) for (int i = 0; i<G.nx; i++) for (int j = 0; j<G.ny; j++) for (int k = 0; k<G.nz; k++)







#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>
//#include <fstream>
#include <complex>
//#include <vector>
//#include <iomanip>



#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include "mathlib/tools/arrays/array1.h"
#include "../vectors/euVec.h"

#include "mathlib/math/LC.h"
#include "mathlib/math/norms.h"



#define gridPars_ Xtype,Ytype,YtypeScalar,vecType,arrayType
#define classgridPars_ class Xtype, class Ytype, class YtypeScalar, class vecType, class arrayType 




template<classgridPars_ >
class grid_virtual
{
public:
	
	virtual int dimension()=0;
	virtual void debug(bool displayData=true)=0;
};





template<class Xtype=double, class Ytype=double, class YtypeScalar=double, class vecType=euVec<1,Xtype > , class arrayType=array1<1,Ytype > >
class grid1D : public virtual grid_virtual<gridPars_ >
{
public: // basics
	
	int dimension() {return 1;}
	void debug(bool displayData=true);
	
	arrayType * array;
	Xtype xa,xb;
	int nx;
	
	Xtype dx() { return (xb-xa)/nx;}
	Xtype dk() { return (Xtype)1.0/(xb-xa);}
	
	Xtype x_(int i)
	{
		return ((xb-xa)/nx)*i+xa;
	}
	
	Xtype k_(int i)
	{
		return ((Xtype)(1.0)/(xb-xa))*i;
	}
	
	Ytype & operator() (int i)
	{
		return (*array)[i];
	}
	
	Ytype & operator[] (int i) { return (*array)[i]; }
	
public: // constructors, etc

	//void copy();
	
	grid1D():array(0),nx(0) {};
	grid1D(const grid1D<gridPars_ > & rhs);
	~grid1D() { if (array) delete array; }
	
	grid1D( int NX, Xtype Xa, Xtype Xb );
	
public: // simple math

	Ytype integrate_trap();
	YtypeScalar L2norm_trap();
	
public: // operators
	
	
	template<class Xt2,class Yt2,class St2,class Vt2,class At2 >
	void operator=( grid1D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

	void operator=( Ytype rhs );

	template<class Xt2, class T2 >
	void operator=( functor<Xt2, T2 > & rhs );
	
	template<class Xt2, class T2 >
	void operator=( T2 (*rhs)(Xt2) );
	
public: // other

	void dispData( bool clean );	
			
};



template<class Xtype=double, class Ytype=double, class YtypeScalar=double, class vecType=euVec<2,Xtype > , class arrayType=array1<2,Ytype > >
class grid2D : public grid_virtual<gridPars_ >
{
public: // basics
	
	int dimension() {return 2;}
	void debug(bool displayData=true);
	
	arrayType * array;
	Xtype xa,xb;
	Xtype ya,yb;
	int nx;
	int ny;
	
	vecType dx() {return vecType((xb-xa)/nx,(yb-ya)/ny);}
	vecType dk() {return vecType( (Xtype)1.0/(xb-xa),(Xtype)1.0/(yb-ya) );}
	
	vecType x_(int i, int j) {	return vecType( ((xb-xa)/nx)*i+xa, ((yb-ya)/ny)*j+ya ); }
	
	Xtype x0_(int i) {	return ((xb-xa)/nx)*i+xa; }
	Xtype x1_(int j) {	return ((yb-ya)/ny)*j+ya; }
	
	vecType k_(int i, int j) {	return vecType( ((Xtype)(1.0)/(xb-xa))*i, ((Xtype)(1.0)/(yb-ya))*j ); }
	
	Ytype & operator() (int i, int j) { return (*array)(i,j); }
	inline Ytype & operator[] (const int i) { return (*array)[i]; }
		
public: // constructors, etc

	//void copy();
	
	grid2D():array(0),nx(0),ny(0) {};
	grid2D(const grid2D<gridPars_ > & rhs);
	~grid2D() { if (array) delete array; }
	
	grid2D( int NX, Xtype Xa, Xtype Xb, int NY, Xtype Ya, Xtype Yb );
	
public: // simple math

	Ytype integrate_trap();
	YtypeScalar L2norm_trap();
	
public: // operators
	
	
	template<class Xt2,class Yt2,class St2,class Vt2,class At2 >
	void operator=( grid2D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

	void operator=( Ytype rhs );

	template<class Xt2, class T2 >
	void operator=( functor2<Xt2, T2 > & rhs );
	
	template< class T2 >
	void operator=( functor<vecType, T2 > & rhs );	
		
};




template<class Xtype=double, class Ytype=double, class YtypeScalar=double, class vecType=euVec<3,Xtype > , class arrayType=array1<3,Ytype > >
class grid3D : public grid_virtual<gridPars_ >
{
public: // basics
	
	int dimension() {return 2;}
	void debug(bool displayData=true);
	
	arrayType * array;
	Xtype xa,xb;
	Xtype ya,yb;
	Xtype za,zb;
	int nx;
	int ny;
	int nz;
	
	vecType dx() {return vecType((xb-xa)/nx,(yb-ya)/ny,(zb-za)/nz);}
	vecType dk() {return vecType( (Xtype)1.0/(xb-xa),(Xtype)1.0/(yb-ya),(Xtype)1.0/(zb-za) );}
	
	vecType x_(int i, int j, int k)
	{
		return vecType( ((xb-xa)/nx)*i+xa, ((yb-ya)/ny)*j+ya, ((zb-za)/nz)*k+za );
	}
	
	vecType k_(int i, int j, int k)
	{
		return vecType( ((Xtype)(1.0)/(xb-xa))*i, ((Xtype)(1.0)/(yb-ya))*j, ((Xtype)(1.0)/(zb-za))*k );
	}
	
	Ytype & operator() (int i, int j, int k)
	{
		return (*array)(i,j,k);
	}
	
	Ytype & operator[] (int i) { return (*array)[i]; }
	
public: // constructors, etc

	//void copy();
	
	grid3D():array(0),nx(0),ny(0),nz(0) {};
	grid3D(const grid3D<gridPars_ > & rhs);
	~grid3D() { if (array) delete array; }
	
	grid3D( int NX, Xtype Xa, Xtype Xb, int NY, Xtype Ya, Xtype Yb, int NZ, Xtype Za, Xtype Zb );
	
public: // operators
	
	//void operator= ();	
	
public: // simple math

	Ytype integrate_trap();
	YtypeScalar L2norm_trap();
	
	
public: // operators
	
	
	template<class Xt2,class Yt2,class St2,class Vt2,class At2 >
	void operator=( grid3D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

	void operator=( Ytype rhs );

	template<class Xt2, class T2 >
	void operator=( functor3<Xt2, T2 > & rhs );
	
	template< class T2 >
	void operator=( functor<vecType, T2 > & rhs );
		
};


// trasform functions

#ifdef GRID1D_DOUBLE_FFT
	
template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void FFT(
	grid1D<Xt1,double,St1,Vt1,At1 > & in,
	grid1D<Xt2,complex<double>,St2,Vt2,At2 > & out );
	
template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void IFFT(
	grid1D<Xt1,complex<double>,St1,Vt1,At1 > & in,
	grid1D<Xt2,double,St2,Vt2,At2 > & out );

#endif

#ifdef GRID1D_DOUBLE_DCT

template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void DCT4(
	grid1D<Xt1,double,St1,Vt1,At1 > & in,
	grid1D<Xt2,double,St2,Vt2,At2 > & out );
	
template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void IDCT4(
	grid1D<Xt1,double,St1,Vt1,At1 > & in,
	grid1D<Xt2,double,St2,Vt2,At2 > & out );

#endif	


// Laplacians

#ifdef GRID1D_DOUBLE_FD_DEL2

#define GRID_PERIODIC_BD 0
#define GRID_NONPER_BD_1 1


template<class Xt1,class Ytype,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void Del2_FD(
	grid1D<Xt1,Ytype,St1,Vt1,At1 > & in,
	grid1D<Xt2,Ytype,St2,Vt2,At2 > & out,
	int order,
	int bdcond=GRID_NONPER_BD_1 );


template<class Xt1,class Ytype,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void Del2_FD(
	grid2D<Xt1,Ytype,St1,Vt1,At1 > & in,
	grid2D<Xt2,Ytype,St2,Vt2,At2 > & out,
	int order,
	int bdcond=GRID_NONPER_BD_1 );

#endif

#ifdef GRID1D_DOUBLE_HDAF_CONVO_DEL2

#define GRID_PERIODIC_BD 0
#define GRID_NONPER_BD_1 1

template<class Xt1,class Ytype,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void Del2_HDAF_CONVO(
	grid1D<Xt1,Ytype,St1,Vt1,At1 > & in,
	grid1D<Xt2,Ytype,St2,Vt2,At2 > & out,
	int order,
	int bdcond=GRID_NONPER_BD_1 );
	
template<class Xt1,class Ytype,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void Del2_HDAF_CONVO(
	grid2D<Xt1,Ytype,St1,Vt1,At1 > & in,
	grid2D<Xt2,Ytype,St2,Vt2,At2 > & out,
	int order,
	int bdcond=GRID_NONPER_BD_1 );

#endif

#ifdef GRID1D_DOUBLE_SPEC_HDAF1_DEL2

#define GRID1D_DOUBLE_FFT

template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void Del2_FD(
	grid1D<Xt1,double,St1,Vt1,At1 > & in,
	grid1D<Xt2,double,St2,Vt2,At2 > & out,
	int m,
	double eps=1E-12,
	double fudge=1.4 );
	
	

#endif


//	other

template<classgridPars_ >
Ytype max(grid1D<gridPars_ > & rhs, int O=0);

template<classgridPars_ >
int max_loc(grid1D<gridPars_ > & rhs, int O=0);

template<classgridPars_ >
Ytype min(grid1D<gridPars_ > & rhs, int O=0);



// norms

template<classgridPars_ >
YtypeScalar maxNorm(grid1D<gridPars_ > & rhs, int O=0);

template<classgridPars_ >
YtypeScalar maxNorm(grid2D<gridPars_ > & rhs, int O=0);

template<classgridPars_ >
YtypeScalar maxNorm(grid3D<gridPars_ > & rhs, int O=0);


template<classgridPars_ >
YtypeScalar norm(grid1D<gridPars_ > & rhs);

template<classgridPars_ >
YtypeScalar norm2(grid1D<gridPars_ > & rhs);

template<classgridPars_ >
YtypeScalar norm(grid2D<gridPars_ > & rhs);

template<classgridPars_ >
YtypeScalar norm2(grid2D<gridPars_ > & rhs);

template<classgridPars_ >
YtypeScalar norm(grid3D<gridPars_ > & rhs);

template<classgridPars_ >
YtypeScalar norm2(grid3D<gridPars_ > & rhs);



// grid1D operator declarations

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator+=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid1D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator-=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid1D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator*=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid1D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator/=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid1D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator+=( grid1D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator-=( grid1D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator*=( grid1D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator/=( grid1D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
void operator+=( grid1D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator-=( grid1D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator*=( grid1D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator/=( grid1D<gridPars_ > & lhs, YtypeScalar rhs );


template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator+=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Xt2, T2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator-=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Xt2, T2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator*=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Xt2, T2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator/=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Xt2, T2 > & rhs );


template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator+=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, T2 (*rhs)(Xt2) );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator-=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, T2 (*rhs)(Xt2) );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator*=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, T2 (*rhs)(Xt2) );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator/=( grid1D<Xt1,Yt1,St1,Vt1,At1 > & lhs, T2 (*rhs)(Xt2) );

// grid2D operator declarations

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator+=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid2D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator-=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid2D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator*=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid2D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator/=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid2D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator+=( grid2D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator-=( grid2D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator*=( grid2D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator/=( grid2D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
void operator+=( grid2D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator-=( grid2D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator*=( grid2D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator/=( grid2D<gridPars_ > & lhs, YtypeScalar rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator+=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator-=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator*=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator/=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator+=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor2<Xt2, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator-=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor2<Xt2, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator*=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor2<Xt2, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator/=( grid2D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor2<Xt2, T2 > & f );


// grid3D operator declarations

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator+=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid3D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator-=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid3D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator*=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid3D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2,class Yt2,class St2,class Vt2,class At2 >
void operator/=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, grid3D<Xt2,Yt2,St2,Vt2,At2 > & rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator+=( grid3D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator-=( grid3D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator*=( grid3D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
typename boost::disable_if<boost::is_same<Ytype,YtypeScalar >, void>::type
operator/=( grid3D<gridPars_ > & lhs, Ytype rhs );

template<classgridPars_ >
void operator+=( grid3D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator-=( grid3D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator*=( grid3D<gridPars_ > & lhs, YtypeScalar rhs );

template<classgridPars_ >
void operator/=( grid3D<gridPars_ > & lhs, YtypeScalar rhs );


template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator+=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator-=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator*=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class T2 >
void operator/=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor<Vt1, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator+=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor3<Xt2, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator-=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor3<Xt2, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator*=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor3<Xt2, T2 > & f );

template<class Xt1,class Yt1,class St1,class Vt1,class At1, class Xt2, class T2 >
void operator/=( grid3D<Xt1,Yt1,St1,Vt1,At1 > & lhs, functor3<Xt2, T2 > & f );


















// CRAP FROM HERE ON OUT -------------------------


// other

template<class Xtype, class YtypeScalar, class vecType, class arrayType >
void dispData(grid1D<Xtype,double,YtypeScalar,vecType,arrayType > & rhs, bool clean=true )
{

	
	if (rhs.nx >= 1)
		for (int i = 0; i<rhs.nx-1; i++)
				if ( clean && fabs(rhs(i)) < 1E-14 )
					std::cout << 0.0 << ", ";
				else
					std::cout << rhs(i) << ", ";
	if (rhs.nx > 0)
		if ( clean && fabs(rhs(rhs.nx-1)) < 1E-14 )
			std::cout << 0.0 << " ";
		else
			std::cout << rhs(rhs.nx-1) << " ";
	
	cout << endl;
}

template<class Xtype, class YtypeScalar, class vecType, class arrayType >
void dispData(grid1D<Xtype,complex<double>,YtypeScalar,vecType,arrayType > & rhs, bool clean=true )
{

	double eps = 1E-14;
	
	if (rhs.nx >= 1)
	for (int i = 0; i<rhs.nx-1; i++)
	{
		if (clean)
		{
			std::cout << '(';
			if ( fabs(real(rhs(i))) < eps )
				std::cout << 0.0 << ',';
			else
				std::cout << real(rhs(i)) << ',';
				
			if ( fabs(imag(rhs(i))) < eps )
				std::cout << 0.0 << "), ";
			else
				std::cout << imag(rhs(i)) << "), ";
		}
		else
		{
			std::cout << '(';
			std::cout << real(rhs(i)) << ',';
			std::cout << imag(rhs(i)) << "), ";
		}		
	}	
	
	if (rhs.nx > 0)
		if (clean)
		{
			std::cout << '(';
			if ( fabs(real(rhs(rhs.nx-1))) < eps )
				std::cout << 0.0 << ',';
			else
				std::cout << real(rhs(rhs.nx-1)) << ',';
				
			if ( fabs(imag(rhs(rhs.nx-1))) < eps )
				std::cout << 0.0 << "), ";
			else
				std::cout << imag(rhs(rhs.nx-1)) << ")";
		}
		else
		{
			std::cout << '(';
			std::cout << real(rhs(rhs.nx-1)) << ',';
			std::cout << imag(rhs(rhs.nx-1)) << ")";
		}		
	
	cout << endl;
}



//plan.execute( &in(0), &in(1)-&in(0), &out(0), &out(1)-&out(0) );

#include "mathlib/tools/graphing/plot2D.h"

template<class Xtype, class Ytype, class YtypeScalar, class vecType, class arrayType >
void plotGrid1D(grid1D<Xtype,Ytype,YtypeScalar,vecType,arrayType > & G,  plot2D & plot,  color3 & C )
{
	for (int i = 0;i<G.nx-1;i++)
		 plot.ptLine(G.x_(i),G(i),G.x_(i+1),G(i+1),C);
}

template<class Xtype, class Ytype, class YtypeScalar, class vecType, class arrayType >
void plotGrid1D_dots(grid1D<Xtype,Ytype,YtypeScalar,vecType,arrayType > & G,  plot2D & plot, int size,  color3 & C )
{
	for (int i = 0;i<G.nx;i++)
		 plot.ptDot(G.x_(i),G(i),size,C);
}










template<class Xtype, class Ytype, class YtypeScalar, class vecType, class arrayType >
void plotGrid2D_1(grid2D<Xtype,Ytype,YtypeScalar,vecType,arrayType > & G,  char * fname, color3 (*cmap)(double I) )
{	
	plot2D plot(G.xa,G.xb,G.ya,G.yb, G.nx,G.ny);
	
	for (int i = 0; i<plot.NX; i++)
	for (int j = 0; j<plot.NY; j++)
	{
		plot.set_px(i,j,cmap(G(i,j)));
	}
		
	plot.png(fname);
}




double pltG2D_f(double x)
{
	return exp(-pow(2.0*x,8));
}

template<class Xtype, class Ytype, class YtypeScalar, class vecType, class arrayType >
void plotGrid2D_2(grid2D<Xtype,Ytype,YtypeScalar,vecType,arrayType > & G,  plot2D & plot )
{
	//for (int i = 0;i<G.nx-1;i++)
	//	 plot.ptLine(G.x_(i),G(i),G.x_(i+1),G(i+1),C);
	
	plot = c3_white;
	
	color3 background = c3_black;
	
	for (int i = 0; i<plot.NX; i++)
	for (int j = 0; j<plot.NY; j++)
	{
		double x = plot.to_x(i);
		double y = plot.to_y(j);
		
		if (x < G.xa) plot.set_px(i,j,background);
		else if (x >= G.xb) plot.set_px(i,j,background);
		else if (y < G.ya) plot.set_px(i,j,background);
		else if (y >= G.yb) plot.set_px(i,j,background);		
		else
		{
			int ig,jg;
			for (ig = 0; G.x0_(ig)<x; ig++ ) {}
			for (jg = 0; G.x1_(jg)<y; jg++ ) {}
			ig --;
			jg --;
			
			double a,b,c,d;
			
			a = norm( G.x_(ig,jg)-euVec<2,double>(x,y) )/norm(G.dx());
			b = norm( G.x_(ig,jg+1)-euVec<2,double>(x,y) )/norm(G.dx());
			c = norm( G.x_(ig+1,jg)-euVec<2,double>(x,y) )/norm(G.dx());
			d = norm( G.x_(ig+1,jg+1)-euVec<2,double>(x,y) )/norm(G.dx());
			
			color3 color = c3_white;
			
			color *= 
			pltG2D_f(a)*G(ig,jg)+
			pltG2D_f(b)*G(ig,jg+1)+
			pltG2D_f(c)*G(ig+1,jg)+
			pltG2D_f(d)*G(ig+1,jg+1);
					
			plot.set_px(i,j,color);
		}
	
		
	//	for (int k = 0; k<= 
		
	}	
}









#endif

