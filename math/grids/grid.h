#ifndef GRID_H
#define GRID_H


// Design objectives:
/*

1. A grid is a, discetized, truncated, representation of a function
	from 1,2, 3 or 4 -dimensional real space (represented by class type Xtype)
	to a vector space, represented by the class type, Ytype. The vector 
	space is over the scalar field YtypeScalar.
	
	dimensions will be refered to as (x0,x1,x2,x3)
	
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
	
		f_i = f(x_i), x_i = i*dx+xa, dx = Lx/nx, Lx = xb-xa. 0<=i<N-1,
		
	so a loop through the grid points will look like,
	
		for (int i = 0;i<nx;i++)
			cout << dx*i << endl; // output x_i
			
	so x_i will not quite reach x_b.
	
	Then, mechanisms will be implemented to map to k-space via fourier transforms.
	
	the nyquist frequency will be 1/(2*dx). let g = FFT(f); g is discretized as,
	
	g_j = g(k_j), k_j = dk*j, dk = 2pi/Lx = 2pi/(xb-xa).
	
4. can specify either row or column major. frist dimension ~ row. second dimension ~ column.
	third dimension ~ row of matrices, fourth dimension ~ matrix of matrices, etc.


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
		
	for column major ordering
		#define GRID_COLUMN_MAJOR
	
	for row major ordering (default)
		#define GRID_COLUMN_MAJOR


*/


// some defs for laplaicians

#define GRID_PERIODIC_BD 0
#define GRID_NONPER_BD_1 1


//-------------------------

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>
//#include <fstream>
#include <complex>

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

//#include "mathlib/tools/arrays/array1.h"
//#include "../vectors/euVec.h"

#include <mathlib/math/std_math.h>
#include <mathlib/tools/std_tools.h>

#include <mathlib/math/LC.h>
#include <mathlib/math/norms.h>

#include "grid_convolution.cpp"

#define GRID1D_LOOP( G ) for (int i = 0; i<(G).n1; i++)
#define GRID2D_LOOP( G ) for (int j = 0; j<(G).n2; j++) for (int i = 0; i<(G).n1; i++)
#define GRID3D_LOOP( G ) for (int k = 0; k<(G).n3; k++) for (int j = 0; j<(G).n2; j++) for (int i = 0; i<(G).n1; i++) 

#define GRID1D_LINLOOP( G ) for (int i = 0; i<(G).n1; i++)
#define GRID2D_LINLOOP( G ) for (int i = 0; i<(G).n1* (G).n2; i++)
#define GRID3D_LINLOOP( G ) for (int i = 0; i<(G).n1* (G).n2* (G).n3; i++) 

#ifndef GRID_COLUMN_MAJOR
	#define GRID_ROW_MAJOR
#endif

#ifdef GRID_ROW_MAJOR
	#define GRID1D_ORDER( i ) i
	#define GRID2D_ORDER( i,j ) i*n2+j
	#define GRID3D_ORDER( i,j,k ) n3*(i*n2+j)+k
#endif

#ifdef GRID_COLUMN_MAJOR // check this
	#error "GRID_COLUMN_MAJOR not checked"
	#define GRID1D_ORDER( i ) i
	#define GRID2D_ORDER( i,j ) j*n1+i
	#define GRID3D_ORDER( i,j,k ) n1*(k*n2+j)+i
#endif

#define GRID1D_ARRAY( i ) array[ GRID1D_ORDER(i) ]
#define GRID2D_ARRAY( i,j ) array[ GRID2D_ORDER(i,j) ]
#define GRID3D_ARRAY( i,j,k ) array[ GRID3D_ORDER(i,j,k) ]

#define gridPars_ Ytype,YtypeScalar,Xtype 
#define classgridPars_ class Ytype, class YtypeScalar, class Xtype 


template<classgridPars_ >
class grid_virtual
{
public:
	
	virtual int dimension()=0;
	virtual void debug(bool displayData=true)=0;
};

#endif

