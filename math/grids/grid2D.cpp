#ifndef GRID2D_CPP
#define GRID2D_CPP

#include <typeinfo>

#ifndef GRID_H
	#include "grid.h"
#endif

// grid2D debug()
template<classgridPars_ >
void grid2D<gridPars_ >
::debug(bool displayData )
{	
	std::cout << "debug info: grid2D<";
	std::cout << typeid(Xtype).name() << ", ";
	std::cout << typeid(Ytype).name() << ", ";
	std::cout << typeid(YtypeScalar).name() << " >\n";
	std::cout << "\taddress:" << this << std::endl;
	std::cout << "\tarray:" << array << std::endl;
	std::cout << "\tn1: " << n1 << "\tn2: " << n2 << endl;
	std::cout << "\ta1,b1: " << a1 << ",\t" << b1 << "\n";
	std::cout << "\ta2,b2: " << a2 << ",\t" << b2 << "\n";
	
//	std::cout << "\t[a1,b1): " << "[ " << a1 << " , " << b1 << " )\n";
//	std::cout << "\t[a2,b2): " << "[ " << a2 << " , " << b2 << " )\n";
//	std::cout << "\t[a3,b3): " << "[ " << a3 << " , " << b3 << " )\n";
	
	//std::cout << "\t[a1,b1)X[a2,b3)X[a2,b3): " << "[" << a1 << ',' << b1 << ")X[" << a2 << ',' << b3 << ")X[" << a2 << ',' << b3 << ")\t" << endl;
	//std::cout << "\ta1: " << a1 << "\tb1: " << b1 << "\ta2: " << a2 << "\tb2: " << b2 << "\ta3: " << a3 << "\tb3: " << b3 << std::endl;		
	
//	if (displayData)
//		dispData();
}

// grid2D dispData()
/*
template<classgridPars_ >
void grid2D<gridPars_ >
::dispData( bool clean, YtypeScalar eps )
{
	std::cout << "start data:\n\n";
	if (clean)
	for (int j = 0; j<n2; j++)
	{
		if (n1 >= 1)
			for (int i = 0; i<n1-1; i++)
				std::cout << zeroRound((*this)(i,j),eps) << ", ";
		if (n1 > 0)
			std::cout << zeroRound((*this)(n1-1,j),eps) << ";\n";
	}
	else
	for (int j = 0; j<n2; j++)
	{
		if (n1 >= 1)
			for (int i = 0; i<n1-1; i++)
				std::cout << (*this)(i,j) << ", ";
		if (n1 > 0)
			std::cout << (*this)(n1-1,j) << ";\n";
	}
	std::cout << "end data.\n";
}
*/

// grid2D copy
template<classgridPars_ >
template<class Yt2,class St2,class Xt2 >
void grid2D<gridPars_ >
::copy(const grid2D<Yt2,St2,Xt2 > & rhs)
{
	a1 = rhs.a1;
	b1 = rhs.b1;
	a2 = rhs.a2;
	b2 = rhs.b2;
	
	if ( n1 != rhs.n1 || n2 != rhs.n2 || array == 0 )
	{
		if (array) delete [] array;
		array = new Ytype [rhs.n1*rhs.n2];
	}
	
	n1 = rhs.n1;
	n2 = rhs.n2;
}

// grid2D copyData
template<classgridPars_ >
template<class Yt2,class St2,class Xt2 >
void grid2D<gridPars_ >
::copyData(const grid2D<Yt2,St2,Xt2 > & rhs)
{
	copy(rhs);
	if ( rhs.array )
	GRID2D_LINLOOP( *this )
		array[i] = rhs.array[i];
}

// grid2D copy constructor
template<classgridPars_ >
grid2D<gridPars_ >
::grid2D(const grid2D<gridPars_ > & rhs)
:array(0),n1(0),a1(0),b1(0),n2(0),a2(0),b2(0)
{	
	copy(rhs);
}

// grid2D basic constructor
template<classgridPars_ >
grid2D<gridPars_ >
::grid2D( int N1, Xtype A1, Xtype B1, int N2, Xtype A2, Xtype B2 )
{
	a1 = A1;
	b1 = B1;
	n1 = N1;
	
	a2 = A2;
	b2 = B2;
	n2 = N2;
	
    array = ml_alloc<Ytype > (n1*n2);
}

template<classgridPars_ >
void grid2D<gridPars_ >
::operator=( grid2D<gridPars_ > const & rhs )
{
	copy(rhs);
	copyData(rhs);
}

template<classgridPars_ >
template<class Yt2,class St2,class Xt2 >
void grid2D<gridPars_ >
::operator=( grid2D<Yt2,St2,Xt2 > const & rhs )
{
	copy(rhs,true);
}

template<classgridPars_ >
void grid2D<gridPars_ >
::operator=( const Ytype & rhs )
{	
	GRID2D_LINLOOP(*this)
		array[i] = rhs;
}

template<classgridPars_ >
template<class Xt2, class T2 >
void grid2D<gridPars_ >
::operator=( const functor2<Xt2, T2 > & rhs )
{
	Xt2 DX1 = dx1();
	Xt2 DX2 = dx2();
	
	GRID2D_LOOP(*this)
		GRID2D_ARRAY(i,j) = rhs(DX1*i+a1, DX2*j+a2 );
}

template<classgridPars_ >
template<class Xt2, class T2 >
void grid2D<gridPars_ >
::operator=( T2 (*rhs)(Xt2,Xt2) )
{
	Xt2 DX1 = dx1();
	Xt2 DX2 = dx2();
	
	GRID2D_LOOP(*this)
		GRID2D_ARRAY(i,j) = rhs(DX1*i+a1, DX2*j+a2 );
}

#endif
