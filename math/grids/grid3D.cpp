#ifndef GRID3D_CPP
#define GRID3D_CPP

#include <typeinfo>

#ifndef GRID_H
	#include "grid.h"
#endif

// grid3D debug()
template<classgridPars_ >
void grid3D<gridPars_ >
::debug(bool displayData )
{	
	std::cout << "debug info: grid3D<";
	std::cout << typeid(Xtype).name() << ", ";
	std::cout << typeid(Ytype).name() << ", ";
	std::cout << typeid(YtypeScalar).name() << " >\n";
	std::cout << "\taddress:" << this << std::endl;
	std::cout << "\tarray:" << array << std::endl;
	std::cout << "\tn1: " << n1 << "\tn2: " << n2 << "\tn3: " << n3 << endl;
	std::cout << "\ta1,b1: " << a1 << ",\t" << b1 << "\n";
	std::cout << "\ta2,b2: " << a2 << ",\t" << b2 << "\n";
	std::cout << "\ta3,b3: " << a3 << ",\t" << b3 << "\n";
	
//	std::cout << "\t[a1,b1): " << "[ " << a1 << " , " << b1 << " )\n";
//	std::cout << "\t[a2,b2): " << "[ " << a2 << " , " << b2 << " )\n";
//	std::cout << "\t[a3,b3): " << "[ " << a3 << " , " << b3 << " )\n";
	
	//std::cout << "\t[a1,b1)X[a2,b3)X[a2,b3): " << "[" << a1 << ',' << b1 << ")X[" << a2 << ',' << b3 << ")X[" << a2 << ',' << b3 << ")\t" << endl;
	//std::cout << "\ta1: " << a1 << "\tb1: " << b1 << "\ta2: " << a2 << "\tb2: " << b2 << "\ta3: " << a3 << "\tb3: " << b3 << std::endl;		
	
	if (displayData)
		dispData();
}

// grid3D dispData()
template<classgridPars_ >
void grid3D<gridPars_ >
::dispData( bool clean, YtypeScalar eps)
{
	
	std::cout << "start data:\n\n";
	if (clean)
	for (int k = 0; k<n3; k++)
	{
		for (int j = 0; j<n2; j++)
		{
			if (n1 >= 1)
				for (int i = 0; i<n1-1; i++)
					std::cout << zeroRound((*this)(i,j,k),eps) << ", ";
			if (n1 > 0)
				std::cout << zeroRound((*this)(n1-1,j,k),eps) << ";\n";
		}
		std::cout << std::endl;
	}
	else
	for (int k = 0; k<n3; k++)
	{
		for (int j = 0; j<n2; j++)
		{
			if (n1 >= 1)
				for (int i = 0; i<n1-1; i++)
					std::cout << (*this)(i,j,k) << ", ";
			if (n1 > 0)
				std::cout << (*this)(n1-1,j,k) << ";\n";
		}
		std::cout << std::endl;
	}
	std::cout << "end data.\n";
}

// grid3D copy
template<classgridPars_ >
template< class Yt2,class St2,class Xt2 >
void grid3D<gridPars_ >
::copy(const grid3D<Yt2,St2,Xt2 > & rhs)
{
	a1 = rhs.a1;
	b1 = rhs.b1;
	a2 = rhs.a2;
	b2 = rhs.b2;
	a3 = rhs.a3;
	b3 = rhs.b3;
	
	if ( n1 != rhs.n1 || n2 != rhs.n2 || n3 != rhs.n3 || array == 0 )
	{
		if (array) delete [] array;
		array = new Ytype [rhs.n1*rhs.n2*rhs.n3];
		
		n1 = rhs.n1;
	    n2 = rhs.n2;
	    n3 = rhs.n3;
	}
	

}

// grid3D copyData
template<classgridPars_ >
template<class Yt2,class St2,class Xt2 >
void grid3D<gridPars_ >
::copyData(const grid3D<Yt2,St2,Xt2 > & rhs)
{
	copy(rhs);
	if ( rhs.array )
	GRID3D_LINLOOP( *this )
		array[i] = rhs.array[i];
}

// grid3D copy constructor
template<classgridPars_ >
grid3D<gridPars_ >
::grid3D(const grid3D<gridPars_ > & rhs)
:array(0),n1(0),a1(0),b1(0),n2(0),a2(0),b2(0),n3(0),a3(0),b3(0)
{	
	copy(rhs);
}

// grid3D basic constructor
template<classgridPars_ >
grid3D<gridPars_ >
::grid3D( int N1, Xtype A1, Xtype B1, int N2, Xtype A2, Xtype B2, int N3, Xtype A3, Xtype B3 )
{
	a1 = A1;
	b1 = B1;
	n1 = N1;
	
	a2 = A2;
	b2 = B2;
	n2 = N2;
	
	a3 = A3;
	b3 = B3;
	n3 = N3;
	
	array = new Ytype [n1*n2*n3];
}

template<classgridPars_ >
void grid3D<gridPars_ >
::operator=( grid3D<gridPars_ > const & rhs )
{
	copyData(rhs);
}

template<classgridPars_ >
template<class Yt2,class St2,class Xt2 >
void grid3D<gridPars_ >
::operator=( grid3D<Yt2,St2,Xt2 > const & rhs )
{
	copy(rhs,true);
}

template<classgridPars_ >
void grid3D<gridPars_ >
::operator=( const Ytype & rhs )
{	
	GRID3D_LINLOOP(*this)
		array[i] = rhs;
}

template<classgridPars_ >
template<class Xt2, class T2 >
void grid3D<gridPars_ >
::operator=( const functor3<Xt2, T2 > & rhs )
{
	Xt2 DX1 = dx1();
	Xt2 DX2 = dx2();
	Xt2 DX3 = dx3();
	
	GRID3D_LOOP(*this)
		GRID3D_ARRAY(i,j,k) = rhs(DX1*i+a1, DX2*j+a2, DX3*k+a3 );
}

template<classgridPars_ >
template<class Xt2, class T2 >
void grid3D<gridPars_ >
::operator=( T2 (*rhs)(Xt2,Xt2,Xt2) )
{
	Xt2 DX1 = dx1();
	Xt2 DX2 = dx2();
	Xt2 DX3 = dx3();
	
	GRID3D_LOOP(*this)
		GRID3D_ARRAY(i,j,k) = rhs(DX1*i+a1, DX2*j+a2, DX3*k+a3 );
}

#endif
