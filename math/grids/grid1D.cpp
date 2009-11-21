#ifndef GRID1D_CPP
#define GRID1D_CPP

#include <typeinfo>

#ifndef GRID_H
	#include "grid.h"
#endif

// grid1D debug()
template<classgridPars_ >
void grid1D<gridPars_ >
::debug(bool displayData )
{	
	std::cout << "debug info: grid1D<";
	std::cout << typeid(Xtype).name() << ", ";
	std::cout << typeid(Ytype).name() << ", ";
	std::cout << typeid(YtypeScalar).name() << " >\n";
	std::cout << "\taddress:" << this << std::endl;
	std::cout << "\tarray:" << array << std::endl;
	std::cout << "\tn1: " << n1 << endl;
	std::cout << "\ta1,b1: " << a1 << ",\t" << b1 << "\n";
	
//	std::cout << "\t[a1,b1): " << "[ " << a1 << " , " << b1 << " )\n";
//	std::cout << "\t[a2,b2): " << "[ " << a2 << " , " << b2 << " )\n";
//	std::cout << "\t[a3,b3): " << "[ " << a3 << " , " << b3 << " )\n";
	
	//std::cout << "\t[a1,b1)X[a2,b3)X[a2,b3): " << "[" << a1 << ',' << b1 << ")X[" << a2 << ',' << b3 << ")X[" << a2 << ',' << b3 << ")\t" << endl;
	//std::cout << "\ta1: " << a1 << "\tb1: " << b1 << "\ta2: " << a2 << "\tb2: " << b2 << "\ta3: " << a3 << "\tb3: " << b3 << std::endl;		
	
	if (displayData)
		dispData();
}

// grid1D dispData()
template<classgridPars_ >
void grid1D<gridPars_ >
::dispData( bool clean, YtypeScalar eps )
{
	std::cout << "start data:\n\n";
	if (clean)
	{
		if (n1 >= 1)
			for (int i = 0; i<n1-1; i++)
				std::cout << zeroRound((*this)(i),eps) << ", ";
		if (n1 > 0)
			std::cout << zeroRound((*this)(n1-1),eps) << ";\n";
	}
	else
	{
		if (n1 >= 1)
			for (int i = 0; i<n1-1; i++)
				std::cout << (*this)(i) << ", ";
		if (n1 > 0)
			std::cout << (*this)(n1-1) << ";\n";
	}
	std::cout << "end data.\n";
}

// linear interpolation
template<classgridPars_ >
inline Ytype grid1D<gridPars_ >::operator() (const double & x)
{
    if ( x > a1 && x <= b1 )
    {
        int i1 = index_down(x);
        int i2 = index_up(x);
        
        if ( fabs(x1(i2)-x1(i1)) > 1E-15 )
            return ((*this)(i2)-(*this)(i1))*(x-x1(i1))/(x1(i2)-x1(i1))+(*this)(i1);
        else
            return (*this)(i1);
    }
    return 0;
}

// grid1D copy
template<classgridPars_ >
template<class Yt2,class St2,class Xt2 >
void grid1D<gridPars_ >
::copy(const grid1D<Yt2,St2,Xt2 > & rhs)
{
	a1 = rhs.a1;
	b1 = rhs.b1;
	
	if ( n1 != rhs.n1 || array == 0 )
	{
		if (array) delete [] array;
		array = new Ytype [rhs.n1];
	}
	
	n1 = rhs.n1;
}

// grid1D copyData
template<classgridPars_ >
template< class Yt2,class St2,class Xt2 >
void grid1D<gridPars_ >
::copyData(const grid1D<Yt2,St2,Xt2 > & rhs)
{
	if ( rhs.array )
		GRID1D_LINLOOP( *this )
			array[i] = rhs.array[i];
}

// grid1D copy constructor
template<classgridPars_ >
grid1D<gridPars_ >
::grid1D(const grid1D<gridPars_ > & rhs)
:array(0),n1(0),a1(0),b1(0)
{	
	copy(rhs);
}

// grid1D basic constructor
template<classgridPars_ >
grid1D<gridPars_ >
::grid1D( int N1, Xtype A1, Xtype B1 )
{
	a1 = A1;
	b1 = B1;
	n1 = N1;
	
	array = new Ytype [n1];
}

template<classgridPars_ >
void grid1D<gridPars_ >
::operator=( grid1D<gridPars_ > const & rhs )
{
	copy(rhs);
	copyData(rhs);
}

template<classgridPars_ >
void grid1D<gridPars_ >
::resize(int n_new)
{
	if ( n_new != n1 || array == 0 )
	{
		if (array) delete [] array;
		array = new Ytype [n_new];
	}
	
	n1 = n_new;
}

template<classgridPars_ >
template<class Yt2,class St2,class Xt2 >
void grid1D<gridPars_ >
::operator=( grid1D<Yt2,St2,Xt2 > const & rhs )
{
	copy(rhs,true);
}

template<classgridPars_ >
void grid1D<gridPars_ >
::operator=( const Ytype & rhs )
{	
	GRID1D_LINLOOP(*this)
		array[i] = rhs;
}

template<classgridPars_ >
template<class Xt2, class T2 >
void grid1D<gridPars_ >
::operator=( const functor<Xt2, T2 > & rhs )
{
	Xt2 DX1 = dx1();
	
	GRID1D_LOOP(*this)
		GRID1D_ARRAY(i) = rhs(DX1*i+a1 );
}

template<classgridPars_ >
template<class Xt2, class T2 >
void grid1D<gridPars_ >
::operator=( T2 (*rhs)(Xt2) )
{
	Xt2 DX1 = dx1();
	
	GRID1D_LOOP(*this)
		GRID1D_ARRAY(i) = rhs(DX1*i+a1 );
}

#endif
