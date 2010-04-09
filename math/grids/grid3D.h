#ifndef GRID3D_H
#define GRID3D_H

#include "grid.h"

template<class Ytype=double, class YtypeScalar=double, class Xtype=double >
class grid3D : public virtual grid_virtual<gridPars_ >
{
public: // members
	
	Ytype * array;
	
	Xtype a1,b1;
	Xtype a2,b2;
	Xtype a3,b3;
	
	int n1;
	int n2;
	int n3;
	
public: // discretization
	
	inline Xtype dx1() { return (b1-a1)/n1; }
	inline Xtype dx2() { return (b2-a2)/n2; }
	inline Xtype dx3() { return (b3-a3)/n3; }
	
	inline Xtype dk1() { return (Xtype)ml_2pi/(b1-a1); }
	inline Xtype dk2() { return (Xtype)ml_2pi/(b2-a2); }
	inline Xtype dk3() { return (Xtype)ml_2pi/(b3-a3); }
	
	inline Xtype x1(const int & i) { return ((b1-a1)/n1)*i+a1; }
	inline Xtype x2(const int & j) { return ((b2-a2)/n2)*j+a2; }
	inline Xtype x3(const int & k) { return ((b3-a3)/n3)*k+a3; }
	
	inline Xtype k1(const int & i) { 
		if (i <= n1/2) return dk1()*i;
		 else return -dk1()*(n1-i);
	}

	inline Xtype k2(const int & j) { 
		if (j <= n2/2) return dk2()*j;
		 else return -dk2()*(n2-j);
	}

	inline Xtype k3(const int & k) { 
		if (k <= n3/2) return dk3()*k;
		 else return -dk3()*(n3-k);
	}	

public: // access
	
	inline Ytype & operator() (const int & i, const int & j, const int & k ) { return GRID3D_ARRAY(i,j,k); }
	inline Ytype & operator[] (const int & i) { return array[ i ]; }
	
	inline Ytype operator() (const double & x1, const double & x2, const double & x3 );
	
public: // constructors, etc
	
	template< class Yt2,class St2,class Xt2 >
	void copy( grid3D<Yt2,St2,Xt2 > const & rhs);
	
	template<class Yt2,class St2,class Xt2 >
	void copyData( grid3D<Yt2,St2,Xt2 > const & rhs);
	
	grid3D():array(0),n1(0),a1(0),b1(0),n2(0),a2(0),b2(0),n3(0),a3(0),b3(0) {};
	
	grid3D(grid3D<gridPars_ > const & rhs);
	
	~grid3D() { if (array) delete array; }
	
	grid3D( int N1, Xtype A1, Xtype B1, int N2, Xtype A2, Xtype B2, int N3, Xtype A3, Xtype B3 );
	
public: // assignment operators	
	
	template<class Yt2,class St2,class Xt2 >
	void operator=( grid3D<Yt2,St2,Xt2 > const & rhs );
	
	void operator=( grid3D<gridPars_ > const & rhs );
	
	void operator=( const Ytype & rhs );
	
	template<class Xt2, class T2 >
	void operator=( const functor3<Xt2, T2 > & rhs );
	
	template<class Xt2, class T2 >
	void operator=( T2 (*rhs)(Xt2,Xt2,Xt2) );
	
public: // other
	
	int dimension() {return 3;}
	
	void debug(bool displayData=false);
	
	void dispData( bool clean=false, YtypeScalar eps=1E-14 );
	
};

#include "grid3D_LC.h"
#include "grid3D_norms.h"
#include "grid3D_mod_assign_op.cpp"
#include "grid3D_misc.h"

#endif

