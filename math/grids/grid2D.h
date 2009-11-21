#ifndef GRID2D_H
#define GRID2D_H

#include "grid.h"

template<class Ytype=double, class YtypeScalar=double, class Xtype=double >
class grid2D : public virtual grid_virtual<gridPars_ >
{
public: // members
	
	Ytype * array;
	
	Xtype a1,b1;
	Xtype a2,b2;
	
	int n1;
	int n2;
	
public: // discretization
	
	inline Xtype dx1() { return (b1-a1)/n1; }
	inline Xtype dx2() { return (b2-a2)/n2; }
	
	inline Xtype dk1() { return (Xtype)_2pi/(b1-a1); }
	inline Xtype dk2() { return (Xtype)_2pi/(b2-a2); }
	
	inline Xtype x1(const int & i) { return ((b1-a1)/n1)*i+a1; }
	inline Xtype x2(const int & j) { return ((b2-a2)/n2)*j+a2; }
	
	inline Xtype k1(const int & i) { 
		if (i <= n1/2) return dk1()*i;
		 else return -dk1()*(n1-i);
	}

	inline Xtype k2(const int & j) { 
		if (j <= n2/2) return dk2()*j;
		 else return -dk2()*(n2-j);
	}
	

public: // access
	
	inline Ytype & operator() (const int & i, const int & j) { return GRID2D_ARRAY(i,j); }
	inline Ytype & operator[] (const int & i) { return array[ i ]; }
	
public: // constructors, etc
	
	template<class Yt2,class St2,class Xt2 >
	void copy( grid2D<Yt2,St2,Xt2 > const & rhs);
	
	template<class Yt2,class St2,class Xt2 >
	void copyData( grid2D<Yt2,St2,Xt2 > const & rhs);
	
	grid2D():array(0),n1(0),a1(0),b1(0),n2(0),a2(0),b2(0) {};
	
	grid2D(grid2D<gridPars_ > const & rhs);
	
	~grid2D() { if (array) delete array; array=0; }
	
	grid2D( int N1, Xtype A1, Xtype B1, int N2, Xtype A2, Xtype B2 );
	
public: // assignment operators	
	
	template<class Yt2,class St2,class Xt2 >
	void operator=( grid2D<Yt2,St2,Xt2 > const & rhs );
	
	void operator=( grid2D<gridPars_ > const & rhs );
	
	void operator=( const Ytype & rhs );
	
	template<class Xt2, class T2 >
	void operator=( const functor2<Xt2, T2 > & rhs );
	
	template<class Xt2, class T2 >
	void operator=( T2 (*rhs)(Xt2,Xt2) );
	
public: // other
	
	int dimension() {return 3;}
	
	void debug(bool displayData=false);
	
	void dispData( bool clean=false, YtypeScalar eps=1E-14 );
	
};


#include "grid2D_LC.h"
#include "grid2D_norms.h"
#include "grid2D_mod_assign_op.h"
#include "grid2D_misc.h"

#endif

