#ifndef GRID1D_H
#define GRID1D_H

#include "grid.h"

template<class Ytype=double, class YtypeScalar=double, class Xtype=double >
class grid1D : public virtual grid_virtual<gridPars_ >
{
public: // members
	
	Ytype * array;
	
	Xtype a1,b1;
	
	int n1;
	
public: // discretization
	
	inline Xtype dx1() { return (b1-a1)/n1; }
	
	inline Xtype dk1() { return (Xtype)_2pi/(b1-a1); }
	
	inline Xtype x1(const int & i) { return ((b1-a1)/n1)*i+a1; }
	
	inline Xtype k1(const int & i) { 
		if (i <= n1/2) return dk1()*i;
		 else return -dk1()*(n1-i);
	}

public: // access
	
	inline Ytype & operator() (const int & i) { return GRID1D_ARRAY(i); }
	inline Ytype & operator[] (const int & i) { return array[ i ]; }
	
public: // interpolation
    
    inline Ytype operator() (const double & x);
	
public: // constructors, etc
	
	template<class Yt2,class St2,class Xt2 >
	void copy( grid1D<Yt2,St2,Xt2 > const & rhs);
	
	template<class Yt2,class St2,class Xt2 >
	void copyData( grid1D<Yt2,St2,Xt2 > const & rhs);
	
	grid1D():array(0),n1(0),a1(0),b1(0) {};
	
	grid1D(grid1D<gridPars_ > const & rhs);
	
	~grid1D() { if (array) delete array; }
	
	grid1D( int N1, Xtype A1, Xtype B1 );
	
	void resize(int n_new);
	
public: // assignment operators	
	
	template<class Yt2,class St2,class Xt2 >
	void operator=( grid1D<Yt2,St2,Xt2 > const & rhs );
	
	void operator=( grid1D<gridPars_ > const & rhs );
	
	void operator=( const Ytype & rhs );
	
	template<class Xt2, class T2 >
	void operator=( const functor<Xt2, T2 > & rhs );
	
	template<class Xt2, class T2 >
	void operator=( T2 (*rhs)(Xt2) );
	
public: // other
	
	int dimension() {return 3;}
	
	void debug(bool displayData=false);
	
	void dispData( bool clean=false, YtypeScalar eps=1E-14 );
	
	inline int index_down( double const & x )
	{
	    return (int)floor(((x-a1)/(b1-a1))*n1);
	}
	
	inline int index_up( double const & x )
	{
	    return (int)ceil(((x-a1)/(b1-a1))*n1);
	}
	
};

template<classgridPars_ >
void rectify_fftw_order(grid1D<gridPars_ > & rhs)
{
    for (int i=0; i<=rhs.n1/2-1; i++)
        swap(rhs(i),rhs(i+rhs.n1/2));
}





#include "grid1D_LC.h"
#include "grid1D_norms.h"
#include "grid1D_mod_assign_op.h"
#include "grid1D_misc.h"

#endif

