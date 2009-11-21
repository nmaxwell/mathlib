#ifndef ML_POLY_H
#define ML_POLY_H

#define ML_POLY_LOOP( P ) for (int i = 0; i<=(P).n; i++)

#ifndef ML_POLY_NOT_INIT_ZEROS
    #define ML_POLY_INIT_ZEROS
#endif

//#ifndef ml_poly_TYPE_NAME
 //   #define ml_poly_TYPE_NAME ml_poly
//endif

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

template<class T=double >
class ml_poly
{
public:
	int n; // degree
	T * c_; // coefficients
	
	void debug();

public:
	ml_poly();
	ml_poly( int n);
    
    ~ml_poly( )
    { if(c_) delete [] c_; }
    
	template<class T2 > ml_poly(const ml_poly<T2 > & rhs);
	
	template<class T2 > inline void copy(const ml_poly<T2 > & rhs);
	inline void copyData(const ml_poly<T > & rhs);

public:
    
    inline T & operator[](int const & i);
    
    template< class Xtype >
	inline T operator() (Xtype const & x) const;
    
	
public:
    
    template<class T2 > inline void operator=(T2 const & rhs);
	template<class T2 > inline void operator*=(T2 const & rhs);
	template<class T2 > inline void operator/=(T2 const & rhs);
	template<class T2 > inline void operator+=(T2 const & rhs);
	template<class T2 > inline void operator-=(T2 const & rhs);
	
	void operator=(ml_poly const & rhs);
	template< class T2 > void operator=(ml_poly<T2 > const & rhs);
	template< class T2 > inline void operator+=(ml_poly<T2 > const & rhs);
	template< class T2 > inline void operator-=(ml_poly<T2 > const & rhs);
	template< class T2 > inline void operator*=(ml_poly<T2 > const & rhs);
	
public:
    
    void resize(int const & n);
    inline int degree();    
    inline int degree(int const & n);
    
public:
    
    void differentiate(int d=1);
    void differentiate(ml_poly<T > & D, int d=1);
    
    
    void integrate();
    void integrate(ml_poly<T > & I);
    template< class Xtype >
	T integrate(Xtype const & a, Xtype const & b);
	
public:

 /*   template<class T2 >
    typename boost::disable_if<boost::is_same<T2,int >, void>::type
    ml_poly( T2 const & c0);*/
	    	
	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1);	
	
	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2);

	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3);	
	
	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4);
	
	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5);	

	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5,
	    T2 const & c6);

	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5,
	    T2 const & c6,
	    T2 const & c7);
	
	template<class T2 > ml_poly<T >(
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5,
	    T2 const & c6,
	    T2 const & c7,
	    T2 const & c8);
	    
};












#endif

