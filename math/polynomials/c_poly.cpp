#ifndef C_POLY_CPP
#define C_POLY_CPP


#ifndef ML_POLY_H
	#include "ml_poly.h"
#endif


template<class T >
void ml_poly<T >::debug()
{
	cout << "ml_poly, n = " << n << endl;
	if (c_)
	for (int i = 0; i<=n; i++)
		cout << "\t" << i  << ":\t"<< c_[i] << endl;
}

template<class T >
ml_poly<T >::ml_poly()
:c_(0),n(0)
{
}

template<class T >
ml_poly<T >::ml_poly(int n)
:n(n),c_(0)
{
    c_ = new T [n+1];
    #ifdef ml_poly_INIT_ZEROS
        ML_POLY_LOOP(*this)
            c_[i] = 0.0;
    #endif
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly(const ml_poly<T2 > & rhs)
{
	copy(rhs);
	//	copyData(rhs);
}

template<class T >
template<class T2 >
inline void ml_poly<T >::copy(const ml_poly<T2 > & rhs)
{
	n = rhs.n;
	if (c_) delete [] c_;
	c_ = new T [n+1];
}

template<class T >
inline void ml_poly<T >::copyData(const ml_poly<T > & rhs)
{
	ML_POLY_LOOP( *this )
	    c_[i] = rhs.c_[i];
}

template< class T >
inline T & ml_poly<T >::operator[](int const & i)
{
	return c_[i];
}

template< class T >
template< class Xtype >
inline T ml_poly<T >::operator()(Xtype const & x) const
{
	T s = 0.0;
	for (int i = n; i>=0; i--)
	{
		s *= x;
		s += c_[i];
	}
	return s;
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator=(T2 const & rhs)
{
	ML_POLY_LOOP( *this )
		c_[i] = rhs;
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator*=(T2 const & rhs)
{
	ML_POLY_LOOP( *this )
		c_[i] *= rhs;
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator/=(T2 const & rhs)
{
	ML_POLY_LOOP( *this )
		c_[i] /= rhs;
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator+=(T2 const & rhs)
{
	ML_POLY_LOOP( *this )
		c_[i] += rhs;
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator-=(T2 const & rhs)
{
	ML_POLY_LOOP( *this )
		c_[i] -= rhs;
}

template<class T >
void ml_poly<T >::operator=(ml_poly const & rhs)
{
    copy(rhs);
	copyData(rhs);
}


template<class T >
template<class T2 >
void ml_poly<T >::operator=(ml_poly<T2 > const & rhs)
{
    copy(rhs);
	copyData(rhs);
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator+=(ml_poly<T2 > const & rhs)
{
	if (n>=rhs.n)
	ML_POLY_LOOP(rhs)
		c_[i] += rhs.c_[i];
	else 
	{
        T * d_ = new T [rhs.n+1];
        
        ML_POLY_LOOP(rhs)
            d_[i] = rhs.c_[i];
        
        ML_POLY_LOOP(*this)
            d_[i] += c_[i];
               
        if (c_)
            delete [] c_;
        
        c_ = d_;
        n = rhs.n;
	}
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator-=(ml_poly<T2 > const & rhs)
{
	if (n>=rhs.n)
	ML_POLY_LOOP(rhs)
		c_[i] -= rhs.c_[i];
	else 
	{
        T * d_ = new T [rhs.n+1];
        
        ML_POLY_LOOP(rhs)
            d_[i] = -rhs.c_[i];
        
        ML_POLY_LOOP(*this)
            d_[i] += c_[i];
               
        if (c_)
            delete [] c_;
        
        c_ = d_;
        n = rhs.n; 
	}
}

template<class T >
template<class T2 >
inline void ml_poly<T >::operator*=(ml_poly<T2 > const & rhs)
{
    int m = rhs.n;
    ml_poly<T > prod(n+m);
    prod = 0.0;
    for (int k = 0;k<=n+m;k++)
    for (int s = max(k-m,0);s<=min(k,n);s++)
        prod.c_[k] += c_[s]*rhs.c_[k-s];
    (*this) = prod;
}

template< class T >
void ml_poly<T >::resize(int const & n_new)
{
    if (!c_) n = 0;
	
	T * d_ = new T [n_new+1];
	
	if (n && c_)
	for (int i=0; i<=min(n,n_new); i++)
		d_[i] = c_[i];
    
    #ifdef ML_POLY_INIT_ZEROS
        for (int i=min(n,n_new)+1; i<=n_new; i++)
            d_[i] = 0.0;
    #endif
    
    if (c_) delete [] c_;
    
    c_ = d_;
    n = n_new;
}

template< class T >
inline int ml_poly<T >::degree()
{
	return n;
}

template< class T >
inline int ml_poly<T >::degree(int const & n_new)
{
    resize(n_new);
    return n;
}

template<class T >
void ml_poly<T >::differentiate(int d)
{
    int n_new = n-d;
    for (int j = 0; j<d; j++)
    if (n)
    {
        for (int i = 0;i<n;i++)
            c_[i] = c_[i+1]*(T)(i+1);
        n--;
    }
    else c_[0] = 0;
    resize(n_new);
}

template<class T >
void ml_poly<T >::differentiate(ml_poly<T > & D, int d)
{
    D = *this;
    D.differentiate(d);
}

template<class T >
void ml_poly<T >::integrate()
{
    if(n)
    {
        resize(n+1);
        for (int i=n; i>=1; i--)
            c_[i] = c_[i-1]/(T)(i);      
        c_[0] = 0;      
    }
}

template<class T >
void ml_poly<T >::integrate(ml_poly<T > & I)
{
    I.resize(n+1);
    for (int i=n+1; i>=1; i--)
       I.c_[i] = c_[i-1]/(T)(i);
}

template<class T >
template< class Xtype >
T ml_poly<T >::integrate(Xtype const & a, Xtype const & b)
{
    ml_poly<T > P;
    integrate(P);
    return P(b)-P(a);
}








/*
template<class T, class T2 >
typename boost::disable_if<boost::is_same<T2,int >, void>::type
ml_poly<T >::ml_poly<T > (
	    T2 const & c0)
:n(0)
{
    c_ = new T [n+1];
	c_[0]= c0;
}
*/

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1)
:n(1)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2)
:n(2)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
	c_[2]= c2;
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3)
:n(3)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
	c_[2]= c2;
	c_[3]= c3;
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4)
:n(4)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
	c_[2]= c2;
	c_[3]= c3;
	c_[4]= c4;
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5)
:n(5)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
	c_[2]= c2;
	c_[3]= c3;
	c_[4]= c4;
    c_[5]= c5;
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5,
	    T2 const & c6)
:n(6)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
	c_[2]= c2;
	c_[3]= c3;
	c_[4]= c4;
    c_[5]= c5;
    c_[6]= c6;
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5,
	    T2 const & c6,
	    T2 const & c7)
:n(7)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
	c_[2]= c2;
	c_[3]= c3;
	c_[4]= c4;
    c_[5]= c5;
    c_[6]= c6;
 	c_[7]= c7;
}

template<class T >
template<class T2 >
ml_poly<T >::ml_poly (
	    T2 const & c0,
	    T2 const & c1,
	    T2 const & c2,
	    T2 const & c3,
	    T2 const & c4,
	    T2 const & c5,
	    T2 const & c6,
	    T2 const & c7,
	    T2 const & c8)
:n(8)
{
    c_ = new T [n+1];
	c_[0]= c0;
	c_[1]= c1;
	c_[2]= c2;
	c_[3]= c3;
	c_[4]= c4;
    c_[5]= c5;
    c_[6]= c6;
 	c_[7]= c7;
    c_[8]= c8;
}


#endif
