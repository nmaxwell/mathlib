#ifndef STD_MATH_H
#define STD_MATH_H

#include <complex>
#include <math.h>
#include <arprec/mp_real.h>

#include "../tools/std_tools.h"
#include "std_math_constants.h"
#include "norms.h"
//#include "LC.h"

/*
 * some commonly used functions, etc. 
 */

template<class T >
inline T norm_sinc(T const & x)
{
    if (x != 0.0) return sin(x*ml_pi)/(x*ml_pi);
    else return 1.0;
}

inline double log_factorial(double const & x)
{
    return gamma(x+1); // this is log(gamma(x+1))
}

inline double factorial(double const & x)
{
    return exp(gamma(x+1)); // this is exp(log(gamma(x+1)))
}

inline double binom(double const & n, double const & k)
{   
    return exp( gamma(n+1.0)- gamma(k+1.0)- gamma(n-k+1.0));
}


template<class T >
T max(T * x, int n)
{
    T M = x[0];
    for (int i=0; i<n; i++)
        if (x[i] > M) M = x[i];
    return M; 
}

template<class T >
T min(T * x, int n)
{
    T M = x[0];
    for (int i=0; i<n; i++)
        if (x[i] < M) M = x[i];
    return M; 
}

template<class T >
T max_norm(T * x, int n)
{
    T M = norm(x[0]);
    for (int i=0; i<n; i++)
        if (norm(x[i]) > M) M = norm(x[i]);
    return M;
}

inline double sign(double const & x)
{
	if (x>=0) return x; else return -x;
}

template<class X = double, class Y = X >
class functor
{
public:
	virtual Y operator() (X const & x) const =0;
};

template<class X = double, class Y = double >
class functor2
{
public:
	virtual Y operator() (X const & x, X const & y) const =0;
};

template<class X = double, class Y = double >
class functor3
{
public:
	virtual Y operator() (X const & x, X const & y, X const & z)const =0;
};

template<class X = double, class Y = double >
class functor4
{
public:
	virtual Y operator() (X const & x, X const & y, X const & z, X const & t) const =0;
};

template<class X = double, class Y = X >
class composition : public functor <X,Y >
{
    // f composed with g, so f(g(x))
public:
    functor<X,Y > const * f;
    functor<X,Y > const * g;
    
    composition(functor<X,Y > const & F, functor<X,Y > const & G )
    :f(&F),g(&G) {}
    
    ~composition () {f=0; g=0; }
    
    composition ()
    :f(0),g(0) {}
    
public:
    Y operator() (X const & x) const
    {
        return (*f)((*g)(x));
    }
};



















#endif











