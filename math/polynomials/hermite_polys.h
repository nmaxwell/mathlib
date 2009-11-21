#ifndef HERMITE_POLYS_H
#define HERMITE_POLYS_H


// if want double hermites coefficients,
// 	#define DOUBLE_HERMITE_COEFF

// if want double hermites coefficients,
// 	#define DOUBLE_HERMITE_COEFF_OLD

// if want long long int hermites coefficients,
// 	#define INT64_HERMITE_COEFF

// if want double hermites, evaluated by boost:
// 	#define DOUBLE_HERMITE_BOOST

// if want double hermites, evaluated by laguerre polys and mp_real:
// 	#define DOUBLE_HERMITE_MP_REAL_LAG

// if want mp_real hermites, evaluated by laguerre polys and mp_real:
// 	#define MP_REAL_HERMITE_LAG


#include "ml_poly.h"

template<class T=double>
class hermitePoly
{
public:
	unsigned n;//degree 
	hermitePoly (int N);
	T & operator[](unsigned i);
	
	T operator()(T x); 
};









#ifdef DOUBLE_HERMITE_BOOST
	
	//#include "poly1.cpp"
	#include "hermite_polys.h"
	#include <boost/math/special_functions/hermite.hpp>
	
	double hermite(int n,double x)
	{
		return boost::math::hermite<double>(n,x);
	}
	
#endif

#ifdef DOUBLE_HERMITE_MP_REAL_LAG
	
	#include <arprec/mp_real.h>
	#include <arprec/mp_complex.h>
	#include <arprec/mp_int.h>
	
	#include "laguerre_polys.h"
	
	double hermite(int n, double x)
	{
		if (!(n%2)) // even
			return dble(laguerre<mp_real>(n/2, -0.5, x*x  )* gamma((mp_real)(n/2+1))*pow(mp_real(-4.0),n/2));
		else // odd
			return dble(laguerre<mp_real>((n-1)/2, +0.5, x*x  )* gamma((mp_real)( ((n-1)/2)+1))*pow(mp_real(-4.0),(n-1)/2)*2.0*x);
	}
	
#endif

#ifdef MP_REAL_HERMITE_LAG
	
	#include <arprec/mp_real.h>
	#include <arprec/mp_complex.h>
	#include <arprec/mp_int.h>
	
	#include "laguerre_polys.h"
	
	mp_real hermite(int n, mp_real x)
	{
		if (!(n%2)) // even
			return (laguerre<mp_real>(n/2, -0.5, x*x  )* gamma((mp_real)(n/2+1))*pow(mp_real(-4.0),n/2));
		else // odd
			return (laguerre<mp_real>((n-1)/2, +0.5, x*x  )* gamma((mp_real)( ((n-1)/2)+1))*pow(mp_real(-4.0),(n-1)/2)*2.0*x);
	}

#endif	




#ifdef DOUBLE_HERMITE_COEFF
	
	#include <arprec/mp_real.h>
	#include <arprec/mp_complex.h>
	#include <arprec/mp_int.h>
	
	#include <vector>
	#include "poly1.cpp"
	#include "hermite_polys.h"
	
	double hermiteCoefficient(int n, int m)
		// mth coefficient of nth hermite
	{
		mp::mp_init(30);
		static mp_real two(2.0);
		static mp_real h;
		static mp_real n_;
		static mp_real m_;
						
		if (!(n%2)) // even n
		{
			if (!(m%2))// even m
			{			
				n_ = n;
				m_ = m;
				h = gamma(n_+1)*pow(two,m_)/(gamma(m_+1)*gamma((n_-m_)/2.0+1));
				if ((((n+m)/2)%2))
					return -dble(h);
				else return dble(h);
			}
			else return 0.0;
		}
		else // odd n
		{
			if (!(m%2))// even m
				return 0.0;
			else // odd m
			{			
				n_ = n;
				m_ = m;
				h = gamma(n_+1)*pow(two,m_)/(gamma(m_+1)*gamma((n_-m_)/2.0+1));
				if ((((n+m)/2)%2))
					return dble(h);
				else return -dble(h);
			}
		}
	}
	
	vector< vector<double> > hermite_c_;
	
	void set_hermite(int N)
	{
		static int N_ = 0;
		if (N_ < N)
			for (int n = N_;n<=N;n++)
			{
				vector<double> p;
				for (int m = 0;m<=n;m++)
					p.push_back( hermiteCoefficient(n,m) );
				hermite_c_.push_back(p);
			}		
	}
	
#endif

















#endif





