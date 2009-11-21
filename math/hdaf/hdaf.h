#ifndef HDAF_H
#define HDAF_H

#include <arprec/mp_real.h>
#include <mathlib/math/std_math.h>
#include <mathlib/tools/std_tools.h>

#include <mathlib/math/polynomials/ml_poly.h>
#include <mathlib/math/convolution/convolution.h>

#include <boost/math/special_functions/erf.hpp>


#include "../polynomials/laguerre_polys.cpp"

//#include "grids/grid1D.cpp"
//#include "hdaf.h"



struct hdaf_coefficients
{
public:
	
	int mmax;
	double * data;
	
	hdaf_coefficients(int mmax):mmax(mmax),data(0) {
		load();				
	}
	
	hdaf_coefficients(const hdaf_coefficients & rhs):mmax(rhs.mmax),data(0) {
		load();
	}
	
	~hdaf_coefficients() { if (data) delete [] data; }
	
	void operator= (const hdaf_coefficients & rhs) {
		mmax = rhs.mmax;
		load();
	}
	
	inline const double operator() (const int & m, const int & n) const 
	{
	    assert(0 <= n && n <= m);
	    assert(m <= mmax);
	    return data[n+m*(m+1)/2 ];}
	
	void load()
	{
		ifstream in;
		sprintf(fname,"%s/hdaf_coefficients_%d.dat",ML_DATA_DIR,mmax);
		in.open(fname,ios::binary);
		
		if ( !(in.good() && in.is_open()))
		{
			hdaf_coefficients::generate(mmax);
			
			sprintf(fname,"%s/hdaf_coefficients_%d.dat",ML_DATA_DIR,mmax);
			in.open(fname,ios::binary);
			assert(in.good() && in.is_open());
		}
		
		if (data) delete [] data;
		data = new double [ (mmax+1)*(mmax+2)/2 +1];
		in.read((char *)data,((mmax+1)*(mmax+2)/2+1)*8);
		in.close();
	}
	
	static void generate(int M )
	{
		mp::mp_init(100);
		
		ofstream out;
		sprintf(fname,"%s/hdaf_coefficients_%d.dat",ML_DATA_DIR,M);
		out.open(fname, fstream::binary);
		assert(out.good() && out.is_open());
		
		double * data_ = new double [(M+1)*(M+2)/2+1];
	    
		
        for (int m=0; m<=M; m++)
        {
            cout << m << endl;
            static mp_real mp_sqrtpi = sqrt( atan(mp_real(1.0))*4.0 );
            mp_real sum = 0.0;
            mp_real Bn = 0.0;
            mp_real P;
            for (int n=0; n<=m; n++)
            {
                P = 1.0;
                Bn = 1.0;
                for (int k=n+1; k<=m; k++)
                {
                    P *= (mp_real)((mp_real)k-0.5)/(k-n);
                    Bn += P;
                }
                Bn /= gamma((mp_real)(n+1));
                if ( n%2 ) Bn = -Bn;
                
                data_[m*(m+1)/2+n] = dble(Bn/mp_sqrtpi);
            }
            
          /*  mp_real sum = 0.0;
            for (int k=0; k<=m; k++)
            {
                sum = 0;
                for (int n=0; n<=m-k; n++)
                    sum += gamma((mp_real)(2*(n+k)+1))/
                    ( gamma((mp_real)(2*k+1))*gamma((mp_real)(n+1))*gamma((mp_real)(n+k+1))*pow(4,-n) );
                
                data_[m*(m+1)/2+k] = dble(sum/sqrt( atan(mp_real(1.0))*4.0 ))*pow(-1,k);
            }*/
            
            
        }
		
		out.write((char *)(data_),8*(M*(M+1)/2+1));
		out.close();
		delete [] data_;
	}
};

template<class T >
void make_hdaf_ml_poly(ml_poly<T > & P, int m )
{
	// sets P such that 
	//
	// delta_m (x) = exp(-x*x) * P(x)
	//
	// then delta_{m,sig} (x) = 1/(sqrt2*sig)  *  exp(-x*x/(sig*sig*2)) * P(x/(sig*sqrt2));
	//
	
    static int m_max = 100;
    static hdaf_coefficients C(m_max);
    
    if (m > m_max)
    {
        m_max *= 1 << ((int)ceil(log2(m/m_max)));
        assert(m <= m_max); // this is crap, dont care, fix later, 100 enough for now
        C = hdaf_coefficients(m_max);
    }
    
    P.degree(m*2);
    P = 0.0;
    for (int k=0; k<=m; k++)
        P[2*k] = C(m,k);
}

template<class T >
void differentiate_hdaf_poly( ml_poly<T > & P, int d=1 )
{
	// takes the P set by make_hdaf_ml_poly,
	// and makes it the result of differentiating 
	//	delta_m (x) = exp(-x*x) * P(x);
	//	d times, so
	//
	//  delta^(d)_m (x) = exp(-x*x) * P(x);
	
    ml_poly<T > Q;    
    ml_poly<T > gd(0.0,-2.0);
    
    // differentiate d times
    for (int k=0; k<d; k++)
    {
        P.differentiate(Q); // Q = d/dx P
        P *= gd;
        P += Q;
    }
}



template<class T >
void square_hdaf_poly( ml_poly<T> & P )
{
	ml_poly<T > Q;  
	
    // square
    Q = P;
    P *= Q;
}

template<class T >
T hdaf_truncate_error(double x, ml_poly<T > & P, double sigma=ml_sqrt2/2 )
{
	// let f(x) = ( delta^(d)_{m,sig} (x) )^2
	//	then F(x) =  sqrt(  \int_x^\infty f(y) dy )
	//
	// so the fractional errror, in truncating a concolution by delta, to x, is F(x)/F(0)
	//
	// this is what this function returns. uses integration by partsand erfc.
	
	
	x /= ml_sqrt2*sigma;
	int n = P.degree();
    T F_0 = 0.0, F_x = 0.0; // cumulative distributions
    T x2 = x*x;
    T ex2 = exp(-x2*2.0);
    T I_k = 0.0;    // integral steps in int by parts, recursion
    T xp = 0.0; // appropriate powers of x
    
    I_k = ml_sqrt2pi*erfc(x*ml_sqrt2)/4;
    F_x += P[0]*I_k;
    xp = x;
    for (int k=2; k<=n; k+=2)
    {
		I_k = ex2*xp/4 + (I_k*(k-1))/4;
		F_x += P[k]*I_k;
		xp *= x2;
	}
	
    I_k = ex2/4;
    if(n>=1) F_x += P[1]*I_k;
    xp = x2;
    for (int k=3; k<=n; k+=2)
    {
		I_k = ex2*xp/4 + (I_k*(k-1))/4;
		F_x += P[k]*I_k;
		xp *= x2;
	}
	
    I_k = ml_sqrtpi/(ml_sqrt2*2);
    F_0 += I_k*P[0];
    for (int k=2; k<=n; k+=2)
    {
		I_k *= 0.25*(k-1);
		F_0 += I_k*P[k];
	}
	
	I_k = ml_sqrt2/8;
    if(n>=1) F_0 += I_k*P[1];
    for (int k=1; k<=n; k+=2)
    {
		I_k *= 0.5*((k+1)*k);
		F_0 += I_k*P[k];
	}
    //cout << "|" ; clean_print(F_x); cout << "   "; clean_print(F_x);  cout << "   |  "; 
    //cout << "|" ;  clean_print(F_x/F_0);  cout << "   |  ";    
    
    return sqrt(norm(F_x/F_0));
}

double hdaf_truncate_error (double x, int m, double sigma, int d=0 )
{
	// returns the error in truncating the convoltion of a function by delta^(d)_{m,sig} to x,
	// in the L2 sense. error(0) = 1, error(\infty) = 0, decreases monotonically.
	//
	// finding the inverse of this function would be great...
	
	static ml_poly<double > P;
	static int m_ = -1;
	
	if (m != m_)
	{
		m_ = m;
		make_hdaf_ml_poly( P, m );
		differentiate_hdaf_poly( P, d );
		square_hdaf_poly( P );
	}
	
	return hdaf_truncate_error (x, P, sigma );
}

template< class T >
T hdaf_truncate_error_ (double x, int m, double sigma, int d=0 )
{
	// this version for aribtrary type
	// returns the error in truncating the convoltion of a function by delta^(d)_{m,sig} to x,
	// in the L2 sense. error(0) = 1, error(\infty) = 0, decreases monotonically.
	//
	// finding the inverse of this function would be great...
	
	static ml_poly<T > P;
	static int m_ = -1;
	
	if (m != m_)
	{
		m_ = m;
		make_hdaf_ml_poly( P, m );
		differentiate_hdaf_poly( P, d );
		square_hdaf_poly( P );
	}
	
	return hdaf_truncate_error (x, P, sigma );
}

double hdaf_truncate_point (double eps_min, double eps_max, int m, double sigma, int d=0 )
{
	// determines the truncate point for the hdaf kernel, delta^(d)_{m,sig} (x).
	//
	//	finds point so that eps_min < error(x) < eps_max
	//						eps_max > error(x) > eps_min
	
	if (eps_min > eps_max){
		double temp = eps_min;
		eps_min = eps_max;
		eps_max = temp;
	}
	
	// eps_max > eps_min
	
	if (eps_min == eps_max) eps_min -= 1E-8;
    
    static double eps_min_save = 0;
    static double eps_max_save = 0;
    static int m_save = 0;
    static double sigma_save = 0;
    static int d_save = 0;
    static double x_0_save = 0;    
    
	if ( eps_min_save == eps_min && eps_max_save == eps_max && m_save == m && sigma_save == sigma && d_save == d )
        return x_0_save;
    
	mp::mp_init(50);	
	int max_iter1 = 130;
	int max_iter2 = 200;
	double search_factor = 1.5;
	
	double lb,ub;
	double x = 1.0;
	double eps = 0.0;
	int k,count = 0;
	
	// first find upper and lower bounds by searching between powers of search_power
	
	x = 1.0;
	eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;	
	if ( eps_min < eps && eps < eps_max ) {x_0_save = x; return x; } // condition met
	
	if (eps > eps_max) // too low
	{
		// look for lower bound
		for (int k=0; k<max_iter1; k++)
		{
			x *= search_factor;
			eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
			if ( eps_min < eps && eps < eps_max ) { x_0_save = x; return x; } // condition met
			
			if (eps >= eps_max) continue; // still too low
			else { lb = x/search_factor; break; } // found lower bound
		}
		
		// now look for upper bound		
		if (eps < eps_min)
			ub = x; // found
		else for (k=0; k<max_iter1; k++)
		{
			x *= search_factor;
			eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
			if ( eps_min < eps && eps < eps_max ) { x_0_save = x; return x; } // condition met
			
			if (eps >= eps_min) continue; // still too low
			else { ub = x; break; } // found upper bound
		}
	}
	else // too high
	{
		// look for upper bound
		for (k=0; k<max_iter1; k++)
		{
			x /= search_factor;
			eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
			if ( eps_min < eps && eps < eps_max ) { x_0_save = x; return x; } // condition met
			
			if (eps <= eps_min) continue; // still too high
			else { ub = x*search_factor; break; } // found upper bound
		}
		
		// now look for lower bound
		if (eps >= eps_max)
			lb = x; // found
		else for (k=0; k<max_iter1; k++)
		{
			x /= search_factor;
			eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
			if ( eps_min < eps && eps < eps_max ) { x_0_save = x; return x; } // condition met
			
			if (eps < eps_max) continue; // still too low
			else { lb = x; break; } // found upper bound
		}
	}
	//cout << k << endl;
	// found upper and lower bounds, check
	x = lb;
	eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
	if ( eps_min < eps && eps < eps_max ) { x_0_save = x; return x; } // condition met
	assert( eps >= eps_max );
	
	x = ub;
	eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
	if ( eps_min < eps && eps < eps_max ) { x_0_save = x; return x; } // condition met
	assert( eps <= eps_min );
	
	// now do a linear type search on interval [lb,ub].
	
	//	finds point so that eps_min < eps < eps_max
	//						eps_max > eps > eps_min
	
	x = (ub+lb)/2;
	eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
	for (k=0; k<max_iter2; k++)
	{
		if ( eps_min < eps && eps < eps_max ) { 
			//cout << "took this many iterations:  " << count << endl;
			 return x; } // condition met
		if (eps >= eps_max) lb = x; // upper half
		if (eps <= eps_min) ub = x; // lower half
		
		x = (ub+lb)/2;
		eps = dble( hdaf_truncate_error_<mp_real> (x, m, sigma, d ) ); count ++;
	}
	if ( eps_min < eps && eps < eps_max ) {x_0_save = x; return x;}
	
	std::cout << "warning: hdaf_truncate_point unsatisfactory, hdaf.h.\n";
	x_0_save = x; return x;
}








double fast_hdaf( double const & x, int const & m, double const & sigma )
{
    static int mmax = 100;
    static hdaf_coefficients C(mmax);
    
    if (m > mmax) 
    {
        mmax *= 2;
        C.mmax = mmax;
        C.load();
        return fast_hdaf(x,m,sigma);
    }
    
    double Z = x*x/(sigma*sigma*2);
    double sum = 0;
    
    for (int i=m; i>=0; i--)
    {
        sum *= Z;
        sum += C(m,i);
    }
    
    return exp(-Z)*sum/(ml_sqrt2*sigma);
}


template<class Xtype, class Ytype=Xtype >
class hdaf_delta :  public functor<Xtype, Ytype>
{
public:
	Ytype sigma;
	int m;
		
	hdaf_delta(int m_, Ytype sig)
	:m(m_),sigma(sig) {};
	
	Ytype operator() (Xtype const & x) const 
	{
	    Ytype S = 0.0;
	    Ytype Z = x*x/(sigma*sigma*2);
	    for ( int n=0; n<=m; n++)
	        S += laguerre(n,(Ytype)(-0.5),Z);
	    S *= exp(-Z);
	    S /= ml_sq2pi*sigma;
	    return S;
	}
};
/*
class hdaf_delta<double, double > :  public functor<double, double >
{
public:
	double sigma;
	int m;
		
	hdaf_delta(int m_, double sig)
	:m(m_),sigma(sig) {};
	
	double operator() (double const & x) const 
	{
	    mp_real S = 0.0;
	    mp_real Z = x*x/(sigma*sigma*2);
	    for ( int n=0; n<=m; n++)
	        S += laguerre(n,(mp_real)(-0.5),Z);
	    S *= exp(-Z);
	    S /= sq2pi*sigma;
	    return dble(S);
	}
};
*/

template<class T >
double hdaf_bp_interpolate( 
    T * & s_in,  T * & s_out, int n, double sampling_period,
    const int & m_low, const double & low_cut, const int & m_high, const double & high_cut,
    double eps = 1E-3 )
{    
    double high_sigma = sqrt(2*m_high+1)/(high_cut*ml_2pi);
    double low_sigma = sqrt(2*m_low+1)/(low_cut*ml_2pi);
       
    double ierfc_eps = boost::math::erfc_inv( eps );
    double t_max = 2.0*low_sigma*ierfc_eps;
    int k_max = (int)ceil(t_max/sampling_period);
    
    T * ker_ = new T [2*k_max+1];
    T * ker = &(ker_[k_max]);
    
    for (int k = 0; k<=k_max; k++)
        ker[k] = (fast_hdaf( k*sampling_period, m_high, high_sigma )
            -fast_hdaf( k*sampling_period, m_low, low_sigma ))*sampling_period;
    
    for (int k = 1; k<=k_max; k++)
        ker[-k] = (fast_hdaf( -k*sampling_period, m_high, high_sigma )
            -fast_hdaf( -k*sampling_period, m_low, low_sigma ))*sampling_period;
    
    non_periodic_convolution( ker, s_in, s_out, k_max, k_max, n );
}

template<class T >
double hdaf_bp_interpolate_periodic( 
    T * & s_in,  T * & s_out, int n, double sampling_period,
    const int & m_low, const double & low_cut, const int & m_high, const double & high_cut)
{
    double high_sigma = sqrt(2*m_high+1)/(high_cut*ml_2pi);
    double low_sigma = sqrt(2*m_low+1)/(low_cut*ml_2pi);
    
    int k_max = n/2-1;
    
    T * ker_ = new T [2*k_max+1];
    T * ker = &(ker_[k_max]);
    
    for (int k = 0; k<=k_max; k++)
        ker[k] = (fast_hdaf( k*sampling_period, m_high, high_sigma )
            -fast_hdaf( k*sampling_period, m_low, low_sigma ))*sampling_period;
        
    for (int k = 1; k<=k_max; k++)
        ker[-k] = (fast_hdaf( -k*sampling_period, m_high, high_sigma )
            -fast_hdaf( -k*sampling_period, m_low, low_sigma ))*sampling_period;
    
    //hdaf_delta<Ttype,Stype > hdaf1(m_low,sigma1)
    
    periodic_convolution( ker, s_in, s_out, k_max, k_max, n );    
}




double hdaf_bp_interpolate( 
    double * & s_in, int n_in, double in_sampling_period, double in_initial_time,
    double * & s_out, int n_out, double out_sampling_period, double out_initial_time, 
    const int & m_low, const double & low_cut, const int & m_high, const double & high_cut,
    double eps = 1E-3 )
{
    double ierfc_eps = boost::math::erfc_inv( eps );
    
    double high_sigma = sqrt(2*m_high+1)/(high_cut*ml_2pi);
    double low_sigma = sqrt(2*m_low+1)/(low_cut*ml_2pi);
    
    double t_max = 2.0*low_sigma*ierfc_eps;
    
    if (!s_in) s_in = new double [n_in];
    if (!s_out) s_out = new double [n_out];
    
    double t_j,t_k;
    int k_min,k_max;
    
    //non_periodic_convolution(T * ker, T * in, T *& out, int n_l, int n_r, int n )
    
    for (int j=0; j<n_out; j++)
    {
        cout << j << "\t" << n_out << endl;
        
        t_j = out_initial_time + out_sampling_period*j;
        s_out[j] = 0;
        
        k_min = max((int)floor((t_j-t_max-in_initial_time)/in_sampling_period), 0);
        k_max = min((int)ceil((t_j+t_max-in_initial_time)/in_sampling_period)+1, n_in);
        
        for (int k = k_min; k<k_max; k++)
        {
            t_k = in_initial_time + in_sampling_period*k;
            s_out[j] += s_in[k]*fast_hdaf( t_j-t_k, m_high, high_sigma );
            s_out[j] -= s_in[k]*fast_hdaf( t_j-t_k, m_low, low_sigma );
        }
        
        s_out[j] *= out_sampling_period;
    }
}




double hdaf_bp_interpolate( 
    float * & s_in, int n_in, double in_sampling_period, double in_initial_time,
    float * & s_out, int n_out, double out_sampling_period, double out_initial_time, 
    const int & m_low, const double & low_cut, const int & m_high, const double & high_cut,
    double eps_min = 1E-3, double eps_max = 1E-4 )
{
    double high_sigma = sqrt(2*m_high+1)/(high_cut*ml_2pi);    
    double low_sigma = sqrt(2*m_low+1)/(low_cut*ml_2pi);  
    double t_max = hdaf_truncate_point (eps_min, eps_max, m_low, low_sigma, 0 );
    
    if (!s_out) s_out = new float [n_out];    
    
    double t_j,t_k;
    int k_min,k_max;
    
    for (int j=0; j<n_out; j++)
    {
        cout << j << "\t" << n_out << endl;
        
        t_j = out_initial_time + out_sampling_period*j;
        s_out[j] = 0;
        
        k_min = max((int)floor((t_j-t_max-in_initial_time)/in_sampling_period), 0);
        k_max = min((int)ceil((t_j+t_max-in_initial_time)/in_sampling_period)+1, n_in);
        
        for (int k = k_min; k<k_max; k++)
        {            
            t_k = in_initial_time + in_sampling_period*k;
            s_out[j] += s_in[k]*fast_hdaf( t_j-t_k, m_high, high_sigma );
            s_out[j] -= s_in[k]*fast_hdaf( t_j-t_k, m_low, low_sigma );
        }
        
        s_out[j] *= out_sampling_period;
    }
}




class hdaf_deltaHat :  public functor<double, double>
{
public:
	double sigma;
	int m;
	dlognfact lognfact;

		
	hdaf_deltaHat(int m_, double sig, int mmax_ = 2000)
	:m(m_),sigma(sig),lognfact(mmax_) {};
	
	double operator() (double const & k) const 
	{
		if (m >= lognfact.nmax)
			*(const_cast<dlognfact *> (&lognfact)) = dlognfact(lognfact.nmax*2);
		
		double s = 0.0;
		double x = k*k*sigma*sigma/2.0;
		double r = log(fabs(k))*2.0+log(sigma)*2.0;
		if (fabs(k)<1E-12) return 1.0;
		for (int n = 0;n<=m;n++)
			s += exp(-x+r*n-ml_log2*n-lognfact(n));
		return s;
	}
};

template<class Xtype, class Ytype>
class hdaf_delta_nosigma :  public functor<Xtype, Ytype>
{
public:
	int m;
		
	hdaf_delta_nosigma (int m_ )
	:m(m_) {};
	
	Ytype operator() (Xtype const & x) const 
	{
	    Ytype S = 0.0;
	    Ytype Z = x*x;
	    for ( int n=0; n<=m; n++)
	        S += laguerre(n,(Ytype)(-0.5),Z);
	    S *= exp(-Z);
	    S /= ml_sqrtpi;
	    return S;
	}
};










double hdaf_integral(int m,double sigma, double a, double b)
{
/*    a /= sigma*sqrt2;
    b /= sigma*sqrt2;
    
    double B[m+1];
    
    for (int n=0; n<=m; n++)
    {
        B[n] = 0;
        int k = n+1;
        double P = 1.0/(n+1);
        if (k<=m)
            B[n] += P;
        k++;
        for (; k<=m; k++)
        {
            P *= (double)(2*k*k-3*k+1)/(2*k*(k-n-1));
            B[n] += P;
        }
        B[n] /= (2*sqrtpi*dfactorial[n]);
        if (n%2)
            B[n] = -B[n];
    }
    
    double sum_a = 0;
    double sum_b = 0;
    
    for (int i=m; i>=0; i--)
    {
        sum_a *= a*a;
        sum_a += B[i];

        sum_b *= b*b;
        sum_b += B[i];
    }
    
    sum_a *= a*exp(-a*a);
    sum_b *= b*exp(-b*b);
    
    return ((erf(b)-erf(a))/2 +sum_b-sum_a);*/
    
    
    static int m_max = 100;
    static hdaf_coefficients C(m_max);
    
    a /= sigma*ml_sqrt2;
    b /= sigma*ml_sqrt2;
    
    double sum = 0.0;
    
    double a2 = a*a;
    double b2 = b*b;
    double pa = 0.5*a*exp(-a*a);
    double pb = 0.5*b*exp(-b*b);
    
    double I_k = (erf(b)-erf(a))*ml_sqrtpi/2;
    sum += I_k*C(m,0);
    
    for (int k=1; k<=m; k++ )
    {
     //   I_k = (pa-pb)+ I_k*(2*k-1)/2;
        
        I_k = -0.5*( pow(b,2*k-1)*exp(-b*b) - pow(a,2*k-1)*exp(-a*a) ) + (I_k*(2*k-1))/2;
        
        pa *= a2;
        pb *= b2;
        
        sum += C(m,k)*I_k;
    }
    
    return sum;//*sqrt2*sigma;
}


double hdaf_integral_accurate(int m,double sigma, double a, double b)
{
    mp::mp_init(1);
    
    a /= sigma*ml_sqrt2;
    b /= sigma*ml_sqrt2;
    
    mp_real B[m+1];
    
    for (int n=0; n<=m; n++)
    {
        //cout << n << endl;
        B[n] = 0;
        int k = n+1;
        mp_real P = (mp_real)(1.0)/(n+1);
        if (k<=m)
            B[n] += P;
        k++;
        for (; k<=m; k++)
        {
            P *= ((mp_real)(2*k*k-3*k+1))/(2*k*(k-n-1));
            B[n] += P;
        }
        B[n] /= (2*ml_sqrtpi*dfactorial[n]);
        B[n] /= (gamma((mp_real)0.5)*gamma(mp_real(n+1)))*2.0;
        if (n%2)
            B[n] = -B[n];
    }
    
    mp_real sum_a = 0.0;
    mp_real sum_b = 0.0;
    
    for (int i=m; i>=0; i--)
    {
        sum_a *= a*a;
        sum_a += B[i];

        sum_b *= b*b;
        sum_b += B[i];
    }
    
    sum_a *= a*exp(-a*a);
    sum_b *= b*exp(-b*b);
    
    return dble(((erf((mp_real)b)-erf((mp_real)a))/2 +sum_b-sum_a));
}




class hdaf_bp_filter_response :  public functor<double, double>
{
public:
	hdaf_deltaHat hdaf_low;
	hdaf_deltaHat hdaf_high;
			
	hdaf_bp_filter_response(int m_low, double low_cut, int m_high, double high_cut )
    :hdaf_low(m_low, sqrt(2*m_low+1)/(low_cut*ml_2pi) ),hdaf_high(m_high, sqrt(2*m_high+1)/(high_cut*ml_2pi) )
    {
        cout << low_cut << "\t" << high_cut << endl;
    };
	
	double operator() (double const & k) const 
	{
		return hdaf_high(k)-hdaf_low(k);
	}
};










/*



template<class D> D hdaf_delta(D x, int M, D sigma);









class HDAFfilter :  public functor<double, double>
{
public:
	double k0;
	double delta;
	double * lognfact;
	int mmax;
	
	HDAFfilter(int K, double D, double * lognfact_, int M) : k0(K),delta(D),lognfact(lognfact_),mmax(M) {};
	
	double operator() (double k)
	{
		double sigma = sqrt2/(delta);
		int m = ceil(k0*k0/(delta*delta)-0.5);
		assert(m <= mmax);
		double s = 0.0;
		double x = k*k*sigma*sigma/2.0;
		double r = log(fabs(k))*2.0+log(sigma)*2.0;
		if (fabs(k)<1E-12) return 1.0;
		for (int n = 0;n<=m;n++)
			s += exp(-x+r*n-_log2*n-lognfact[n]);
		return s;
	}
};





template<class D> D hdaf_delta(D x, int M, D sigma);

template<> double hdaf_delta(double x, int M, double sigma);

template<> mp_real hdaf_delta(mp_real x, int M, mp_real sigma);

template<class D> D hdaf_delta_2(D x, int M, D sigma);

template<> double hdaf_delta_2(double x, int M, double sigma);

template<> mp_real hdaf_delta_2(mp_real x, int M, mp_real sigma);



*/


/*
void HDAF1(grid1D<double> in, grid1D<double> out, double M, double sig);


double hdafDelta(double x, double M, double sig);
double hdafDelta_2(double x, double M, double sig);
double hdafDeltaHat_M(double xi, double M, double sig);
*/




#endif
