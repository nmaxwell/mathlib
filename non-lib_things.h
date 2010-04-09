
/*
 * Some usefull things which mathlib does not depend on. Include after link.cpp if you use it at all.
 * 
 */

#ifndef iu
    #define iu      (complex<double>(0.0,1.0))
#endif
#ifndef pi
    #define pi      3.141592653589793238463
#endif
#ifndef _2pi
    #define _2pi    6.2831853071795865
#endif
#ifndef _4pi
    #define _4pi    6.2831853071795865*2.0
#endif
#ifndef _4pi2
    #define _4pi2   3.9478417604357434e+1
#endif
#ifndef pi2
    #define pi2     9.8696044010893586
#endif
#ifndef sq2pi
    #define sq2pi   2.5066282746310005
#endif
#ifndef sqrt2pi
    #define sqrt2pi 2.5066282746310005
#endif
#ifndef sqrt2
    #define sqrt2   1.414213562373095
#endif
#ifndef sqrtpi
    #define sqrtpi  1.772453850905516
#endif
#ifndef _log2
    #define _log2   6.931471805599453e-1
#endif
#ifndef logsqrtpi
    #define logsqrtpi   0.572364942924700081938738
#endif
#ifndef logpi
    #define logpi   1.144729885849400174143427
#endif




template< class T >
void output( T * data, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    char str[40];    
    for ( int k=0; k<n; k++ )
    {
        sprintf(str, "%1.18e", data[k] );
        out << str << delim1;
    }    
    
    out.close();
}


template< class T >
void output( vector<T > data, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    int n = data.size();
    
    char str[40];
    
    for ( int k=0; k<n; k++ )
    {
    //    sprintf(str, "%1.18e", data[k][j] );
        out << data[k] << delim2;
    }
    
    out.close();
}

template< class T, int n_lines, int n_cols >
void output( T data[n_lines][n_cols], const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    // not working yet, figure out...
    
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n_lines; k++ )
    {
        for ( int j=0; j<n_cols-1; j++ )
            out << data[k][j] << delim1;
            out << data[k][n_cols-1] << delim2;
    }
    
    out.close();
}



template< class T >
void output( T **data, int n_lines, int n_cols, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    char str[40];
    
    for ( int k=0; k<n_lines; k++ )
    {
        for ( int j=0; j<n_cols-1; j++ )
        {
            sprintf(str, "%1.18e", data[k][j] );
            out << str << delim1;
        }
            sprintf(str, "%1.18e", data[k][n_cols-1] );
            out << data[k][n_cols-1] << delim2;
    }
    
    out.close();
}






template< class T1, class T2 >
void output( T1 * data1, T2 * data2, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data1[k] << delim1 << data2[k] << delim2;
    
    out.close();
}

template< class T1, class T2, class T3 >
void output( T1 * data1, T2 * data2, T3 * data3, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data1[k] << delim1 << data2[k] << delim1 << data3[k] << delim2;
    
    out.close();
}

#define reg_discr regular_discretization

double * regular_discretization(double a, double b, int n)
{
    double * p = ml_alloc<double > ( n );
    for (int k=0; k<n; k++ )
        p[k] = ((b-a)/n)*k + a;
    
    return p;
}

double ml_mean ( double * data, int n )
{
    double mean =0;
    
    for ( int k=0; k<n; k++ )
        mean += data[k];
    
    return mean/n;
}

double std_dev ( double * data, int n )
{
    double mu = 0;
    double sig = 0;
    
    for ( int k=0; k<n; k++ )
        mu += data[k];
    
    for ( int k=0; k<n; k++ )
        sig += data[k]*data[k];
    
    return sqrt((sig-mu*mu/n)/n);
}

# define mean_stdev mean_std_dev
void mean_std_dev ( double * data, int n, double & mu, double & sig )
{
    mu = 0;
    sig = 0;
    
    for ( int k=0; k<n; k++ )
        mu += data[k];
    
    mu /= n;
    
    for ( int k=0; k<n; k++ )
        sig += data[k]*data[k];
    
    sig = sqrt(sig/n-mu*mu);
}

template<class T1, class T2 >
void ml_copy( T1 * x, T2 * y, int n)
{
    // x = y
    for (int k=0; k<n; k++)
        x[k] = y[k];
}


double ml_slope( double * x, double * y, int n, double & b, double & a )
{
    double mx = ml_mean( x, n );
    double my = ml_mean( y, n );
    
    double Sxy = 0.0;
    for ( int k=0; k<n; k++ )
        Sxy += (x[k] - mx )*(y[k] - my );
    
    double Sx = 0.0;
    for ( int k=0; k<n; k++ )
        Sx += (x[k] - mx )*(x[k] - mx );
    
    double Sy = 0.0;
    for ( int k=0; k<n; k++ )
        Sy += (y[k] - my )*(y[k] - my );
    
    b = Sxy/Sx;
    a = my - b*mx;
}





#ifdef include_integration

#include <gsl/gsl_integration.h>



double integrate_R (double (*func)(double,void *), void * arg_ptr=0,  double epsabs=1E-6, double epsrel=1E-10)
{
    double result,real_abserr,imag_result,imag_abserr;
    static gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000000);
    
    gsl_function F;
    F.function = func;
    F.params = arg_ptr;
    
    gsl_integration_qagi (
        &F,
        epsabs,
        epsrel,
        200000,
        workspace,
        &result,
        &real_abserr );
    
    return result;
}


double integrate(double (*func)(double,void *), double a, double b, void * arg_ptr=0,  double epsabs=1E-6, double epsrel=1E-10)
{
    double result,real_abserr,imag_result,imag_abserr;
    static gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000000);
    
    gsl_function F;
    F.function = func;
    F.params = arg_ptr;
    
    gsl_integration_qag (
        &F,
        a,
        b,
        epsabs,
        epsrel,
        200000,
        2,
        workspace,
        &result,
        &real_abserr );
    
    return result;
}



#endif
