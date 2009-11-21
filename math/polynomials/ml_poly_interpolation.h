#ifndef ML_POLY_INTERPOLATION
#define ML_POLY_INTERPOLATION

#include "ml_poly.h"

template<class T>
inline T nev_Pij(T const & x,int const & i,int const & j, T const * X, T const * Y )
{
    // function used for neville's algorithm for polynomial interpolation
    //
    
	//assert(i>=0); assert(i<N); assert(j<N);assert(i<=j);
	if (i == j)
		return Y[i];
	else
		return ((x-X[j])*nev_Pij(x,i,j-1,X,Y) +(X[i]-x)*nev_Pij(x,i+1,j,X,Y))/(X[i]-X[j]);
}

template<class T>
T nev_interp(T x, T * X, T * Y, int N)
{
    // neville's algorithm for polynomial interpolation    
    // interpolates at the point x, between points (x_i,y_i) for x_i in X and y_i in Y,
    // and 0 <= i < N
    
	return nev_Pij(x,0,N-1,X,Y);
}




template< class T >
inline const T divided_difference(const int & k, const T * x_, const T * y_ )
{
    // computes the forward divided difference of the k+1 number of 2-tuples (x_i,y_i)
    if (!k)
        return y_[0];
    else
        return (divided_difference(k-1,x_+1,y_+1)-divided_difference(k-1,x_,y_))/(x_[k]-x_[0]);
}

template< class T >
inline const T divided_difference_save(const int & k, const T * x_, const T * y_, T * a_ )
{
    // computes the forward divided difference of the k+1 number of 2-tuples (x_i,y_i),
    // saving the results a_j = [y_0,y_1,..,y_j]
    //
    // usage: if there are n 2-tuples, then call 
    //    divided_difference_save(n-1,x_,y_,a_+n-1); 
    //
    //    where a_ is a pointer to an array of size n, it will be written to in reverse order
    //
    
    if (!k)
    {
        *a_ = y_[0];
        return y_[0];
        
    }
    else
    {
        *a_ = (divided_difference(k-1,x_+1,y_+1)-divided_difference_save(k-1,x_,y_,a_-1))/(x_[k]-x_[0]);
        return *a_;
    }
}

template< class T >
void expand_product_binomial_roots(int degree, T * roots, ml_poly<T > & E)
{
    // takes the polynomial
    //    product( (x-roots_i) , 0 <= i < degree )
    //
    // and expands it into the monomial basis
    //
    
    E.resize(1);    
    ml_poly<T > q(1);
    
    E[0] = -roots[0];
    E[1] = 1.0;
    
    for (int k=1; k<degree; k++)
    {
        q[0] = -roots[k];
        q[1] = 1.0;
        E *= q;
    }    
}

template< class T >
ml_poly<T > * gen_newton_basis(int n, T * roots)
{
    // generatres the newtom basis,
    //    { N_j(x) }, 1 <= j < n
    //
    // where N_j(x) = product( (x-roots_i) , 0 <= i < j )
    //
    // N_0(x) = 1.0, N_1(x) = (x-x_0), ...
    //
    // delete N when done
    
    ml_poly<T > * N = new ml_poly<T > [n];
    
    N[0].resize(0);
    N[0][0] = 1.0;
    
    for (int j=1; j<n; j++)
    {
        N[j].resize(1);
        N[j][0] = -roots[j-1];
        N[j][1] = 1.0;
        
        N[j] *= N[j-1];
    }
    
    return N;
}

template<class T >
void gen_newton_interp_poly(int n, T * x_, T * y_, ml_poly<T > & P )
{
    // generatres the newton interpolation polynomial,
    //  n-1 th degree, for n data points (x_i,y_i)
    //
        
    ml_poly<T > * N = gen_newton_basis(n,x_ );
    T * a_ = new T [n];
    divided_difference_save(n-1,x_,y_,a_+n-1);
    
    P.resize(n-1);
    P = 0;
    
    for (int j=0; j<n; j++)
    {        
        N[j] *= a_[j];
        P += N[j];
    }
    
    delete [] N;
    delete [] a_;
}




template<class D >
ml_poly<D > * gen_orthoNorm_polys( int N, double a, double b )
{
    // generate a set of N, N-1 th degree, polynomials,  
    // using modified gram-schmidt, starting with monomials,
    // using the L2 inner product over the interval (a,b)
    //
    
    ml_poly<D > * V = new ml_poly<D >  [N];
    for (int i = 0; i<N; i++)
    {
        V[i].resize(N-1);
        V[i] = 0.0;
        V[i][i] = 1.0;
    }
    ml_poly<D > p,q;
    
    p = V[0];
    p *= V[0]; 
    V[0] /= mp_real(sqrt(p.integrate(a,b)));
    
    for (int k=1; k<N; k++)
    for (int j = 0; j<k; j++)
    {            
        p = V[j];
        p *= V[k];
        q = V[j];
        q *= p.integrate(a,b);
        V[k] -= q;
        
        p = V[k];
        p *= V[k];
        V[k] /= mp_real(sqrt(p.integrate(a,b)));
    }
    
    return V;
}


#ifdef ARPREC_MPREAL_H

ml_poly<double > * gen_orthoNorm_polys_mp_real_double( int N, double a, double b, int digits=70 )
{
    // generate a set of N, N-1 th degree, polynomials,  
    // using modified gram-schmidt, starting with monomials,
    // using the L2 inner product over the interval (a,b)
    //
    // using mp_reals, converted to double.
    // 
    
    mp::mp_init(digits);
    ml_poly<mp_real > * B_ = gen_orthoNorm_polys<mp_real > (N,a,b);
    ml_poly<double > * B = new ml_poly<double> [N];
    for (int j = 0; j<N; j++)
    {
        B[j].copy(B_[j]);
        ML_POLY_LOOP( B[j] )
            B[j][i] = dble(B_[j][i]);
    }
    delete [] B_;
    return B;    
}

#endif









#endif
