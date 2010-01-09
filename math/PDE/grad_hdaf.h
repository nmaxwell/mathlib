#ifndef LAPLACIAN_HDAF_H
#define LAPLACIAN_HDAF_H


#include <mathlib/math/std_math.h>

#include "grad.h"

#include "../hdaf/hdaf.h"

#define mp_0 ((mp_real)(0.0))
#define mp_1 ((mp_real)(1.0))
#define mp_2 ((mp_real)(2.0))
#define mp_pi (atan(mp_1)*4.0)
#define mp_2pi (mp_pi*2.0)
#define mp_4pi2 (mp_2pi*mp_2pi)
#define mp_sqrt2 (sqrt(mp_2))
#define mp_iu ( complex<mp_real > (mp_0, mp_1 ) )





class grad_2d_hdaf
{
public:
    complex<double > * ft_ker_1;
    complex<double > * ft_ker_2;
    
    complex<double > * workspace_1;
    complex<double > * workspace_2;
    int n1,n2;
    
    grad_2d_hdaf():ft_ker_1(0),ft_ker_2(0),workspace_1(0),workspace_2(0),n1(0),n2(0) {}
    
    ~grad_2d_hdaf() {
        if (ft_ker_1) fftw_free(ft_ker_1);
        ft_ker_1 = 0;
        if (ft_ker_2) fftw_free(ft_ker_2);
        ft_ker_2 = 0;
        if (workspace_1) fftw_free(workspace_1);
        workspace_1 = 0;
        if (workspace_2) fftw_free(workspace_2);
        workspace_2 = 0;  }
    
public:
    
    void init( int n1_, int n2_, double L1, double L2, int m1, int m2, double gamma1, double gamma2 )
    {
        n1 = n1_;
        n2 = n2_;
        
        L1 = fabs(L1);
        L2 = fabs(L2);
        
        mp::mp_init(50);
        
        mp_real sigma1 = (((mp_real)L1)/n1)*sqrt( mp_2*m1+ mp_1)/(mp_pi*gamma1);
        mp_real sigma2 = (((mp_real)L2)/n2)*sqrt( mp_2*m2+ mp_1)/(mp_pi*gamma2);
        
        ml_poly<mp_real > P1;
        make_hdaf_ml_poly(P1,m1 );
        differentiate_hdaf_poly( P1, 1 );
        
        ml_poly<mp_real > P2;
        make_hdaf_ml_poly(P2,m2 );
        differentiate_hdaf_poly( P2, 1 );
        
        int w1 = (int)ceil(hdaf_truncate_point (1E-16, 1E-17, m1, dble(sigma1), 2 )/(L1/n1));
        int w2 = (int)ceil(hdaf_truncate_point (1E-16, 1E-17, m2, dble(sigma2), 2 )/(L2/n2));
        
        double * kernel1_ = ml_alloc<double > ( 2*w1+1 );
        double * kernel2_ = ml_alloc<double > ( 2*w2+1 );
        double * kernel1 = &(kernel1_[w1]);
        double * kernel2 = &(kernel2_[w2]);
        
        mp_real h1 = ((mp_real)L1)/n1;
        mp_real s1 = h1/(mp_sqrt2*sigma1);
        mp_real ss1 = s1*s1;
        mp_real h2 = ((mp_real)L2)/n2;
        mp_real s2 = h2/(mp_sqrt2*sigma2);
        mp_real ss2 = s2*s2;
        mp_real f1 = pow(mp_sqrt2*sigma1,-2)*(h1/(n1*n2));
        mp_real f2 = pow(mp_sqrt2*sigma2,-2)*(h2/(n1*n2));
        
        for (int k=-w1; k<=w1; k++ )
            kernel1[k] = dble( exp(-ss1*(k*k)) *P1(dble(s1*k)) *f1 );
        
        for (int k=-w2; k<=w2; k++ )
            kernel2[k] = dble( exp(-ss2*(k*k)) *P2(dble(s2*k)) *f2 );
        
        double * ker_ext_1=0;
        grad_2d_kernel_extension ( ker_ext_1, n1, n2, kernel1, w1, w1, 1 );
        
        double * ker_ext_2=0;
        grad_2d_kernel_extension ( ker_ext_2, n1, n2, kernel2, w2, w2, 2 );
        
        if (ft_ker_1 != 0) fftw_free( ft_ker_1 );
        ft_ker_1 = 0;
        fft_2d( ker_ext_1, ft_ker_1, n1, n2, 1, 1 );
        
        if (ft_ker_2 != 0) fftw_free( ft_ker_2 );
        ft_ker_2 = 0;
        fft_2d( ker_ext_2, ft_ker_2, n1, n2, 1, 1 );
        
        if (workspace_1 != 0) fftw_free( workspace_1 );
        workspace_1 = (complex<double> *)fftw_malloc( sizeof(complex<double>)*n1*(n2/2+1) );
        
        if (workspace_2 != 0) fftw_free( workspace_2 );
        workspace_2 = (complex<double> *)fftw_malloc( sizeof(complex<double>)*n1*(n2/2+1) );
        
        ml_free (kernel1_ );
        ml_free (kernel2_ );
        ml_free (ker_ext_1 );
        ml_free (ker_ext_2 );
    }
    
    void execute( double * in, double *& out1, double *& out2 )
    {
        if ( out1 == 0 ) out1 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        if ( out2 == 0 ) out2 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        
        int N = n1*(n2/2+1);
        
        fft_2d( in, workspace_1, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] = ft_ker_1[k]*workspace_1[k];
        
        ifft_2d( workspace_2, out1, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] = ft_ker_2[k]*workspace_1[k];
        
        ifft_2d( workspace_2, out2, n1, n2 );
    }
    
    void execute( double * in1, double * in2, double *& out11, double *& out12, double *& out21, double *& out22 )
    {
        if ( out11 == 0 ) out11 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        if ( out12 == 0 ) out12 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        if ( out21 == 0 ) out21 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        if ( out22 == 0 ) out22 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        
        int N = n1*(n2/2+1);
        
        fft_2d( in1, workspace_1, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] = ft_ker_1[k]*workspace_1[k];
        
        ifft_2d( workspace_2, out11, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] = ft_ker_2[k]*workspace_1[k];
        
        ifft_2d( workspace_2, out12, n1, n2 );
        
        
        fft_2d( in2, workspace_1, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] = ft_ker_1[k]*workspace_1[k];
        
        ifft_2d( workspace_2, out21, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] = ft_ker_2[k]*workspace_1[k];
        
        ifft_2d( workspace_2, out22, n1, n2 );
    }
    
    
    void grad_dot( double * in1, double * in2, double *& out )
    {
        if ( out == 0 ) out = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        
        int N = n1*(n2/2+1);
        
        fft_2d( in1, workspace_1, n1, n2 );
        
        fft_2d( in2, workspace_2, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_1[k] = ft_ker_1[k]*workspace_1[k] + ft_ker_2[k]*workspace_2[k];
        
        ifft_2d( workspace_1, out, n1, n2 );
    }
    
    void grad_dot_sym( double * in11, double * in12, double * in22, double *& out1, double *& out2 )
    {
        if ( out1 == 0 ) out1 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        if ( out2 == 0 ) out2 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        
        int N = n1*(n2/2+1);
        
        fft_2d( in11, workspace_1, n1, n2 );
        
        fft_2d( in12, workspace_2, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_1[k] = ft_ker_1[k]*workspace_1[k]  + ft_ker_2[k]*workspace_2[k];
        
        ifft_2d( workspace_1, out1, n1, n2 );
        
        fft_2d( in22, workspace_1, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_1[k] = ft_ker_1[k]*workspace_2[k] + ft_ker_2[k]*workspace_1[k];
        
        ifft_2d( workspace_1, out2, n1, n2 );
    }
    
    void linear_strain_tensor( double * in1, double * in2, double * out11, double * out12, double * out22 )
    {
        if ( out11 == 0 ) out11 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        if ( out12 == 0 ) out12 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        if ( out22 == 0 ) out22 = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        
        int N = n1*(n2/2+1);
        
        fft_2d( in1, workspace_1, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] = ft_ker_1[k]*workspace_1[k];
        
        ifft_2d( workspace_2, out11, n1, n2 );
        
        fft_2d( in2, workspace_2, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_1[k] = (ft_ker_2[k]*workspace_1[k] + ft_ker_1[k]*workspace_2[k])*0.5;
        
        ifft_2d( workspace_1, out12, n1, n2 );
        
        for (int k=0; k<N; k++)
            workspace_2[k] *= ft_ker_2[k];
        
        ifft_2d( workspace_2, out22, n1, n2 );
    }
    
};



#endif
