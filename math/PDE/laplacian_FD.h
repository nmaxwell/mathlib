#ifndef LAPLACIAN_FD_H
#define LAPLACIAN_FD_H

#include <mathlib/math/std_math.h>

#include "laplacian.h"

#define INCLUDE_FD_D2_LR

#include "../kernels/FD_DEL2_CENTERED.h"


class laplacian_2d_fd
{
public:
    complex<double > * ft_ker_ext;
    complex<double > * workspace;
    int n1,n2;
    
    laplacian_2d_fd():ft_ker_ext(0),workspace(0),n1(0),n2(0) {}
    
    ~laplacian_2d_fd() {
        if (ft_ker_ext) fftw_free(ft_ker_ext);
        ft_ker_ext = 0;
        if (workspace) fftw_free(workspace);
        workspace = 0; }
    
public:
    
    void init( int n1_, int n2_, double L1, double L2, int m1, int m2 )
    {
        n1 = n1_;
        n2 = n2_;
        
        double * kernel1_ = ml_alloc<double > ( 2*m1+1 );
        double * kernel2_ = ml_alloc<double > ( 2*m2+1 );
        double * kernel1 = &(kernel1_[m1]);
        double * kernel2 = &(kernel2_[m2]);
        
        double s1 = ((double)(n1))/(L1*L1*n2);
        double s2 = ((double)(n2))/(L2*L2*n1);
        
        for (int k=-m1; k<0; k++ )
            kernel1[k] = s1*FD_DEL2_CENTERED_[m1][-k];
        for (int k=0; k<=m1; k++ )
            kernel1[k] = s1*FD_DEL2_CENTERED_[m1][k];
        
        for (int k=-m2; k<0; k++ )
            kernel2[k] = s1*FD_DEL2_CENTERED_[m2][-k];
        for (int k=0; k<=m2; k++ )
            kernel2[k] = s1*FD_DEL2_CENTERED_[m2][k];
        
        double * ker_ext=0;
        
        laplacian_2d_kernel_extension ( ker_ext, n1, n2, kernel1, m1, m1,  kernel2, m2, m2 );
        
        if (ft_ker_ext != 0) fftw_free( ft_ker_ext );
        ft_ker_ext = 0;
        fft_2d( ker_ext, ft_ker_ext, n1, n2 );
        
        if (workspace != 0) fftw_free( workspace );
        workspace = (complex<double> *)fftw_malloc( sizeof(complex<double>)*n1*(n2/2+1) );
        
        ml_free (kernel1_ );
        ml_free (kernel2_ );
        ml_free (ker_ext );
    }
    
    void execute( double * in, double *& out )
    {
        if ( out == 0 ) out = (double *)fftw_malloc( sizeof(double)*n1*n2 );
        
        fft_2d( in, workspace, n1, n2 );
        
        int N = n1*(n2/2+1);
        for (int k=0; k<N; k++)
            workspace[k] *= ft_ker_ext[k];
        
        ifft_2d( workspace, out, n1, n2 );
    }
};



#endif
