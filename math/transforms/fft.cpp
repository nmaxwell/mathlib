#ifndef FFT_CPP
#define FFT_CPP

#include "fft.h"



int fftw_mz_1d::plan( int req_size, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size[pos] == req_size )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size.push_back( req_size );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in  = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size };
        
        *(plans[pos]) = fftw_plan_many_dft(
            1, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_FORWARD,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_mz_1d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_mz_1d::~fftw_mz_1d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}

int fftw_mz_2d::plan( int req_size_1, int req_size_2, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in  = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2 };
        
        *(plans[pos]) = fftw_plan_many_dft(
            2, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_FORWARD,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_mz_2d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_mz_2d::~fftw_mz_2d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}


int fftw_mz_3d::plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 and size_3[pos] == req_size_3 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        size_3.push_back( req_size_3 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in  = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_size_3*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_size_3*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2, req_size_3 };
        
        *(plans[pos]) = fftw_plan_many_dft(
            3, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_FORWARD,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_mz_3d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_mz_3d::~fftw_mz_3d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}







int ifftw_mz_1d::plan( int req_size, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size[pos] == req_size )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size.push_back( req_size );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in  = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size };
        
        *(plans[pos]) = fftw_plan_many_dft(
            1, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_BACKWARD,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: ifftw_mz_1d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

ifftw_mz_1d::~ifftw_mz_1d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}

int ifftw_mz_2d::plan( int req_size_1, int req_size_2, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in  = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2 };
        
        *(plans[pos]) = fftw_plan_many_dft(
            2, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_BACKWARD,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: ifftw_mz_2d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

ifftw_mz_2d::~ifftw_mz_2d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}


int ifftw_mz_3d::plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 and size_3[pos] == req_size_3 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        size_3.push_back( req_size_3 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in  = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_size_3*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_size_3*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2, req_size_3 };
        
        *(plans[pos]) = fftw_plan_many_dft(
            3, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_BACKWARD,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: ifftw_mz_3d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

ifftw_mz_3d::~ifftw_mz_3d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}








int fftw_r2c_mz_1d::plan( int req_size, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size[pos] == req_size )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size.push_back( req_size );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        double * t_in = (double *)fftw_malloc( sizeof(double)*req_size*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*(req_size/2+1)*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size };
        
        *(plans[pos]) = fftw_plan_many_dft_r2c(
            1, // rank 
            &(sizes[0]),
            1,
            t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_r2c_mz_1d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_r2c_mz_1d::~fftw_r2c_mz_1d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}

int fftw_r2c_mz_2d::plan( int req_size_1, int req_size_2, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        double * t_in = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*(req_size_2/2+1)*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2 };
        
        *(plans[pos]) = fftw_plan_many_dft_r2c(
            2, // rank 
            &(sizes[0]),
            1,
            t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_r2c_mz_2d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_r2c_mz_2d::~fftw_r2c_mz_2d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}

int fftw_r2c_mz_3d::plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 and size_3[pos] == req_size_3 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        size_3.push_back( req_size_3 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        double * t_in = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_size_3*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*(req_size_3/2+1)*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2, req_size_3 };
        
        *(plans[pos]) = fftw_plan_many_dft_r2c(
            3, // rank 
            &(sizes[0]),
            1,
            t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_r2c_mz_3d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_r2c_mz_3d::~fftw_r2c_mz_3d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}







int fftw_c2r_mz_1d::plan( int req_size, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size[pos] == req_size )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size.push_back( req_size );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*(req_size/2+1)*req_in_stride );
        double * t_out = (double *)fftw_malloc( sizeof(double)*(req_size/2+1)*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size };
        
        *(plans[pos]) = fftw_plan_many_dft_c2r(
            1, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            t_out,
            0,
            req_out_stride,
            0,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_c2r_mz_1d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_c2r_mz_1d::~fftw_c2r_mz_1d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}

int fftw_c2r_mz_2d::plan( int req_size_1, int req_size_2, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*(req_size_2/2+1)*req_in_stride );
        double * t_out = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2 };
        
        *(plans[pos]) = fftw_plan_many_dft_c2r(
            2, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            t_out,
            0,
            req_out_stride,
            0,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_c2r_mz_2d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_c2r_mz_2d::~fftw_c2r_mz_2d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}

int fftw_c2r_mz_3d::plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 and size_3[pos] == req_size_3 )
            if ( in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        size_3.push_back( req_size_3 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*(req_size_3/2+1)*req_in_stride );
        double * t_out = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_size_3*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2, req_size_3 };
        
        *(plans[pos]) = fftw_plan_many_dft_c2r(
            3, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            t_out,
            0,
            req_out_stride,
            0,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_c2r_mz_3d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_c2r_mz_3d::~fftw_c2r_mz_3d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}
























int fftw_r2r_mz_1d::plan( int req_size, fftw_r2r_kind req_kind, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size[pos] == req_size )
            if ( kind[pos] == req_kind and in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size.push_back( req_size );
        kind.push_back( req_kind );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        double * t_in  = (double *)fftw_malloc( sizeof(double)*req_size*req_in_stride );
        double * t_out = (double *)fftw_malloc( sizeof(double)*req_size*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        
        *(plans[pos]) = fftw_plan_many_r2r(
            1, // rank
            &(req_size),
            1,
            t_in,
            0,
            req_in_stride,
            0,
            t_out,
            0,
            req_out_stride,
            0,
            &(req_kind),
            FFTW_PLAN_MODE );
        
        if (plans[pos] == NULL )
            std::cerr << "error: fftw_r2r_mz_1d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_r2r_mz_1d::~fftw_r2r_mz_1d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}
    


int fftw_r2r_mz_2d::plan( int req_size_1,int req_size_2, fftw_r2r_kind req_kind, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 )
            if ( kind[pos] == req_kind and in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        kind.push_back( req_kind );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        double * t_in = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_in_stride );
        double * t_out = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2 };
        fftw_r2r_kind kinds[] = { req_kind, req_kind };
        
        *(plans[pos]) = fftw_plan_many_r2r(
            2, // rank
            &(sizes[0]),
            1,
            t_in,
            0,
            req_in_stride,
            0,
            t_out,
            0,
            req_out_stride,
            0,
            &(kinds[0]),
            FFTW_PLAN_MODE );
        
        if (plans[pos] == NULL )
            std::cerr << "error: fftw_r2r_mz_2d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_r2r_mz_2d::~fftw_r2r_mz_2d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}






int fftw_r2r_mz_3d::plan( int req_size_1, int req_size_2, int req_size_3, fftw_r2r_kind req_kind, int req_in_stride, int req_out_stride )
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 and size_3[pos] == req_size_3  )
            if ( kind[pos] == req_kind and in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        size_3.push_back( req_size_2 );
        kind.push_back( req_kind );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        double * t_in = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_size_3*req_in_stride );
        double * t_out = (double *)fftw_malloc( sizeof(double)*req_size_1*req_size_2*req_size_3*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2, req_size_3 };
        fftw_r2r_kind kinds[] = { req_kind, req_kind };
        
        *(plans[pos]) = fftw_plan_many_r2r(
            3, // rank
            &(sizes[0]),
            1,
            t_in,
            0,
            req_in_stride,
            0,
            t_out,
            0,
            req_out_stride,
            0,
            &(kinds[0]),
            FFTW_PLAN_MODE );
        
        if (plans[pos] == NULL )
            std::cerr << "error: fftw_r2r_mz_3d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_r2r_mz_3d::~fftw_r2r_mz_3d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}









int fftw_mz_2d_vec::plan( int req_vec_size, int req_size_1, int req_size_2, int req_in_stride, int req_out_stride );
{
    bool found = false;
    int pos=0;
    
    for ( pos=0; pos<plans.size(); pos++ )
    {
        if ( size_1[pos] == req_size_1 and size_2[pos] == req_size_2 )
            if ( vec_size[pos] == req_vec_size and in_stride[pos] == req_in_stride and out_stride[pos] == req_out_stride )
            {
                found = true;
                break;
            }
    }
    
    if (!found)
    {
        vec_size.push_back( req_vec_size );
        size_1.push_back( req_size_1 );
        size_2.push_back( req_size_2 );
        in_stride.push_back( req_in_stride );
        out_stride.push_back( req_out_stride );
        plans.push_back( new fftw_plan );
        
        pos = plans.size()-1;
        
        complex<double> * t_in  = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_in_stride );
        complex<double> * t_out = (complex<double> *)fftw_malloc( sizeof(fftw_complex)*req_size_1*req_size_2*req_out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        int sizes[] = { req_size_1, req_size_2 };
        
        *(plans[pos]) = fftw_plan_many_dft(
            2, // rank 
            &(sizes[0]),
            1,
            (fftw_complex *)t_in,
            0,
            req_in_stride,
            0,
            (fftw_complex *)t_out,
            0,
            req_out_stride,
            0,
            FFTW_FORWARD,
            FFTW_PLAN_MODE );
        
        if ( plans[pos] == NULL )
            std::cerr << "error: fftw_mz_2d returned null\n";
        
        fftw_free(t_in);
        fftw_free(t_out);
    }
    
    return pos;
}

fftw_mz_2d_vec::~fftw_mz_2d()
{
    for (int pos=0; pos<plans.size(); pos++)
    {
        //  fftw_destroy_plan( *(plans[pos]) );
        //  delete plans[pos];
    }
}











#endif
