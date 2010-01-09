#ifndef FFT_H
#define FFT_H

#include <fftw3.h>

#include <mathlib/math/std_math.h>

#ifndef N_FFT_THREADS
	#define N_FFT_THREADS  N_NATIVE_THREADS
#endif

#ifndef FFTW_PLAN_MODE
	#define FFTW_PLAN_MODE  FFTW_ESTIMATE
	// #define FFTW_PLAN_MODE FFTW_PATIENT
#endif


/*
http://www.fftw.org/fftw3_doc/Real_002dto_002dReal-Transform-Kinds.html#Real_002dto_002dReal-Transform-Kinds
* 
* 
# FFTW_R2HC computes a real-input DFT with output in “halfcomplex” format, i.e. real and imaginary parts for a transform of size n stored as:

r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
(Logical N=n, inverse is FFTW_HC2R.)
# FFTW_HC2R computes the reverse of FFTW_R2HC, above. (Logical N=n, inverse is FFTW_R2HC.)
# FFTW_DHT computes a discrete Hartley transform. (Logical N=n, inverse is FFTW_DHT.)
# FFTW_REDFT00 computes an REDFT00 transform, i.e. a DCT-I. (Logical N=2*(n-1), inverse is FFTW_REDFT00.)
# FFTW_REDFT10 computes an REDFT10 transform, i.e. a DCT-II (sometimes called “the” DCT). (Logical N=2*n, inverse is FFTW_REDFT01.)
# FFTW_REDFT01 computes an REDFT01 transform, i.e. a DCT-III (sometimes called “the” IDCT, being the inverse of DCT-II). (Logical N=2*n, inverse is FFTW_REDFT=10.)
# FFTW_REDFT11 computes an REDFT11 transform, i.e. a DCT-IV. (Logical N=2*n, inverse is FFTW_REDFT11.)
# FFTW_RODFT00 computes an RODFT00 transform, i.e. a DST-I. (Logical N=2*(n+1), inverse is FFTW_RODFT00.)
# FFTW_RODFT10 computes an RODFT10 transform, i.e. a DST-II. (Logical N=2*n, inverse is FFTW_RODFT01.)
# FFTW_RODFT01 computes an RODFT01 transform, i.e. a DST-III. (Logical N=2*n, inverse is FFTW_RODFT=10.)
# FFTW_RODFT11 computes an RODFT11 transform, i.e. a DST-IV. (Logical N=2*n, inverse is FFTW_RODFT11.) 

*/





    /* 
     * 
     * memoizes fftw plans,
     * pretty self explanatory
     * 
     * fftw_mz_1d
     * ifftw_mz_2d
     * fftw_mz_2d
     * ifftw_mz_2d
     * fftw_mz_3d
     * ifftw_mz_3d
     * 
     * fftw_r2c_mz_1d
     * fftw_c2r_mz_1d
     * fftw_r2c_mz_2d
     * fftw_c2r_mz_2d
     * fftw_r2c_mz_3d
     * fftw_c2r_mz_3d
     * 
     * fftw_r2r_mz_1d
     * fftw_r2r_mz_2d
     * fftw_r2r_mz_3d
     * 
     * 
     * 
     */




struct fftw_mz_1d
{
    
    vector<int > size;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_mz_1d():size(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_mz_1d();
    
    int plan( int req_size, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft( *(plans[pos]), (fftw_complex *)in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};

struct ifftw_mz_1d
{
    
    vector<int > size;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    ifftw_mz_1d():size(),in_stride(),out_stride(),plans() {  }
    
    ~ifftw_mz_1d();
    
    int plan( int req_size, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft( *(plans[pos]), (fftw_complex *)in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};




struct fftw_mz_2d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_mz_2d():size_1(),size_2(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_mz_2d();
    
    int plan( int req_size_1, int req_size_2, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft( *(plans[pos]), (fftw_complex *)in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};

struct ifftw_mz_2d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    ifftw_mz_2d():size_1(),size_2(),in_stride(),out_stride(),plans() {  }
    
    ~ifftw_mz_2d();
    
    int plan( int req_size_1, int req_size_2, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft( *(plans[pos]), (fftw_complex *)in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};


struct fftw_mz_3d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > size_3;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_mz_3d():size_1(),size_2(),size_3(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_mz_3d();
    
    int plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft( *(plans[pos]), (fftw_complex *)in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};

struct ifftw_mz_3d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > size_3;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    ifftw_mz_3d():size_1(),size_2(),size_3(),in_stride(),out_stride(),plans() {  }
    
    ~ifftw_mz_3d();
    
    int plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft( *(plans[pos]), (fftw_complex *)in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};





struct fftw_r2c_mz_1d
{
    
    vector<int > size;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_r2c_mz_1d():size(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_r2c_mz_1d();
    
    int plan( int req_size, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, double *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft_r2c( *(plans[pos]), in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};


struct fftw_c2r_mz_1d
{
    
    vector<int > size;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_c2r_mz_1d():size(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_c2r_mz_1d();
    
    int plan( int req_size, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, double *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft_c2r( *(plans[pos]), (fftw_complex *)in, out );
        else return 1;
        
        return 0;
    }
};


struct fftw_r2c_mz_2d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_r2c_mz_2d():size_1(),size_2(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_r2c_mz_2d();
    
    int plan( int req_size_1, int req_size_2, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, double *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft_r2c( *(plans[pos]), in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};


struct fftw_c2r_mz_2d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_c2r_mz_2d():size_1(),size_2(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_c2r_mz_2d();
    
    int plan( int req_size_1, int req_size_2, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, double *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft_c2r( *(plans[pos]), (fftw_complex *)in, out );
        else return 1;
        
        return 0;
    }
};


struct fftw_r2c_mz_3d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > size_3;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_r2c_mz_3d():size_1(),size_2(),size_3(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_r2c_mz_3d();
    
    int plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, double *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft_r2c( *(plans[pos]), in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};


struct fftw_c2r_mz_3d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > size_3;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_c2r_mz_3d():size_1(),size_2(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_c2r_mz_3d();
    
    int plan( int req_size_1, int req_size_2, int req_size_3, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, double *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft_c2r( *(plans[pos]), (fftw_complex *)in, out );
        else return 1;
        
        return 0;
    }
};














struct fftw_r2r_mz_1d
{
    
    vector<int > size;
    vector<fftw_r2r_kind > kind;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_r2r_mz_1d():size(),kind(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_r2r_mz_1d();
    
    int plan( int req_size, fftw_r2r_kind req_kind, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, double *in, double*out )
    {
        if ( pos <= plans.size() )
            fftw_execute_r2r( *(plans[pos]), in, out );
        else return 1;
        
        return 0;
    }
};


struct fftw_r2r_mz_2d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<fftw_r2r_kind > kind;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_r2r_mz_2d():size_1(),size_2(),kind(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_r2r_mz_2d();
    
    int plan( int req_size_1, int req_size_2, fftw_r2r_kind req_kind, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, double *in, double*out )
    {
        if ( pos <= plans.size() )
            fftw_execute_r2r( *(plans[pos]), in, out );
        else return 1;
        
        return 0;
    }
};


struct fftw_r2r_mz_3d
{
    
    vector<int > size_1;
    vector<int > size_2;
    vector<int > size_3;
    vector<fftw_r2r_kind > kind;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_r2r_mz_3d():size_1(),size_2(),size_3(),kind(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_r2r_mz_3d();
    
    int plan( int req_size_1, int req_size_2, int req_size_3, fftw_r2r_kind req_kind, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, double *in, double*out )
    {
        if ( pos <= plans.size() )
            fftw_execute_r2r( *(plans[pos]), in, out );
        else return 1;
        
        return 0;
    }
};






/*

struct fftw_mz_2d_vec
{
    vector<int > vec_size;
    vector<int > size_1;
    vector<int > size_2;
    vector<int > in_stride;
    vector<int > out_stride;
    vector<fftw_plan * > plans;
    
    fftw_mz_2d():vec_size(),size_1(),size_2(),in_stride(),out_stride(),plans() {  }
    
    ~fftw_mz_2d();
    
    int plan( int req_vec_size, int req_size_1, int req_size_2, int req_in_stride=1, int req_out_stride=1 );
    
    int execute( int pos, complex<double> *in, complex<double> *out )
    {
        if ( pos <= plans.size() )
            fftw_execute_dft( *(plans[pos]), (fftw_complex *)in, (fftw_complex *)out );
        else return 1;
        
        return 0;
    }
};



*/






#include "wrappers.h"


// house keeping stuff
static bool fftw_initialized = false;

void initialize_fftw()
{
	if (!fftw_initialized)
		assert(fftw_init_threads());
	
	fftw_initialized = true;
}

void release_fftw()
{
	//if (fftw_initialized)	
	//	fftw_cleanup_threads();
    
	fftw_initialized = false;
}


























#endif












