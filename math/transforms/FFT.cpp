#ifndef FFT_CPP
#define FFT_CPP

#include "FFT.h"

void FFT(
	double * in,
    complex<double > * & out,
    int n,
    int in_stride,
    int out_stride )
{
	if (out == 0) out = ml_alloc<complex<double> > ( (int)(n/2+1)*out_stride );
	
	static vector<int > dims;
	static vector<fftw_plan * > plans;
	static vector<int > in_strides;
	static vector<int > out_strides;
	
	static int I = 0;
	bool found = 0;
	
	if (I)
        if (dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride)
            found = true;
	
	if (!found)
    {
        for ( I = 0; I<plans.size(); I++)
        {
            if ( dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride )
            {
                found = true;
                break;
            }
        }
    }
	
	if (!found)
	{
		I = plans.size();
		dims.push_back( n );
		in_strides.push_back( in_stride );
		out_strides.push_back( out_stride );
		
		plans.push_back( new fftw_plan );
		
		double * t_in = ml_alloc<double> ( dims[I]*in_stride );
		complex<double > * t_out = ml_alloc<complex<double> > ( (int)(dims[I]/2+1)*out_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        
        *plans[I] = 
            fftw_plan_many_dft_r2c(
                1,
                &(dims[I]),
                1,
                t_in,
                0,
                in_strides[I],
                0,
                (fftw_complex *)t_out,
                0,
                out_strides[I],
                0,
                FFTW_ESTIMATE_MODE | FFTW_PRESERVE_INPUT );
                //FFTW_ESTIMATE | FFTW_PRESERVE_INPUT ); // FFTW_ESTIMATE
                //FFTW_PATIENT | FFTW_PRESERVE_INPUT ); // FFTW_ESTIMATE
          
          ml_free(t_in);
          ml_free(t_out);
	}
		
	fftw_execute_dft_r2c(
		*plans[I],
		in,
		(fftw_complex *)out );

}


void IFFT(
    complex<double > * in,
	double * & out,    
    int n,
    int in_stride,
    int out_stride )
{
	if (out == 0) out = ml_alloc<double> (n*out_stride);
	
	static vector<int > dims;
	static vector<fftw_plan * > plans;
	static vector<int > in_strides;
	static vector<int > out_strides;
	
	int I = 0;
	bool found = 0;
	
	if (I)
	if (dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride)
	    found = true;
	
	if (!found)
	for ( I = 0; I<plans.size(); I++)
	{
		if ( dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride )
		{
			found = true;
			break;
		}
	}	
	
	if (!found)
	{
		I = plans.size();
		dims.push_back( n );
		in_strides.push_back( in_stride );
        out_strides.push_back( out_stride );
		
		plans.push_back( new fftw_plan );		

		double * t_out = ml_alloc<double> ( dims[I]*out_stride );
		complex<double > * t_in = ml_alloc<complex<double> > ( (int)(dims[I]/2+1)*in_stride );
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        
        *plans[I] = 
            fftw_plan_many_dft_c2r(
                1,
                &(dims[I]),
                1,
                (fftw_complex *)t_in,
                0,
                in_strides[I],
                0,
                t_out,
                0,
                out_strides[I],
                0,
                FFTW_ESTIMATE_MODE | FFTW_PRESERVE_INPUT );
                //FFTW_PATIENT ); // FFTW_ESTIMATE
                //FFTW_ESTIMATE ); // FFTW_ESTIMATE
                
          
          ml_free(t_in);
          ml_free(t_out);
	}
		
	fftw_execute_dft_c2r(
		*plans[I],
		(fftw_complex *)in,
		out );

}





void FFT(
	float * in,
    complex<float > * & out,
    int n,
    int in_stride,
    int out_stride )
{
	//if (out == 0) out = new complex<float > [(int)(n/2+2)*out_stride];
	if (out == 0) out = (complex<float > *    )fftwf_malloc(sizeof(fftwf_complex) * (int)(n/2+2));
	
	static vector<int > dims;
	static vector<fftwf_plan * > plans;
	
	static vector<int > in_strides;
	static vector<int > out_strides;
	
	int I = 0;
	bool found = 0;
	
	//if (I)
	//if (dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride)
	//    found = true;
	
	//if (!found)
	for ( I = 0; I<plans.size(); I++)
	{
		if ( dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride )
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
		//cout << "new" << endl;
		I = plans.size();
		dims.push_back( n );
		in_strides.push_back( in_stride );
        out_strides.push_back( out_stride );
		
		plans.push_back( new fftwf_plan );
		
		float * t_in = new float [dims[I]*in_stride ];
		complex<float > * t_out = new complex<float > [ (int)(dims[I]/2+1)*out_stride ];
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        
        *plans[I] = 
            fftwf_plan_many_dft_r2c(
                1,
                &(dims[I]),
                1,
                t_in,
                0,
                in_strides[I],
                0,
                (fftwf_complex *)t_out,
                0,
                out_strides[I],
                0,
                FFTW_ESTIMATE | FFTW_PRESERVE_INPUT ); // FFTW_ESTIMATE
                //FFTW_PATIENT | FFTW_PRESERVE_INPUT ); // FFTW_ESTIMATE
          
          delete [] t_in;
          delete [] t_out;
	}
    
	fftwf_execute_dft_r2c(
		*(plans[I]),
		in,
		(fftwf_complex *)out );
}


void IFFT(
    complex<float > * in,
	float * & out,    
    int n,
    int in_stride,
    int out_stride )
{
	if (out == 0) out = new float [n*out_stride];
	
	static vector<int > dims;
	static vector<fftwf_plan * > plans;
	static vector<int > in_strides;
	static vector<int > out_strides;
	
	int I = 0;
	bool found = 0;
	
	if (I)
	if (dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride)
	    found = true;
	
	if (!found)
	for ( I = 0; I<plans.size(); I++)
	{
		if ( dims[I] == n && in_strides[I] == in_stride && out_strides[I] == out_stride )
		{
			found = true;
			break;
		}
	}	
	
	if (!found)
	{
		//cout << "new" << endl;
		I = plans.size();
		dims.push_back( n );
		in_strides.push_back( in_stride );
        out_strides.push_back( out_stride );
		
		plans.push_back( new fftwf_plan );
		
		float * t_out = new float [dims[I]*out_stride ];
		complex<float > * t_in = new complex<float > [ (int)(dims[I]/2+1)*in_stride ];
        
        fftw_plan_with_nthreads( N_FFT_THREADS );
        
        *plans[I] = 
            fftwf_plan_many_dft_c2r(
                1,
                &(dims[I]),
                1,
                (fftwf_complex *)t_in,
                0,
                in_strides[I],
                0,
                t_out,
                0,
                out_strides[I],
                0,
                //FFTW_PATIENT ); // FFTW_ESTIMATE
                FFTW_ESTIMATE ); // FFTW_ESTIMATE
                
          
          delete [] t_in;
          delete [] t_out;
	}
		
	fftwf_execute_dft_c2r(
		*plans[I],
		(fftwf_complex *)in,
		out );

}    



#endif
