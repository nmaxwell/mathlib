#ifndef GRID_TRANSFORMS_CPP
#define GRID_TRANSFORMS_CPP








#ifdef GRID1D_DOUBLE_FFT
#ifndef GRID1D_DOUBLE_FFT__
#define GRID1D_DOUBLE_FFT__

#include <complex>
#include <fftw3.h>


template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void FFT(
	grid1D<Xt1,double,St1,Vt1,At1 > & in,
	grid1D<Xt2,complex<double>,St2,Vt2,At2 > & out )
{
	if (out.nx != in.nx || out.array == 0)
	{
		out.nx = in.nx;
		if (!out.array)
			out.array = new At2;
		out.array->resize(out.nx);
	}
		
	out.xa = in.xa;
	out.xb = in.xb;
	
	int nx = in.nx;
	static int nx_ = 0;
	
	if ( (int)&in(1) - (int)&in(0) == 8 && (int)&out(1) - (int)&out(0) == 16 )
	{	//contiguous data 
		
		static fftw_plan plan;
		static fftw_complex * complex_in = 0;
				
		if (nx != nx_)
		{
			//if (!nx_)
				fftw_destroy_plan(plan);
				
			if (complex_in)
				fftw_free(complex_in);
			
			complex_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);		
			assert (complex_in);
			
			plan = fftw_plan_dft_1d(
				nx,
				complex_in,
				(fftw_complex* )&(out(0)),
				FFTW_FORWARD,
				FFTW_MEASURE );
			
			nx_ = nx;	
		}
		
		for (int i = 0; i<nx; i++)
		{
			complex_in[i][0] = (*in.array)[i];
			complex_in[i][1] = 0.0;
		}
		
		fftw_execute_dft(
			plan,
			complex_in,
			(fftw_complex* )&(out(0)) );

		
		out /= (double)nx;
	}
	else
	{	//non-contiguous data
		std::cerr << "non-contiguous data FFT not yet implemented.\n"; 
	}
}
	
template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void IFFT(
	grid1D<Xt1,complex<double>,St1,Vt1,At1 > & in,
	grid1D<Xt2,double,St2,Vt2,At2 > & out )
{
	if (in.nx != out.nx)
	{
		out.nx = in.nx;
		if (!out.array)
			out.array = new At2;
		out.array->resize(in.nx);
	}
	
	out.xa = in.xa;
	out.xb = in.xb;
	
	int nx = in.nx;
	static int nx_ = 0;
	
	if ( (int)&in(1) - (int)&in(0) == 16 && (int)&out(1) - (int)&out(0) == 8 )
	{	//contiguous data 
		
		static fftw_plan plan;
		static fftw_complex * complex_out = 0;
		
		if (nx != nx_)
		{
			//if (!nx_)
				fftw_destroy_plan(plan);
			
			if (complex_out)
				fftw_free(complex_out);
			
			complex_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);		
			
			plan = fftw_plan_dft_1d(
				nx,
				(fftw_complex* )&(in(0)),
				complex_out,
				FFTW_FORWARD,// i give up... wont work on backwards
				FFTW_MEASURE );
			
			nx_ = nx;
		}
		
		fftw_execute_dft(
			plan,
			(fftw_complex* )&(in(0)),
			complex_out );
		
		for (int i = 0; i<nx; i++)
			(*out.array)[i] = complex_out[i][0];
		
	}
	else
	{	//non-contiguous data
		std::cerr << "non-contiguous data IFFT not yet implemented.\n"; 
	}
}



#endif
#endif




#ifdef GRID1D_DOUBLE_DCT
#ifndef GRID1D_DOUBLE_DCT__
#define GRID1D_DOUBLE_DCT__

#include <complex>
#include <fftw3.h>


template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void DCT2(
	grid1D<Xt1,double,St1,Vt1,At1 > & in,
	grid1D<Xt2,double,St2,Vt2,At2 > & out )
{
	
}
	
template<class Xt1,class St1,class Vt1,class At1, class Xt2,class St2,class Vt2,class At2 >
void IDCT4(
	grid1D<Xt1,double,St1,Vt1,At1 > & in,
	grid1D<Xt2,double,St2,Vt2,At2 > & out )
{
	
}

	



#endif
#endif



















#endif
