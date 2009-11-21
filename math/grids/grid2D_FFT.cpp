#ifndef GRID2D_FFT_CPP
#define GRID2D_FFT_CPP

#ifndef GRID2D_FFT_H
	#include "grid2D_FFT.h"
#endif


#ifndef GRID_COLUMN_MAJOR
template<class St1,class Xt1, class St2,class Xt2 >
void FFT(
	grid2D<double,St1,Xt1 > & in,
	grid2D<complex<double>,St2,Xt2 > & out )
{
	out.copy(in);
	//assert( (int)&in(1) - (int)&in(0) == 8 && (int)&out(1) - (int)&out(0) == 16 );	
	
	int n1 = in.n1;
	int n2 = in.n2;
	int N = n1*n2;
	
	static vector<int > dim1;
	static vector<int > dim2;
	static vector<fftw_plan * > plans;
	static vector<fftw_complex * > complex_in;
	static vector<fftw_complex * > complex_out;
	
	int I = 0;
	bool found = 0;
	
	for ( I = 0; I<plans.size(); I++)
	{
		if ( dim1[I] == n1 && dim2[I] == n2 )
		{
			found = 1;
			break;
		}
	}
	
	if (!found)
	{
		//cout << "new" << endl;
		I = plans.size();
		dim1.push_back( n1 );
		dim2.push_back( n2 );
		
		complex_in.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		complex_out.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		
		plans.push_back( new fftw_plan );
		*plans[I] = 
			fftw_plan_dft_2d(
				n1,n2,
				complex_in[I],
				complex_out[I],
				FFTW_FORWARD,
				FFTW_MEASURE );
	}
	
	GRID2D_LINLOOP( in )
	{
		complex_in[I][i][0] = in.array[i];
		complex_in[I][i][1] = 0.0;
	}
	
	fftw_execute_dft(
		*plans[I],
		complex_in[I],
		complex_out[I] );
	
	GRID2D_LINLOOP( out )
	{
		real(out.array[i]) = complex_out[I][i][0]/N;
		imag(out.array[i]) = complex_out[I][i][1]/N;
	}
}
#endif
	
#ifndef GRID_COLUMN_MAJOR
template<class St1,class Xt1, class St2,class Xt2 >
void IFFT(
	grid2D<complex<double>,St1,Xt1 > & in,
	grid2D<double,St2,Xt2 > & out )
{
	out.copy(in);
	//assert( (int)&in(1) - (int)&in(0) == 8 && (int)&out(1) - (int)&out(0) == 16 );	
	
	int n1 = in.n1;
	int n2 = in.n2;
	int N = n1*n2;
	
	static vector<int > dim1;
	static vector<int > dim2;
	static vector<fftw_plan * > plans;
	static vector<fftw_complex * > complex_in;
	static vector<fftw_complex * > complex_out;
	
	int I = 0;
	bool found = 0;
	
	for ( I = 0; I<plans.size(); I++)
	{
		if ( dim1[I] == n1 && dim2[I] == n2 )
		{
			found = 1;
			break;
		}
	}
	
	if (!found)
	{
		I = plans.size();
		dim1.push_back( n1 );
		dim2.push_back( n2 );
		
		complex_in.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		complex_out.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		
		plans.push_back( new fftw_plan );
		*plans[I] = 
			fftw_plan_dft_2d(
				n1,n2,
				complex_in[I],
				complex_out[I],
				FFTW_BACKWARD,
				FFTW_MEASURE );
	}
	
	GRID2D_LINLOOP( in )
	{
		complex_in[I][i][0] = real(in.array[i]);
		complex_in[I][i][1] = imag(in.array[i]);
	}
	
	fftw_execute_dft(
		*plans[I],
		complex_in[I],
		complex_out[I] );
	
	GRID2D_LINLOOP( out )
	{
		out.array[i] = complex_out[I][i][0];
	}
}
#endif


#endif
