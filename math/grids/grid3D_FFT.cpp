#ifndef GRID3D_FFT_CPP
#define GRID3D_FFT_CPP

#ifndef GRID3D_FFT_H
	#include "grid3D_FFT.h"
#endif


#ifndef GRID_COLUMN_MAJOR
template<class St1,class Xt1, class St2,class Xt2 >
void FFT(
	grid3D<double,St1,Xt1 > & in,
	grid3D<complex<double>,St2,Xt2 > & out )
{
    //double t0 = getRealTime();
    
	out.copy(in);
	//assert( (int)&in(1) - (int)&in(0) == 8 && (int)&out(1) - (int)&out(0) == 16 );	
	
	int n1 = in.n1;
	int n2 = in.n2;
	int n3 = in.n3;
	int N = n1*n2*n3;
	
	static vector<int > dim1;
	static vector<int > dim2;
	static vector<int > dim3;
	static vector<fftw_plan * > plans;
	static vector<fftw_complex * > complex_in;
	static vector<fftw_complex * > complex_out;
	
	static int I = 0;
	bool found = 0;
	
    if (I)
    if ( dim1[I] == n1 && dim2[I] == n2 && dim3[I] == n3 )
        found = true;
	
	if (!found)
	for ( I = 0; I<plans.size(); I++)
	{
		if ( dim1[I] == n1 && dim2[I] == n2 && dim3[I] == n3 )
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
		//cout << "new" << endl;
		I = plans.size();
		dim1.push_back( n1 );
		dim2.push_back( n2 );
		dim3.push_back( n3 );
		
		complex_in.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		complex_out.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		
		plans.push_back( new fftw_plan );
		*plans[I] = 
			fftw_plan_dft_3d(
				n1,n2,n3,
				complex_in[I],
				complex_out[I],
				FFTW_FORWARD,
				FFTW_MEASURE );
	}
	
	//double t1 = getRealTime();
	
	GRID3D_LINLOOP( in )
	{
		complex_in[I][i][0] = in.array[i];
		complex_in[I][i][1] = 0.0;
	}
	
	//double t2 = getRealTime();
	
	fftw_execute_dft(
		*plans[I],
		complex_in[I],
		complex_out[I] );
		
	//double t3 = getRealTime();
	
	GRID3D_LINLOOP( out )
	{
		real(out.array[i]) = complex_out[I][i][0]/N;
		imag(out.array[i]) = complex_out[I][i][1]/N;
	}
	
	//double t4 = getRealTime();
	
	//cout << t4-t0 << "\t" << (t1-t0)/(t4-t0) << "\t" << (t2-t1)/(t4-t0) << "\t" << (t3-t2)/(t4-t0) << "\t" << (t4-t3)/(t4-t0) << "\t" << endl2;
}
#endif
	
#ifndef GRID_COLUMN_MAJOR
template<class St1,class Xt1, class St2,class Xt2 >
void IFFT(
	grid3D<complex<double>,St1,Xt1 > & in,
	grid3D<double,St2,Xt2 > & out )
{
	out.copy(in);
	//assert( (int)&in(1) - (int)&in(0) == 8 && (int)&out(1) - (int)&out(0) == 16 );	
	
	int n1 = in.n1;
	int n2 = in.n2;
	int n3 = in.n3;
	int N = n1*n2*n3;
	
	static vector<int > dim1;
	static vector<int > dim2;
	static vector<int > dim3;
	static vector<fftw_plan * > plans;
	static vector<fftw_complex * > complex_in;
	static vector<fftw_complex * > complex_out;
	
	int I = 0;
	bool found = 0;
	
	for ( I = 0; I<plans.size(); I++)
	{
		if ( dim1[I] == n1 && dim2[I] == n2 && dim3[I] == n3 )
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
		dim3.push_back( n3 );
		
		complex_in.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		complex_out.push_back( (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N) );
		
		plans.push_back( new fftw_plan );
		*plans[I] = 
			fftw_plan_dft_3d(
				n1,n2,n3,
				complex_in[I],
				complex_out[I],
				FFTW_BACKWARD,
				FFTW_MEASURE );
	}
	
	GRID3D_LINLOOP( in )
	{
		complex_in[I][i][0] = real(in.array[i]);
		complex_in[I][i][1] = imag(in.array[i]);
	}
	
	fftw_execute_dft(
		*plans[I],
		complex_in[I],
		complex_out[I] );
	
	GRID3D_LINLOOP( out )
	{
		out.array[i] = complex_out[I][i][0];
	}
}
#endif


#endif
