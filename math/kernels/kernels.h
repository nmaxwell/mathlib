#ifndef KERNELS_H
#define KERNELS_H


#include <mathlib/math/transforms/FFT.h>

// design notes : convolution_plan
/*
	
	convolution_plan is a plan to apply a kernel on a strip 
	of N data points, of type X, with stride x_tride, x_ is 
	address of frist point.
	
	The kernel may change from point to point, mostly to take
	care of boundary conditions. K will be the type of data
	of the kernel, and K ** ker_ will be an array, of K 
	pointers, the kernels to be applied. T will be the type
	used in the sum, so extended precision might be used.
	
	y_ will be a pointer to the first element in the return
	array, which will have stride y_stride.
	
	The arrays ker_, n_l, n_r, will be contiguous, and and be of length N.
	
	nr_i will be the points to the right of the current x_i.
	nl_i will be the points to the left of the current x_i.
	
	So the ith kernel must be of length ( nr_i + nl_i + 1 ).
	
	kernels[i][0] will be the element at the center of the kernel.
	So, the loop will look like		
		
		for (int i = 0; i<N; i++)
		for (int j = nl_[i]; j <= nr_[j]; j++)
			sum += ker_[i][j]*x_[i+j];
			
	If j goes of the grid, it will wrap around.
		
	
	execute() will run the convolution.
	
	this only adds the convolution to y. so, y += sum*scale, so remember to set y=0 before hand.
	
	dont call copy constructor...
*/

template<class X=double, class K=X, class T=X >
class convolution_plan
{
public:
	
	int N;
	
	K ** ker_;
	
	int * nl_;
	int * nr_;
	
public:
	
	convolution_plan(int N):N(N)
	{		
		ker_ = new K* [N];
		nl_ = new int [N];
		nr_ = new int [N];
	}
	
	~convolution_plan() {
		if (ker_) delete [] ker_;
		if (nl_) delete [] nl_;
		if (nr_) delete [] nr_; }
		
	void resize(int N_)
	{
		N = N_;
		
		if (ker_) delete [] ker_;
		if (nl_) delete [] nl_;
		if (nr_) delete [] nr_;
		
		ker_ = new K* [N];
		nl_ = new int [N];
		nr_ = new int [N];
	}
	
	void set_ker(int i, K* ker) { ker_[i] = ker; }
	void set_nl(int i, int nl) { nl_[i] = nl; }
	void set_nr(int i, int nr) { nr_[i] = nr; }
	
public:

	int x_stride;
	X * x_;
	
	int y_stride;
	X * y_;
	
	void execute( X * x_, int x_stride, X * y_, int y_stride, T scale = 1.0 )
	{
	/*	T sum;
		int i,j;
		
		for ( i = 0; i<N; i++)
		{
			sum = 0.0;
						
			if ( i-nl_[i] >= 0)
			{
				if ( nr_[i]+i < N )
				{
					// no wrap on left and no wrap on right
					for ( j = -nl_[i]; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
				}
				else
				{
					// no wrap on left, but is wrap on right
				//	cout << "no wrap on left, but is wrap on right " << endl;
					for ( j = -nl_[i]; i+j < N; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
					for ( ; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j-N)*x_stride];
				}	
			}
			else
			{
				if ( nr_[i]+i < N )
				{
					// is wrap on left, but no wrap on right
				//	cout << "is wrap on left, but no wrap on right " << endl;
					for ( j = -nl_[i]; j+i < 0 ; j++)
						sum += ker_[i][j]*x_[(i+j+N)*x_stride];
					for ( ; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
				}
				else
				{
					// is wrap on left, and is wrap on right
				//	cout << "is wrap on left, and is wrap on right " << endl;
					for ( j = -nl_[i]; j+i < 0 ; j++)
						sum += ker_[i][j]*x_[(i+j+N)*x_stride];
					for ( ;  i+j < N; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
					for ( ; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j-N)*x_stride];
				}
			}
			
			y_[i*y_stride] += sum*scale;
		}*/
	}
	
};


// design notes : hybrid_convolution_plan
/*
    Does much the same as convolution_plan, except that it applies a kernel , cen_ker_, 
    to the data first, using the DFT convolution theorem, then overites that result at the ends with
    the kernels, ** bd_ker_, applied the same way as in convolution_plan.
    
    hybrid_convolution_plan replaces y_ with the result (i.e. not =, not += )
    
    dont call copy constructor...
*/

template<class X=double, class K=X, class T=X >
class hybrid_convolution_plan
{
public:
	
	int N;
	
	complex<K > * cen_ker_fft;
	
	int l_N;	
	K ** l_ker_;	
	int * l_nl_;
	int * l_nr_;	
	
	int r_N;
	K ** r_ker_;	
	int * r_nl_;
	int * r_nr_;
	
	complex<X > * x_fft_;
	
public:
	
	hybrid_convolution_plan():N(0),l_N(0),r_N(0)
	{
	    cen_ker_fft = 0;		
		l_ker_ = 0;
		l_nl_ = 0;
		l_nr_ = 0;		
		r_ker_ = 0;
		r_nl_ = 0;
		r_nr_ = 0;		
		x_fft_ = 0;
	}	
	
	hybrid_convolution_plan(int N, int l_N, int r_N):N(N),l_N(l_N),r_N(r_N)
	{		
		cen_ker_fft = 0;
		
		l_ker_ = new K* [l_N];
		l_nl_ = new int [l_N];
		l_nr_ = new int [l_N];
		
		r_ker_ = new K* [r_N];
		r_nl_ = new int [r_N];
		r_nr_ = new int [r_N];
		
		x_fft_ = new complex<X > [N/2+1];
	}
	
	void resize(int _N, int _l_N, int _r_N)
	{
	    N = _N;
	    l_N = _l_N;
        r_N = _r_N;
        
		if (cen_ker_fft) delete [] cen_ker_fft; cen_ker_fft = 0;
		if (l_ker_) delete [] l_ker_; l_ker_ = 0;
		if (l_nl_) delete [] l_nl_; l_nl_ = 0;
		if (l_nr_) delete [] l_nr_; l_nr_ = 0;
		if (r_ker_) delete [] r_ker_; r_ker_ = 0;
		if (r_nl_) delete [] r_nl_; r_nl_ = 0;
		if (r_nr_) delete [] r_nr_; r_nr_ = 0;
		if (x_fft_) delete [] x_fft_; x_fft_ = 0;
	    
		cen_ker_fft = 0;
		
		l_ker_ = new K* [l_N];
		l_nl_ = new int [l_N];
		l_nr_ = new int [l_N];
		
		r_ker_ = new K* [r_N];
		r_nl_ = new int [r_N];
		r_nr_ = new int [r_N];
		
		x_fft_ = new complex<X > [N/2+1];
	}
	
	~hybrid_convolution_plan() {
		if (cen_ker_fft) delete [] cen_ker_fft; cen_ker_fft = 0;
		if (l_ker_) delete [] l_ker_; l_ker_ = 0;
		if (l_nl_) delete [] l_nl_; l_nl_ = 0;
		if (l_nr_) delete [] l_nr_; l_nr_ = 0;
		if (r_ker_) delete [] r_ker_; r_ker_ = 0;
		if (r_nl_) delete [] r_nl_; r_nl_ = 0;
		if (r_nr_) delete [] r_nr_; r_nr_ = 0;
		if (x_fft_) delete [] x_fft_; x_fft_ = 0; }
		
public:

	void set_cen_ker( K* ker,int nl, int nr ) {

	        K * ker_ext_ = new K [N];
	        
	        for (int i=0; i<N; i++)
	            ker_ext_[i] = 0;
	        
	        for (int i=0; i<=nr; i++)
	            ker_ext_[i] = ker[i];
	            
	        for (int i=1; i<=nl; i++)
	            ker_ext_[N-i] = ker[-i];
	        
	        FFT(ker_ext_,cen_ker_fft,N,1,1); 
	        
	        delete [] ker_ext_;
	        
	    }
	
	void set_l_ker(int i, K* ker) { l_ker_[i] = ker; }
	void set_l_nl(int i, int nl) { l_nl_[i] = nl; }
	void set_l_nr(int i, int nr) { l_nr_[i] = nr; }
	
	void set_r_ker(int i, K* ker) { r_ker_[i] = ker; }
	void set_r_nl(int i, int nl) { r_nl_[i] = nl; }
	void set_r_nr(int i, int nr) { r_nr_[i] = nr; }
	
public:
	
	void execute( X * x_, int x_stride, X * y_, int y_stride, T scale = 1.0 )
	{
	    FFT(x_,x_fft_,N, x_stride,1 );
	    
	    int N2 = N/2+1;
	    for (int i=0; i<N2; i++)
	        x_fft_[i] *= cen_ker_fft[i]/(X)N;
	    	    
	    IFFT(x_fft_,y_,N, 1,y_stride );
	    
	    if (scale != 1.0)
	    for ( int i = l_N; i<N-r_N; i++)
	        y_[i*y_stride] *= scale;
	    	    
    	T sum;
	    int i,j,k;
	    
	    for ( i = 0; i<l_N; i++)
		{
			sum = 0.0;
			
			for ( j = -l_nl_[i]; j <= l_nr_[i]; j++)
			{
			    if ( (i+j) >= N ) k = -N;
			    else if ( (i+j) < 0 ) k = +N;
			    else k = 0;
			    
				sum += l_ker_[i][j]*x_[(i+j+k)*x_stride];
			}
			
			y_[i*y_stride] = sum*scale;
	    }
	    
	    for ( i = 0; i<r_N; i++)
		{
			sum = 0.0;
			
			for ( j = -r_nl_[i]; j <= r_nr_[i]; j++)
			{
			    if ( (N-1-i+j) >= N ) k = -N;
			    else if ( (N-1-i+j) < 0 ) k = +N;
			    else k = 0;
			    
				sum += r_ker_[i][j]*x_[(N-1-i+j+k)*x_stride];
			}
			
			y_[(N-1-i)*y_stride] = sum*scale;
	    }
	}
	
};

template<class X, class K >
void nonperiodic_convolve(
    X * in,
    X * out,
    int n,
    int in_stride,
    int out_stride,
    K * kernel,
    int nl,
    int nr,
    int ker_stride,
    X scale,
    int bdl,
    int bdr );








#endif
