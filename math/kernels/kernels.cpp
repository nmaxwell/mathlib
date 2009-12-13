#ifndef KERNELS_CPP
#define KERNELS_CPP

#include "kernels.h"



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
    int bdr )
{
    
 //   cout << "here1\n";
    
    
    /*
    
    double * ker_ext = new double [n ];	    
    complex<double > * ker_fft = new complex<double > [n/2+2];
    complex<double > * Q = new complex<double > [n/2+2] ;
	
	
    for (int i=0; i<n; i++)
        ker_ext[i] = 0.0;
    for (int i=0; i<=nr; i++)
        ker_ext[i] = kernel[i];
    for (int i=1; i<=nl; i++)
        ker_ext[n-i] = kernel[-i];	    
    
    FFT(ker_ext,ker_fft,n,1,1 );    
    FFT(in,Q,n,1,1);
    
    for (int i=0; i<=(n/2+1); i++) 
        Q[i] *= ker_fft[i]*scale/((double)n);
    
    IFFT(Q,out,n,1,1);  
    
    
    delete [] ker_ext;	    
    delete []  ker_fft;
    delete []  Q;
    
    
    */
    
    
    
    /*  Convolves the data 'in' with the kernel 'kernel', and outputs into 'out'. 
        kernel[0] is the center point of the kernel, it has 'nl' points to the left,
        and 'nr' to the right, kernel is accessed as kernel[i] for -L<=i<=R.
        'n' is the length of 'in,out'.
        
        First a periodic convulution is applied, which, assuming kernel is shorter than 'in,out',
        gives the right values in the center region.
        
        At the edges, defined by 'bdl' points on the left, and 'bdr' points on the right, 
        the convolution is applied from the definition, and the points that get sampled off the 
        grid are assumed to be zero.
        
        The formula used is out[i] = sum[j=-R to L]{ x[i-j]*ker[j] }
    */
    
 /*   if (out == 0) out = new X [ n*out_stride ];
	
	static vector<int > ns;
	static vector<K * > ker_ptrs;
	static vector<int > bdls;
	static vector<int > bdrs;
	
	static vector<complex<K > * > ker_ffts;
	static vector<K * > ker_ext;
    static vector<complex<X > * > Q;
	
	int I = 0;
	bool found = 0;
	
	if (I)
	if ( ns[I] == n && ker_ptrs[I] == kernel && bdls[I] == bdl && bdrs[I] == bdr )
	    found = true;
	
	if (!found)
	for ( I = 0; I<ns.size(); I++)
	{
		if ( ns[I] == n && ker_ptrs[I] == kernel && bdls[I] == bdl && bdrs[I] == bdr )
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
	    ns.push_back(n);
	    bdls.push_back(bdl);
        bdrs.push_back(bdr);
	    ker_ptrs.push_back(kernel);
	    
	    ker_ext.push_back ( new K [n ] );	    
	    ker_ffts.push_back (new complex<K > [n/2+2]);
	    Q.push_back( new complex<X > [n/2+2] );
	    
	    for (int i=0; i<n; i++)
	        ker_ext[I][i] = 0.0;
	    for (int i=0; i<=nr; i++)
	        ker_ext[I][i] = kernel[i*ker_stride];
	    for (int i=1; i<=nl; i++)
	        ker_ext[I][n-i] = kernel[-i*ker_stride];	    
	    
	    FFT(ker_ext[I],ker_ffts[I],n,ker_stride,1 );
	    
	}*/
    
  /*  FFT(in,Q[I],n,in_stride,1);
    for (int i=0; i<=(n/2+1); i++) 
        Q[I][i] *= ker_ffts[I][i];
    IFFT(Q[I],out,n,1,out_stride); */
    
 /*   for (int i=0; i<=bdl; i++)
        out[i*out_stride] = 0;
    for (int i=0; i<=bdr; i++)
        out[(n-1-i)*out_stride] = 0;
    
    
    for (int i=0; i<=bdl; i++)
    {
     //   for (int j=-R; j<=L; j++)
     //       out[i] += ker[j]*in[i-j];
     
     //   for (int j=-bdl; j<=i; j++)
     //       out[i] += ker[j]*in[i-j];     
        
        for (int j=-bdl; j<=-1; j++)
            out[i*out_stride] += ker_ext[I][(n+j)*ker_stride]*in[(i-j)*in_stride];
        
        for (int j=0; j<=i; j++)
            out[i*out_stride] += ker_ext[I][j*ker_stride]*in[(i-j)*in_stride];        
    }
    
    for (int i=0; i<=bdr; i++)
    {
     //   for (int j=-R; j<=L; j++)
     //       out[i] += ker[j]*in[i-j];
     
     //   for (int j=-i; j<=bdr; j++)
     //       out[n-1-i] += ker[j]*in[n-1-i -j]; 
     
        for (int j=-i; j<=-1; j++)
            out[(n-1-i)*out_stride] += ker_ext[(n+j)*ker_stride]*in[(n-1-i -j)*in_stride];
         
        for (int j=0; j<=bdr; j++)
            out[(n-1-i)*out_stride] += ker_ext[j*ker_stride]*in[(n-1-i -j)*in_stride];       
    }
    
     */   
};





#endif


