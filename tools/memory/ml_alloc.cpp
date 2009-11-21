#ifndef ML_ALLOC_CPP
#define ML_ALLOC_CPP



// general (default to new)

template<class T >
T * ml_alloc( int n )
{
	return new T [ n ];
}

template<class T >
void ml_free( T* &p )
{
	if (p) delete [] p;
	p=0;
}

// double, using fftw_alloc

template<>
double * ml_alloc( int n )
{
	return (double *)fftw_malloc( sizeof(double)*n );
}

template< >
void ml_free( double* &p )
{
	if (p) fftw_free(p);
	p=0;
}

// complex<double>, using fftw_alloc

template<>
complex<double> * ml_alloc( int n )
{
	return (complex<double> *)fftw_malloc( sizeof(fftw_complex)*n );
}

template< >
void ml_free( complex<double>* &p )
{
	if (p) fftw_free(p);
	p=0;
}






#endif









