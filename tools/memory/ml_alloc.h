#ifndef ML_ALLOC_H
#define ML_ALLOC_H

// uses fftw's memory allocator
#include <fftw3.h>

template<class T >
T * ml_alloc( int n );


template<class T >
void ml_free( T* & );



template<>
double * ml_alloc( int n );

template< >
void ml_free( double* &p );


#endif



