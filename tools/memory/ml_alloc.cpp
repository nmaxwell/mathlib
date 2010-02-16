#ifndef ML_ALLOC_CPP
#define ML_ALLOC_CPP


#include <map>

map< uint64, uint64 > ml_memory_usage;
map< uint64, uint64 >::iterator ml_memory_iterator;



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

template<class T >
T ** ml_alloc(  int m, int n )
{
    T ** p = ml_alloc< T* > (m);
    
    for (int k=0; k<m; k++)
        p[k] = ml_alloc<T > (n);
    
    return p;
}

template<class T >
void ml_free( T** &p, int m )
{
	if (p)
    {
        for (int k=0; k<m; k++)
            ml_free ( p[k] );
    }
    delete [] p;
	p=0;
}


#ifndef ML_NOT_USE_FFTW_MALLOC

// double, using fftw_alloc

template<>
double * ml_alloc( int n )
{
    double * p = (double *)fftw_malloc( sizeof(double)*n );
    
    #ifdef ML_TRACK_MEMORY
    
        ml_memory_usage[ (int64)p ] = sizeof(double)*n;
        
        int64 sum = 0;
        
        for ( ml_memory_iterator=ml_memory_usage.begin() ; ml_memory_iterator != ml_memory_usage.end(); ml_memory_iterator++ )
            sum += (*ml_memory_iterator).second;
        
        std::cout << "ml_alloc called; memory usage (MiB):\t" << sum/(1024*1024) << std::endl;
        
    #endif
    
    return p;
}

template< >
void ml_free( double* &p )
{
    #ifdef ML_TRACK_MEMORY
    
        ml_memory_usage.erase( (int64)p );
        
        int64 sum = 0;
        
        for ( ml_memory_iterator=ml_memory_usage.begin() ; ml_memory_iterator != ml_memory_usage.end(); ml_memory_iterator++ )
            sum += (*ml_memory_iterator).second;
        
        std::cout << "ml_free called; memory usage (MiB):\t" << sum/(1024*1024) << std::endl;
    
    #endif
    
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




#endif









