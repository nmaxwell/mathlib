#ifndef ML_RAND_WALK_CPP
#define ML_RAND_WALK_CPP


#include "random_walk.h"

inline int ml_bitsum( uint32 const & x )
{
    /*
     * sum all the set bits in x
     * 
     * see http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
     * 
     */
    
    /*int sum = 0;
    for (int k=0; k<32; k++)
        sum += (x & (1 << k)) >> k;
    return sum; */
    
    uint32 sum;

    sum = x - ((x >> 1) & 0x55555555);
    sum = ((sum >> 2) & 0x33333333) + (sum & 0x33333333);
    sum = ((sum >> 4) + sum) & 0x0F0F0F0F;
    sum = ((sum >> 8) + sum) & 0x00FF00FF;
    sum = ((sum >> 16) + sum) & 0x0000FFFF;
    
    return sum;
}

/*
void binom_to_random_walk( uint32 * pool, int n, int32 * walk )
{
    walk[0] = 0;
    
    for ( int k=1; k<n; k++ )
        walk[k] = walk[k-1] + 2*ml_bitsum( pool[k-1] ) -32;
    for ( int k=1; k<n; k++ )
        walk[k] /= 2;
}*/

void binom_to_random_walk( uint32 * pool, int n, int32 * walk )
{
    walk[0] = 0;
    
    int i=0, j=0, k=0;
    int R = pool[0];
    
    while (k<n)
    {
        //cout << (R&1);
        if (R & 1) walk[++k] = walk[k-1]+1;
        else walk[++k] = walk[k-1]-1;
        
        if ( i++ < 32 )
            R >>= 1;
        else { R = pool[++j]; i=0; }
    }
    
}



void gen_BM( double h, int n_steps, double *& W, int32 * walk, int32 * pool )
{
    static ml_random rng;
    
    for (int k=0; k<div_up(n_steps,32); k++)
        pool[k] = rng.gen_int();
    
    binom_to_random_walk ( pool, n_steps, walk );
    
    for ( int k=0; k<n_steps; k++ )
        W[k] = h*walk[k];
}

void gen_BM( double h, int n_steps, float *& W, int32 * walk, int32 * pool )
{
    static ml_random rng;
    
    for (int k=0; k<div_up(n_steps,32); k++)
        pool[k] = rng.gen_int();
    
    binom_to_random_walk ( pool, n_steps, walk );
    
    for ( int k=0; k<n_steps; k++ )
        W[k] = h*walk[k];
}





void gen_BM( double h, int n_steps, double *& W )
{
    if ( W == 0 ) W = ml_alloc<double > ( n_steps );
    
    int32 *walk = ml_alloc<int32 > (n_steps);
    int32 *pool = ml_alloc<int32 > ( div_up(n_steps,32) );
    
    gen_BM( h, n_steps, W, walk, pool);
    
    ml_free (walk);
    ml_free (pool);
}

void gen_BM( double h, int n_steps, double **& W, int n_runs )
{
    if ( W == 0 ) W = ml_alloc<double > ( n_steps, n_runs );
    
    int32 *walk = ml_alloc<int32 > (n_steps);
    int32 *pool = ml_alloc<int32 > ( div_up(n_steps,32) );
    
    for (int j=0; j<n_runs; j++)
        gen_BM( h, n_steps, W[j], walk, pool);
    
    ml_free (walk);
    ml_free (pool);
}

void gen_BM( double h, int n_steps, float *& W )
{
    if ( W == 0 ) W = ml_alloc<float > ( n_steps );
    
    int32 *walk = ml_alloc<int32 > (n_steps);
    int32 *pool = ml_alloc<int32 > ( div_up(n_steps,32) );
    
    gen_BM( h, n_steps, W, walk, pool);
    
    ml_free (walk);
    ml_free (pool);
}

void gen_BM( double h, int n_steps, float **& W, int n_runs )
{
    if ( W == 0 ) W = ml_alloc<float > ( n_steps, n_runs );
    
    int32 *walk = ml_alloc<int32 > (n_steps);
    int32 *pool = ml_alloc<int32 > ( div_up(n_steps,32) );
    
    for (int j=0; j<n_runs; j++)
        gen_BM( h, n_steps, W[j], walk, pool);
    
    ml_free (walk);
    ml_free (pool);
}


#endif
