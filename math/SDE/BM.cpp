#ifndef ML_BM_CPP
#define ML_BM_CPP

#include "BM.h"

#ifndef BM_mode_1
    #ifndef BM_mode_3
        #define BM_mode_2
    #endif
#endif

inline int ml_bitsum( uint32 const & x )
{
    /*
     * sum all the set bits in x
     * 
     * see http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
     * 
     */
    
    int sum = 0;
    for (int k=0; k<32; k++)
        sum += (x & (1 << k)) >> k;
    return sum;
    
    /*
    uint32 sum;
    sum = x - ((x >> 1) & 0x55555555);
    sum = ((sum >> 2) & 0x33333333) + (sum & 0x33333333);
    sum = ((sum >> 4) + sum) & 0x0F0F0F0F;
    sum = ((sum >> 8) + sum) & 0x00FF00FF;
    sum = ((sum >> 16) + sum) & 0x0000FFFF;
    */
    
    return sum;
}


//----------------------------



void gen_BM( double h, int n_steps, double *& W, int mode )
{
    if ( W == 0 ) W = ml_alloc<double > ( n_steps );
    static ml_random rng;
    
    if ( mode == BM_mode_1 )
    {
        int n_32 = div_up(n_steps,32);
        
        int32 *walk = ml_alloc<int32 > ( n_steps );
        int32 *pool = ml_alloc<int32 > ( n_32 );
        
            for (int k=0; k<n_32; k++)
                pool[k] = rng.gen_int();
            
            walk[0] = 0;
                
            int i=0, j=0;
                
            for ( int k=1; k<n_steps; k++ )
            {
                if ((pool[j] >> i) & 1) walk[k] = walk[k-1]+1;
                else walk[k] = walk[k-1]-1;
                    
                i++;
                    
                if ( i >= 32 )
                {
                    i=0;
                    j++;
                }
            }
            
            for ( int k=0; k<n_steps; k++ )
                W[k] = sqrt(h)*walk[k];
        
        ml_free (walk);
        ml_free (pool);
    }
    
    if ( mode == BM_mode_2 )
    {
        int32 *walk = ml_alloc<int32 > ( n_steps );
        int32 *pool = ml_alloc<int32 > ( n_steps );
        
            for (int k=0; k<n_steps; k++)
                pool[k] = rng.gen_int();
            
            walk[0] = 0;
            
            for ( int k=1; k<n_steps; k++ )
                walk[k] = walk[k-1] + 2*ml_bitsum( pool[k-1] ) -32;
            
            for ( int k=0; k<n_steps; k++ )
                W[k] = sqrt(h/32)*walk[k];
        
        ml_free (walk);
        ml_free (pool);
    }
    
    if ( mode == BM_mode_3 )
    {
        double * dW = ml_alloc<double > ( n_steps );
        
            rng.std_normal_rv( dW, n_steps);
            
            W[0] = 0;
            
            for ( int j=1; j<n_steps; j++ )
                W[j] = W[j-1] + dW[j-1]*sqrt(h);
        
        ml_free ( dW );
    }
}




void gen_BM( double h, int n_steps, double **& W, int n_runs, int mode )
{
    if ( W == 0 ) W = ml_alloc<double > ( n_runs, n_steps );
    static ml_random rng;
    
    if ( mode == BM_mode_1 )
    {
        int n_32 = div_up(n_steps,32);
        
        int32 *walk = ml_alloc<int32 > ( n_steps );
        int32 *pool = ml_alloc<int32 > ( n_32 );
        
        for (int run=0; run<n_runs; run++)
        {
            for (int k=0; k<n_32; k++)
                pool[k] = rng.gen_int();
            
            walk[0] = 0;
                
            int i=0, j=0;
                
            for ( int k=1; k<n_steps; k++ )
            {
                if ((pool[j] >> i) & 1) walk[k] = walk[k-1]+1;
                else walk[k] = walk[k-1]-1;
                    
                i++;
                    
                if ( i >= 32 )
                {
                    i=0;
                    j++;
                }
            }
            
            for ( int k=0; k<n_steps; k++ )
                W[run][k] = sqrt(h)*walk[k];
        }
        
        ml_free (walk);
        ml_free (pool);
    }
    
    if ( mode == BM_mode_2 )
    {
        int32 *walk = ml_alloc<int32 > ( n_steps );
        int32 *pool = ml_alloc<int32 > ( n_steps );
        
        for (int run=0; run<n_runs; run++)
        {
            for (int k=0; k<n_steps; k++)
                pool[k] = rng.gen_int();
            
            walk[0] = 0;
            
            for ( int k=1; k<n_steps; k++ )
                walk[k] = walk[k-1] + 2*ml_bitsum( pool[k-1] ) -32;
            
            for ( int k=0; k<n_steps; k++ )
                W[run][k] = sqrt(h/32)*walk[k];
        }
        
        ml_free (walk);
        ml_free (pool);
    }
    
    if ( mode == BM_mode_3 )
    {
        double * dW = ml_alloc<double > ( n_steps );
        
        for (int run=0; run<n_runs; run++)
        {
            rng.std_normal_rv( dW, n_steps);
            
            W[run][0] = 0;
            
            for ( int j=1; j<n_steps; j++ )
                W[run][j] = W[run][j-1] + dW[j-1]*sqrt(h);
        }
        
        ml_free ( dW );
    }
}



#endif
