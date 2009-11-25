#ifndef ML_RANDOM_H
#define ML_RANDOM_H

#include <mathlib/math/std_math.h>
#include <mathlib/tools/std_tools.h>

#include <gsl/gsl_rng.h>

void int32_to_bool_decomp(int x, bool * decomp )
{
    for (int k=0; k<32; k++)
        decomp[k] = ((x >> k) & 1);
}

void int32_to_bool_decomp(int n, int * x, bool * decomp )
{
    for (int j=0; j<n; j++)
    for (int k=0; k<32; k++)
        decomp[32*j+k] = ((x[j] >> k) & 1);
}


#define box_muller( u1, u2, n1, n2 ) \
    n1 = sqrt(-log(u1)*2)*cos(ml_2pi*u2); \
    n2 = sqrt(-log(u1)*2)*sin(ml_2pi*u2);


class ml_random
{
public:
    char * pool;
    int pool_size; // in bytes
    gsl_rng * gen; // generator
    unsigned pos; // byte position in pool 
    
    ml_random();
    ~ml_random();
    ml_random(const ml_random & rhs );
    ml_random(int pool_size_ );
    
    bool setup();
    
    void * get_pool() { return pool; }
    
    void set_pool(const int pool_size);
    void refresh_pool();
    
    double next_double();
    float next_float();
    int next_int();
    uint next_uint();
    bool next_bool();
    
    int gen_int();
    double gen_double();
    double gen_double_nonzero();
    
    void std_normal_rv( double * x, int n );
};




      
        
#endif
