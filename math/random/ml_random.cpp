#ifndef ML_RANDOM_CPP
#define ML_RANDOM_CPP


#include "ml_random.h"


ml_random::ml_random()
:pool(0),pool_size(0),gen(0)
{
    assert(setup());
}

ml_random::ml_random( const ml_random & rhs )
:pool(0),pool_size(0),gen(0)
{
    assert(setup());
    set_pool( rhs.pool_size );
    refresh_pool();
}

ml_random::ml_random(int pool_size_ )
:pool(0),pool_size(0),gen(0)
{
    assert(setup());
    set_pool(pool_size_);
    refresh_pool();
}

bool ml_random::setup()
{
    bool OK=true;
    
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gen = gsl_rng_alloc(T);
    
    timespec currentTime;
    clock_gettime(CLOCK_REALTIME, &currentTime);
    gsl_rng_set(gen,currentTime.tv_nsec);
    
    /*if ( gsl_rng_max(gen) != (uint32)(-1) )
        OK = false;
    if ( gsl_rng_min(gen) != 0 )
        OK = false;
    */
    
    return OK;
}

ml_random::~ml_random() 
{
    gsl_rng_free (gen);
    
    if (pool) delete [] pool;
    pool = 0;
}


void ml_random::set_pool(int pool_size_ )
{
    if (!pool || pool_size != pool_size_ )
    {
        if (pool) delete [] pool;
        pool_size = pool_size_;
        pool = new char [ pool_size ];
    }
    
}

void ml_random::refresh_pool()
{
    if (pool_size < 800) set_pool(800);
    if (pool_size%8) set_pool(8-pool_size%8+pool_size);
    assert(!(pool_size%8));
    
    for (pos=0; pos<pool_size; pos+=4)
        *((int *)(pool+pos)) = gen_int();
    
    pos=0;
}    
    
double ml_random::next_double()
{
    if (pos+8 >= pool_size) refresh_pool();
    double val;
    
    uint64 & intval = *((uint64 *)(pool+pos));
    
    if ( (uint64)intval & (uint64)4503599627370496 )
        *((uint64 *)(pool+pos)) &= (uint64)4602678819172646911;
    else
        *((uint64 *)(pool+pos)) &= (uint64)4607182418800017407;
    
    val = *((double *)(pool+pos));
    pos += 8;
    
    return val;
}

float ml_random::next_float()
{
    if (pos+4 >= pool_size) refresh_pool();
    pos += 4;
    return *((float *)(pool+pos-4));
}

int ml_random::next_int()
{
    if (pos+4 >= pool_size) refresh_pool();
    pos += 4;
    return *((int *)(pool+pos-4));
}

uint ml_random::next_uint()
{
    if (pos+4 >= pool_size) refresh_pool();
    pos += 4;
    return *((uint *)(pool+pos-4));
}

bool ml_random::next_bool()
{
    return 0;
}

int ml_random::gen_int()
{
    return gsl_rng_get(gen);
}

double ml_random::gen_double()
{
    return gsl_rng_uniform(gen);
}

double ml_random::gen_double_nozero()
{
    return gsl_rng_uniform_pos (gen);
}
          		
#endif
