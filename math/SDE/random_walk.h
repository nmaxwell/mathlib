#ifndef ML_RAND_WALK_H
#define ML_RAND_WALK_H

#include "../random/ml_random.h"

inline int ml_bitsum( uint32 const & x );

void int_to_random_walk( uint32 * pool, int n, int32 * walk );

void int_to_random_walk( int32 * pool, int n, int32 * walk )
{
    int_to_random_walk( (uint32 *)pool, n, walk );
}


void gen_BM( double time_step, int n_steps, double *& W, int32 *walk, int32 *pool );

void gen_BM( double time_step, int n_steps, double *& W );

void gen_BM( double time_step, int n_steps, double **& W, int n_runs );



#endif
