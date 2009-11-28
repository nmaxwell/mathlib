#ifndef ML_BM_H
#define ML_BM_H

#include "../random/ml_random.h"


#define BM_mode_1 1   // straight random walk
#define BM_mode_2 2  // modified random walk, so that each step is 32 steps of a regular random walk.
#define BM_mode_3 3  // W is a cumulative sum of dW, dW  is series of standard normal random variables, gotten from box-muller




void gen_BM( double h, int n_steps, double **& W, int n_runs, int mode = BM_mode_2 );

void gen_BM( double h, int n_steps, double *& W, int32 *walk, int32 *pool, int mode = BM_mode_2 );

//     gen_BM( double h, int n_steps, float **& W, int n_runs, int mode )
//     gen_BM( double h, int n_steps, float *& W, int mode )


//-----------

inline int ml_bitsum( uint32 const & x );

#endif
