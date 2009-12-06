#ifndef GRID2D_DEL2_fft_H
#define GRID2D_DEL2_fft_H

#include "grid2D.h"

//#include <mathlib/math/hdaf/hdaf.h>





template<class St1,class Xt1,class St2,class Xt2 >
void Del2_fft(
	grid2D<double,St1,Xt1 > & in,
	grid2D<double,St2,Xt2 > & out,
    double eps = 5E-15 )
 )
{
    out.copy(in);
    
    grid2D<complex<double>, St1, Xt1 > Q;
    
    FFT(in, Q );
    
    
}




#endif
