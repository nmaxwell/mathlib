
/*
 * this file is in liu of precompiling the library and linking to it, as it is still in developement
 * 
 * include this file last when using mathlib
 * 
 */

#ifdef ML_RAND_WALK_H
    #include "math/SDE/random_walk.cpp"
#endif

#ifdef ML_ALLOC_H
	#include "tools/memory/ml_alloc.cpp"
#endif

#ifdef ML_ALLOC_H
	#include "tools/memory/ml_alloc.cpp"
#endif

#ifdef APPLY_HDAF_H
	#include "math/hdaf/apply_hdaf.cpp"
#endif

#ifdef STD_MATH_H
	#include "math/std_math.cpp"
#endif

#ifdef STD_TOOLS_H
	#include "tools/std_tools.cpp"
#endif

#ifdef ML_CONVOLUTION_H
	#include "math/convolution/convolution.cpp"
#endif

#ifdef ML_RANDOM_H
	#include "math/random/ml_random.cpp"
#endif

#ifdef KERNELS_H
	#include "math/kernels/kernels.cpp"
#endif

#ifdef ML_SIGNAL_REG_H
	#include "math/signals/ml_signal_reg.cpp"
#endif

#ifdef DIFFERENTIATION_H
	#include "math/differentiation/differentiation.cpp"
#endif

#ifdef FFT_H
	#include "math/transforms/FFT.cpp"
#endif

#ifdef EUVEC_H
	#include "math/vectors/euVec.cpp"
#endif

#ifdef FD_KERNELS_H
	#include "math/kernels/FDkernels.cpp"
#endif

#ifdef SIGNALS_H
	#include "math/signals/signals.cpp"
#endif

#ifdef SIGNAL_1_2_FFT_H
	#include "math/signals/fft.cpp"
#endif

#ifdef HERMITE_POLYS_H
	#include "math/polynomials/hermite_polys.cpp"
#endif

#ifdef LAGUERRE_POLYS_H
	#include "math/polynomials/laguerre_polys.cpp"
#endif

#ifdef ML_POLY_H
	#include "math/polynomials/ml_poly.cpp"
#endif

#ifdef MATRIX1_H
	#include "math/matrices/matrix1.cpp"
#endif

#ifdef HDAF_H
	#include "math/hdaf/hdaf.cpp"
#endif

#ifdef GRID_H
	#include "math/grids/grid.cpp"
#endif

#ifdef FD_KERNELS_H
	#include "math/kernels/FDkernels.cpp"
#endif


#ifdef PLOT2D_H
	#include "tools/graphing/plot2D.cpp"
#endif

#ifdef ARRAY1_H
	#include "tools/arrays/array1.cpp"
#endif

/*
#ifndef EUVEC_H
	#include "euVec.h"
#endif
*/











