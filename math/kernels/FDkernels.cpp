#ifndef FD_KERNELS_CPP
#define FD_KERNELS_CPP

#ifndef FD_KERNELS_H
	#include "FDkernels.h"
#endif

#include "FD_DEL2_CENTERED.h"
#include "FD_DEL2_LR.h"




inline double FD_Del2_cen::operator()(int const & m,int const & i)
{
	if (i < 0)
	{
	    if (i <= m) return FD_DEL2_CENTERED_[m][i];
	    else return 0.0;
	}
	else
	{
	    if (-i <= m) return FD_DEL2_CENTERED_[m][-i];
	    else return 0.0;
	}
}


inline double FD_Del2_LR::operator()(int const & L,int const & R,int const & i)
{
/*	if (L<R)
		return FD_DEL2_LR_[L][R][i+L];
	else
		return FD_DEL2_LR_[R][L][i+R-1];
	*/	
	return FD_DEL2_LR_[L][R][i+L];
}




#ifdef INCLUDE_FD_D1_LR

	#include "FD_D1_LR.h"
	
	double FD_D1_LR::operator()(int L,int R,int i)
	{
		return FD_D1_LR_[L][R][i+L];
	}

#endif

#ifdef INCLUDE_FD_D2_LR

	#include "FD_D2_LR.h"
	
	double FD_D2_LR::operator()(int L,int R,int i)
	{
		return FD_D2_LR_[L][R][i+L];
	}

#endif

#ifdef INCLUDE_FD_D3_LR

	#include "FD_D3_LR.h"
	
	double FD_D3_LR::operator()(int L,int R,int i)
	{
		return FD_D3_LR_[L][R][i+L];
	}

#endif

#ifdef INCLUDE_FD_D4_LR

	#include "FD_D4_LR.h"
	
	double FD_D4_LR::operator()(int L,int R,int i)
	{
		return FD_D4_LR_[L][R][i+L];
	}

#endif

#ifdef INCLUDE_FD_D5_LR

	#include "FD_D5_LR.h"
	
	double FD_D5_LR::operator()(int L,int R,int i)
	{
		return FD_D5_LR_[L][R][i+L];
	}

#endif

#ifdef INCLUDE_FD_D6_LR

	#include "FD_D6_LR.h"
	
	double FD_D6_LR::operator()(int L,int R,int i)
	{
		return FD_D6_LR_[L][R][i+L];
	}

#endif










/*
template <class D, class Dint> 
Dint applyCenStencil( D * d, D * f, D M, int dstride, int fstride)
{
	// d is the stencil, applied to f
	// uses d[0,1,2,..,M]
	// applies d symmetircally, so accesses f[-M,-M+1,...,-1,0,1,...,M]
	
	Dint S = 0.0;
	for (int i = 0;i<=M;i++)
		S += (Dint(d[i]))*Dint(f[i]);
	for (int i = 1;i<=M;i++)
		S += (Dint(d[i]))*Dint(f[-i]);
	return S;
}

template <class D, class Dint> 
Dint applyStencil( D * d, D * f, D M, int dstride, int fstride)
{
	// d is the stencil, applied to f
	// uses d[0,1,2,..,M]
	// accesses f[0,1,2,...,M]
	
	Dint S = 0.0;
	for (int i = 0;i<=M;i++)
		S += ((Dint)d[i])*f[i];
	return S;
}*/


#endif


