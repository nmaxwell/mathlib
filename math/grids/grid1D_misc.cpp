#ifndef GRID1D_MISC_CPP
#define GRID1D_MISC_CPP

#ifndef GRID1D_MISC_H
	#include "grid1D_MISC.h"
#endif

template<classgridPars_ >
Ytype max(grid1D<gridPars_ > & rhs, int O=0)
{
	Ytype M = rhs(O);
	for (int i = O; i<rhs.n1-O; i++)
		if (rhs(i) > M) M = rhs(i);
	return M;
}


#endif
