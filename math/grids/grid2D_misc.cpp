#ifndef GRID2D_MISC_CPP
#define GRID2D_MISC_CPP

#ifndef GRID2D_MISC_H
	#include "grid2D_MISC.h"
#endif

template<classgridPars_ >
Ytype max(grid2D<gridPars_ > & rhs, int O)
{
	Ytype M = rhs(O,O);
	for (int i = O; i<rhs.n1-O; i++)
	for (int j = O; j<rhs.n2-O; j++)
		if (rhs(i,j) > M) M = rhs(i,j);
	return M;
}


#endif
