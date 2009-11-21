#ifndef GRID3D_MISC_CPP
#define GRID3D_MISC_CPP

#ifndef GRID3D_MISC_H
	#include "grid3D_MISC.h"
#endif

template<classgridPars_ >
Ytype max(grid3D<gridPars_ > & rhs, int O=0)
{
	Ytype M = rhs(O,O,O);
	for (int i1 = O; i1<rhs.n1-O; i1++)
	for (int i2 = O; i2<rhs.n2-O; i2++)
	for (int i3 = O; i3<rhs.n3-O; i3++)
		if (rhs(i1,i2,i3) > M) M = rhs(i1,i2,i3);
	return M;
}


#endif
