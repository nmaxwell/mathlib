#ifndef GRID3D_LC_H
#define GRID3D_LC_H

#include "mathlib/math/LC.h"

#include "grid3D.h"

template<classgridPars_ >
void LC( grid3D<gridPars_ > & R, grid3D<gridPars_ > * L, YtypeScalar * C, int N );

template<classgridPars_ >
void LC ( 
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0
	    );

template<classgridPars_ >
void LC ( 
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0,
	    grid3D<gridPars_ > & L1, YtypeScalar C1
	    );

template<classgridPars_ >
void LC ( 
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0,
	    grid3D<gridPars_ > & L1, YtypeScalar C1,
	    grid3D<gridPars_ > & L2, YtypeScalar C2
	    );

template<classgridPars_ >
void LC ( 
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0,
	    grid3D<gridPars_ > & L1, YtypeScalar C1,
	    grid3D<gridPars_ > & L2, YtypeScalar C2,
	    grid3D<gridPars_ > & L3, YtypeScalar C3
	    );

template<classgridPars_ >
void LC ( 
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0,
	    grid3D<gridPars_ > & L1, YtypeScalar C1,
	    grid3D<gridPars_ > & L2, YtypeScalar C2,
	    grid3D<gridPars_ > & L3, YtypeScalar C3,
	    grid3D<gridPars_ > & L4, YtypeScalar C4
	    );

template<classgridPars_ >
void LC ( 
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0,
	    grid3D<gridPars_ > & L1, YtypeScalar C1,
	    grid3D<gridPars_ > & L2, YtypeScalar C2,
	    grid3D<gridPars_ > & L3, YtypeScalar C3,
	    grid3D<gridPars_ > & L4, YtypeScalar C4,
	    grid3D<gridPars_ > & L5, YtypeScalar C5
	    );
	    
template<classgridPars_ >
void LC ( 
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0,
	    grid3D<gridPars_ > & L1, YtypeScalar C1,
	    grid3D<gridPars_ > & L2, YtypeScalar C2,
	    grid3D<gridPars_ > & L3, YtypeScalar C3,
	    grid3D<gridPars_ > & L4, YtypeScalar C4,
	    grid3D<gridPars_ > & L5, YtypeScalar C5,
	    grid3D<gridPars_ > & L6, YtypeScalar C6
	    );

template<classgridPars_ >
void LC (
		 grid3D<gridPars_ > & R,
		 grid3D<gridPars_ > & L0, YtypeScalar C0,
	    grid3D<gridPars_ > & L1, YtypeScalar C1,
	    grid3D<gridPars_ > & L2, YtypeScalar C2,
	    grid3D<gridPars_ > & L3, YtypeScalar C3,
	    grid3D<gridPars_ > & L4, YtypeScalar C4,
	    grid3D<gridPars_ > & L5, YtypeScalar C5,
	    grid3D<gridPars_ > & L6, YtypeScalar C6,
	    grid3D<gridPars_ > & L7, YtypeScalar C7
	    );




#endif
