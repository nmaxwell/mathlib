#ifndef GRID2D_LC_CPP
#define GRID2D_LC_CPP

#ifndef GRID2D_LC_H
	#include "grid2D_LC.h"
#endif


template<classgridPars_ >
void LC( grid2D<gridPars_ > & R, grid2D<gridPars_ > * L, YtypeScalar * C, int N )
{
	GRID2D_LINLOOP( R )
		R.array[i] = L[0].array[i]*C[0];
	for (int n = 1; n<N; n++)
	GRID2D_LINLOOP( R )
		 R.array[i] = L[n].array[i]*C[n];
}


template<classgridPars_ >
void LC ( 
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0
	    )
{
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0;
}

template<classgridPars_ >
void LC ( 
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0,
	    grid2D<gridPars_ > & L1, YtypeScalar C1
	    )
{
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0 +
				L1.array[i]*C1;
}

template<classgridPars_ >
void LC ( 
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0,
	    grid2D<gridPars_ > & L1, YtypeScalar C1,
	    grid2D<gridPars_ > & L2, YtypeScalar C2
	    )
{
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0 +
				L1.array[i]*C1 +
				L2.array[i]*C2;
}

template<classgridPars_ >
void LC ( 
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0,
	    grid2D<gridPars_ > & L1, YtypeScalar C1,
	    grid2D<gridPars_ > & L2, YtypeScalar C2,
	    grid2D<gridPars_ > & L3, YtypeScalar C3
	    )
{
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0 +
				L1.array[i]*C1 +
				L2.array[i]*C2 +
				L3.array[i]*C3;
}

template<classgridPars_ >
void LC ( 
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0,
	    grid2D<gridPars_ > & L1, YtypeScalar C1,
	    grid2D<gridPars_ > & L2, YtypeScalar C2,
	    grid2D<gridPars_ > & L3, YtypeScalar C3,
	    grid2D<gridPars_ > & L4, YtypeScalar C4
	    )
{
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0 +
				L1.array[i]*C1 +
				L2.array[i]*C2 +
				L3.array[i]*C3 +
				L4.array[i]*C4;
}

template<classgridPars_ >
void LC ( 
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0,
	    grid2D<gridPars_ > & L1, YtypeScalar C1,
	    grid2D<gridPars_ > & L2, YtypeScalar C2,
	    grid2D<gridPars_ > & L3, YtypeScalar C3,
	    grid2D<gridPars_ > & L4, YtypeScalar C4,
	    grid2D<gridPars_ > & L5, YtypeScalar C5
	    )
{
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0 +
				L1.array[i]*C1 +
				L2.array[i]*C2 +
				L3.array[i]*C3 +
				L4.array[i]*C4 +
				L5.array[i]*C5;
}

template<classgridPars_ >
void LC ( 
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0,
	    grid2D<gridPars_ > & L1, YtypeScalar C1,
	    grid2D<gridPars_ > & L2, YtypeScalar C2,
	    grid2D<gridPars_ > & L3, YtypeScalar C3,
	    grid2D<gridPars_ > & L4, YtypeScalar C4,
	    grid2D<gridPars_ > & L5, YtypeScalar C5,
	    grid2D<gridPars_ > & L6, YtypeScalar C6
	    )
{
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0 +
				L1.array[i]*C1 +
				L2.array[i]*C2 +
				L3.array[i]*C3 +
				L4.array[i]*C4 +
				L5.array[i]*C5 +
				L6.array[i]*C6;
}

template<classgridPars_ >
void LC (
		 grid2D<gridPars_ > & R,
		 grid2D<gridPars_ > & L0, YtypeScalar C0,
	    grid2D<gridPars_ > & L1, YtypeScalar C1,
	    grid2D<gridPars_ > & L2, YtypeScalar C2,
	    grid2D<gridPars_ > & L3, YtypeScalar C3,
	    grid2D<gridPars_ > & L4, YtypeScalar C4,
	    grid2D<gridPars_ > & L5, YtypeScalar C5,
	    grid2D<gridPars_ > & L6, YtypeScalar C6,
	    grid2D<gridPars_ > & L7, YtypeScalar C7
	    )
{				
	GRID2D_LINLOOP( L0 )
		R.array[i] =
				L0.array[i]*C0 +
				L1.array[i]*C1 +
				L2.array[i]*C2 +
				L3.array[i]*C3 +
				L4.array[i]*C4 +
				L5.array[i]*C5 +
				L6.array[i]*C6 +
				L7.array[i]*C7;
}




#endif
