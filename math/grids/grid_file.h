#ifndef GRID_FILE_OPS_H
#define GRID_FILE_OPS_H

#include "grid.h"

#ifdef GRID1D_H
	
	template<classgridPars_ >
	int writeFile( grid1D<gridPars_ > & G, const char * fname );
	
	template<classgridPars_ >
	int readFile( grid1D<gridPars_ > & G, const char * fname );
	
#endif

#ifdef GRID2D_H
	
	template<classgridPars_ >
	int writeFile( grid2D<gridPars_ > & G, const char * fname );
	
	template<classgridPars_ >
	int readFile( grid2D<gridPars_ > & G, const char * fname );
	
#endif

#ifdef GRID3D_H
	
	template<classgridPars_ >
	int writeFile( grid3D<gridPars_ > & G, const char * fname );
	
	template<classgridPars_ >
	int readFile( grid3D<gridPars_ > & G, const char * fname );
	
#endif

#endif
