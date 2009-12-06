#ifndef GRID_CPP
#define GRID_CPP

#ifndef GRID_H
	#include "grid.h"
#endif

#ifdef GRID1D_H
	#include "grid1D.cpp"
	#include "grid1D_LC.cpp"
	#include "grid1D_norms.cpp"
	#include "grid1D_mod_assign_op.cpp"
	#include "grid1D_misc.cpp"
	
	#ifdef GRID1D_FFT_H
		#include "grid1D_FFT.cpp"
	#endif

	#ifdef GRID1D_DEL2_FD_H
		#include "grid1D_Del2_FD.cpp"
	#endif
	
	#ifdef GRID1D_DEL2_SPEC_H
		#include "grid1D_Del2_SPEC.cpp"
	#endif
	
	#ifdef GRID1D_Dn_H
		#include "grid1D_Dn.cpp"
	#endif
	
#endif

#ifdef GRID2D_H
	#include "grid2D.cpp"
	#include "grid2D_LC.cpp"
	#include "grid2D_norms.cpp"
	#include "grid2D_mod_assign_op.cpp"
	#include "grid2D_misc.cpp"	
	
	#ifdef GRID2D_FFT_H
		#include "grid2D_FFT.cpp"
	#endif

	#ifdef GRID2D_DEL2_FD_H
		#include "grid2D_Del2_FD.cpp"
	#endif
    
    #ifdef GRID2D_DEL2_fft_H
		#include "grid2D_Del2_fft.cpp"
	#endif
	
	#ifdef GRID2D_HDAF_H
		#include "grid2D_hdaf.cpp"
	#endif
	
#endif

#ifdef GRID3D_H
	#include "grid3D.cpp"
	#include "grid3D_LC.cpp"
	#include "grid3D_norms.cpp"
	#include "grid3D_mod_assign_op.cpp"
	#include "grid3D_misc.cpp"
	
	#ifdef GRID3D_FFT_H
		#include "grid3D_FFT.cpp"
	#endif

	#ifdef GRID3D_DEL2_FD_H
		#include "grid3D_Del2_FD.cpp"
	#endif
	
#endif

#ifdef GRID_FILE_OPS_H
	#include "grid_file.cpp"
#endif



#endif
