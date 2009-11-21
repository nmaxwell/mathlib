#ifndef FD_KERNELS_H
#define FD_KERNELS_H

#include "kernels.h"
 
/*

#define INCLUDE_FD_D1_LR
#define INCLUDE_FD_D2_LR
#define INCLUDE_FD_D3_LR
#define INCLUDE_FD_D4_LR
#define INCLUDE_FD_D5_LR
#define INCLUDE_FD_D6_LR

*/



struct FD_Del2_cen {
	inline double operator()(int const & m,int const & i);
};

struct FD_Del2_LR {
	inline double operator()(int const & L,int const & R,int const & i);
};


struct FD_D1_LR {
	double operator()(int L, int R,int i);
};

struct FD_D2_LR {
	double operator()(int L, int R,int i);
};

struct FD_D3_LR {
	double operator()(int L, int R,int i);
};

struct FD_D4_LR {
	double operator()(int L, int R,int i);
};

struct FD_D5_LR {
	double operator()(int L, int R,int i);
};

struct FD_D6_LR {
	double operator()(int L, int R,int i);
};

double FD_Dn_LR(int n, int L, int R, int i)
{
    static FD_D1_LR D1;
    static FD_D2_LR D2;
    static FD_D3_LR D3;
    static FD_D4_LR D4;
    static FD_D5_LR D5;
    static FD_D6_LR D6;
    
    if (n==0)
		if(i==0) return 1.0;
		else return 0.0;
  
    #ifdef INCLUDE_FD_D1_LR
    if (n == 1)
        return -D1(L,R,i);
    #endif
    
    #ifdef INCLUDE_FD_D2_LR
    if (n == 2)
        return D2(L,R,i);
    #endif
    
    #ifdef INCLUDE_FD_D3_LR
    if (n == 3)
        return -D3(L,R,i);
    #endif
    
    #ifdef INCLUDE_FD_D4_LR
    if (n == 4)
        return D4(L,R,i);
    #endif
    
    #ifdef INCLUDE_FD_D5_LR
    if (n == 5)
        return -D5(L,R,i);
    #endif
    
    #ifdef INCLUDE_FD_D6_LR
    if (n == 6)
        return D6(L,R,i);
    #endif
    
    
    _debug_here(999);
    cout << "you havent included the diff kernel: #define INCLUDE_FD_D?_LR " << endl;
    
    
    return 0;
}





#endif
