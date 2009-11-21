#ifndef LAGUERRE_POLYS_CPP
#define LAGUERRE_POLYS_CPP

#ifndef LAGUERRE_POLYS_H
	#include "laguerre_polys.h"
#endif

template<class D>
inline D laguerre(int n, D const & alpha, D const & x)
{
	// scheme copied form wikipedia page
	// http://en.wikipedia.org/wiki/Laguerre_polynomials
	
	D R = 1.0;
	D bin = 1.0;
	for (int i = n;i>=1;i--)
	{
		bin *= (alpha+i)/(n+1-i);
		R = bin - x*R/D(i);
	}
	return R;
} 



#endif


