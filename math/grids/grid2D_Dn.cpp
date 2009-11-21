#ifndef GRID2D_Dn_CPP
#define GRID2D_Dn_CPP

#ifndef GRID2D_Dn_H
	#include "grid2D_Dn.h"
#endif




#ifdef GRID2D_Del1_FD

template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del1_FD(
	grid2D<Ytype,St1,Xt1 > & in,
	grid2D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond )
{
	out.copy(in,0);
	
	
	if ( bdcond == GRID_PERIODIC_BD )
	{		
		std::cerr << "GRID_PERIODIC_BD not yet implemented for 2D\n";
	}

	if ( bdcond == GRID_NONPER_BD_1 )
	{
		static int order_ = 0;
		static int N1_ = 0;
		static int N2_ = 0;
		
		int& M = order;
		int& N1 = in.n1;
		int& N2 = in.n2;
		
		static convolutionPlan<Ytype,double, double > plan1(N1);
		static convolutionPlan<Ytype,double, double > plan2(N2);
		
		static FD_D1_LR LR;
		
		if (N1 != N1_ || order != order_)
		{
			N1_ = N1;
			order_ = order;
						
			plan1.resize(N1);
			
			for (int i=0; i<N1; i++)
			{
				int L = i;
				int R = N1-i-1;
				
				if (L>=M && R>=M) { L = M; R = M; }
				if (L>24) { L = 24; }
				if (R>24) { R = 24; }
				
				double * ker = new double [L+R+1];
				
				for (int j = -L; j<=R; j++)
					ker[j+L] = LR(L,R,j);
				
				plan1.set_nl(i,L);
				plan1.set_nr(i,R);
				plan1.set_ker(i,&ker[L]);
			}
		}
		
		if (N2 != N2_ || order != order_)
		{
			N2_ = N2;
			order_ = order;
						
			plan2.resize(N2);
			
			for (int i=0; i<N2; i++)
			{
				int L = i;
				int R = N2-i-1;
				
				if (L>=M && R>=M) { L = M; R = M; }
				if (L>24) { L = 24; }
				if (R>24) { R = 24; }
				
				double * ker = new double [L+R+1];
				
				for (int j = -L; j<=R; j++)
					ker[j+L] = LR(L,R,j);
				
				plan2.set_nl(i,L);
				plan2.set_nr(i,R);
				plan2.set_ker(i,&ker[L]);
			}
		}		
		
		out = 0.0;
		
		for (int j = 0; j<N2; j++)		
			plan1.execute( &in(0,j), &in(1,j)-&in(0,j), &out(0,j), &out(1,j)-&out(0,j), 1.0/(out.dx1()) );
		
		for (int i = 0; i<N1; i++)
			plan2.execute( &in(i,0), &in(i,1)-&in(i,0), &out(i,0), &out(i,1)-&out(i,0), 1.0/(out.dx2()) );
	}
}

#endif







#endif
