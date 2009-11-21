#ifndef GRID1D_Dn_CPP
#define GRID1D_Dn_CPP

#ifndef GRID1D_Dn_H
	#include "grid1D_Dn.h"
#endif




#ifdef GRID1D_Del1_FD

template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del1_FD(
	grid1D<Ytype,St1,Xt1 > & in,
	grid1D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond )
{
	out.copy(in);
	
	if ( bdcond == GRID_PERIODIC_BD )
	{		
		static int order_ = 0;
		static int N_ = 0;		
		
		int& M = order;
		int& N = in.n1;
		
		static convolutionPlan<Ytype,double, double > plan(in.n1);		
		
		if (N != N_ || order != order_)
		{
			N_ = N;
			order_ = order;
			static double * ker = 0;
			static FD_D1_LR LR;
			
			if (ker) delete [] ker;
			ker = new Ytype [2*M+1];				
			for (int i = -M; i<=M; i++)
				ker[i+M] = LR(M,M,i);				
			plan.resize(N);
					
			for (int i = 0; i<N; i++)
			{
				plan.set_nl(i,M);
				plan.set_ker(i,&ker[M]);
				plan.set_nr(i,M);
			}
		}
		
		out = 0.0;
				
		plan.execute( &in(0), &in(1)-&in(0), &out(0), &out(1)-&out(0), 1.0/(out.dx1()) );
	}

	if ( bdcond == GRID_NONPER_BD_1 )
	{
		static int order_ = 0;
		static int N_ = 0;
		
		int& M = order;
		int& N = in.n1;
		
		static convolutionPlan<Ytype,double, double > plan(in.n1);
		
		if (N != N_ || order != order_)
		{
			N_ = N;
			order_ = order;
			static FD_D1_LR LR;
			
			plan.resize(N);
			
			for (int i=0; i<N; i++)
			{
				int L = i;
				int R = N-i-1;
				
				if (L>=M && R>=M) { L = M; R = M; } // center
				if (L>24) { L = 24; }
				if (R>24) { R = 24; }
				
				//cout << i << "\t" << L << "\t" << R << endl;
				
				double * ker = new double [L+R+1];
				
				for (int j = -L; j<=R; j++)
					ker[j+L] = LR(L,R,j);
				
				plan.set_nl(i,L);
				plan.set_nr(i,R);
				plan.set_ker(i,&ker[L]);
			}
		}
		
		out = 0.0;
		
		plan.execute( &in(0), &in(1)-&in(0), &out(0), &out(1)-&out(0), 1.0/(out.dx1()) );
	}
}

#endif







#endif
