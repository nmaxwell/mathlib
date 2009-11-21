#ifndef GRID1D_DEL2_SPEC_CPP
#define GRID1D_DEL2_SPEC_CPP

#ifndef GRID1D_DEL2_SPEC_H
	#include "grid1D_del2_SPEC.h"
#endif


int cutoff_freq_tool(complex<double> * X, int N,int stride, double eps)
{
	double S = 0.0;
	for (int i = 0; i<N-1; i++)
		S += fabs(X[i*stride]);
	double s = 0.0;
	int i;
	for (i = 0; i<N-1; i++)
	{
		s += abs(X[i*stride]);
		if (1.0 - s/S < eps) break;
	}
	return i;
}

template<class Xtype >
void Del2_SPEC1(
	grid1D<complex<double>,double,Xtype > & in,
	double eps )
{
	grid1D<complex<double>,double,Xtype > Q;
	FFT(in,Q);
	GRID1D_LOOP( Q )
		Q(i) *= -Q.k1(i)*Q.k(i);
		
	double kc = Q.k1(cutoff_freq_tool(Q.array[Q.n1/2],Q.n1/2-1,-1,eps));
	
	if (isnan(kc) || isinf(kc))
		
	
	// too narrow
	if (c <= N*0.1)
		c = N*0.25;
	
	if (c >= N*0.1 and c <= N*0.9)
	{
		//cut off frequency is in a nice area
		
		double sig = 1.0/(in.k1(c)*fudge);		
		return hdaf_deltaHat(m,sig*sqrt(2.0*m+1.0));
	}
	else
	{
		//cut off frequency is in a nice area
		
		double sig = 1.0/(in.k1(c)*fudge);		
		return hdaf_deltaHat(100,sqrt(2.0*100+1.0)*(in.k1(N)*10.0));
	}	
	
}








template<class Xtype >
hdaf_deltaHat
	hdaf_autofilter_fftw_complex(
		grid1D<complex<double>,double,Xtype > & in,
		int m,
		double eps,
		double fudge )
{
	int c = 0;
	double S = 0.0;
	double s = 0.0;
	int N = in.n1/2-1;
	
	//S = sum of entire grid, excuding DC term
	for (int i = 1;i<=N;i++)
		S += norm(in[i]);
	
	assert(!isnan(S));
	assert(!isinf(S));
	
	//loop c from outside inwards, untill s, 
	//the cumulative sum along c, is the fraction eps of S 
	for (c=N;c>0;c--)
	{
		s += norm(in(c));
		if (s > eps*S) break;
	}
	
	if (isnan(s) || isinf(s))
		//return constfunctor<double,double> (1.0);
		return hdaf_deltaHat(100,sqrt(2.0*100+1.0)*(in.k1(N)*10.0));
	
	// too narrow
	if (c <= N*0.1)
		c = N*0.25;
	
	if (c >= N*0.1 and c <= N*0.9)
	{
		//cut off frequency is in a nice area
		
		double sig = 1.0/(in.k1(c)*fudge);		
		return hdaf_deltaHat(m,sig*sqrt(2.0*m+1.0));
	}
	else
	{
		//cut off frequency is in a nice area
		
		double sig = 1.0/(in.k1(c)*fudge);		
		return hdaf_deltaHat(100,sqrt(2.0*100+1.0)*(in.k1(N)*10.0));
	}
}





/*
template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
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
		
		static FD_Del2_LR LR;
		
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
			plan1.execute( &in(0,j), &in(1,j)-&in(0,j), &out(0,j), &out(1,j)-&out(0,j), 1.0/(out.dx1()*out.dx1()) );
		
		for (int i = 0; i<N1; i++)
			plan2.execute( &in(i,0), &in(i,1)-&in(i,0), &out(i,0), &out(i,1)-&out(i,0), 1.0/(out.dx2()*out.dx2()) );
	}	
}

*/


#endif
