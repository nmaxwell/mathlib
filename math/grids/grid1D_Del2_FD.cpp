#ifndef GRID1D_DEL2_FD_CPP
#define GRID1D_DEL2_FD_CPP

#ifndef GRID1D_DEL2_FD_H
	#include "grid1D_del2_FD.h"
#endif

template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid1D<Ytype,St1,Xt1 > & in,
	grid1D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond )
{
	out.copy(in);
	
	if ( bdcond == GRID_PERIODIC_BD )
	{        
        static vector<grid1D<complex<double>,St1,Xt1 > * > kernels;
        static vector<int > orders;
        static vector<int > dims;
        
        static int I = 0;
        bool found = 0;
        
        if (I)
            if ( dims[I] == in.n1 && orders[I] == order )
                found = true;
        
        if (!found)
        for ( I = 0; I<kernels.size(); I++)
        {
        	if ( dims[I] == in.n1 && orders[I] == order )
        	{
        		found = 1;
        		break;
        	}
        }
        
        if (!found)
        {
        	//cout << "new" << endl;
        	I = kernels.size();
        	dims.push_back( in.n1 );
        	orders.push_back( order );
        	kernels.push_back( new grid1D<complex<double>,St1,Xt1 > );
        	
        	*(kernels[I]) = grid1D<complex<double> > ( dims[I], 0, 0 );
        	*kernels[I] = 0.0;
        	
        	for( int i=0; i<=orders[I]; i++)
        	    kernels[I][i] += FD_DEL2_CENTERED_[orders[I] ][i];
        	for( int i=1; i<=orders[I]; i++)
        	    kernels[I][dims[I]-i] += FD_DEL2_CENTERED_[orders[I] ][i];
        }
        
        static grid1D<complex<double>,St1,Xt1 > Q;
        FFT(in,Q);
        Q *= *kernels[I];
        Q *= (double)in.n1;
        IFFT(Q,out);
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
			static FD_Del2_LR LR;
			
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
		
		plan.execute( &in(0), &in(1)-&in(0), &out(0), &out(1)-&out(0), 1.0/(out.dx1()*out.dx1()) );
	}	
}


template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid1D<Ytype,St1,Xt1 > & in,
	grid1D<Ytype,St2,Xt2 > & out,
	int order,
	int L_BD,
	int R_BD )
{
	out.copy(in);
	
	static int order_ = 0;
	static int N_ = 0;
	static int L_BD_ = 0;
	static int R_BD_ = 0;
	
	int& M = order;
	int& N = in.n1;
	
	static convolutionPlan<Ytype,double, double > plan(in.n1);
	
	if (N != N_ || L_BD != L_BD_ || R_BD != R_BD_ || order != order_)
	{
		N_ = N;
		L_BD_ = L_BD;
		R_BD_ = R_BD;
		
		order_ = order;
		static FD_Del2_LR LR;
		
		plan.resize(N-L_BD-R_BD);
		
		for (int i=L_BD; i<N-R_BD; i++)
		{
			int L = i;
			int R = N-i-1;
			
			if (L>=M && R>=M) { L = M; R = M; } // center
			if (L>24) { L = 24; }
			if (R>24) { R = 24; }
			
			
			double * ker = new double [L+R+1];
			
			for (int j = -L; j<=R; j++)
				ker[j+L] = LR(L,R,j);
			
			plan.set_nl(i-L_BD,L);
			plan.set_nr(i-L_BD,R);
			plan.set_ker(i-L_BD,&ker[L]);
		}
	}
	
	out = 0.0;
	
	plan.execute( &in(L_BD), &in(1)-&in(0), &out(L_BD), &out(1)-&out(0), 1.0/(out.dx1()*out.dx1()) );
}

#endif
