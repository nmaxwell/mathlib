#ifndef GRID3D_DEL2_FD_CPP
#define GRID3D_DEL2_FD_CPP

#ifndef GRID3D_DEL2_FD_H
	#include "grid3D_del2_FD.h"
#endif


template<class Ytype,class St1,class Xt1,class St2,class Xt2 >
void Del2_FD(
	grid3D<Ytype,St1,Xt1 > & in,
	grid3D<Ytype,St2,Xt2 > & out,
	int order,
	int bdcond,
	int BD_order )
{
	out.copy(in);
	
	if ( bdcond == GRID_PERIODIC_BD )
	{	    
		static vector<grid3D<complex<double>,St1,Xt1 > * > kernels;
        static vector<int > orders;
        static vector<int > N1,N2,N3;
        
        static int I = 0;
        bool found = 0;
        
        if (I)
            if ( N1[I] == in.n1 && N2[I] == in.n2 && N3[I] == in.n3 && orders[I] == order )
                found = true;
        
        if (!found)
        for ( I = 0; I<kernels.size(); I++)
        {
        	if ( N1[I] == in.n1 && N2[I] == in.n2 && N3[I] == in.n3 && orders[I] == order )
        	{
        		found = 1;
        		break;
        	}
        }
        
        if (!found)
        {
        	//cout << "new" << endl;
        	I = kernels.size();
        	N1.push_back( in.n1 );
        	N2.push_back( in.n2 );
            N3.push_back( in.n3 );
        	orders.push_back( order );
        	kernels.push_back( new grid3D<complex<double>,St1,Xt1 > );
        	
            *(kernels[I]) = grid3D<complex<double>,St1,Xt1 > ( N1[I], 0.0, 0.0, N2[I], 0.0, 0.0, N3[I], 0.0, 0.0  );
        	*kernels[I] = 0.0;
        	
        	for( int i=0; i<=orders[I]; i++)
        	{
        	    (*kernels[I])(i,0,0) += FD_DEL2_CENTERED_[orders[I] ][i];
        	    (*kernels[I])(0,i,0) += FD_DEL2_CENTERED_[orders[I] ][i];
        	    (*kernels[I])(0,0,i) += FD_DEL2_CENTERED_[orders[I] ][i];
        	}
        	
        	for( int i=1; i<=orders[I]; i++)
        	{
        	    (*kernels[I])(N1[I]-i,0,0) += FD_DEL2_CENTERED_[orders[I] ][i];
        	    (*kernels[I])(0,N2[I]-i,0) += FD_DEL2_CENTERED_[orders[I] ][i];
        	    (*kernels[I])(0,0,N3[I]-i) += FD_DEL2_CENTERED_[orders[I] ][i];
        	}        	
        	
        }
        
        
        static grid3D<complex<double>,St1,Xt1 > Q;
        
    //    double t0 = getRealTime();
        FFT(in,Q);
        
    //    double t1 = getRealTime();
        Q *= *kernels[I];
        
    //    double t2 = getRealTime();
        Q *= (double)(N1[I]*N2[I]*N3[I]);
        
    //    double t3 = getRealTime();
        IFFT(Q,out);
        
    //    double t4 = getRealTime();
        
    //    cout << t4-t0 << "\t" << (t1-t0)/(t4-t0) << "\t" << (t2-t1)/(t4-t0) << "\t" << (t3-t2)/(t4-t0) << "\t" << (t4-t3)/(t4-t0) << "\t" << endl2;	
        
	   // double t3 = getRealTime();
	   // cout << t3-t1 << "\t" << t3-t2 << "\t" << (t2-t1)/(t3-t1) << "\t"  << endl;	
	}

	if ( bdcond == GRID_NONPER_BD_1 )
	{
		static int order_ = 0;
		static int N1_ = 0;
		static int N2_ = 0;
		static int N3_ = 0;
		
		const int& M = order;
		const int& BDM = BD_order;
		const int& N1 = in.n1;
		const int& N2 = in.n2;
		const int& N3 = in.n3;
		
		static convolutionPlan<Ytype,double, double > plan1(N1);
		static convolutionPlan<Ytype,double, double > plan2(N2);
		static convolutionPlan<Ytype,double, double > plan3(N3);
		
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
				if (L>BDM) { L = BDM; }
				if (R>BDM) { R = BDM; }
				
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
				if (L>BDM) { L = BDM; }
				if (R>BDM) { R = BDM; }
				
				double * ker = new double [L+R+1];
				
				for (int j = -L; j<=R; j++)
					ker[j+L] = LR(L,R,j);
				
				plan2.set_nl(i,L);
				plan2.set_nr(i,R);
				plan2.set_ker(i,&ker[L]);
			}
		}	
		
		if (N3 != N3_ || order != order_)
		{
			N3_ = N3;
			order_ = order;
						
			plan3.resize(N3);
			
			for (int i=0; i<N3; i++)
			{
				int L = i;
				int R = N3-i-1;
				
				if (L>=M && R>=M) { L = M; R = M; }
				if (L>BDM) { L = BDM; }
				if (R>BDM) { R = BDM; }
				
				double * ker = new double [L+R+1];
				
				for (int j = -L; j<=R; j++)
					ker[j+L] = LR(L,R,j);
				
				plan3.set_nl(i,L);
				plan3.set_nr(i,R);
				plan3.set_ker(i,&ker[L]);
			}
		}		
		
		out = 0.0;

        double t1,t2;
		t1 = getRealTime();
        
		for (int j = 0; j<N2; j++)
		for (int k = 0; k<N3; k++)
			plan1.execute( &in(0,j,k), &in(1,j,k)-&in(0,j,k), &out(0,j,k), &out(1,j,k)-&out(0,j,k), 1.0/(out.dx1()*out.dx1()) );
		
		
		for (int i = 0; i<N1; i++)
		for (int k = 0; k<N3; k++)
			plan2.execute( &in(i,0,k), &in(i,1,k)-&in(i,0,k), &out(i,0,k), &out(i,1,k)-&out(i,0,k), 1.0/(out.dx2()*out.dx2()) );
		
				
		for (int i = 0; i<N1; i++)
		for (int j = 0; j<N2; j++)
			plan3.execute( &in(i,j,0), &in(i,j,1)-&in(i,j,0), &out(i,j,0), &out(i,j,1)-&out(i,j,0), 1.0/(out.dx3()*out.dx3()) );

			t2 = getRealTime();
			cout << t2-t1 <<  endl;
		
	}	
}




#endif
