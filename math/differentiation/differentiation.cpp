#ifndef DIFFERENTIATION_CPP
#define DIFFERENTIATION_CPP

#include "differentiation.h"



template<class X, class Y >
Y cheap_diff(X x, functor<X,Y > const & f, int p, X h, int M)
{
	// returns pth derivative of f around x, 
	// using finite differences with a step size of h, and using m points
	
//	FD_Dn_LR(int n, int L, int R, int i);
	
	Y sum = 0.0;
	for (int k=-M; k<=M; k++)
		sum += f(x-h*k)*FD_Dn_LR(p,M,M, k);
		
	return sum/(h*h);
}

double cheap_diff(double x, double (&f)(double), int p, double h, int M)
{
	// returns pth derivative of f around x, 
	// using finite differences with a step size of h, and using m points
	
//	FD_Dn_LR(int n, int L, int R, int i);
	
	double sum = 0.0;
	for (int k=-M; k<=M; k++)
		sum += f(x-h*k)*FD_Dn_LR(p,M,M, k);
		
	return sum*pow(h,-p);
}






template<class T >
void differentiate(
	T * in,
    T * & out,
    int n,
    T h,
    int diff_order,
    int in_stride,
    int out_stride )
{
	static bool do_once = 1;
    
	if(do_once)
	{
		cout << "warning: may be problem with odd derivatives" << endl;
		do_once = false;
	}
	
	if (out == 0) out = new T [ n*out_stride ];
	
	static vector<int > ns;
	static vector<int > diff_orders;
	
	static vector<hybrid_convolution_plan<T > > plans;
	
	static int I = 0;
	bool found = 0;
	
	if (I)
	if ( ns[I] == n && diff_orders[I] == diff_order )
	    found = true;
	
	if (!found)
	for ( I = 0; I<plans.size(); I++)
	{
		if ( ns[I] == n && diff_orders[I] == diff_order )
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
		I = plans.size();
		ns.push_back( n );
        diff_orders.push_back( diff_order );
        
        if (n <= 100)
        {
           /* plans.push_back ( 
                hybrid_convolution_plan<T >( n,  )
             );*/
        }
        else
        {
            int M = 24;
            
            plans.push_back ( hybrid_convolution_plan<T > () );
            plans[I].resize( n,M,M );
            
            T * ker_cen = new T [2*M+1];
            
            for (int i=0; i<=M; i++ )
            {
                ker_cen[M+i] = FD_Dn_LR(diff_order,M,M, i)*pow(-1,diff_order);
                ker_cen[M-i] = FD_Dn_LR(diff_order,M,M, -i)*pow(-1,diff_order);
                //cout << FD_Dn_LR(diff_order,M,M, i) << "\t" << ker_cen[M+i] << endl;
            }
            
            plans[I].set_cen_ker( &(ker_cen[M]), M,M );
            
            delete [] ker_cen;
            
            for (int i=0; i<M; i++)
            {
                int L = i;
				int R = M;
				
				if (L>M) { L = M; }
				
				T * ker = new T [L+R+1];
				
				for (int j = -L; j<=R; j++)
					ker[j+L] = FD_Dn_LR(diff_order,L,R,j);
				
				plans[I].set_l_nl(i,L);
				plans[I].set_l_nr(i,R);
				plans[I].set_l_ker(i,&ker[L]);
            }
            
            for (int i=0; i<M; i++)
            {                
                int L = M;
				int R = i;
				
				if (R>M) { R = M; }
				
				T * ker = new T [L+R+1];
				
				for (int j = -L; j<=R; j++)
					ker[j+L] = FD_Dn_LR(diff_order,L,R,j);
				
				plans[I].set_r_nl(i,L);
				plans[I].set_r_nr(i,R);
				plans[I].set_r_ker(i,&ker[L]);
            }
            
        }
        
	}	
	
	plans[I].execute(in,in_stride,out,out_stride,pow(h,-diff_order));	
}


#endif
