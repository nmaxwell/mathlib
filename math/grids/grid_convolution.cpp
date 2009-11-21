#ifndef GRID_CONVOLUTION_CPP
#define GRID_CONVOLUTION_CPP

#include <typeinfo>


// design notes : convolution_plan
/*
	
	convolution_plan is a plan to apply a kernel on a strip 
	of N data points, of type X, with stride x_tride, x_ is 
	address of frist point.
	
	The kernel may change from point to point, mostly to take
	care of boundary conditions. K will be the type of data
	of the kernel, and K ** ker_ will be an array, of K 
	pointers, the kernels to be applied. T will be the type
	used in the sum, so extended precision might be used.
	
	y_ will be a pointer to the first element in the return
	array, which will have stride y_stride.
	
	The arrays ker_, n_l, n_r, will be contiguous, and and be of length N.
	
	nr_i will be the points to the right of the current x_i.
	nl_i will be the points to the left of the current x_i.
	
	So the ith kernel must be of length ( nr_i + nl_i + 1 ).
	
	kernels[i][0] will be the element at the center of the kernel.
	So, the loop will look like		
		
		for (int i = 0; i<N; i++)
		for (int j = nl_[i]; j <= nr_[j]; j++)
			sum += ker_[i][j]*x_[i+j];
			
	If j goes of the grid, it will wrap around.
		
	
	execute() will run the convolution.
	
	this only adds the convolution to y. so, y += sum*scale, so remember to set y=0 before hand.
*/

template<class X=double, class K=X, class T=X >
class convolutionPlan
{
public:
	
	int N;
	
	K ** ker_;
	
	int * nl_;
	int * nr_;
	
public:
	
	convolutionPlan(int N):N(N)
	{		
		ker_ = new K* [N];
		nl_ = new int [N];
		nr_ = new int [N];
	}
	
	~convolutionPlan() {
		if (ker_) delete [] ker_;
		if (nl_) delete [] nl_;
		if (nr_) delete [] nr_; }
		
	void resize(int N_)
	{
		N = N_;
		
		if (ker_) delete [] ker_;
		if (nl_) delete [] nl_;
		if (nr_) delete [] nr_;
		
		ker_ = new K* [N];
		nl_ = new int [N];
		nr_ = new int [N];
	}
	
	void set_ker(int i, K* ker) { ker_[i] = ker; }
	void set_nl(int i, int nl) { nl_[i] = nl; }
	void set_nr(int i, int nr) { nr_[i] = nr; }	
	
public:

	int x_stride;
	X * x_;
	
	int y_stride;
	X * y_;
	
	void execute( X * x_, int x_stride, X * y_, int y_stride, T scale = 1.0 )
	{
	//	X x;
	//	K y;
	//	T t;
	//	cout << typeid(x).name() << "\t" << typeid(y).name() << "\t" << typeid(t).name() <<  endl;
		
		T sum;
		int i,j;
		
		for ( i = 0; i<N; i++)
		{
			sum = 0.0;
						
			if ( i-nl_[i] >= 0)
			{
				if ( nr_[i]+i < N )
				{
					// no wrap on left and no wrap on right
					for ( j = -nl_[i]; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
				}
				else
				{
					// no wrap on left, but is wrap on right
				//	cout << "no wrap on left, but is wrap on right " << endl;
					for ( j = -nl_[i]; i+j < N; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
					for ( ; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j-N)*x_stride];
				}	
			}
			else
			{
				if ( nr_[i]+i < N )
				{
					// is wrap on left, but no wrap on right
				//	cout << "is wrap on left, but no wrap on right " << endl;
					for ( j = -nl_[i]; j+i < 0 ; j++)
						sum += ker_[i][j]*x_[(i+j+N)*x_stride];
					for ( ; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
				}
				else
				{
					// is wrap on left, and is wrap on right
				//	cout << "is wrap on left, and is wrap on right " << endl;
					for ( j = -nl_[i]; j+i < 0 ; j++)
						sum += ker_[i][j]*x_[(i+j+N)*x_stride];
					for ( ;  i+j < N; j++)
						sum += ker_[i][j]*x_[(i+j)*x_stride];
					for ( ; j <= nr_[i]; j++)
						sum += ker_[i][j]*x_[(i+j-N)*x_stride];
				}
			}
			//cout << sum << " " ;
			y_[i*y_stride] += sum*scale;
		}
	}
	
};


//for ( j = -nl_[i]; j <= nr_[i]; j++)
//	sum += ker_[i][j]*x_[(i+j)*x_stride];

template<class X, class T >
void convolve_periodic( int n, X * x_, X * y_, X * z_,T scale = 1.0, int x_stride=1, int y_stride=1, int z_stride=1 )
{
    T sum = 0.0;
    int i,j,k;
    
    for (i=0; i<n; i++)
    {
        sum = 0.0;
        
        for (j=0; j<n; j++)
        {
            k = 0;
            if ( (i-j) < 0 ) k = n;
            if ( (i-j) >= n ) k = -n;
            
            sum += (T)x_[j*x_stride]*(T)y_[(i-j+k)*y_stride];
        }
        
        z_[i*z_stride] += scale*sum;
    }
};




/*
// not sure why this wont compile... 
#ifdef GRID3D_H    

    template<classgridPars_ >
    void convolve_periodic(grid3D<gridPars_ > & x, grid3D<gridPars_ > & y, grid3D<gridPars_ > & z, YtypeScalar scale)
    {
        z.copy(x);
        Ytype sum;
        
        int i1,i2,i3, j1,j2,j3, k1,k2,k3;
        
        for ( i1=0; i1<x.n1; i1++)
        for ( i2=0; i2<x.n2; i2++)
        for ( i3=0; i3<x.n3; i3++)
        {
            sum = 0.0;
            
            for ( j1=0; j1<x.n1; j1++)
            for ( j2=0; j2<x.n2; j2++)
            for ( j3=0; j3<x.n3; j3++)
            {
                k1=0;
                k2=0;
                k3=0;
                
                if ( (i1-j1) < 0 ) k1 = x.n1;
                if ( (i2-j2) < 0 ) k2 = x.n2;
                if ( (i3-j3) < 0 ) k3 = x.n3;
                
                if ( (i1-j1) >= x.n1 ) k1 = -x.n1;
                if ( (i2-j2) >= x.n2 ) k2 = -x.n2;
                if ( (i3-j3) >= x.n3 ) k3 = -x.n3;
                
                sum += x(i1-j1+k1, i2-j2+k2, i3-j3+k3 )*y(j1, j2, j3 );
            }
            
            z(i1,i2,i3) = sum*scale;
        }
    }

#endif*/


#endif

