#ifndef ML_HDAF_CPP
#define ML_HDAF_CPP

#ifndef HDAF_H
	#include "hdaf.h"
#endif



class HDAFfilter :  public functor<double, double>
{
public:
	double k0;
	double delta;
	double * lognfact;
	int mmax;
	
	HDAFfilter(int K, double D, double * lognfact_, int M) : k0(K),delta(D),lognfact(lognfact_),mmax(M) {};
	
	double operator() (double k)
	{
		double sigma = ml_sqrt2/(delta);
		int m = ceil(k0*k0/(delta*delta)-0.5);
		assert(m <= mmax);
		double s = 0.0;
		double x = k*k*sigma*sigma/2.0;
		double r = log(fabs(k))*2.0+log(sigma)*2.0;
		if (fabs(k)<1E-12) return 1.0;
		for (int n = 0;n<=m;n++)
			s += exp(-x+r*n-ml_log2*n-lognfact[n]);
		return s;
	}
};




















/*


template<class D> D hdaf_delta(D x, int M, D sigma)
{
	D y = x/(sigma*sqrt2);
	D d = 0.0;
	for (int n = 0;n<=M/2;n++)
		d += laguerre<D>(n, 0.5, x);
	d /= sigma*sqrt2pi;
	d *= exp(-y*y);
	return d;
}

	double y = x/(sigma*sqrt2);
	double d = 0.0;
	for (int n = 0;n<=M/2;n++)
		d += laguerre<double>(n, 0.5, x);
	d /= sigma*sqrt2pi;
	d *= exp(-y*y);
	return d;
}

template<> mp_real hdaf_delta(mp_real x, int M, mp_real sigma)
{
	static mp_real sqrt2 = exp(x._log2*0.5);	
	mp_real y = x/(sigma*sqrt2);
	mp_real d = 0.0;
	for (int n = 0;n<=M/2;n++)
		d += laguerre<mp_real>(n, 0.5, x);
	d /= sigma*sqrt(x._pi*2.0);
	d *= exp(-y*y);
	return d;
}


template<class D> D hdaf_delta_2(D x, int M, D sigma)
{

	D y = x*x/(sigma*sigma*2.0);
	D d = 0.0;
	for (int n = 0;n<=M/2;n++)
		d += laguerre<D>(n+1, -0.5, y);
	d /= -pow(sigma,3)*sqrt2pi*0.25;
	d *= exp(-y);
	return d;
}

template<> double hdaf_delta_2(double x, int M, double sigma)
{
	double y = x/(sigma*sqrt2);
	double d = 0.0;
	for (int n = 0;n<=M/2;n++)
	// d += hermite(2*n+2,y)*pow(-4.0,-n)/dfactorial[n];
		d += laguerre<double>(n+1,-0.5, y*y)*(double)(-4*n-4);
	d /= sigma*sigma*sigma*sqrt2pi*2.0;
	d *= exp(-y*y);
	return d;
}

template<> mp_real hdaf_delta_2(mp_real x, int M, mp_real sigma)
{
	static mp_real sqrt2 = exp(x._log2*0.5);
	mp_real y = x/(sigma*sqrt2);
	mp_real y2 = x*x/(sigma*sigma*2.0);
	mp_real d = 0.0;
	for (int n = 0;n<=M/2;n++)
	// d += hermite(2*n+2,y)*pow(-4.0,-n)/dfactorial[n];
		d += laguerre<mp_real>(n+1,-0.5, y2)*(mp_real)(-2*n-2);
	d /=pow(sigma,3)*sqrt(x._pi*2.0);
	d *= exp(-y2);
	return d;
}


*/






/*

template<class D> D partialExpProduct(D x, int n)
{
	D sum = 1.0;
	D s = 1.0;
	for (int j = 1;j<=n;j++)
	{
		s *= x/j;
		sum += s;
	}
	return sum*exp(-x);
}

template<class D> D partialExpProduct_d1(D x, int n);
template<class D> D partialExpProduct_d2(D x, int n);

template<>
double partialExpProduct_d1(double x, int n)
{
	return -pow(x,n)*exp(-x)/dfactorial[n];
}

template<>
double partialExpProduct_d2(double x, int n)
{
	return pow(x,n-1)*exp(-x)*(x-n)/dfactorial[n];
}



template<class D> D hdaf_window(D k, int n, D sigma)
{
	return partialExpProduct<D> ( k*k*sigma*sigma/2.0, n );
}

template<class D> D hdaf_window_d1(D k, int n, D sigma);
template<class D> D hdaf_window_d2(D k, int n, D sigma);

template<>
double hdaf_window_d1(double k, int n, double sigma)
{
	return partialExpProduct_d1<double> ( k*k*sigma*sigma/2.0,n )*k*sigma*sigma;
}

template<>
double hdaf_window_d2(double k, int n, double sigma)
{
	return (partialExpProduct_d2<double> ( k*k*sigma*sigma/2.0,n )*k*k*sigma*sigma+
		partialExpProduct_d1<double> ( k*k*sigma*sigma/2.0,n ))*sigma*sigma;
}


*/


















/*
void HDAF_interpolate(vector<double>& X, vector<double> & Y, grid1D<double> & R, double M, double sig)
{
	double S;
	for (int j = 0;j<R.NX;j++)
	{
		S = 0;
		for (int i = 1;i<X.size()-1;i++)
			S += hdaf_delta( R.X(j)-X[i], M,sig) * Y[i] * (X[i]-X[i-1]);
		R(j) = S;
	}
}
*/




/*
void Del2_5(grid1D<double> & G, grid1D<double> & Del, int M, double sig, int W)
{
	static double * kernel;
	static int M_,W_;
	static double sig_;
	static double DX_;
	
	
		
	if ( (int)(kernel) == 0 || W != W_  || M != M_ || sig != sig_ || DX_ != G.DX)
	{
		mp::mp_init(100);
		
		if (kernel) delete [] kernel;
		kernel = new double [W/2+1];
		M_ = M;
		W_ = W;
		sig_ = sig;
		DX_ = G.DX;
		
		for (int i = 0;i<=W/2;i++)	
			kernel[i] =dble( hdaf_delta_2<mp_real>( (double)i *  G.DX  ,M,sig));
	}
	
	int i,k;
		
	double s;
	for (i = 0;i<G.NX;i++)
	{
		s = 0.0;
		for (k = 1;k<=W/2;k++) {
			if ( i-k < 0 ) s += kernel[k] *  G(i-k+ G.NX);
			else s += kernel[k] *  G(i-k); }
		for (k = 0;k<=W/2;k++) {
			if ( i+k >=  G.NX ) s += kernel[k] * G(i+k- G.NX);
			else s += kernel[k] * G(i+k); }
		Del(i) = s* G.DX;
	} 
	
}
*/










/*



double HDAF0(double x, int M, double sig)
{
	double  X = x/(sqrt2*sig);
	double D = 0.0;
	for (int i = 0;i<=M;i++)
		D += hermite(2*i, X)/(dfactorial[i]*pow(-1.0,i)*pow(2.0,2*i));
	D *= exp(-X*X)/(sq2pi*sig);
	return D;
}


mp_real HDAF1(double x, double M, double sig)
{
	mp_real X(x);
	X /= sqrt(mp_real(2.0))*sig;
	mp_real D(0.0);
	for (int i = 0;i<=M;i++)
		D += hermite(2*i, X)/(gamma((mp_real)(i+1))*pow(-1.0,i)*pow(2.0,2*i));
	D *= exp(-X*X)/(sqrt(X._pi*2.0)*sig);
	return D;
}

mp_real HDAF2(double x, int M, double sig)
{
	mp_real X(x);
	X /= sqrt(mp_real(2.0))*sig;
	mp_real D(0.0);
	mp_real Un(0.0);
	mp_real Un1(0.0);
	mp_real Un2(0.0);	
	int n = 2;
	Un2 = 1.0;
	Un1 = X*X*(-2.0) + 1.0;
	if (M == 0)
		D = Un2;
	if (M >= 1)
		D = Un2+Un1;
	while (n <= M)
	{
		Un  = Un1 * (X*X*(-2.0) + 4.0*n-3.0)/(double)n + Un2 * (double)(-4*n*n+n-6)/(double)(n*n-n);
		D += Un;
		n++;
		Un2 = Un1;
		Un1 = Un; 
	}	
	D *= exp(-X*X)/(sqrt(X._pi*2.0)*sig);
	return D;
}


double HDAF3(double x, int M, double sig)
{
	double X = x;
	X /= sqrt2*sig;
	double D,Un,Un1,Un2;
	int n = 2;
	Un2 = 1.0;
	Un1 = X*X*(-2.0) + 1.0;
	if (M == 0)
		D = Un2;
	if (M >= 1)
		D = Un2+Un1;
	while (n <= M)
	{
		Un  = Un1 * (X*X*(-2.0) + 4.0*n-3.0)/(double)n + Un2 * (double)(-4*n*n+n-6)/(double)(n*n-n);
		D += Un;
		n++;
		Un2 = Un1;
		Un1 = Un; 
	}	
	D *= exp(-X*X)/(sq2pi*sig);
	return D;
}











void HDAF_1(vector<double>& X, vector<double> & Y, grid1D<double> & R, double M, double sig)
{
	double S;
	for (int j = 0;j<R.NX;j++)
	{
		S = 0;
		for (int i = 1;i<X.size()-1;i++)
			S += HDAF0(R.X(j)-X[i],M,sig)*Y[i]*(X[i]-X[i-1]);
		R(j) = S;
	}		
}


void HDAF_1(grid1D<double> & F, grid1D<double> & R, double M, double sig)
{
	double S;
	for (int j = 0;j<R.NX;j++)
	{
		S = 0;
		for (int i = 0;i<F.NX;i++)
			S += HDAF0(R.X(j)-F.X(i),M,sig)*F(i);
		R(j) = S*F.DX;
	}		
}

double HDAF1(vector<double>& X, vector<double> & Y, double x, double M, double sig)
{
	double S = 0;
		for (int i = 1;i<X.size()-1;i++)
			S += HDAF0(x-X[i],M,sig)*Y[i]*(X[i]-X[i-1]);
	return S;
}










double hdafDelta(double x, double M, double sig)
{
	x /= sqrt2*sig;
	double D = 0.0;
	for (int n = 0; n<=M/2; n++)
		D += hermite(2*n,x)*pow(-4.0,-n)/dfactorial[n];
	D *= exp(-x*x)/(sq2pi*sig);
	return D;
}


double hdafDelta_2(double x, double M, double sig)
{
	double y = x /( sig * sqrt(2.0));
	double D = 0.0;
	for (int n = 0; n<=M/2; n++)
		D += hermite(2*n+2,y)*pow(-4.0,-n)/dfactorial[n];
	D *= exp(-x*x/(2.0*sig*sig))/(2.0*sqrt(2.0*pi)*sig*sig*sig);
	return D;
}

double hdafDeltaHat_M(double xi, double M, double sig)
{
	double xisig2 = xi*xi*sig*sig;
	double D = 0.0;
	for (int n = 0; n<=M; n++)
		D += pow(xisig2,n)/((double)pow(2,n)*dfactorial[n]);
	D *= exp(-xisig2/2.0);	
	return D;
}

*/
/*

void Del2_5(grid1D<double> & G, grid1D<double> & Del, int M, double sig, int W)
{
	static double * kernel;
	static int M_,W_;
	static double sig_;
	static double DX_;
	
	if ( (int)(kernel) == 0 || W != W_  || M != M_ || sig != sig_ || DX_ != G.DX)
	{
		if (kernel) delete [] kernel;
		kernel = new double [W/2+1];
		M_ = M;
		W_ = W;
		sig_ = sig;
		DX_ = G.DX;
		
		for (int i = 0;i<=W/2;i++)
			kernel[i] = hdafDelta_2( (double)i *  G.DX  ,M,sig);
	}
	
	int i,k;
		
	double s;
	for (i = 0;i<G.NX;i++)
	{
		s = 0.0;
		for (k = 1;k<=W/2;k++) {
			if ( i-k < 0 ) s += kernel[k] *  G(i-k+ G.NX);
			else s += kernel[k] *  G(i-k); }
		for (k = 0;k<=W/2;k++) {
			if ( i+k >=  G.NX ) s += kernel[k] * G(i+k- G.NX);
			else s += kernel[k] * G(i+k); }
		Del(i) = s* G.DX;
	} 
	
}

*/



//----- Laplacian operators 

/*
template<>
void grid1D<double>::Del2_3(grid1D<double> & Del, int M, double sig)
{
	static grid1D<complex<double> > q;
	static grid1D<complex<double> > delHat;
	if (q.NX != NX)
		q = grid1D<complex<double> >(XMIN,XMAX,NX);	
	
	if (delHat.NX != NX)
	for (int i = 0;i<NX;i++)
	{
		delHat = grid1D<complex<double> >(XMIN,XMAX,NX);
		double f = fftMap(i,NX,XMAX-XMIN);
		delHat(i) = hdafDeltaHat_M(f, M, sig);
		delHat(i) *= -_4pi2*f*f;
	}
	
	FFT(*this,q);
	
	q *= delHat;
	
	IFFT(q,Del);
}

*/















#endif


