#include "mathlib/tools/graphing/plot2D.h"


#ifdef GRID1D_H

template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid1D(grid1D<Ytype,YtypeScalar,Xtype > & G,  plot2D & plot,  color3 & C )
{
	for (int i = 0;i<G.n1-1;i++)
		 plot.ptLine(G.x1(i),G(i),G.x1(i+1),G(i+1),C);
}

template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid1D(grid1D<Ytype,YtypeScalar,Xtype > & G, double y0, double y1,int nx, int ny, char * fname, color3 const & C=c3_black,color3 const & C_back=c3_creme)
{
	plot2D plot(G.a1,G.b1,y0,y1, nx,ny);
	plot= C_back;
	for (int i = 0;i<nx-1;i++)
		 plot.ptLine(G.x1(i),G(i),G.x1(i+1),G(i+1),C);
    plot.png(fname);
}

template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid1D(grid1D<Ytype,YtypeScalar,Xtype > & G, double x0, double x1, double y0, double y1,int nx, int ny, char * fname, color3 const & C=c3_black,color3 const & C_back=c3_creme)
{
	plot2D plot(x0,x1,y0,y1, nx,ny);
	plot= C_back;
	plot.axes_1(0,0, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 2,5,1,5  );
	for (int i = 0;i<nx-1;i++)
		 plot.ptLine(G.x1(i),G(i),G.x1(i+1),G(i+1),C);
    plot.png(fname);
}



template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid1D_fftw(grid1D<complex<Ytype>,YtypeScalar,Xtype > & G, double Kmax, double y0, double y1,int nx, int ny, char * fname, color3 const & C_real=c3_blue,color3 const & C_imag=c3_red,color3 const & C_back=c3_creme)
{
	plot2D plot(0,Kmax,y0,y1, nx,ny);
	plot = C_back;
	for (int i = 0;i<G.n1/2;i++)
	{
	    plot.ptLine((G.k1(i)),real(G(i)),(G.k1(i+1)),real(G(i+1)),C_real);
        plot.ptLine((G.k1(i)),imag(G(i)),(G.k1(i+1)),imag(G(i+1)),C_imag);
    }
    plot.png(fname);
}







template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid1D(grid1D<complex<Ytype>,YtypeScalar,Xtype > & G, double y0, double y1,int nx, int ny, char * fname, color3 const & C_real=c3_blue,color3 const & C_imag=c3_red,color3 const & C_back=c3_creme)
{
	plot2D plot(G.a1,G.b1,y0,y1, nx,ny);
	plot= C_back;
	for (int i = 0;i<nx-1;i++)
	{
		 plot.ptLine((G.x1(i)),real(G(i)),(G.x1(i+1)),real(G(i+1)),C_real);
		 plot.ptLine((G.x1(i)),imag(G(i)),(G.x1(i+1)),imag(G(i+1)),C_imag);
    }
    plot.png(fname);
}


template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid1D(grid1D<complex<Ytype>,YtypeScalar,Xtype > & G, double x0, double x1, double y0, double y1,int nx, int ny, char * fname, color3 const & C_real=c3_blue,color3 const & C_imag=c3_red,color3 const & C_back=c3_creme)
{
	plot2D plot(x0,x1,y0,y1, nx,ny);
	plot= C_back;
	for (int i = 0;i<nx-1;i++)
	{
		 plot.ptLine((G.x1(i)),real(G(i)),(G.x1(i+1)),real(G(i+1)),C_real);
		 plot.ptLine((G.x1(i)),imag(G(i)),(G.x1(i+1)),imag(G(i+1)),C_imag);
    }
    plot.png(fname);
}

template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid1D_dots(grid1D<Ytype,YtypeScalar,Xtype > & G,  plot2D & plot, int size,  color3 & C )
{
	for (int i = 0;i<G.n1;i++)
		 plot.ptDot(G.x1(i),G(i),size,C);
}


#endif 






#ifdef GRID2D_H

template<class Ytype, class YtypeScalar,class Xtype >
void plotGrid2D_1(grid2D<Ytype,YtypeScalar,Xtype > & G,  const char * fname, ml_color (*cmap)(Ytype I) )
{	
	plot2D plot(G.a1,G.b1,G.a2,G.b2, G.n1,G.n2);
	
	for (int i = 0; i<plot.NX; i++)
	for (int j = 0; j<plot.NY; j++)
	{
		plot.set_px(i,j,cmap(G(i,j)));
	}
	
	plot.png(fname);
}

#endif













/*
double pltG2D_f(double x)
{
	return exp(-pow(2.0*x,8));
}

template<class Xtype, class Ytype, class YtypeScalar, class vecType, class arrayType >
void plotGrid2D_2(grid2D<Xtype,Ytype,YtypeScalar,vecType,arrayType > & G,  plot2D & plot )
{
	//for (int i = 0;i<G.nx-1;i++)
	//	 plot.ptLine(G.x_(i),G(i),G.x_(i+1),G(i+1),C);
	
	plot = c3_white;
	
	color3 background = c3_black;
	
	for (int i = 0; i<plot.NX; i++)
	for (int j = 0; j<plot.NY; j++)
	{
		double x = plot.to_x(i);
		double y = plot.to_y(j);
		
		if (x < G.xa) plot.set_px(i,j,background);
		else if (x >= G.xb) plot.set_px(i,j,background);
		else if (y < G.ya) plot.set_px(i,j,background);
		else if (y >= G.yb) plot.set_px(i,j,background);		
		else
		{
			int ig,jg;
			for (ig = 0; G.x0_(ig)<x; ig++ ) {}
			for (jg = 0; G.x1_(jg)<y; jg++ ) {}
			ig --;
			jg --;
			
			double a,b,c,d;
			
			a = norm( G.x_(ig,jg)-euVec<2,double>(x,y) )/norm(G.dx());
			b = norm( G.x_(ig,jg+1)-euVec<2,double>(x,y) )/norm(G.dx());
			c = norm( G.x_(ig+1,jg)-euVec<2,double>(x,y) )/norm(G.dx());
			d = norm( G.x_(ig+1,jg+1)-euVec<2,double>(x,y) )/norm(G.dx());
			
			color3 color = c3_white;
			
			color *= 
			pltG2D_f(a)*G(ig,jg)+
			pltG2D_f(b)*G(ig,jg+1)+
			pltG2D_f(c)*G(ig+1,jg)+
			pltG2D_f(d)*G(ig+1,jg+1);
					
			plot.set_px(i,j,color);
		}
	
		
	//	for (int k = 0; k<= 
		
	}	
}*/


