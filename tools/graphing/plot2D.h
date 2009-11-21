#ifndef PLOT2D_H
#define PLOT2D_H



/*
 *  this is NOT meant to be good code in any way shape or form
 * 
 *  this is terrible terrible code.
 * 
 *  it's meant to be sloppy so I don't spend time writing a 
 *  good graphics library which is done already and would take months
 * 
 *  it's meant to output a lot of data quickly, without much hassle.
 * 
 *  use python + matplotlib for nice graphics.
 * 
 */

#include <mathlib/tools/std_tools.h>

#include "pngwriter.h"




//icpc .cpp -lrt  -I /workSpace/math `freetype-config --cflags` -I/usr/local/include  -L/usr/local/lib -lpng -lpngwriter -lz -lfreetype

// `freetype-config --cflags` 



void plot2D_stdFormat (double xy,bool xory,char * S)
{
	if ( xory )
		sprintf(S,"%2.2f",xy);
	else
		sprintf(S,"%2.2f",xy);
}



class plot2D
{
public:// members
	int NX,NY,NXY;
	ml_color * U;
	
public:// canonicals
	plot2D(double rhsXMIN,double rhsXMAX,double rhsYMIN,double rhsYMAX,unsigned rhsNX,unsigned rhsNY);
	plot2D();
	plot2D(const plot2D & rhs);
	void operator=(const plot2D & rhs);
	~plot2D() {	looseMem(); }
	
public:// operators
	ml_color & operator[](int);
	ml_color & operator()(int,int);
	void operator=(ml_color (*func)(double x,double y));
	
	void operator=(const ml_color & rhs);
	
	
public:// other funcs
	double to_x(int i);
	double to_y(int j);
	int to_i(double x);
	int to_j(double y);	
	
public://memory stuff
	void getMem();
	void looseMem();
	
public: // output
	void png(const char * fname);
	void debug();
	
public: // basic graphics

	void set_px(int,int,ml_color); // with protection
	
	void pxLine( int x1, int y1, int x2, int y2, ml_color const  & C);
	void ptLine(double x1,double y1,double x2,double y2, ml_color const & C);
	
	template<class D> void plotFunc_1(D (*F)(D x),ml_color & C,int n);
	
	void plotDot_1(double x,double y, double r, ml_color & C);
	void plotDot_2(double x,double y, int size, ml_color & C);
	
public:
	
	void text_1(int x, int y, const char * text, ml_color bkgndColor = ml_white, ml_color faceColor=ml_black, int height=12, int width=0,const char * font_fname="/usr/share/fonts/truetype/latex-xft-fonts/cmr10.ttf"); 

public: // math graphics
	
	double XMIN,XMAX,YMIN,YMAX;
	
	void axes_1( double x0=0, double y0=0, bool drawtext=0, ml_color bkgndColor = ml_white, ml_color faceColor=ml_black, void (*format)(double xy,bool xory,char * S) = plot2D_stdFormat, ml_color C=ml_black,double xMajorTicSpace=1.0, int NxMinorTics=5, double yMajorTicSpace=1.0, double NyMinorTics=5,double scale = 0.015,double text_scale=2.0);
	
	void plot_1(vector<double> & X, vector<double> & Y,  ml_color & C);
	
	//void Xaxis_labels_1( vector<double> majorTicks, vector<string> labels);
	
	
	void ptLine2(double x1,double y1,double x2,double y2,int n,ml_color & C);
	
	void pxDot(int i,int j, int size, ml_color const & C);
	void ptDot(double x1,double y1, int size, ml_color & C);
	void pxBresLine1(int x0, int y0, int x1, int y1,  ml_color const & C);
	void pxLineThick(int x0, int y0, int x1, int y1, int size, ml_color const & C);
	void ptLineThick(double x1,double y1,double x2,double y2, int thick,ml_color const  & C);
	
	int text_2(int x, int y, const unsigned int * text, ml_color bkgndColor, ml_color faceColor, int height, int width, const char * font_fname);
	
	void hline(double y,ml_color & c);
	void vline(double x,ml_color & c);
	
	template< class X, class Y >
	void plot( functor<X,Y > const & f, ml_color & c );
	
	template< class plotType >
	void plot( plotType const & f, ml_color & c );
	
	template< class plotType >
	void plot_abs( plotType const & f, ml_color & c );
	
	template< class plotType >
	void plot_real( plotType const & f, ml_color & c );
	
	template< class plotType >
	void plot_imag( plotType const & f, ml_color & c );
	
	template< class Xtype, class Ytype >
	void plot( Xtype * x, Xtype * y, int n, ml_color & c );
};













#endif
