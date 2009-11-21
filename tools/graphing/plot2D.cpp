#ifndef PLOT2D_CPP
#define PLOT2D_CPP

#ifndef PLOT2D_H
	#include "plot2D.h"
#endif






//#include "freetype2/freetype.h"
//#include "freetype/ft2build.h"
 



//icpc .cpp -lrt  -I /workSpace/math `freetype-config --cflags` -I/usr/local/include  -L/usr/local/lib -lpng -lpngwriter -lz -lfreetype

// `freetype-config --cflags` 




plot2D::plot2D(double rhsXMIN,double rhsXMAX,double rhsYMIN,double rhsYMAX,unsigned rhsNX,unsigned rhsNY)
:XMIN(rhsXMIN),XMAX(rhsXMAX),YMIN(rhsYMIN),YMAX(rhsYMAX),NX(rhsNX),NY(rhsNY),U(0)
{
	NXY = NX*NY;
	getMem();
}


plot2D::plot2D():XMIN(0),XMAX(0),YMIN(0),YMAX(0),NX(0),NY(0),U(0)
{
}

plot2D::plot2D(const plot2D & rhs):XMIN(0),XMAX(0),YMIN(0),YMAX(0),NX(0),NY(0),U(0)
{
	this->operator=(rhs);
}

void plot2D::operator=(const plot2D & rhs)
{	
	this->looseMem();
	XMIN = rhs.XMIN;
	XMAX = rhs.XMAX;
	YMIN = rhs.YMIN;
	YMAX = rhs.YMAX;
	NX = rhs.NX;
	NY = rhs.NY;
	NXY = rhs.NXY;
	
	if(rhs.U) {
		this->getMem();		
		for (unsigned i = 0;i<NXY;i++)
			(*this)[i] = (const_cast<plot2D &>(rhs))[i];		
	}
}

void plot2D::getMem()
{
	assert(NXY && !U);
	U = new ml_color [NXY];
	assert(U);
}

void plot2D::looseMem()
{
	if (U!=0) {
	delete [] U;
	U = 0; }
	assert(U==0);
}

ml_color & plot2D::operator[](int i)
{
	return U[i];
}

ml_color & plot2D::operator()(int X,int Y)
{
	return U[Y*NX+X];
}

void plot2D::operator=(ml_color (*func)(double x,double y))
{
	for (int i=0;i<NX;i++)
	for (int j=0;j<NY;j++)
	{
		(*this)(i,j) = func(to_x(i),to_y(j));
	}
}

void plot2D::operator=(const ml_color & rhs)
{
	for (int i=0;i<NX;i++)
	for (int j=0;j<NY;j++)
	{
		(*this)(i,j) = rhs;
	}
}

double plot2D::to_x(int i)
{
	return ((XMAX-XMIN)*(double)i/(double)(NX) + XMIN);
}

double plot2D::to_y(int j)
{
	return ((YMAX-YMIN)*(double)j/(double)(NY) + YMIN);
}

int plot2D::to_i(double x)
{
	return (x-XMIN)/(XMAX-XMIN)*NX;	
}

int  plot2D::to_j(double y)
{
	return (y-YMIN)/(YMAX-YMIN)*NY;
}






void plot2D::png(const char * fname)
{
	pngwriter png(NX,NY,1.0,fname);
	
 	for (int i = 0; i<NX;i++)
 	for (int j = 0; j<NY;j++)
 	{
 		png.plot(i+1,j+1, operator()(i,j).red*256,operator()(i,j).green*256,operator()(i,j).blue*256);
 	}
	
	png.close();
}

void plot2D::pxLine( int x1, int y1, int x2, int y2, ml_color const  & C)
{
	pxBresLine1(x1, y1, x2, y2, C);
}

void plot2D::set_px(int i,int j,ml_color C)
{
	if (i >= NX) i = NX-1;
	if (i < 0) i = 0;
	if (j >= NY) j = NY-1;
	if (j < 0) j = 0;
	operator() (i,j) = C;
}


void plot2D::ptLine(double x1,double y1,double x2,double y2,ml_color const & C)
{
	pxLine(	
	((x1-XMIN)/(XMAX-XMIN))*NX,
	((y1-YMIN)/(YMAX-YMIN))*NY,
	((x2-XMIN)/(XMAX-XMIN))*NX,
	((y2-YMIN)/(YMAX-YMIN))*NY,
	C);
}

void plot2D::ptLineThick(double x1,double y1,double x2,double y2, int thick,ml_color const & C)
{
	pxLineThick(	
	((x1-XMIN)/(XMAX-XMIN))*NX,
	((y1-YMIN)/(YMAX-YMIN))*NY,
	((x2-XMIN)/(XMAX-XMIN))*NX,
	((y2-YMIN)/(YMAX-YMIN))*NY,
	thick,
	C);
}



void plot2D::axes_1(double x0, double y0, bool drawtex,ml_color bkgndColor, ml_color faceColor, void (*format)(double xy,bool xory,char * S), ml_color C, double xMajorTicSpace, int NxMinorTics, double yMajorTicSpace, double NyMinorTics, double scale, double text_scale )
{
	ptLine(XMIN,y0,XMAX,y0,C);
	ptLine(x0,YMIN,x0,YMAX,C);
	double xMinorTicSpace = xMajorTicSpace/(double)NxMinorTics;
	double yMinorTicSpace = yMajorTicSpace/(double)NyMinorTics;	
	double sizey = scale*(YMAX-YMIN);
	double sizex = scale*(XMAX-XMIN);
	int pxsizex = scale*NX*text_scale;
	int pxsizey = scale*NY*text_scale;
	double x = 0,y=0;
	char * S = new char [30];
	
	y = y0;
	x = floor(XMIN/xMinorTicSpace)*xMinorTicSpace;
	while (y0 >= YMIN && y0 <= YMAX && x < XMAX)
	{		
		ptLine(x,y0-sizey/2.0,x,y0,C);
		x +=  xMinorTicSpace;
	}
	x = x0;
	y = floor(YMIN/yMinorTicSpace)*yMinorTicSpace;
	while (x0 >= XMIN && x0 <= XMAX && y < YMAX)
	{		
		ptLine(x0-sizex/2.0,y,x0,y,C);
		y +=  yMinorTicSpace;
	}
	
	y = y0;
	x = floor(XMIN/xMajorTicSpace)*xMajorTicSpace;
	while (y0 >= YMIN && y0 <= YMAX && x < XMAX)
	{		
		ptLine(x,y0-sizey,x,y0+sizey,C);		
		format(x,true,S);
		if (fabs(x -x0)>1E-14 && drawtex) 
			text_1(to_i(x), to_j(y-sizey)-pxsizey,S, bkgndColor, faceColor, pxsizey, 0);
		x +=  xMajorTicSpace;
	}
	x = x0;
	y = ceil(YMIN/yMajorTicSpace)*yMajorTicSpace;
	while (x0 >= XMIN && x0 <= XMAX && y < YMAX)
	{		
		ptLine(x0-sizex,y,x0+sizex,y,C);
		format(y,false,S);
		if (fabs(y -y0)>1E-14 && drawtex)
			text_1(to_i(x+sizex)+pxsizex/2, to_j(y)-pxsizey/2,S, bkgndColor, faceColor, pxsizey, 0);		
		y +=  yMajorTicSpace;
	}
}




template<class D> void plot2D::plotFunc_1(D (*F)(D x),ml_color & C, int n)
{
	D x1,x2;
	int i;
	for (i = 0;i<n;i++)
	{
		x1 = ((XMAX-XMIN)*(double)i/(double)(	n) + XMIN);
		x2 = ((XMAX-XMIN)*(double)(i+1)/(double)(	n) + XMIN);
		ptLine(x1,F(x1),x2,F(x2),C);
	}
}

void plot2D::plotDot_1(double x,double y, double r, ml_color & C)
{
	if (x >= XMIN && x <= XMAX && y >= YMIN && y <= YMAX)
	for (int i = to_i(x-r); i<= to_i(x+r); i++)
	for (int j = to_j(y-r); j<= to_i(y+r); j++)
		if (  (to_x(i)-x)*(to_x(i)-x) + (to_y(j)-y)*(to_y(j)-y)  <= r*r )
			set_px(i,j,C);
}

void plot2D::plotDot_2(double x,double y, int size, ml_color & C)
{
	if (x >= XMIN && x <= XMAX && y >= YMIN && y <= YMAX)
{	
	switch(size)
	{
	case 1:
		(*this).set_px(to_i(x),to_j(y),C);
		break;
	case 2:
		(*this).set_px(to_i(x),to_j(y),C);
		(*this).set_px(to_i(x)+1 ,to_j(y)+0 ,C);
		(*this).set_px(to_i(x)-1 ,to_j(y)+0 ,C);
		(*this).set_px(to_i(x)+0 ,to_j(y)+1 ,C);
		(*this).set_px(to_i(x)+0 ,to_j(y)-1 ,C);
		break;
	case 3:
		(*this).plotDot_2(x,y,2,C);
		(*this).set_px(to_i(x)+1 ,to_j(y)+1 ,C);
		(*this).set_px(to_i(x)-1 ,to_j(y)-1 ,C);
		(*this).set_px(to_i(x)+1 ,to_j(y)-1 ,C);
		(*this).set_px(to_i(x)-1 ,to_j(y)+1 ,C);
		break;
	case 4:
		(*this).plotDot_2(x,y,3,C);
		(*this).set_px(to_i(x)+0 ,to_j(y)+2 ,C);
		(*this).set_px(to_i(x)-0 ,to_j(y)-2 ,C);
		(*this).set_px(to_i(x)+2 ,to_j(y) ,C);
		(*this).set_px(to_i(x)-2 ,to_j(y) ,C);
		break;	
	case 5:
		(*this).plotDot_2(x,y,4,C);
		(*this).set_px(to_i(x)+1 ,to_j(y)+2 ,C);
		(*this).set_px(to_i(x)+2 ,to_j(y)+1 ,C);
		(*this).set_px(to_i(x)+1 ,to_j(y)-2 ,C);
		(*this).set_px(to_i(x)+2 ,to_j(y)-1 ,C);
		(*this).set_px(to_i(x)-1 ,to_j(y)+2 ,C);
		(*this).set_px(to_i(x)-2 ,to_j(y)+1 ,C);
		(*this).set_px(to_i(x)-1 ,to_j(y)-2 ,C);
		(*this).set_px(to_i(x)-2 ,to_j(y)-1 ,C);
	break;			
	default:
		//size not done yet
		break;
	}
	
}
}


void plot2D::pxDot(int i,int j, int size, ml_color const & C)
{
	if (i >= 0 && i < NX && j >= 0 && j < NY)
{
	
	switch(size)
	{
	case 1:
		(*this).set_px(i,j,C);
		break;
	case 2:
		(*this).set_px(i,j,C);
		(*this).set_px(i+1 ,j+0 ,C);
		(*this).set_px(i-1 ,j+0 ,C);
		(*this).set_px(i+0 ,j+1 ,C);
		(*this).set_px(i+0 ,j-1 ,C);
		break;
	case 3:
		(*this).pxDot(i,j,2,C);
		(*this).set_px(i+1 ,j+1 ,C);
		(*this).set_px(i-1 ,j-1 ,C);
		(*this).set_px(i+1 ,j-1 ,C);
		(*this).set_px(i-1 ,j+1 ,C);
		break;
	case 4:
		(*this).pxDot(i,j,3,C);
		(*this).set_px(i+0 ,j+2 ,C);
		(*this).set_px(i-0 ,j-2 ,C);
		(*this).set_px(i+2 ,j ,C);
		(*this).set_px(i-2 ,j ,C);
		break;	
	case 5:
		(*this).pxDot(i,j,4,C);
		(*this).set_px(i+1 ,j+2 ,C);
		(*this).set_px(i+2 ,j+1 ,C);
		(*this).set_px(i+1 ,j-2 ,C);
		(*this).set_px(i+2 ,j-1 ,C);
		(*this).set_px(i-1 ,j+2 ,C);
		(*this).set_px(i-2 ,j+1 ,C);
		(*this).set_px(i-1 ,j-2 ,C);
		(*this).set_px(i-2 ,j-1 ,C);
	break;			
	default:
		//size not done yet
		assert(0);
		break;
	}
	
}
}

void plot2D::text_1(int x, int y, const char * text, ml_color bkgndColor, ml_color faceColor, int height, int width, const char * font_fname)
{	
	float c0r = (float)bkgndColor.red/255.0;
	float c0g = (float)bkgndColor.green/255.0;
	float c0b = (float)bkgndColor.blue/255.0;

	float c1r = (float)faceColor.red/255.0;
	float c1g = (float)faceColor.green/255.0;
	float c1b = (float)faceColor.blue/255.0;
	
	FT_Library library;
	FT_Face face;
	int error = FT_Init_FreeType( &library ); 
	if ( error )
		cout << "free type error: " << error << endl;
	error = FT_New_Face( library,  font_fname, 0, &face );
	if ( error == FT_Err_Unknown_File_Format )
		cout << "free type error: " << error << "; unrecognized font.\n";
	else if ( error )
		cout << "free type error: " << error << "; probably could'nt open font.\n";	
	error = FT_Set_Pixel_Sizes ( face, width,height );
	if ( error )
		cout << "free type error: FT_Set_Char_Size: " << error << endl;

	for (int n = 0;text[n] != 0;n++)
	{
		int glyph_index = FT_Get_Char_Index( face,  text[n] );
		if (!glyph_index) cout << "undefined character code" << endl;
		
		error = FT_Load_Glyph( face, glyph_index, 0 );
		if ( error )
			cout << "free type error: FT_Load_Glyph: " << error << endl;
		if ( face->glyph->format != FT_GLYPH_FORMAT_BITMAP ) {
			error = FT_Render_Glyph( face->glyph, FT_RENDER_MODE_NORMAL );
			if ( error )
				cout << "free type error: FFT_Render_Glyph: " << error << endl;
		}
		FT_Bitmap * bmp =&(face->glyph->bitmap);
		if (
			bmp->pixel_mode == FT_PIXEL_MODE_GRAY  &&
			bmp->num_grays == 256
		 )
		{
			ubyte g = 0; double s;
			 for (int i = 0; i<bmp->rows;i++ )
			 for (int j = 0; j<bmp->width;j++ )
			 {
			 	g = bmp->buffer[ i*bmp->pitch + j ];
			 	s = (float)g/255.0;
				ml_color G(
				(float)s*(c1r-c0r)+c0r,
				(float)s*(c1g-c0g)+c0g,
				(float)s*(c1b-c0b)+c0b );
				set_px(x+j,y+bmp->rows-i-1,G);
			 }
		}
		x += bmp->width;
	}
}


int plot2D::text_2(int x, int y, const unsigned int * text, ml_color bkgndColor, ml_color faceColor, int height, int width, const char * font_fname)
{	
	float c0r = (float)bkgndColor.red/255.0;
	float c0g = (float)bkgndColor.green/255.0;
	float c0b = (float)bkgndColor.blue/255.0;

	float c1r = (float)faceColor.red/255.0;
	float c1g = (float)faceColor.green/255.0;
	float c1b = (float)faceColor.blue/255.0;
	
	FT_Library library;
	FT_Face face;
	int error = FT_Init_FreeType( &library ); 
	if ( error )
		cout << "free type error: " << error << endl;
	error = FT_New_Face( library,  font_fname, 0, &face );
	if ( error == FT_Err_Unknown_File_Format )
		cout << "free type error: " << error << "; unrecognized font.\n";
	else if ( error )
		cout << "free type error: " << error << "; probably could'nt open font.\n";	
	error = FT_Set_Pixel_Sizes ( face, width,height );
	if ( error )
		cout << "free type error: FT_Set_Char_Size: " << error << endl;

	for (int n = 0;text[n] != 0;n++)
	{
		int glyph_index = FT_Get_Char_Index( face,  text[n] );
		//if (!glyph_index) cout << "undefined character code" << endl;
		if (glyph_index)
		{
			error = FT_Load_Glyph( face, glyph_index, 0 );
			if ( error )
				cout << "free type error: FT_Load_Glyph: " << error << endl;
			if ( face->glyph->format != FT_GLYPH_FORMAT_BITMAP ) {
				error = FT_Render_Glyph( face->glyph, FT_RENDER_MODE_NORMAL );
				if ( error )
					cout << "free type error: FFT_Render_Glyph: " << error << endl;
			}
			FT_Bitmap * bmp =&(face->glyph->bitmap);
			if (
				bmp->pixel_mode == FT_PIXEL_MODE_GRAY  &&
				bmp->num_grays == 256
			 )
			{
				ubyte g = 0; double s;
				 for (int i = 0; i<bmp->rows;i++ )
				 for (int j = 0; j<bmp->width;j++ )
				 {
				 	g = bmp->buffer[ i*bmp->pitch + j ];
				 	s = (float)g/255.0;
					ml_color G(
					(float)s*(c1r-c0r)+c0r,
					(float)s*(c1g-c0g)+c0g,
					(float)s*(c1b-c0b)+c0b );
					set_px(x+j,y+bmp->rows-i-1,G);
				 }
			}
		x += bmp->width;
		}
	}
	return x;
}






















void plot2D::plot_1(vector<double> & X, vector<double> & Y,  ml_color & C)
{
	for( int i = 0;i<X.size()-1;i++)
		ptLine(X[i],Y[i],X[i+1],Y[i+1],C);
}


void plot2D::pxBresLine1(int x0, int y0, int x1, int y1,  ml_color const & C)
{	
	if (x0 <= 0 && x1 <= 0) {x0 = 0; x1 = 0; }
	if (x0 >= NX && x1 >= NX) {x0 = NX-1; x1 = NX-1; }
	if (y0 <= 0 && y1 <= 0) {y0 = 0; y1 = 0; }
	if (y0 >= NY && y1 >= NY) {y0 = NY-1; y1 = NY-1; }
	
	if (x1 >= NX)
	{
		if ( x1!=x0 )
			y1 = (float)(y1-y0)*(float)(NX-x0)/(float)(x1-x0)+y0;
		x1 = NX-1;
	}

	if (x1 <= 0)
	{
		if ( x1!=x0 )
			y1 = (float)(y1-y0)*(float)(-x0)/(float)(x1-x0)+y0;
		x1 = 0;
	}		
	
	if (x0 >= NX)
	{
		if ( x1!=x0 )
			y0 = (float)(y1-y0)*(float)(NX-x0)/(float)(x1-x0)+y0;
		x0 = NX-1;
	}

	if (x0 <= 0)
	{
		if ( x1!=x0 )
			y0 = (float)(y1-y0)*(float)(-x0)/(float)(x1-x0)+y0;
		x0 = 0;
	}	
	
	if (y1 >= NY)
	{
		if ( y1!=y0 )
			x1 = (float)(x1-x0)*(float)(NY-y0)/(float)(y1-y0)+x0;
		y1 = NY-1;
	}

	if (y1 <= 0)
	{
		if ( y1!=y0 )
			x1 = (float)(x1-x0)*(float)(-y0)/(float)(y1-y0)+x0;
		y1 = 0;
	}		
	
	if (y0 >= NY)
	{
		if ( y1!=y0 )
			x0 = (float)(x1-x0)*(float)(NY-y0)/(float)(y1-y0)+x0;
		y0 = NY-1;
	}

	if (y0 <= 0)
	{
		if ( y1!=y0 )
			x0 = (float)(x1-x0)*(float)(-y0)/(float)(y1-y0)+x0;
		y0 = 0;
	}
	

	if (x1 >= NX)
	{
		if ( x1!=x0 )
			y1 = (float)(y1-y0)*(float)(NX-x0)/(float)(x1-x0)+y0;
		x1 = NX-1;
	}

	if (x1 <= 0)
	{
		if ( x1!=x0 )
			y1 = (float)(y1-y0)*(float)(-x0)/(float)(x1-x0)+y0;
		x1 = 0;
	}		
	
	if (x0 >= NX)
	{
		if ( x1!=x0 )
			y0 = (float)(y1-y0)*(float)(NX-x0)/(float)(x1-x0)+y0;
		x0 = NX-1;
	}

	if (x0 <= 0)
	{
		if ( x1!=x0 )
			y0 = (float)(y1-y0)*(float)(-x0)/(float)(x1-x0)+y0;
		x0 = 0;
	}	
	
	if (y1 >= NY)
	{
		if ( y1!=y0 )
			x1 = (float)(x1-x0)*(float)(NY-y0)/(float)(y1-y0)+x0;
		y1 = NY-1;
	}

	if (y1 <= 0)
	{
		if ( y1!=y0 )
			x1 = (float)(x1-x0)*(float)(-y0)/(float)(y1-y0)+x0;
		y1 = 0;
	}		
	
	if (y0 >= NY)
	{
		if ( y1!=y0 )
			x0 = (float)(x1-x0)*(float)(NY-y0)/(float)(y1-y0)+x0;
		y0 = NY-1;
	}

	if (y0 <= 0)
	{
		if ( y1!=y0 )
			x0 = (float)(x1-x0)*(float)(-y0)/(float)(y1-y0)+x0;
		y0 = 0;
	}

	if (x0 < 0) cout << "bounds error in pxBresLine1 1\n";
	if (x1 < 0) cout << "bounds error in pxBresLine1 2\n";
	if (y0 < 0) cout << "bounds error in pxBresLine1 3\n";
	if (y1 < 0) cout << "bounds error in pxBresLine1 4\n";
	if (x0 >=NX) cout << "bounds error in pxBresLine1 5\n";
	if (x1 >=NX) cout << "bounds error in pxBresLine1 6\n";
	if (y0 >=NY) cout << "bounds error in pxBresLine1 7\n";
	if (y1 >=NY) cout << "bounds error in pxBresLine1 8\n";

       int Dx = x1 - x0; 
       int Dy = y1 - y0;
       int steep = (abs(Dy) >= abs(Dx));
       if (steep) {
	   swap(x0, y0);
	   swap(x1, y1);
	   // recompute Dx, Dy after swap
	   Dx = x1 - x0;
	   Dy = y1 - y0;
       }
       int xstep = 1;
       if (Dx < 0) {
	   xstep = -1;
	   Dx = -Dx;
       }
       int ystep = 1;
       if (Dy < 0) {
	   ystep = -1;		
	   Dy = -Dy; 
       }
       int TwoDy = 2*Dy; 
       int TwoDyTwoDx = TwoDy - 2*Dx; // 2*Dy - 2*Dx
       int E = TwoDy - Dx; //2*Dy - Dx
       int y = y0;
       int xDraw, yDraw;	
       for (int x = x0; x != x1; x += xstep) {		
	   if (steep) {			
	       xDraw = y;
	       yDraw = x;
	   } else {			
	       xDraw = x;
	       yDraw = y;
	   }
		    set_px(xDraw, yDraw,C);
	   
	   if (E > 0) {
	       E += TwoDyTwoDx; //E += 2*Dy - 2*Dx;
	       y = y + ystep;
	   } else {
	       E += TwoDy; //E += 2*Dy;
	   }
       }
}




void plot2D::pxLineThick(int x0, int y0, int x1, int y1,int size,  ml_color const & C)
{
   int Dx = x1 - x0; 
   int Dy = y1 - y0;
   int steep = (abs(Dy) >= abs(Dx));
   if (steep) {
       swap(x0, y0);
       swap(x1, y1);
       // recompute Dx, Dy after swap
       Dx = x1 - x0;
       Dy = y1 - y0;
   }
   int xstep = 1;
   if (Dx < 0) {
       xstep = -1;
       Dx = -Dx;
   }
   int ystep = 1;
   if (Dy < 0) {
       ystep = -1;		
       Dy = -Dy; 
   }
   int TwoDy = 2*Dy; 
   int TwoDyTwoDx = TwoDy - 2*Dx; // 2*Dy - 2*Dx
   int E = TwoDy - Dx; //2*Dy - Dx
   int y = y0;
   int xDraw, yDraw;	
   for (int x = x0; x != x1; x += xstep) {		
       if (steep) {			
           xDraw = y;
           yDraw = x;
       } else {			
           xDraw = x;
           yDraw = y;
       }
       pxDot(xDraw, yDraw,size,  C);
       
       if (E > 0) {
           E += TwoDyTwoDx; //E += 2*Dy - 2*Dx;
           y = y + ystep;
       } else {
           E += TwoDy; //E += 2*Dy;
       }
   }
}

void plot2D::ptDot(double x,double y, int size, ml_color & C)
{
	pxDot(to_i(x),to_j(y), size, C);
}

	
	
void plot2D::hline(double y,ml_color & c)
{
	ptLine(XMIN,y,XMAX,y,c);
}

void plot2D::vline(double x,ml_color & c)
{
	ptLine(x,YMIN,x,YMAX,c);
}



template< class X, class Y >
void plot2D::plot( functor<X,Y > const & f, ml_color & c )
{
    for (int i = 0; i<NX-1; i++)
        ptLine(to_x(i),f(to_x(i)),to_x(i+1),f(to_x(i+1)),c);
}


template< class plotType >
void plot2D::plot( plotType const & f, ml_color & c )
{
    for (int i = 0; i<NX-1; i++)
        ptLine(
            to_x(i),
            ( const_cast<plotType & > (f)(to_x(i)) ),
            to_x(i+1),
            ( const_cast<plotType & > (f)(to_x(i+1)) ),
            c);
}

template< class plotType >
void plot2D::plot_abs( plotType const & f, ml_color & c )
{
    for (int i = 0; i<NX-1; i++)
        ptLine(
            to_x(i),
            abs( const_cast<plotType & > (f)(to_x(i)) ),
            to_x(i+1),
            abs( const_cast<plotType & > (f)(to_x(i+1)) ),
            c);
}

template< class plotType >
void plot2D::plot_real( plotType const & f, ml_color & c )
{
    for (int i = 0; i<NX-1; i++)
        ptLine(
            to_x(i),
            real( const_cast<plotType & > (f)(to_x(i)) ),
            to_x(i+1),
            real( const_cast<plotType & > (f)(to_x(i+1)) ),
            c);
}

template< class plotType >
void plot2D::plot_imag( plotType const & f, ml_color & c )
{
    for (int i = 0; i<NX-1; i++)
        ptLine(
            to_x(i),
            imag( const_cast<plotType & > (f)(to_x(i)) ),
            to_x(i+1),
            imag( const_cast<plotType & > (f)(to_x(i+1)) ),
            c);
}



template< class Xtype, class Ytype >
void plot2D::plot( Xtype * x, Xtype * y, int n, ml_color & c )
{
    for (int i = 1; i<n; i++)
        ptLine(
            x[i-1],
            y[i-1],
            x[i],
            y[i],
            c);
}






void plot2D::debug()
{
	cout << "NX: " << NX << endl;
	cout << "NY: " << NY << endl;
	cout << "XMIN: " << XMIN << endl;
	cout << "XMAX: " << XMAX << endl;
	cout << "YMIN: " << YMIN << endl;
	cout << "YMAX: " << YMAX << endl;
}

















#endif
