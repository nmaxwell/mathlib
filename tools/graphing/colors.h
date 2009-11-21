#ifndef MATHLIB_COLORS_H
#define MATHLIB_COLORS_H

/*
 * this is pretty straight-forward...RGB, RGBA
 * 
 */


struct ml_color
{
public:
	ubyte red;
	ubyte green;
	ubyte blue;
	
	ml_color(float R,float G,float B) // [0,1]_float -> [0,255]_int
	{
		if (R<0) R = 0; if (R>1.0) R = 1.0;
		if (G<0) G = 0; if (G>1.0) G = 1.0;
		if (B<0) B = 0; if (B>1.0) B = 1.0;	
		red = ubyte(R*255);
		green = ubyte(G*255);
		blue = ubyte(B*255);
	}
	
	ml_color(ubyte R,ubyte G,ubyte B) 
	:red( R ), green( G ), blue( B )	{}
	
public:
	void operator -= (ml_color & rhs)
	{
		red -= rhs.red;
		green -= rhs.green;
		blue -= rhs.blue;
	}
	void operator /= (ml_color & rhs)
	{
		red /= rhs.red;
		green /= rhs.green;
		blue /= rhs.blue;
	}
	
	void operator *= (double rhs)
	{
		red *= rhs;
		green *= rhs;
		blue *= rhs;
	}	
	
public:
	ml_color() {};
	~ml_color() {};
	void operator= (const ml_color & rhs) {red = rhs.red; green = rhs.green; blue = rhs.blue;}
	ml_color(const ml_color & rhs) {operator=(rhs);}
};





struct ml_color4
{
public:
	ubyte red;
	ubyte green;
	ubyte blue;
	ubyte alpha;	
	
	ml_color4(float R,float G,float B,float A) // [0,1]_float -> [0,255]_int
	{
		if (R<0) R = 0; if (R>1.0) R = 1.0;
		if (G<0) G = 0; if (G>1.0) G = 1.0;
		if (B<0) B = 0; if (B>1.0) B = 1.0;	
		if (A<0) A = 0; if (A>1.0) A = 1.0;
		red = ubyte(R*255);
		green = ubyte(G*255);
		blue = ubyte(B*255);
		alpha = ubyte(A*255); 
	}
	
	ml_color4(ubyte R,ubyte G,ubyte B,ubyte A) 
	:red( R ), green( G ), blue( B ), alpha( A )	{}
	
public:
	ml_color4() {};
	~ml_color4() {};
	void operator= (const ml_color4 & rhs) {red = rhs.red; green = rhs.green; blue = rhs.blue; alpha = rhs.alpha;}
	ml_color4(const ml_color4 & rhs) {operator=(rhs);}
};



ml_color ml_smoothblue((ubyte)39,(ubyte)108,(ubyte)166);

ml_color ml_black(0.f,0.f,0.f);
ml_color ml_red(1.f,0.f,0.f);
ml_color ml_green(0.f,1.f,0.f);
ml_color ml_blue(0.f,0.f,1.f);
ml_color ml_white(1.f,1.f,1.f);
ml_color ml_yellow(1.f,1.f,0.f);
ml_color ml_magenta(1.f,0.f,1.f);
ml_color ml_cyan(0.f,1.f,1.f);
ml_color ml_grey1(0.1f,0.1f,0.1f);
ml_color ml_grey2(0.2f,0.2f,0.2f);
ml_color ml_grey3(0.3f,0.3f,0.3f);
ml_color ml_grey4(0.4f,0.4f,0.4f);
ml_color ml_grey5(0.5f,0.5f,0.5f);
ml_color ml_grey6(0.6f,0.6f,0.6f);
ml_color ml_grey7(0.7f,0.7f,0.7f);
ml_color ml_grey8(0.8f,0.8f,0.8f);
ml_color ml_grey9(0.9f,0.9f,0.9f);

ml_color ml_orange(1.f,0.5f,0.f);

ml_color ml_creme((ubyte)238,(ubyte)220,(ubyte)130);







#endif
