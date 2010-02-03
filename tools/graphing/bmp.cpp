#ifndef ML_BMP_CPP
#define ML_BMP_CPP


/*
 * this isn't mean to be good code by any means. dont use it.
 */

#include <mathlib/math/std_math.h>



struct bmp_data
{
	unsigned nx;
	unsigned ny;
	unsigned bpp;
	ubyte * data;
	
	bmp_data():data(0),nx(0),ny(0),bpp(0) {};
	
	~bmp_data() { if(data) ml_free(data); data=0; }
};


int read_bmp( const char * fname, bmp_data & bmp)
{
	if (bmp.data) ml_free(bmp.data);
	bmp.data = 0;
	
	ifstream in;
	in.open(fname, fstream::in);
	if (!in.good() || !in.is_open())
		return 1;
	
	char * temp = ml_alloc<char > (8);
	
	in.read(temp, 8);
	
	int size = *((unsigned * )(temp+2));
	
	ml_free(temp);
	temp = ml_alloc<char > (size);
	
	in.read(temp, size);
	
	bmp.nx = *((unsigned * )(temp+10));
	bmp.ny = *((unsigned * )(temp+14));
	bmp.bpp = *((unsigned * )(temp+20));
	unsigned offset = *((unsigned * )(temp+2));	
	
	bmp.data = ml_alloc<ubyte > ( bmp.nx*bmp.ny*(bmp.bpp/8) );
	
	if ( bmp.bpp == 8 and *((unsigned * )(temp+2)) == 1078 )
	{
		for (int i=0; i<bmp.nx; i++ )
		for (int j=0; j<bmp.ny; j++ )
		{
			int k = i*bmp.ny + j + offset-8;
			bmp.data[ i + j*bmp.nx ] = (int)*((ubyte *)(temp+  k));
		
	}
	else  return 2;
	
	
	ml_free(temp);
	return 0;
}


int read_bmp( const char * fname, bmp_data & bmp)
{
	if (bmp.data) ml_free(bmp.data);
	bmp.data = 0;
	
	ifstream in;
	in.open(fname, fstream::in);
	if (!in.good() || !in.is_open())
		return 1;
	
	char * temp = ml_alloc<char > (8);
	
	in.read(temp, 8);
	
	int size = *((unsigned * )(temp+2));
	
	ml_free(temp);
	temp = ml_alloc<char > (size);
	
	in.read(temp, size);
	
	bmp.nx = *((unsigned * )(temp+10));
	bmp.ny = *((unsigned * )(temp+14));
	bmp.bpp = *((unsigned * )(temp+20));
	unsigned offset = *((unsigned * )(temp+2));	
	
	bmp.data = ml_alloc<ubyte > ( bmp.nx*bmp.ny*(bmp.bpp/8) );
	
	if ( bmp.bpp == 8 and *((unsigned * )(temp+2)) == 1078 )
	{
		for (int i=0; i<bmp.nx; i++ )
		for (int j=0; j<bmp.ny; j++ )
		{
			int k = i*bmp.ny + j + offset-8;
			bmp.data[ i + j*bmp.nx ] = (int)*((ubyte *)(temp+  k));
		
	}
	else  return 2;
	
	
	ml_free(temp);
	return 0;
}






#endif
