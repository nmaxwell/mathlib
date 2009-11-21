#ifndef STD_TOOLS_H
#define STD_TOOLS_H

/*
 * some commonly used things. 
 */

/*
    ifstream in;
    sprintf(fname,"pwd.txt");
    in.open(fname,ios::in);
    assert(in.good() && in.is_open());
    assert( in >> pwd );
    in.close();

    ofstream out;
    sprintf(fname,"%s/out/out.txt");
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    out.close();
*/

#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>
#include <stddef.h>
#include <sys/types.h>
#include <assert.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

#define filename_string_size 200
char * fname = new char [filename_string_size];
char * pwd = new char [filename_string_size];

#include <mathlib/machine_info.h>
#include "memory/ml_alloc.h"
#include "graphing/colors.h"
#include "../math/std_math.h"

#define endl1 cout << endl;
#define endl2 endl << endl
#define endl3 endl << endl << endl

#define _debug_here( pos ) ( printf( "debug: file %s; \n line %d; \t code %d\n", __FILE__ , __LINE__, pos ) );

#define div_up(a,b) ((a%b)?(a/b+1):(a/b))
    
int mkdir(const char * fname )
{
    return mkdir(fname, S_IRUSR | S_IREAD | S_IWUSR | S_IWRITE | S_IXUSR | S_IEXEC | S_IRWXU | S_IRGRP | S_IWGRP  | S_IXGRP | S_IRWXG | S_IROTH | S_IWOTH | S_IXOTH | S_IRWXO | S_ISUID | S_ISGID);
}

void clean_print(double x)
{    
    if (fabs(log10(fabs(x))) <= 2 ) // normal format
        if (x >= 0) printf(" %010.6f",x);
        else printf("%010.6f",x);
    else  // exponential format
        if (x >=0 ) printf(" %08.4E",x);
        else printf("%08.4E",x);
}

double get_real_time();

inline void swap(double & x,double & y)
{
    double z = x;
	x=y;
	y=z;
}

inline void swap(float & x,float & y)
{
    float z = x;
	x=y;
	y=z;
}

#ifdef ARPREC_MPREAL_H
inline void swap(mp_real & x,mp_real & y)
{
    mp_real z = x;
    
    x=y;
    y=z;
}
#endif

inline void swap(int & x, int & y)
{
    x ^= y;
    y ^= x;
    x ^= y;
    
   // int z = x;
   // x=y;
   // y=z;
}



void std_setup();

void std_exit();















#endif









