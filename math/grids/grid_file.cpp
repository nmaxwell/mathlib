#ifndef GRID_FILE_OPS_CPP
#define GRID_FILE_OPS_CPP

#ifndef GRID_FILE_OPS_H
	#include "grid_file.h"
#endif







#ifdef GRID1D_H

template<classgridPars_ >
int writeFile( grid1D<gridPars_ > & G, const char * fname )
{
	int N_data = G.n1*sizeof(Ytype );	
	int N_header = 2*sizeof(Xtype )+1*sizeof(int );
		
	static char * H = (char *) new char [N_header];
		
	ofstream out;
	out.open(fname, fstream::out);
	if (!out.good() || !out.is_open())
		return 1;
	
	((int *)H)[0] = G.n1;
	((Xtype *)&(H[4]))[0] = G.a1;
	((Xtype *)&(H[4]))[1] = G.b1;
	
	out.write(H,N_header);
	out.write((char *)(G.array),N_data);
	out.close();
	
	return 0;
}


template<classgridPars_ >
int readFile( grid1D<gridPars_ > & G, const char * fname )
{
	int N_header = 2*sizeof(Xtype )+1*sizeof(int );
	
	static char * H = (char *) new char [N_header];
	
	ifstream in;
	in.open(fname, fstream::in);
	if (!in.good() || !in.is_open())
		return 1;
	
	in.read(H,N_header);
	
	int n1;
	
	n1 = ((int *)H)[0];
	G.a1 = ((Xtype *)&(H[4]))[0];
	G.b1 = ((Xtype *)&(H[4]))[1];
	
	if ( !G.array || n1 != G.n1 )
	{
		if (G.array) delete [] G.array;
		G.n1 = n1;
		G.array = new Ytype [n1];
	}
	
	int N_data = G.n1*sizeof(Ytype );
	
	in.read((char *)G.array,N_data);
	in.close();
	
	return 0;
}

#endif 



#ifdef GRID2D_H

template<classgridPars_ >
int writeFile( grid2D<gridPars_ > & G, const char * fname )
{
	int N_data = G.n1*G.n2*sizeof(Ytype );	
	int N_header = 4*sizeof(Xtype )+2*sizeof(int );
		
	static char * H = (char *) new char [N_header];
    
	ofstream out;
	out.open(fname, fstream::out);
	if (!out.good() || !out.is_open())
		return 1;
	
	((int *)H)[0] = G.n1;
	((int *)H)[1] = G.n2;
	((Xtype *)&(H[8]))[0] = G.a1;
	((Xtype *)&(H[8]))[1] = G.b1;
	((Xtype *)&(H[8]))[2] = G.a2;
	((Xtype *)&(H[8]))[3] = G.b2;
	
	out.write(H,N_header);
	out.write((char *)(G.array),N_data);
	out.close();
	
	return 0;
}


template<classgridPars_ >
int readFile( grid2D<gridPars_ > & G, const char * fname )
{
	int N_header = 4*sizeof(Xtype )+2*sizeof(int );
	
	static char * H = (char *) new char [N_header];
	
	ifstream in;
	in.open(fname, fstream::in);
	if (!in.good() || !in.is_open())
		return 1;
	
	in.read(H,N_header);
	
	int n1,n2;
	
	n1 = ((int *)H)[0];
	n2 = ((int *)H)[1];
	G.a1 = ((Xtype *)&(H[8]))[0];
	G.b1 = ((Xtype *)&(H[8]))[1];
	G.a2 = ((Xtype *)&(H[8]))[2];
	G.b2 = ((Xtype *)&(H[8]))[3];
	
	if ( !G.array || n1 != G.n1 || n2 != G.n2 )
	{
		if (G.array) delete [] G.array;
		G.n1 = n1;
		G.n2 = n2;
		G.array = new Ytype [n1*n2];
	}
	
	int N_data = G.n1*G.n2*sizeof(Ytype );
	
	in.read((char *)G.array,N_data);
	in.close();
	
	return 0;
}

#endif 


#ifdef GRID3D_H

template<classgridPars_ >
int writeFile( grid3D<gridPars_ > & G, const char * fname )
{
	int N_data = G.n1*G.n2*G.n3*sizeof(Ytype );	
	int N_header = 6*sizeof(Xtype )+3*sizeof(int );
		
	static char * H = (char *) new char [N_header];
		
	ofstream out;
	out.open(fname, fstream::out);
	if (!out.good() || !out.is_open())
		return 1;
	
	((int *)H)[0] = G.n1;
	((int *)H)[1] = G.n2;
	((int *)H)[2] = G.n3; 
	((Xtype *)&(H[12]))[0] = G.a1;
	((Xtype *)&(H[12]))[1] = G.b1;
	((Xtype *)&(H[12]))[2] = G.a2;
	((Xtype *)&(H[12]))[3] = G.b2;
	((Xtype *)&(H[12]))[4] = G.a3;
	((Xtype *)&(H[12]))[5] = G.b3;	
	
	out.write(H,N_header);
	out.write((char *)(G.array),N_data);
	out.close();	
	
	return 0;
}


template<classgridPars_ >
int readFile( grid3D<gridPars_ > & G, const char * fname )
{
	int N_header = 6*sizeof(Xtype )+3*sizeof(int );
	
	static char * H = (char *) new char [N_header];
	
	ifstream in;
	in.open(fname, fstream::in);
	if (!in.good() || !in.is_open())
		return 1;
	
	in.read(H,N_header);
	
	int n1,n2,n3;
	
	n1 = ((int *)H)[0];
	n2 = ((int *)H)[1];
	n3 = ((int *)H)[2];
	G.a1 = ((Xtype *)&(H[12]))[0];
	G.b1 = ((Xtype *)&(H[12]))[1];
	G.a2 = ((Xtype *)&(H[12]))[2];
	G.b2 = ((Xtype *)&(H[12]))[3];
	G.a3 = ((Xtype *)&(H[12]))[4];
	G.b3 = ((Xtype *)&(H[12]))[5];
	
	if ( !G.array || n1 != G.n1 || n2 != G.n2 || n3 != G.n3 )
	{
		if (G.array) delete [] G.array;
		G.n1 = n1;
		G.n2 = n2;
		G.n3 = n3;
		G.array = new Ytype [n1*n2*n3];
	}
	
	int N_data = G.n1*G.n2*G.n3*sizeof(Ytype );
	
	in.read((char *)G.array,N_data);
	in.close();
	
	return 0;
}

#endif

#endif
