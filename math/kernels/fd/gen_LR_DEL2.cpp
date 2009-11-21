
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <fstream>
#include <complex>
#include <vector>

#include "matrices/matrix1.cpp"

#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include <arprec/mp_int.h>

using namespace std;

void setup();
char * fname = new char [200];
char * pwd = new char [200];


int m;
int L,R;
int k = 2;

mp_real A_ij(int n, int i)
{
	n -= 1;
	i -= L+1;
	
	if (n)
		return mp_real( pow((mp_real)i,n) / gamma((mp_real)(n+1))  );
	else return mp_real(1.0);
}








int main()
{setup();
	mp::mp_init(700);
	//mp::mpsetprecwords(300);
	
	vector<int> Lv;
	vector<int> Rv;
	vector< vector<mp_real> > Cv;
	
	int ML = 25;
	int MR = 25;	
	int mm = 4;
	
	for (L = 0;L<=ML;L++)
	for (R = 0;R<=MR;R++)
	if (R+L>=mm)
	{
		cout << L << " " << R << endl;
		int size =L+R+1;
		matrix1<mp_real> A(size,size);
		A = A_ij;
		matrix1<mp_real> x(size,size);
		matrix1<mp_real> b(size,size);
		b(k+1,1) = 1.0;
		matrix1_LU<mp_real> lu(A);
		lu.solve(x,b);

		vector<mp_real> c;
		for (int i = 1;i<=size;i++)
			c.push_back(x(i,1));
		
		Lv.push_back(L);
		Rv.push_back(R);
		Cv.push_back(c);
	}
	else
	{
		Lv.push_back(L);
		Rv.push_back(R);
		vector<mp_real> c;
		Cv.push_back(c);
	} 
	
	ofstream out;
	sprintf(fname,"FD_DEL2_LR.h");
	out.open(fname,ios::out);	
	assert(out.good() && out.is_open());
	
	int ndigits = 18;
	int exponent;
	char * digits = new char[ndigits+5];

	int u,v;
	
	for (u=0;u<Cv.size();u++)
	{
		L = Lv[u];
		R = Rv[u];
		if (L+R>=mm) 
		{
			out <<  "const double FD_DEL2_LR_" << Lv[u] << "_" << Rv[u] << "_"  << "[] = {\n";
			for (v=0;v<Cv[u].size();v++)
			{
				Cv[u][v].to_digits(digits, exponent, ndigits);
				if (exponent < 0 && exponent < -306)
					{ out << "0";}
				else
				{
					if ( Cv[u][v]< 0.0 ) out << "-";
					else out << "+";
					out << digits[0] << "." << &(digits[1]) << "e";
					if (exponent >= 0) out << "+";
					out << exponent;
				}
				
				if (v < Cv[u].size()-1 ) out << "," << endl;
					else out << "};	" << endl << endl;		
			}
			out << endl;
		}
		
	}
	
	out <<  "const double * FD_DEL2_LR_[" << ML+1 << "][" << MR+1 << "] = {\n";
	for (u=0;u<Cv.size();u++)
	{
		L = Lv[u];
		R = Rv[u];
		if (L+R>=mm) 
		{
			out << "&FD_DEL2_LR_" << L << "_" << R << "_[0]";
		}
		if (L+R<mm)		out << "0";
		
		if (u<Cv.size()-1) out << ",\n";
		else out << "};\n\n";
	}
	out.close();
}


void setup()
{
	ifstream in;
	sprintf(fname,"pwd.txt");
	in.open(fname,ios::in);
	assert(in.good() && in.is_open());
	assert( in >> pwd );
	in.close();	
	setRealTime_0();
}











