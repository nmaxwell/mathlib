
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

#define type mp_real

int m;
int L,R;
int k = 2;

type A_ij(int i, int j)
{
	return mp_real(pow(j,2*(i-1)));
}
	
mp_real Dp_i(int p, int i)
{
	if (i == 0)
		return mp_real(-2.0)/(mp_real(p*p));
	else if (i==p)
		return mp_real(1.0)/(mp_real(p*p));
	else if (i==-p)
		return mp_real(1.0)/(mp_real(p*p));
	else return mp_real(0.0);
}




int main()
{
	srand(time(0));
	mp::mp_init(700);
	
	int M = 35;
	
	mp_real ** C_ = new mp_real * [M+1];
	
	
	for (int m = 1;m<=M;m++)
	{
		cout << m << endl;
		matrix1<mp_real> A(m,m);
		matrix1<mp_real> x(m,1);
		matrix1<mp_real> b(m,1);
		b *= 0.0;
		b(1,1) = 1.0;	
		A = A_ij;
		matrix1_LU<mp_real> lu(A);
		lu.solve(x,b);
		
		mp_real * c_ = new mp_real[m+1];
				
		for (int i = 0;i<=m;i++)
		{
			c_[i] = mp_real(0.0);
			for (int p = 1;p<=m;p++)
				c_[i] += x(p,1)*Dp_i(p,i);
		}
		
		C_[m] = c_;
	}

	ofstream out;
	sprintf(fname,"FD_DEL2_CENTERED.h");
	out.open(fname,ios::out);	
	assert(out.good() && out.is_open());
	
	int ndigits = 18;
	int exponent;
	char * digits = new char[ndigits+5];

	for (int m = 1;m<=M;m++)
	{	
		out <<  "double FD_DEL2_CENTERED_C_" << m << "[] = {\n";
		for (int i = 0;i<=m;i++)
		{
			C_[m][i].to_digits(digits, exponent, ndigits);
			if (C_[m][i] < 0.0) out << "-";
			else out << "+";
			out << digits[0] << "." << &(digits[1]) << "e";
			if (exponent >= 0) out << "+";
			out << exponent;
			if (i != m) out << "," << endl;
				else out << "};	" << endl << endl;
		}
	}
	
	out << endl << "double * FD_DEL2_CENTERED_[] = {\n";
	out << "0,\n";
	for (int m = 1;m<=M;m++)
	{
			out << "&FD_DEL2_CENTERED_C_" << m << "[0]";
			if (m != M) out << "," << endl;
				else out << "};	" << endl << endl;
	}
	
}















