#ifndef LC_H
#define LC_H

/*
 * linear combination of common types, this is a useful tool 
 * that can speed things up a lot, yet maintain abstraction.
 * 
 * 
 * 
 */

template<class T,class D>
void LC( T & R,
		 T * L,D * C, int n )
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 0\n";
		err = 0; }	
		
	R = L[0]*C[0];
	for (int i = 1;i<n;i++)
		R += L[i]*C[i];
}


template<class T,class D>
void LC( T & R,
		 T & L0,D C0 )
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }

	R = L0*C0;	
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1 )
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1;
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2;		
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3;		
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4;		
}

template<class T,class D>
void LC( T & R,
  		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5;		
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5,
	    T & L6,D C6)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5
		+L6*C6;		
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5,
	    T & L6,D C6,
	    T & L7,D C7)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5
		+L6*C6
		+L7*C7;		
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5,
	    T & L6,D C6,
	    T & L7,D C7,
	    T & L8,D C8)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5
		+L6*C6
		+L7*C7
		+L8*C8;		
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5,
	    T & L6,D C6,
	    T & L7,D C7,
	    T & L8,D C8,
	    T & L9,D C9)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5
		+L6*C6
		+L7*C7
		+L8*C8
		+L9*C9;	
}


template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5,
	    T & L6,D C6,
	    T & L7,D C7,
	    T & L8,D C8,
	    T & L9,D C9,
	    T & L10,D C10)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5
		+L6*C6
		+L7*C7
		+L8*C8
		+L9*C9
		+L10*C10;	
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5,
	    T & L6,D C6,
	    T & L7,D C7,
	    T & L8,D C8,
	    T & L9,D C9,
	    T & L10,D C10,
	    T & L11,D C11)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5
		+L6*C6
		+L7*C7
		+L8*C8
		+L9*C9
		+L10*C10
		+L11*C11;	
}

template<class T,class D>
void LC( T & R,
		 T & L0,D C0,
	    T & L1,D C1,
	    T & L2,D C2,
	    T & L3,D C3,
	    T & L4,D C4,
	    T & L5,D C5,
	    T & L6,D C6,
	    T & L7,D C7,
	    T & L8,D C8,
	    T & L9,D C9,
	    T & L10,D C10,
	    T & L11,D C11,
	    T & L12,D C12)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize LC 1\n";
		err = 0; }


	R = L0*C0
		+L1*C1
		+L2*C2
		+L3*C3
		+L4*C4
		+L5*C5
		+L6*C6
		+L7*C7
		+L8*C8
		+L9*C9
		+L10*C10
		+L11*C11
		+L12*C12;	
}










#endif


