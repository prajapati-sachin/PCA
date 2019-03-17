#include <bits/stdc++.h>
#include "matrix.h"
#include <omp.h>


using namespace std;

int main(){
	#pragma omp parallel
	{
		cout << "Welcome\n";
	}

	matrix m1(3,3,0);
	matrix m2(3,3,1);
	matrix m3(3,3,2);
	m1.print();
	cout << endl;
	m2.print();
	cout << endl;
	m3.print();
	cout << endl;

	matrix m4 = m2+m3;
	m4.print();
	cout << endl;
	
	matrix m5 = m3-m2;
	
	m5.print();
	cout << endl;
	
	matrix m6 = m5*m4; 

	m6(1,2)= 9;
	m6(0,2)= 87;
	m6(2,1)= 29;
	m6(2,0)= 39;

	m6.print();
	cout << endl;

	m6.transpose().scalardiv(2).print();
	cout << endl;

	m6.print();

	matrix m7 = m6;
	m7.print();

	return 0;
}