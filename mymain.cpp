#include <bits/stdc++.h>
#include "matrix.h"
#include <omp.h>


using namespace std;

int main(){
	// #pragma omp parallel
	// {
	// 	cout << "Welcome\n";
	// }

	matrix m1(3,3,0);
	matrix m2(3,3,1);
	matrix m3(3,3,2);
	// m1.print();
	// cout << endl;
	// m2.print();
	// cout << endl;
	// m3.print();
	// cout << endl;

	matrix m4 = m2+m3;
	// m4.print();
	// cout << endl;
	
	matrix m5 = m3-m2;
	
	// m5.print();
	// cout << endl;
	
	matrix m6 = m5*m4; 

	m6(1,2)= 9;
	m6(0,2)= 87;
	m6(2,1)= 29;
	m6(2,0)= 39;

	// m6.print();
	// cout << endl;

	// m6.transpose().scalardiv(2).print();
	// cout << endl;

	// m6.print();

	// matrix m7 = m6;
	// m7.print();
	// cout << endl;
	// pair<matrix, matrix> p = m7.qr();
	
	// p.first.print();
	// cout << endl;
	// p.second.print();

	// matrix final = p.first*p.second;
	// cout << endl;
	// final.print();


	// matrix test(3,3 ,0);
	// test(0,0) = 1;
	// test(1,0) = 3;
	// test(2,0) = 6;
	
	// test(0,1) = -3;
	// test(1,1) = -5;
	// test(2,1) = -6;

	// test(0,2) = 3;
	// test(1,2) = 3;
	// test(2,2) = 4;


	matrix test(3,3 ,0);
	test(0,0) = 3;
	test(1,0) = 2;
	test(2,0) = 4;
	
	test(0,1) = 2;
	test(1,1) = 0;
	test(2,1) = 2;

	test(0,2) = 4;
	test(1,2) = 2;
	test(2,2) = 3;

	// matrix test(2,2 ,0);
	// test(0,0) = 1;
	// test(1,0) = 0;
	
	// test(0,1) = 0;
	// test(1,1) = 1;



	test.print();
	cout << endl;

	pair<matrix, matrix> yay = test.eigen();
	yay.first.print();
	cout << endl;
	yay.second.print();


	// m7(0,0)=999;
	// m7.print();
	// p.first.print();
	// p = m7.qr();

	// p.first.print();

	return 0;
}