#ifndef MATRIX
#define MATRIX

#include <bits/stdc++.h>

using namespace std;

class matrix{
public:
	int rows;
	int columns;
	vector<vector<double> > mat;

	//constructor
	//type = 0 (all entries zero)
	//type = 1 (identity matrix)
	//type = 2 (all entries one)
	matrix(int rows1, int columns1, int type);
	//copy constructor
	matrix(const matrix &);


	matrix operator+(matrix &);
	matrix operator-(matrix &);
	matrix operator*(matrix &);

	//Transpose of the matrix
	matrix transpose();
	
	//scalar operations
	matrix scalaradd(double a);
	matrix scalarsub(double a);
	matrix scalarmul(double a);
	matrix scalardiv(double a);

	//access element by ()
 	double& operator()(const unsigned &, const unsigned &);

	void print();
	

};

#endif