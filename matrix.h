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

	//print the matrix
	void print();

	//Return the upper triangular matrix of this matrix
	matrix upper();

	// matrix rotation(double cos, double sin);
	// // Gives a submatrix of 2 rows  from i1 to i2 and columns column
	// matrix submat(int i1, int i2);

	//QR Decomposition of the matrix
	pair<matrix, matrix> qr();

	//Return matrix of eigenvalues and eignevector
	pair<matrix, matrix> eigen();

	//access element by ()
 	double& operator()(const unsigned &, const unsigned &);


	

};

namespace fgivens{
	//function to find [cos, sin, rho]
	vector<double> givens(double a, double b);
	//fuction to find [cos, sin] using rho
	vector<double> givensinv(double rho);
	//fuction return the (2*2) rotation matrix for given cos and sin
	matrix rotation(double cos, double sin);

}


#endif