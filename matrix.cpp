#include "matrix.h"
#include <omp.h>


#define pb push_back
using namespace std;


matrix::matrix(int rows1, int columns1, int type){
	rows = rows1;
	columns = columns1;
	int tofill=-1;
	if(type==2){
		tofill=1.0;
	}
	else{
		tofill=0.0;
	}
	for(int i=0;i<rows1;i++){
		vector<double> temp;
		for(int j=0; j<columns1; j++){
			temp.pb(tofill);
		}
		mat.pb(temp);
	}

	//make it identity
	if(type==1){
		for(int i=0;i<rows1;i++){
			mat[i][i] = 1;
		}
	}
}



// Addition of Two Matrices
matrix matrix::operator+(matrix &A){
    matrix sum(rows, columns, 0);
    #pragma omp parallel for
	    for(int i=0; i<rows; i++){
	        for(int j=0; j<columns; j++){
	            sum(i,j) = mat[i][j] + A(i,j);
	        }
	    }
    return sum;
}

// Subtraction of Two Matrices
matrix matrix::operator-(matrix &A){
    matrix sub(rows, columns, 0);
    #pragma omp parallel for
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            sub(i,j) = mat[i][j] - A(i,j);
        }
    }
    return sub;
}

// Multiplication of Two Matrices
matrix matrix::operator*(matrix &A){
    int r1 = rows;
    int c1 = columns;
    int r2 = A.rows;
    int c2 = A.columns;
	matrix product(r1, c2, 0);
    if(c1==r2){
    	#pragma omp parallel for
     	for(int i=0; i<r1; i++){
    		for(int j=0; j<c2; j++){
 		   		double temp = 0.0;
    			for(int k=0; k<c1; k++){
    				temp += mat[i][k] + A(k,j);
    			}
    			product(i, j) = temp;
    		}
    	}
    }
    else{
    	cout << "Error in muliplication\n";
    }
    return product;

}

// Transpose of the matrix
matrix matrix::transpose(){
	matrix transpose(columns, rows, 0);
    #pragma omp parallel for	
	for(int i=0; i<columns; i++){
		for(int j=0; j<rows; j++){
			transpose(i,j) = mat[j][i];
		}
	}
	return transpose;
}

//scalar methods
matrix matrix::scalaradd(double a){
	matrix scalaradd(rows, columns, 0);
    #pragma omp parallel for	
	for(int i=0; i<rows; i++){
		for(int j=0; j<columns; j++){
			scalaradd(i,j) = mat[i][j] + a;
		}
	}
	return scalaradd;
}


matrix matrix::scalarsub(double a){
	matrix scalarsub(rows, columns, 0);
    #pragma omp parallel for	
	for(int i=0; i<rows; i++){
		for(int j=0; j<columns; j++){
			scalarsub(i,j) = mat[i][j] - a;
		}
	}
	return scalarsub;
}

matrix matrix::scalarmul(double a){
	matrix scalarmul(rows, columns, 0);
    #pragma omp parallel for
	for(int i=0; i<rows; i++){
		for(int j=0; j<columns; j++){
			scalarmul(i,j) = mat[i][j]*a;
		}
	}
	return scalarmul;
}

matrix matrix::scalardiv(double a){
	if(a==0){
		cout << "Div by zero\n";
	}
	matrix scalardiv(rows, columns, 0);
    #pragma omp parallel for
	for(int i=0; i<rows; i++){
		for(int j=0; j<columns; j++){
			scalardiv(i,j) = mat[i][j]/a;
		}
	}
	return scalardiv;
}

// Access elements of matrix using M(row, col)
double& matrix::operator()(const unsigned &row, const unsigned & col)
{
    return mat[row][col];
}

void matrix::print(){
	for(int i=0;i<rows;i++){
		for(int j=0;j<columns;j++){
			cout << mat[i][j] << " | " ;
		}
		cout << endl;
	}
}