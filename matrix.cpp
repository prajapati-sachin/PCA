#include "matrix.h"
#include <omp.h>


#define pb push_back
using namespace std;


namespace fgivens{
	//function to find [cos, sin, rho]
	vector<double> givens(double a, double b){
		vector<double> result;
		if(b==0){
			result.pb(1);
			result.pb(0);
			result.pb(0);
		}
		else if(a==0){
			result.pb(0);
			result.pb(1);
			result.pb(1);
		}
		else{
			if(abs(b)>abs(a)){
				double tau=	(-1)*(a/b);
				double s = 1/sqrt(1+tau*tau);
				double c = s*tau;
				double rho = 2/c; 
				result.pb(c);
				result.pb(s);
				result.pb(rho);
			}
			else{
				double tau=	(-1)*(b/a);
				double c = 1/sqrt(1+tau*tau);
				double s = c*tau;
				double rho = s/2; 
				result.pb(c);
				result.pb(s);
				result.pb(rho);
			}
		}
		return result;
	}
	//fuction to find [cos, sin] using rho
	vector<double> givensinv(double rho){
		vector<double> result; 
		if(rho==0){
			result.pb(1);
			result.pb(0);
		}
		else if(rho==1){
			result.pb(0);
			result.pb(1);
		}
		else if(abs(rho)>2){
			double c = 2/rho;
			double s = sqrt(1-c*c);
			result.pb(c);
			result.pb(s);
		}
		else{
			double s = 2*rho;
			double c = sqrt(1-s*s);
			result.pb(c);
			result.pb(s);
		}
		return result;
	}
	//function to return rotation matrix corresponding to cos and sin
	matrix rotation(double cos, double sin){
		matrix rotation(2, 2, 0);
		rotation(0,0) = cos;
		rotation(1,1) = cos;

		rotation(0,1) = (-1)*sin;
		rotation(1,0) = sin;
		return rotation;
	}
}

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


matrix::matrix(const matrix &A){
	rows = A.rows;
	columns = A.columns;
	for(int i=0;i<A.rows;i++){
		vector<double> temp;
		for(int j=0; j<A.columns; j++){
			temp.pb(A.mat[i][j]);
		}
		mat.pb(temp);
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
    				temp += mat[i][k]*A(k,j);
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

// Return the diagonal of a matrix
vector<double> matrix::diagonal(){
	vector<double> result;
	for(int i=0;i<rows;i++){
		result.pb(mat[i][i]);
	} 
	return result;
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

void matrix::shape(){
	cout << "( "<<rows << ", " << columns << ")" << endl;
}

// matrix matrix::submat(int i1, int i2){
// 	matrix submat(2, columns, 0);
// 	for(int j=0;j<columns;j++){
// 		submat(0,j) = mat[i1][j];
// 		submat(1,j) = mat[i2][j];
// 	}
// 	return submat;
// }

matrix matrix::upper(){
	matrix result(rows, columns, 0);
	for(int i=0;i<rows;i++){
		for(int j=i;j<columns;j++){
			result(i,j)= mat[i][j];
		}
	} 
	return result;
}

pair<matrix, matrix> matrix::qr(){
	if(rows!=columns){
		cout << "Error in QR" << endl;
	}
	matrix copy = *this;
	// copy.print();
	//Q matrix is identity in the beginining
	matrix Q(rows, columns, 1);
	int m = rows;
	// Making the matrix R
	for(int j=0; j<m-1; j++){
		for(int i=m-1; i>=j+1; i--){
			// cout << "Processing: " << copy(i-1,j) <<", " << copy(i,j) << endl;
			vector<double> csp = fgivens::givens(copy(i-1,j), copy(i,j));
			double cos = csp[0];
			double sin = csp[1];			
			double rho = csp[2];
			// cout << "cos, sin, rho: " << cos << ", " << sin << ", "<< rho <<  endl; 
			matrix temp = fgivens::rotation(cos, sin);
			// temp.print();
			// cout << endl;
			matrix submat(2, columns-j, 0);
			int iter=0;
			for(int k=j;k<columns;k++){
				submat(0,iter) = copy(i-1,k);
				submat(1,iter) = copy(i,k);
				iter++;
			}
			// submat.print();
			// cout << endl;			
			matrix product = temp*submat;
			// product.print();
			// cout << endl;
			iter=0;
			for(int k=j;k<columns;k++){
				copy(i-1,k) = product(0,iter);
				copy(i,k) = product(1,iter);
				iter++;
			}
			copy(i,j) = rho;
			// cout << "After Processing: " << i <<", " <<j << endl;
			// copy.print();
			// cout << endl;
		}
	}

	//Making the matrix Q
	for(int j=m-2;j>=0;j--){
		for(int i=j+1; i<=m-1;i++){
			vector<double> cs = fgivens::givensinv(copy(i,j));
			double cos = cs[0];
			double sin = cs[1];

			matrix temp = (fgivens::rotation(cos, sin)).transpose();
			matrix submat(2, columns-j, 0);
			int iter=0;
			for(int k=j;k<columns;k++){
				submat(0,iter) = Q(i-1,k);
				submat(1,iter) = Q(i,k);
				iter++;
			}
			// submat.print();
			// cout << endl;			
			matrix product = temp*submat;
			// product.print();
			// cout << endl;
			iter=0;
			for(int k=j;k<columns;k++){
				Q(i-1,k) = product(0,iter);
				Q(i,k) = product(1,iter);
				iter++;
			}
		}
	}

	return make_pair(Q, copy.upper());
}

pair<matrix, matrix> matrix::eigen(){
	matrix D = *this;
	matrix E(rows, columns, 1);
	// D.print();
	// cout << endl;
	// E.print();
	int count=0;
	while(1){
		count++;
		cout << count << endl;
		//Repeat until convergence
		pair<matrix,matrix> qr = D.qr();
		matrix Q = qr.first;
		matrix R = qr.second;
		// cout << "-----------------------" << endl;		
		// (Q.transpose()).print();
		// matrix temp = (Q.transpose());
		// (Q*temp).print();	
		// cout << endl;
		// (R).print();
		// // cout << endl;
		// cout << "-----------------------" << endl;
		// D = R*Q;
		E = E*Q;
		matrix check = R*Q;
		//check for convergence
		double max= INT_MIN;
		for(int i=0;i<rows;i++){
			if(abs(check(i,i)-D(i,i))>max) max = abs(check(i,i)-D(i,i)); 
		}
		D = check;	
		//if maximum change in the values of eigenvalues is less then epsilon we are converged
		cout << "Max: " << max << endl;
		if(max<1e-5){
			// cout << count << endl;
			break;
		}
	}
	return make_pair(D, E);
}