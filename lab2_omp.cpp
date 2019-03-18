#pragma GCC optimize ("O3")

#include <malloc.h>
#include <omp.h>
#include <bits/stdc++.h>
// #include "matrix.h"

#define pb push_back

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
	
	//Return the diagonal of a matrix in vector
	vector<double> diagonal();

	//scalar operations
	matrix scalaradd(double a);
	matrix scalarsub(double a);
	matrix scalarmul(double a);
	matrix scalardiv(double a);

	//print the matrix
	void print();

	//print the shape of matrix
	void shape();

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
    matrix A_t = A.transpose();
    if(c1==r2){
    	#pragma omp parallel for
     	for(int i=0; i<r1; i++){
    		for(int j=0; j<c2; j++){
 		   		double temp = 0.0;
    			for(int k=0; k<c1; k++){
    				temp += mat[i][k]*A_t(j,k);
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
			// matrix temp = fgivens::rotation(cos, sin);
			// // temp.print();
			// // cout << endl;
			// matrix submat(2, columns-j, 0);
			// int iter=0;
			// for(int k=j;k<columns;k++){
			// 	submat(0,iter) = copy(i-1,k);
			// 	submat(1,iter) = copy(i,k);
			// 	iter++;
			// }
			// // submat.print();
			// // cout << endl;			
			// matrix product = temp*submat;
			// // product.print();
			// // cout << endl;
			// iter=0;
			// for(int k=j;k<columns;k++){
			// 	copy(i-1,k) = product(0,iter);
			// 	copy(i,k) = product(1,iter);
			// 	iter++;
			// }
			for(int k=j;k<columns;k++){
				double num1 = cos*copy(i-1,k) - sin*copy(i,k);
				double num2 = sin*copy(i-1,k) + cos*copy(i,k);
				copy(i-1,k) = num1;
				copy(i,k) = num2;
				// iter++;
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

			// matrix temp = (fgivens::rotation(cos, sin)).transpose();
			// matrix submat(2, columns-j, 0);
			// int iter=0;
			// for(int k=j;k<columns;k++){
			// 	submat(0,iter) = Q(i-1,k);
			// 	submat(1,iter) = Q(i,k);
			// 	iter++;
			// }
			// // submat.print();
			// // cout << endl;			
			// matrix product = temp*submat;
			// // product.print();
			// // cout << endl;
			// iter=0;
			// for(int k=j;k<columns;k++){
			// 	Q(i-1,k) = product(0,iter);
			// 	Q(i,k) = product(1,iter);
			// 	iter++;
			// }
			for(int k=j;k<columns;k++){
				double num1 = cos*Q(i-1,k) + sin*Q(i,k);
				double num2 = -sin*Q(i-1,k) + cos*Q(i,k);
				Q(i-1,k) = num1;
				Q(i,k) = num2;
				// iter++;
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


// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
	matrix d(M, N, 0);
	// printf("%d\n", M);		
	// printf("%d\n", N);		
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			d(i,j) = D[i*N+j];
			// printf("%f |", D[i*N+j]);
		}
		// printf("\n");
	}
	// d.print();
	// d.shape();
	matrix d_t = d.transpose();
	// d_t.shape();
	matrix svd = d_t*d;
	// svd.shape();
	pair<matrix, matrix> eigens = svd.eigen();
	vector<double> eigenvalues = (eigens.first).diagonal();
	matrix eigenvectors = (eigens.second);
	// eigenvectors.print();
	vector<pair<double, int> > eigenv_index;
	for(int i=0;i<eigenvalues.size();i++){
		eigenv_index.pb(make_pair(eigenvalues[i],i));
		// eigenv_index.pb(make_pair(0,i));
		// cout << eigenvalues[i] << endl;
	}
	sort(eigenv_index.begin(), eigenv_index.end());
	// reverse(eigenv_index.begin(), eigenv_index.end());
	// for(int i=0;i<eigenv_index.size();i++){
	// 	cout << eigenv_index[i].first << "| ";
	// 	cout << eigenv_index[i].second << endl;
	// }
	
	// M*N beacuse of transpose and we take svd od D_T
	matrix sigma(N, M, 0);
	int e = eigenv_index.size()-1;
	for(int i=0; i<N;i++){
		// if(eigenv_index[e].first<1e-5){
			// sigma(i,i)= 0;
		// }
		// else{
			sigma(i,i)= sqrt(eigenv_index[e].first);
		// }
		e--;	
	}
	
	// sigma.print();

	matrix sigma_inv(M, N, 0);
	e = eigenv_index.size()-1;
	for(int i=0; i<N;i++){
		if(eigenv_index[e].first<1e-5){
			sigma_inv(i,i)= 0;
		}
		else{
			sigma_inv(i,i)= 1/sqrt(eigenv_index[e].first);
		}
		e--;	
	}
	// sigma_inv.print();
	matrix sigma_invT = sigma_inv.transpose();
	// sigma_invT.shape();

	
	matrix v(M, M, 0);

	v = (d*eigenvectors)*sigma_invT;

	// e = eigenv_index.size()-1;	
	// for(int j=0;j<M;j++){
	// 	int index = eigenv_index[e].second;
	// 	for(int i=0;i<M;i++){
	// 		v(i,j) = eigenvectors(i, index);
	// 	}
	// 	e--;
	// }

	// v.print();
	// matrix u = ((d_t)*(v))*(sigma_inv);
	// u.print();

	// writing the matrices
	//Writing U
	// #pragma omp parallel for
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			(*U)[i*N+j] = eigenvectors(i,j);
		}
	}

	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		cout << (*U)[i*N+j] << " ";
	// 	}
	// 	cout << endl;
	// }

	//Writing sigma
	// #pragma omp parallel for	
 	for(int i=0;i<N;i++){
		// for(int j=0;j<M;j++){
			(*SIGMA)[i] = sigma(i,i);
		// }
	}
	
	//Writing V_T
	matrix vtranspose = v.transpose();
	// #pragma omp parallel for	
 	for(int i=0;i<M;i++){
		for(int j=0;j<M;j++){
			(*V_T)[i*M+j] = vtranspose(i,j);
		}
	}
	// cout << "------------------------------------------------------------------------";
	// vtranspose.print();
	// for(int i=0;i<N;i++){
	// 	for(int j=0;j<N;j++){
	// 		printf("%f |", *U[i*N+j]);
	// 	}
	// 	printf("\n");
	// }

}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
	double num=0;
	int k=0;
	double sigmasqsum=0;
	for(k=0; k<N; k++){
		sigmasqsum+=SIGMA[k]*SIGMA[k];
	}
	// #pragma omp parallel for
	for(k=0; k<N; k++){
		num+=(SIGMA[k]*SIGMA[k])/sigmasqsum;
		if(num>=retention/100.0){
			break;
		}
	}
    
    *K = k+1;
    matrix d(M, N, 0);
	// #pragma omp parallel for    
    for(int i=0;i<M;i++){
    	for(int j=0;j<N;j++){
    		d(i,j)=D[i*N + j];
    	}
    }

    // d.print();
    // d.shape();
    // cout << "\nK: \n" << k+1 << endl;

    matrix newU(N, k+1, 0);
	// #pragma omp parallel for        
    for(int i=0; i<N; i++){
    	for(int j=0;j<k+1;j++){
    		newU(i,j)=U[i*N+ j];
    	}
    }

    // newU.print();
    matrix d_hat = d*newU;

   	// d_hat.shape();

	*D_HAT = (float*) malloc(sizeof(float) * M*(k+1));

	// #pragma omp parallel for    
	for(int i=0; i<M; i++){
    	for(int j=0;j<k+1;j++){
    		(*D_HAT)[i*(k+1)+j] = d_hat(i,j);
    		// newU(i,j)=U[i*(k+1)+ j];
    	}
    }




}
