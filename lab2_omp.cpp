#include <malloc.h>
#include <omp.h>
#include "matrix.h"

#define pb push_back

using namespace std;

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
	matrix d(M, N, 0);
	printf("%d\n", M);		
	printf("%d\n", N);		
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			d(i,j) = D[i*N+j];
			// printf("%f |", D[i*N+j]);
		}
		// printf("\n");
	}
	// d.print();
	d.shape();
	matrix d_t = d.transpose();
	d_t.shape();
	matrix svd = d*d_t;
	svd.shape();
	pair<matrix, matrix> eigens = svd.eigen();
	vector<double> eigenvalues = (eigens.first).diagonal();
	matrix eigenvectors = (eigens.second);
	eigenvectors.shape();
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
	matrix sigma_inv(M, N, 0);
	int e = eigenv_index.size()-1;
	for(int i=0; i<N;i++){
		if(eigenv_index[e].first<1e-5){
			sigma_inv(i,i)= 0;
		}
		else{
			sigma_inv(i,i)= 1/sqrt(eigenv_index[e].first);
		}
		e--;	
	}
	sigma_inv.print();
	
	matrix v(M, M, 0);
	e = eigenv_index.size()-1;	
	for(int j=0;j<M;j++){
		int index = eigenv_index[e].second;
		for(int i=0;i<M;i++){
			v(i,j) = eigenvectors(i, index);
		}
		e--;
	}

	// v.print();
	matrix u = ((d_t)*(v))*(sigma_inv);
	u.print();
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}
