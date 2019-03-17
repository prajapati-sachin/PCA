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
		if(abs(eigenvalues[i])>1e-5) eigenv_index.pb(make_pair(eigenvalues[i],i));
		else eigenv_index.pb(make_pair(0,i));
		// cout << eigenvalues[i] << endl;
	}
	sort(eigenv_index.begin(), eigenv_index.end());
	// reverse(eigenv_index.begin(), eigenv_index.end());
	// for(int i=0;i<eigenv_index.size();i++){
	// 	cout << eigenv_index[i].first << "| ";
	// 	cout << eigenv_index[i].second << endl;
	// }
	matrix sigma_inv(M, N, 0);
	int e = eigenv_index.size()-1;
	for(int i=0; i<N;i++){
		sigma_inv(i,i)= 1/(eigenv_index[e].first);
		e--;	
	}
	sigma_inv.print();
	matrix v(N, N, 0);
	e = eigenv_index.size()-1;	
	// for(int i=0;i<N;i++){
		// int index = eigenv_index[e].second;
		// for()
	// }

}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}
