// Klaudia C - 12.05.2016

const int N = 2000;
#include <algorithm>
#include <cmath>
using namespace std;

void solveEquations(double A[N][N], int n, double b[N], double eps, double xx[N]){
	long double (*LU)[N] = new long double [n][N];
	long double s[N];
	long double r[N];
	long double e[N];
	long double x[N];
	int p[N];

	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			LU[i][j] = A[i][j];

	for(int i = 0; i < n; i++){
		p[i] = i;
		s[i] = 0;
		for(int j = 0; j < n; j++)
			s[i] = max(s[i], fabs(LU[i][j]));
	}
	for(int k = 0; k < n - 1; k++){
		pair<long double, int> J(-1.0, -1);
		for(int j = k; j < n; j++)
			J = max(J, make_pair( fabs(LU[p[j]][k]) / s[p[j]], j ));
		swap(p[k], p[J.second]);

		for(int i = k + 1; i < n; i++){
			long double z = LU[p[i]][k]/LU[p[k]][k];
			LU[p[i]][k] = z;
			for(int j = k + 1; j < n; j++)
				LU[p[i]][j] -= z * LU[p[k]][j];
		}
	}
	for(int i = 0; i < n; i++){
		r[i] = b[i];
		x[i] = 0;
	}

	long double rmax;
	do{
		for(int k = 0; k < n - 1; k++)
			for(int i = k + 1; i < n; i++)
				r[p[i]] -= LU[p[i]][k] * r[p[k]];

		for(int i = n - 1; i >= 0; i--){
			e[i] = r[p[i]];
			for(int j = i + 1; j < n; j++)
				e[i] -= LU[p[i]][j] * e[j];
			e[i] /= LU[p[i]][i];
		}

		for(int i = 0; i < n; i++){
			x[i] += e[i];
			r[i] = b[i];
		}

		rmax = 0;

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++)
				r[i] -= A[i][j] * x[j];
			rmax = max(rmax, fabs(r[i]));
		}
	}
	while(rmax >= eps/10);

	delete LU;
	for(int i=0; i<N; i++)
        xx[i]=x[i];
}
