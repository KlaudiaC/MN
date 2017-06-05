// Klaudia C - 10.06.2016

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>
using namespace std;

double EPS = 1e-15;

void QR(double (*A)[1000], int N){
	double a[1000], b[1000];

	for(int k = 0; k < N - 1; k++){
		a[k]= A[k][k] / sqrt(A[k][k] * A[k][k] + A[k + 1][k] * A[k + 1][k]);
		b[k]= A[k + 1][k] / sqrt(A[k][k] * A[k][k] + A[k + 1][k] * A[k + 1][k]);

		for(int j = k; j < N; j++){
			double c = a[k] * A[k][j] + b[k] * A[k + 1][j];
			double d = b[k] * A[k][j] - a[k] * A[k + 1][j];
			A[k][j] = c;
			A[k + 1][j] = d;
		}
	}
	for(int k = 0; k < N - 1; k++){
		for(int i = 0; i <= k + 1; i++){
			double c = a[k] * A[i][k] + b[k] * A[i][k + 1];
			double d = b[k] * A[i][k] - a[k] * A[i][k + 1];
			A[i][k] = c;
			A[i][k + 1] = d;
		}
	}
}

int main(){
	ios_base::sync_with_stdio(0);
	int N;

	cin >> N;
	double (*A)[1000] = new double[N][1000];

	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			cin >> A[i][j];

	for(int k = 0; k < N - 2; k++){
		pair <double, int> J(-1, -1);
		for(int w = k + 1; w < N; w++)
			J = max(J, make_pair(fabs(A[w][k]), w));
		for(int i = k; i < N; i++)
			swap(A[J.second][i], A[k + 1][i]);
		for(int i = 0; i < N; i++)
			swap(A[i][J.second], A[i][k + 1]);
		if(fabs(A[k + 1][k]) < EPS)
            continue;

		double C[N];
		for(int w = k + 2; w < N; w++){
			C[w] = A[w][k]/A[k + 1][k];
			A[w][k] = 0;
			for(int i = k + 1; i < N; i++)
				A[w][i] -= C[w] * A[k + 1][i];
		}

		for(int w = 0; w < N; w++)
			for(int i = k + 2; i < N; i++)
				A[w][k+1] += C[i] * A[w][i];
	}

	vector<double> res;
	int N1 = N;

	while(N1 > 1){
		if(fabs(A[N1 - 1][N1 - 2]) < EPS){ // deflacja
			res.push_back(A[N1 - 1][N1 - 1]);
			N1--;
			continue;
		}

		double z = A[N1 - 1][N1 - 1];
		for(int i = 0; i < N1; i++)
            A[i][i] -= z;
		QR(A, N1);
		for(int i = 0; i < N1; i++)
            A[i][i] += z;
	}
	res.push_back(A[0][0]);
	sort(res.begin(), res.end());
	cout << setprecision(17);
	for(int i = 0; i < N; i++)
        cout << res[i] << endl;

	delete A;
	return 0;
}
