#include <iostream>

#include "source.cpp"
using namespace std;

int main(){
    double (*a)[N] = new double [N][N];
    double b[N];
    double x[N];
    int n;
    double eps;

    cin >> n;
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j) cin >> a[i][j];
    for(int i = 0; i < n; ++i)
        cin >> b[i];
    cin >> eps;
    solveEquations(a, n, b, eps, x);

    delete [] a;
    return 0;
}
