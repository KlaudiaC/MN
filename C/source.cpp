// Klaudia C - 29.04.2016

#include <cstdio>
#include <cmath>
#include <algorithm>
using namespace std;

void printVector(const double* x, unsigned N){
    for(unsigned i = 0; i < N; ++i) printf("%17.17f ",x[i]);
    printf("\n");
}

typedef void (*FuncPointer)(const double* x, double* y, double* Df);

const double E = 0.5e-14;

struct myStruct{
	FuncPointer f;
	double c;

	myStruct(FuncPointer f, double c): f(f), c(c){}

	void operator()(const double* x, double* y, double* Df){
		double df[6];
		double xx[3];
		xx[0] = x[0];
		xx[1] = x[1];
		xx[2] = c;
		f(xx, y, df);
		Df[0] = df[0];
		Df[1] = df[1];
		Df[2] = df[3];
		Df[3] = df[4];
	}
};

struct myStruct2{
	FuncPointer f;
	double c, d;

	myStruct2(FuncPointer f, double c, double d): f(f), c(c), d(d){}

	void operator()(const double* x, double* y, double* Df){
		double df[3];
		double xx[3];
		xx[0] = x[0];
		xx[1] = c;
		xx[2] = d;
		f(xx, y, df);
		Df[0] = df[0];
	}
};

struct myStruct3{
	FuncPointer f;
	double c, d;

	myStruct3(FuncPointer f, double c, double d): f(f), c(c), d(d){}

	void operator()(const double* x, double* y, double* Df){
		double df[8];
		double xx[4];
		xx[0] = x[0];
		xx[1] = x[1];
		xx[2] = c;
		xx[3] = d;
		f(xx, y, df);
		Df[0] = df[0] - 1.0;
		Df[1] = df[1];
		Df[2] = df[4];
		Df[3] = df[5] - 1.0;
		y[0] -= x[0];
		y[1] -= x[1];
	}
};


template <class T>
int newton2D(T f, double* x){
	for(int i = 0; i < 1000; i++){
		double y[2];
		double df[4];

		f(x, y, df);
		if(max(fabs(y[0]), fabs(y[1])) < E) return 0;

		double det = df[0] * df[3] - df[1] * df[2];
		if(fabs(det) < E) return -1;

		double temp = df[0];
		df[0] = df[3];
		df[3] = temp;
		df[1] *= -1.0;
		df[2] *= -1.0;

		double h[2];
		h[0] = -(y[0] * df[0] + y[1] * df[1])/det;
		h[1] = -(y[0] * df[2] + y[1] * df[3])/det;

		x[0] += h[0];
		x[1] += h[1];
	}
	return -1;
}

template <class T>
int newton1D(T f, double* x){
	for(int i = 0; i < 1000; i++){
		double y;
		double df;

		f(x, &y, &df);
		if(fabs(y) < E) return 0;
		if(fabs(df) < E) return -1;

		double h;
		h = -y/df;
		x[0] += h;
	}
	return -1;
}

int findCurve(FuncPointer f, double* x, unsigned k, double h){
	for(unsigned int i = 1; i <= k; i++){
		double res[2];
		res[0] = x[0];
		res[1] = x[1];
		int zle = newton2D(myStruct(f, x[2] + i * h), res);
		if(zle) return i;

		printf("%17.17f %17.17f %17.17f \n", res[0], res[1], x[2] + i * h);
	}
	return 0;
}

int findSurface(FuncPointer f, double* x, unsigned k1, unsigned k2, double h1, double h2){
	for(unsigned int i = 1; i <= k1; i++){
		for(unsigned int j = 1; j <= k2; j++){
			double res = x[0];
			int zle = newton1D(myStruct2(f, x[1] + i * h1, x[2] + j * h2), &res);
			if(zle) return i * k1 + j;

			printf("%17.17f %17.17f %17.17f \n", res, x[1] + i * h1, x[2] + j * h2);
		}
		printf("\n");
	}
	return 0;
}

int findFixedPoints(FuncPointer f, double* x, unsigned k1, unsigned k2, double h1, double h2){
	for(unsigned int i = 1; i <= k1; i++){
		for(unsigned int j = 1; j <= k2; j++){
			double res[2];
			res[0] = x[0];
			res[1] = x[1];
			int zle = newton2D(myStruct3(f, x[2] + i * h1, x[3] + j * h2), res);
			if(zle) return i * k1 + j;

			printf("%17.17f %17.17f %17.17f %17.17f \n", res[0], res[1], x[2] + i * h1, x[3] + j * h2);
		}
		printf("\n");
	}
	return 0;
}
