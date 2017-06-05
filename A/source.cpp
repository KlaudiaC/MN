// Klaudia C - 14.04.2016

#include <cmath>

double secant(double (*f)(double), double a, double b, int M, double eps, double delta){
    double fa = f(a);
    M--;
	double fb = f(b);
	M--;

    while(fabs(fa - fb) > delta){
        double c = a - fa * ((a - b)/(fa - fb));
        double fc = f(c);
        M--;
        if(fabs(fc) <= eps) return c;

        b = a;
        a = c;
        fb = fa;
        fa = fc;
    }
    return a;
}

double bisection(double (*f)(double), double a, double b, int M, double eps, double delta){
    double fa = f(a);
	M--;
	if(fabs(fa) <= eps) return a;
	double fb = f(b);
	M--;
	if(fabs(fb) <= eps) return b;

    while(M--){
        double c = (a + b)/2;
        double fc = f(c);
        M--;

        if(fabs(fc) <= eps) return c;
        if(fabs(a - b) <= delta) break;
        if(fabs(fc) <= eps) break;
        if(fa * fb > 0){
            a = b;
            fa = fb;
            b = c;
            fb = fc;
        }
        else{
            if(fa * fc < 0){
                b = c;
                fb = fc;
            }
            else{
                a = c;
                fa = fc;
            }
        }
    }
    return secant(f, a, b, M, eps, delta);
}

double findZero(
    double (*f)(double),  // funkcja, ktorej zera szukamy w [a, b]
    double a,             // lewy koniec przedzialu
    double b,             // prawy koniec przedzialu
    int M,                // maksymalna dozwolona liczba wywolan funkcji f
    double eps,           // spodziewana dokladnosc zera
    double delta){        // wystarczajacy blad bezwzgledny wyniku
    return bisection(f, a, b, M, eps, delta);
}
