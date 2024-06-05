#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>

double fun(const double x)
{
    return sin(x);
    // return x*sin(5*x);
}

double trapezoidal(const gsl_matrix *M)
{
    // (a+b)*h/2
    double sum=0;
    for (size_t i = 0; i < M->size1-1; i++)
    {
        const double a = gsl_matrix_get(M,i,1);
        const double b = gsl_matrix_get(M,i+1,1);
        const double h = gsl_matrix_get(M,i+1,0)-gsl_matrix_get(M,i,0);

        sum+= (a+b)*h/2.0;
    }

    return sum;
}


double simson(const gsl_matrix *M)
{
    double sum=0;
    const double h = gsl_matrix_get(M,1,0)-gsl_matrix_get(M,0,0);

    for (size_t i = 1; i < M->size1-1; i+=2)
    {
        const double s = h/3.0*(gsl_matrix_get(M,i-1,1)+4*gsl_matrix_get(M,i,1)+gsl_matrix_get(M,i+1,1));
        sum+=s;
    }

    return sum;
}

int main()
{
    const int n = 201; //liczba węzłów
    const double start_x = 0;
    const double end_x = M_PI;
    
    gsl_matrix *M = gsl_matrix_calloc(n,2);

    // wypełienie macierzy
    for(int i = 0; i < n; ++i)
    {
        const double tmp = start_x + (end_x - start_x) * i / (n - 1);

        gsl_matrix_set(M, i, 0,tmp);
        gsl_matrix_set(M, i, 1, fun(tmp));
    }


    printf("%f\n", trapezoidal(M));
    printf("%f\n", simson(M));



    return 0;
}
