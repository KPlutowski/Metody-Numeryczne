#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

int save_vector_to_file(const gsl_vector *b, const size_t x, const char *name)
{
    FILE *f;
    f = fopen(name, "a");
    if (f == NULL)
    {
        perror("Nie udalo sie otworzyc pliku do zapisu");
        return 0;
    }

    for (int i = 0; i < x; i++)
    {
        fprintf(f, "%5g\t", gsl_vector_get(b, i));
    }
    fprintf(f, "\n");
    fflush(f);
    fclose(f);
    return 1;
}
int save_matrix_to_file(const gsl_matrix *A, const size_t x, const size_t y, const char *name)
{
    gsl_vector *vec = gsl_vector_calloc(x);
    for (int i = 0; i < y; i++)
    {
        gsl_matrix_get_row(vec, A, i);
        save_vector_to_file(vec, y, name);
    }
    return 1;
}
int save_value_to_file(const double value, const char *name)
{
    FILE *f;
    f = fopen(name, "a");
    if (f == NULL)
    {
        perror("Nie udalo sie otworzyc pliku do zapisu");
        return 0;
    }

    fprintf(f, "%5g\n", value);

    fflush(f);
    fclose(f);

    return 1;
}

gsl_matrix *calloc_unit_matrix(const size_t n)
{
    gsl_matrix *A = gsl_matrix_calloc(n, n);
    for (int i = 0; i < n; i++)
        gsl_matrix_set(A, i, i, 0);
    return A;
}
double deviation(const gsl_vector *c, const gsl_vector *b, const size_t n)
{
    // odchylenie
    double deviation = 0;
    for (int i = 0; i < n; i++)
        deviation += pow(gsl_vector_get(c, i) - gsl_vector_get(b, i), 2);
    deviation = sqrt(deviation) / n;
    return deviation;
}
gsl_vector *vec_minus_vec(const gsl_vector *a, const gsl_vector *b, const size_t n)
{
    gsl_vector *c = gsl_vector_calloc(n);
    for (int i = 0; i < n; i++)
    {
        gsl_vector_set(c, i, gsl_vector_get(a, i) - gsl_vector_get(b, i));
    }
    return c;
}

gsl_vector *matrix_vector_multiplication(const gsl_matrix *A, const gsl_vector *x, const size_t n)
{
    gsl_vector *c = gsl_vector_calloc(n);
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum = 0;
        for (int j = 0; j < n; j++)
            sum += gsl_matrix_get(A, i, j) * gsl_vector_get(x, j);
        gsl_vector_set(c, i, sum);
    }
    return c;
}
double vector_multiplication(const gsl_vector *a, const gsl_vector *b, const size_t n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += gsl_vector_get(a, i) * gsl_vector_get(b, i);
    }
    return sum;
}
gsl_matrix *matrix_multiplication(const gsl_matrix *A, const gsl_matrix *B, const size_t n)
{
    // C = A*B
    gsl_matrix *C = gsl_matrix_calloc(n, n);

    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sum = 0;
            for (int k = 0; k < n; k++)
            {
                sum += gsl_matrix_get(A, i, k) * gsl_matrix_get(B, k, j);
            }
            gsl_matrix_set(C, i, j, sum);
        }
    }
    return C;
}

gsl_matrix *tensor_product(const gsl_vector *a, const size_t n)
{
    gsl_matrix *B = gsl_matrix_calloc(n, n);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
        {
            gsl_matrix_set(B, i, j, gsl_vector_get(a, i) * gsl_vector_get(a, j));
        }
    return B;
}

double f(const double D,const double A,const double B,const double K)
{
    return (D / ((A - D) * (B - D))) - K;
}
double pochodna(const double x,const double A,const double B,const double K)
{
    double h = 1e-7;
    return (f(x + h, A, B, K) - f(x, A, B, K)) / h;
}

double miejsce_zerowe_newton(const double D_initial, const double A,const double B, const double K)
{
    double D=D_initial;
    double delta;
    do
    {
        delta = f(D, A, B, K) / pochodna(D, A, B, K);
        D -= delta;
    } while (fabs(delta) > 1e-6);

    return D;
}
void miejsce_zerowe_newton_z_wypisywaniem(const double D_initial,const double A,const double B,const double K)
{
    double D = D_initial;
    double delta;
    printf("Przybl.\tKrok\tWart.\tPoch.\n");
    do
    {
        delta = f(D, A, B, K) / pochodna(D, A, B, K);
        printf("%f\t%f\t%f\t%f\n", D, delta, f(D, A, B, K), pochodna(D, A, B, K));
        D -= delta;
    } while (fabs(delta) > 1e-6);
}


bool compare_two_doubles(const double A, const double B)
{
    return fabs(A - B)<1e-6;
}

int main()
{
    double A = 2.0f;
    double B = 2.0f;
    double K = 0.25f;

    // pkt 1
    // printf("%f\n",miejsce_zerowe_newton(1,A,B,K));

    // pkt 2
    // double d0;
    // double best_c0;
    // for (double c0 = 2.0; c0 < 10.0; c0 += 0.01)
    // {
    //     d0=c0*0.5;
    //     d0 = miejsce_zerowe_newton(d0, c0, c0, K);
    //     printf("%.4f\t%.4f\n", d0, c0 - d0);
    //     if (compare_two_doubles(d0, c0 - d0))
    //         best_c0 = c0;
    // }
    // printf("\n%f\n", best_c0); // pkt 3

    // pkt 4
     miejsce_zerowe_newton_z_wypisywaniem(1,A,B,K);
     miejsce_zerowe_newton_z_wypisywaniem(1.5,A,B,K);
     miejsce_zerowe_newton_z_wypisywaniem(1.99,A,B,K);

    return 0;
}
