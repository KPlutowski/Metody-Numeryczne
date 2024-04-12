#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>

int save_vector_to_file(gsl_vector *b, int x, char *name)
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
        fprintf(f, "%5g\n", gsl_vector_get(b, i));
    }
    fprintf(f, "\n");
    fflush(f);
    fclose(f);
    return 1;
}

int save_matrix_to_file(gsl_matrix *A, int x, int y, char *name)
{
    gsl_vector *vec = gsl_vector_calloc(x);
    for (int i = 0; i < y; i++)
    {
        gsl_matrix_get_row(vec, A, i);
        save_vector_to_file(vec, y, name);
    }
    return 1;
}

gsl_matrix *calloc_unit_matrix(int n)
{
    gsl_matrix *A = gsl_matrix_calloc(n, n);
    for (int i = 0; i < n; i++)
        gsl_matrix_set(A, i, i, 0);
    return A;
}

gsl_vector *gauss_jordan(gsl_matrix *A_org, gsl_vector *b_org, const int n)
{
    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_vector *b = gsl_vector_calloc(n);
    gsl_matrix_memcpy(A, A_org);
    gsl_vector_memcpy(b, b_org);

    for (int i = 0; i < n; i++)
    {
        // dzielenie przez a11
        double a = gsl_matrix_get(A, i, i);

        gsl_vector_set(b, i, gsl_vector_get(b, i) / a);
        for (int j = i; j < n; j++)
            gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) / a);

        for (int y = i + 1; y < n; y++)
        {
            gsl_vector_set(b, y, gsl_vector_get(b, y) - gsl_matrix_get(A, y, i) * gsl_vector_get(b, i));
            for (int x = n - 1; x >= 0; x--)
                gsl_matrix_set(A, y, x, gsl_matrix_get(A, y, x) - gsl_matrix_get(A, y, i) * gsl_matrix_get(A, i, x));
        }
    }
    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = 0; j < i; j++)
        {
            gsl_vector_set(b, j, gsl_vector_get(b, j) - gsl_matrix_get(A, j, i) * gsl_vector_get(b, i));
            gsl_matrix_set(A, j, i, 0);
        }
    }
    return b;
}

double deviation(gsl_vector *c, gsl_vector *b, const int n)
{
    // odchylenie
    double deviation = 0;
    for (int i = 0; i < n; i++)
        deviation += pow(gsl_vector_get(c, i) - gsl_vector_get(b, i), 2);
    deviation = sqrt(deviation) / n;
    return deviation;
}

gsl_vector *mat_vec_multi(const gsl_matrix *A, const gsl_vector *x, const int n)
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

double vec_vec_multi(const gsl_vector *a, const gsl_vector *b, const int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += gsl_vector_get(a, i) * gsl_vector_get(b, i);
    }
    return sum;
}

const double lambda(const double x)
{
    const double d1 = 40;
    const double d2 = 30;
    const double d3 = 30;

    const double l1 = 0.3;
    const double l2 = 0.2;
    const double l3 = 0.1;

    if (x <= d1)
    {
        return l1;
    }
    else if (x <= d1 + d2)
    {
        return l2;
    }
    else
    {
        return l3;
    }
}

void fill_start_matrix(gsl_matrix *A, const size_t n)
{
    const double STEP = (100. / (n +1));
    const double dx05 = STEP / 2;

    gsl_matrix_set(A, 0, 0, -lambda(dx05) - lambda(STEP + dx05));
    gsl_matrix_set(A, 0, 1, lambda(STEP + dx05));

    for (size_t i = 1; i < n - 1; i++)
    {
        gsl_matrix_set(A, i, i - 1, lambda(STEP * i + dx05));
        gsl_matrix_set(A, i, i, -lambda(STEP * i + dx05) - lambda(STEP * (i + 1) + dx05));
        gsl_matrix_set(A, i, i + 1, lambda(STEP * (i +1) + dx05));
    }
    gsl_matrix_set(A, n - 1, n - 2, lambda(STEP * n + dx05));
    gsl_matrix_set(A, n - 1, n - 1, -lambda(STEP * n + dx05) - lambda(STEP * (n + 1) + dx05));
}

void fill_start_wyniki(gsl_vector *b, const size_t n)
{
    const double TL = 1000;
    const double TR = 100;
    const double STEP = (100. / (n + 1));
    const double dx05 = STEP / 2;
    gsl_vector_set(b, 0, -lambda(dx05) * TL);
    for (size_t i = 1; i < n - 1; i++)
    {
        gsl_vector_set(b, i, 0);
    }
    gsl_vector_set(b, n - 1, -lambda(STEP * n + dx05) * TR);
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

void calculate(gsl_matrix *A, gsl_vector *t, gsl_vector *b, const size_t n)
{
    const double STEP = (100. / (n + 1));
    gsl_vector *rk = gsl_vector_calloc(n);

    gsl_vector *At = gsl_vector_calloc(n);
    double alfa = 0;

    double iloczyn_skalarny = 0;
    double mianownik = 0;
    double norma = 0;
    for (int k = 0; k < 60000; k++)
    {
        At = mat_vec_multi(A, t, n);
        rk = vec_minus_vec(b, At, n);
        iloczyn_skalarny = vec_vec_multi(rk, rk, n); // licznik
        if (sqrt(iloczyn_skalarny) < 0.000001)
        {
            return;
        }
        alfa = iloczyn_skalarny / vec_vec_multi(mat_vec_multi(A, rk, n), rk, n);
        for (int i = 0; i < n; i++)
        {
            gsl_vector_set(t, i, gsl_vector_get(t, i) + alfa * gsl_vector_get(rk, i));
        }
        //
        norma = 0;
        for (int i = 0; i < n; i++)
            norma += gsl_vector_get(t, i) * gsl_vector_get(t, i);
        // printf("%d\t%e\t%e\n", x, iloczyn_skal, sqrt(norma));

        FILE *f;
        f = fopen("liczby do sprawozdania.dat", "a");

        fprintf(f, "%d\t%e\t%e\n", k, log(sqrt(iloczyn_skalarny)), sqrt(norma));

        fflush(f);
        fclose(f);
    }
}

int main()
{
    const int n = 99;

    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_vector *b = gsl_vector_calloc(n);
    gsl_vector *t = gsl_vector_calloc(n);

    fill_start_wyniki(b, n);

    fill_start_matrix(A, n);

    calculate(A, t, b, n);

    save_vector_to_file(b, n, "b.dat");
    save_matrix_to_file(A,n,n,"A.dat");
    save_vector_to_file(t, n, "vector t.dat");

    return 0;
}
