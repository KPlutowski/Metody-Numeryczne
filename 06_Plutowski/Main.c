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

double miejsce_zerowe_newton(const double D_initial, const double A, const double B, const double K)
{
    // double D=D_initial;
    // double delta;
    // do
    // {
    //     delta = f(D, A, B, K) / pochodna(D, A, B, K);
    //     D -= delta;
    // } while (fabs(delta) > 1e-6);

    // return D;
}
bool compare_two_doubles(const double A, const double B)
{
    return fabs(A - B) < 1e-6;
}

double fun(const double x, const size_t choice)
{
    switch (choice)
    {
    case 1:
        return exp(-1 * x * x);
        break;
    case 2:
        return x < 0 ? -1 : 1;
        break;
    case 3:
        return cos(2 * x);
        break;
    default:
        return 0;
        break;
    }
}

double lagrange_interpolation(const gsl_vector *x, const gsl_vector *y, const double val)
{
    double result = 0;
    for (size_t j = 0; j < x->size; ++j)
    {
        double term = gsl_vector_get(y, j);
        for (size_t k = 0; k < x->size; ++k)
        {
            if (j != k)
            {
                term *= (val - gsl_vector_get(x, k)) / (gsl_vector_get(x, j) - gsl_vector_get(x, k));
            }
        }
        result += term;
    }

    return result;
}

int save_3_value_to_file(const double value1,const double value2,const double value3, const char *name)
{
    FILE *f;
    f = fopen(name, "a");
    if (f == NULL)
    {
        perror("Nie udalo sie otworzyc pliku do zapisu");
        return 0;
    }

    fprintf(f, "%f\t%f\t%f\n", value1,value2,value3);

    fflush(f);
    fclose(f);

    return 1;
}

int main()
{
    const double start_x = -5.0;                            // Początek przedziału interpolacji
    const double end_x = 5.0;                               // Koniec przedziału interpolacji
    const double step = 0.1;                                // Krok dla generowania punktów na wykresie
    const double step_count = (end_x - start_x) / step + 1; // Ilość kroków

    const size_t num_points = 20; // Liczba węzłów interpolacji
    const size_t choice = 3;
        // 1 = exp(-1 * x * x);
        // 2 = x < 0 ? -1 : 1;
        // 3 = cos(2 * x);

    gsl_vector *x_values = gsl_vector_alloc(num_points); // Wektor położeń węzłów x
    gsl_vector *y_values = gsl_vector_alloc(num_points); // Wektor wartości funkcji w węzłach y

    gsl_vector *x_values_opti = gsl_vector_alloc(num_points); // Wektor położeń węzłów x
    gsl_vector *y_values_opti = gsl_vector_alloc(num_points); // Wektor wartości funkcji w węzłach y

    // Wypełnianie wektorów x i y węzłami interpolacji
    for (double i = 0; i < num_points; ++i)
    {
        double x = start_x + (end_x - start_x) * i / (num_points - 1);
        double x_opti = ((end_x - start_x) * cos((M_PI * (2 * (i) + 1)) / (2 * (num_points - 1) + 2)) + (end_x + start_x)) / 2;
        // printf("%f\n",x);

        gsl_vector_set(x_values, i, x);
        gsl_vector_set(y_values, i, fun(x, choice));

        gsl_vector_set(x_values_opti, i, x_opti);
        gsl_vector_set(y_values_opti, i, fun(x_opti, choice));
    }

    for (size_t i = 0; i < step_count; ++i)
    {
        save_3_value_to_file(fun(start_x + i * step, choice), lagrange_interpolation(x_values, y_values, start_x + i * step), lagrange_interpolation(x_values_opti, y_values_opti, start_x + i * step),"dane.txt");
    }

    return 0;
}
