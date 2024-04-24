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
        fprintf(f, "%g\n", gsl_vector_get(b, i));
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
double fun_a(const double x)
{
    return exp(-1 * x * x);
}
double fun_b(const double x)
{
    return x < 0 ? -1 : 1;
}
double fun_c(const double x)
{
    return cos(2 * x);
}

double (*getFun(const size_t choice))(double)
{
    switch (choice)
    {
    case 1:
        return fun_a;
        break;
    case 2:
        return fun_b;
        break;
    case 3:
        return fun_c;
        break;
    default:
        return 0;
        break;
    }
}

double spline_interpolation(const double val, const gsl_vector *m, const gsl_vector *x, const gsl_vector *y, const gsl_vector *h)
{

    for (int i = 0; i < x->size - 1; ++i)
    {
        if (val >= gsl_vector_get(x, i) && val <= gsl_vector_get(x, i + 1))
        {
            double A1 = (gsl_vector_get(y, i + 1) - gsl_vector_get(y, i)) / gsl_vector_get(h, i + 1);
            double A2 = gsl_vector_get(h, i + 1) / 6.0 * (gsl_vector_get(m, i + 1) - gsl_vector_get(m, i));
            double A = A1 - A2;

            double B = gsl_vector_get(y, i) - gsl_vector_get(m, i) * pow(gsl_vector_get(h, i + 1), 2) / 6.0;

            double S1 = gsl_vector_get(m, i) * pow(gsl_vector_get(x, i + 1) - val, 3) / (6.0 * gsl_vector_get(h, i + 1));
            double S2 = gsl_vector_get(m, i + 1) * pow(val - gsl_vector_get(x, i), 3) / (6.0 * gsl_vector_get(h, i + 1));
            double S3 = A * (val - gsl_vector_get(x, i)) + B;
            return (S1 + S2 + S3);
        }
    }
    return 0.0;
}


double pochodna2(double val, double (*f)(double))
{
    const double hval = 0.01f;

    double wynik = (f(val - hval) - 2 * f(val) + f(val + hval)) / pow(hval, 2);
    return wynik;
}

gsl_vector *gauss_jordan(const gsl_matrix *A_org, const gsl_vector *b_org, const int n)
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

void fill_start_vectors_xyh(double (*f)(double), gsl_vector *x_values, gsl_vector *y_values, gsl_vector *h, double start_x, double end_x, double num_points)
{
    for (double i = 0; i < num_points; ++i)
    {
        double x = start_x + (end_x - start_x) * i / (num_points - 1);
        gsl_vector_set(x_values, i, x);
        // printf("%f\n",x);

        gsl_vector_set(y_values, i, f(x));
        // printf("%g\n",f(x));

        gsl_vector_set(h, i, (end_x - start_x) / (num_points - 1));
        // printf("%f\n",(end_x - start_x) / (num_points - 1));
    }
}
gsl_vector *calc_second_derivative(gsl_vector *y_values,gsl_vector *h, double num_points)
{
    gsl_matrix *A = gsl_matrix_alloc(num_points, num_points); // Macierz A

    gsl_matrix_set(A, 0, 0, 1);
    double lambda;
    for (size_t i = 1; i < num_points - 1; i++)
    {
        // lambda
        lambda=gsl_vector_get(h,i+1)/(gsl_vector_get(h,i)+gsl_vector_get(h,i+1));
        gsl_matrix_set(A, i, i + 1, lambda);

        // u
        gsl_matrix_set(A, i, i - 1, 1-lambda);

        // 2
        gsl_matrix_set(A, i, i, 2);

    }
    gsl_matrix_set(A, num_points - 1, num_points - 1, 1);
    // save_matrix_to_file(A,num_points,num_points,"A.txt");

    gsl_vector *d = gsl_vector_alloc(num_points); // Wektor wyrazów wolnych
    gsl_vector_set(d, 0, 0);
    for (size_t i = 1; i < num_points - 1; i++)
    {

        double _d = 6.0 / (gsl_vector_get(h,i)+gsl_vector_get(h,i+1)) * ((gsl_vector_get(y_values, i + 1) - gsl_vector_get(y_values, i)) / gsl_vector_get(h,i+1) - (gsl_vector_get(y_values, i) - gsl_vector_get(y_values, i - 1)) / gsl_vector_get(h,i));
        gsl_vector_set(d, i, _d);
    }
    gsl_vector_set(d, num_points - 1, 0);
    // save_vector_to_file(d,num_points,"d.txt");

    gsl_vector *m = gsl_vector_alloc(num_points); // drugie pochodne
    m = gauss_jordan(A, d, num_points);
    // save_vector_to_file(m,num_points,"m.txt");

    return m;
}
int save_3_value_to_file(const double value1, const double value2, const double value3, const char *name)
{
    FILE *f;
    f = fopen(name, "a");
    if (f == NULL)
    {
        perror("Nie udalo sie otworzyc pliku do zapisu");
        return 0;
    }

    fprintf(f, "%f\t%g\t%g\n", value1, value2, value3);

    fflush(f);
    fclose(f);

    return 1;
}

int main()
{
    const double start_x = -5.0;  // Początek przedziału interpolacji
    const double end_x = 5.0;     // Koniec przedziału interpolacji
    const size_t num_points = 20; // Liczba węzłów interpolacji
    const double step = 0.1;      // Krok dla generowania punktów na wykresie
    const size_t choice = 3;      // funkcja do wybrania
                // 1 = exp(-1 * x * x);
                // 2 = x < 0 ? -1 : 1;
                // 3 = cos(2 * x);

    const double step_count = (end_x - start_x) / step + 1; // Ilość kroków

    gsl_vector *x_values = gsl_vector_alloc(num_points); // Wektor położeń węzłów x
    gsl_vector *y_values = gsl_vector_alloc(num_points); // Wektor wartości funkcji w węzłach x
    gsl_vector *h = gsl_vector_alloc(num_points);        // odleglosci (stale)

    fill_start_vectors_xyh(getFun(choice), x_values, y_values, h, start_x, end_x, num_points);

    gsl_vector *m = gsl_vector_alloc(num_points);
    m = calc_second_derivative(y_values,h, num_points);
    save_vector_to_file(m, num_points, "pochodna z rownania macierzowego.txt");
    // POCHODNA z iloczynu różnicowego
    for (double i = 0; i < num_points; ++i)
    {
        double x = start_x + (end_x - start_x) * i / (num_points - 1);

        gsl_vector_set(m, i, pochodna2(x, getFun(choice)));
        // printf("%f\n", x);
    }
    save_vector_to_file(m, num_points, "pochodna z wzoru.txt");

    for (size_t i = 0; i < step_count; ++i)
    {
        save_3_value_to_file(start_x + i * step, getFun(choice)(start_x + i * step), spline_interpolation(start_x + i * step, m, x_values, y_values, h),"wyniki.txt");
        // printf("%.2f\t%g\t%g\n", start_x + i * step, getFun(choice)(start_x + i * step), spline_interpolation(start_x + i * step, m, x_values, y_values, h));
    }
    return 0;
}
