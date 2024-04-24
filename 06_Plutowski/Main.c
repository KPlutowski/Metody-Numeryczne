#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

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
