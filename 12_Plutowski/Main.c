#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>

double trapezoidal(const gsl_matrix *M)
{
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


double simpson(const gsl_matrix *M)
{
    double sum=0;
    const double h = gsl_matrix_get(M,1,0)-gsl_matrix_get(M,0,0);

    for (size_t i = 1; i < M->size1-1; i+=2)
    {
        const double s = h/3.0*(gsl_matrix_get(M,i-1,1)+4.0*gsl_matrix_get(M,i,1)+gsl_matrix_get(M,i+1,1));
        sum+=s;
    }

    return sum;
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
    const double start_x = 0;
    const double end_x = M_PI;

    for(int n = 11;n<=201;n+=2)
    {
        gsl_matrix *M1 = gsl_matrix_calloc(n,2);

        for(int i = 0; i < M1->size1; ++i)
        {
            const double x = start_x + (end_x - start_x) * i / n;

            gsl_matrix_set(M1, i, 0,x);
            gsl_matrix_set(M1, i, 1, sin(x));
        }

        double bladTrapez= log(fabs(2-trapezoidal(M1)));
        double bladSimpson= log(fabs(2-simpson(M1)));

        save_3_value_to_file(n,bladTrapez,bladSimpson,"sin(x).csv");
        gsl_matrix_free(M1);
    }


    for(int n = 11;n<=201;n+=2)
    {
        gsl_matrix *M1 = gsl_matrix_calloc(n,2);

        for(int i = 0; i < M1->size1; ++i)
        {
            const double x = start_x + (end_x - start_x) * i / n;

            gsl_matrix_set(M1, i, 0,x);
            gsl_matrix_set(M1, i, 1, x*sin(x*5.0));
        }
        double bladTrapez= log(fabs(M_PI/5.0-trapezoidal(M1)));
        double bladSimpson= log(fabs(M_PI/5.0-simpson(M1)));

        save_3_value_to_file(n,bladTrapez,bladSimpson,"xsin(x5).csv");
        gsl_matrix_free(M1);
    }




    return 0;
}
