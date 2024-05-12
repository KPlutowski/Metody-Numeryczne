#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>
#include <time.h>
int save_vector_to_file(const gsl_vector *b, const char *name)
{
    FILE *f;
    f = fopen(name, "a");
    if (f == NULL)
    {
        perror("Nie udalo sie otworzyc pliku do zapisu");
        return 0;
    }

    for (int i = 0; i < b->size; i++)
    {
        fprintf(f, "%5g\t", gsl_vector_get(b, i));
    }
    fprintf(f, "\n");
    fflush(f);
    fclose(f);
    return 1;
}
int save_matrix_to_file(const gsl_matrix *A, const char *name)
{
    gsl_vector *vec = gsl_vector_calloc(A->size2);
    
    for (int i = 0; i < A->size1; i++)
    {
        gsl_matrix_get_row(vec, A, i);
        save_vector_to_file(vec, name);
    }
    return 1;
}


gsl_vector *matrix_vector_multiplication(const gsl_matrix *A, const gsl_vector *x)
{

    gsl_vector *c = gsl_vector_calloc(A->size1);
    
    for (int i = 0; i < A->size1; i++)
    {
        double sum = 0;
        for (int j = 0; j < A->size2; j++)
            sum += gsl_matrix_get(A, i, j) * gsl_vector_get(x, j);
        gsl_vector_set(c, i, sum);
    }
    return c;
}

double vector_multiplication(const gsl_vector *a, const gsl_vector *b)
{
    double sum = 0;
    for (int i = 0; i < a->size; i++)
    {
        sum += gsl_vector_get(a, i) * gsl_vector_get(b, i);
    }
    return sum;
}

gsl_matrix *matrix_multiplication(const gsl_matrix *A, const gsl_matrix *B)
{
    //  A*B = C

    gsl_matrix *C = gsl_matrix_calloc(A->size1, B->size2);
    
    for (int i = 0; i < A->size1; i++)
    {
        for (int j = 0; j < B->size2; j++)
        {
            double sum = 0;
            for (int k = 0; k < A->size2; k++)
            {
                sum += gsl_matrix_get(A, i, k) * gsl_matrix_get(B, k, j);
            }
            gsl_matrix_set(C, i, j, sum);
        }
    }
    return C;
}

gsl_matrix *matrix_transpose( const gsl_matrix *A)
{
    gsl_matrix *B1 = gsl_matrix_alloc(A->size2, A->size1);

    for (int i = 0; i < A->size1; i++)
    {
        for (int j = 0; j < A->size2; j++)
        {
            gsl_matrix_set(B1, j, i, gsl_matrix_get(A, i, j));
        }
    }
    return B1;
}

double fun(const double x){
    const double delta= 5.0*(((double)rand() / (double)RAND_MAX)-0.5);
    printf("%f\n",delta);
    return (-(0.25*x*x*x)-(0.5*x*x)+(5*x)+5.0+delta);
}


int main()
{
    srand( time( NULL ) );
    

    const size_t n = 101;
    const size_t m = 4;
    const double start_x = -5;
    const double end_x = 5;

    gsl_vector *y = gsl_vector_alloc(n); 
    gsl_vector *x = gsl_vector_alloc(n); 
    gsl_matrix *X = gsl_matrix_alloc(n,m); 


    for(int i = 0; i < n; ++i)
    {
        const double tmp = start_x + (end_x - start_x) * i / (n - 1);

        for (int j = 0 ; j < m ; j++)
            gsl_matrix_set(X, i, j, pow(tmp, j));
        gsl_vector_set(y, i, fun(tmp));
    
        gsl_vector_set(x, i, tmp);
    }

    gsl_matrix *Xtrans=matrix_transpose(X);
    gsl_matrix *D = matrix_multiplication(Xtrans,X);

    save_vector_to_file(x,"x");
    save_vector_to_file(y,"y");
    save_matrix_to_file(X,"X");
    save_matrix_to_file(Xtrans,"Xtrans");
    save_matrix_to_file(D,"D");


    gsl_permutation *p = gsl_permutation_calloc(D->size1);
    gsl_vector *r = matrix_vector_multiplication(Xtrans,y);
    gsl_vector *b = gsl_vector_alloc(m);
    int signum = 1;
    gsl_linalg_LU_decomp(D, p, &signum);
    gsl_linalg_LU_solve(D, p, r, b);

    save_vector_to_file(b,"b");
    save_vector_to_file(r,"r");


    gsl_matrix *wyniki = gsl_matrix_alloc(n,3);
    const double a_ = gsl_vector_get(b,3);
    const double b_ = gsl_vector_get(b,2);
    const double c_ = gsl_vector_get(b,1);
    const double d_ = gsl_vector_get(b,0);
    for(int i = 0; i < n; ++i)
    {
        // y = a*x^3 + b*x^2 +c*x + d
        const double x_= gsl_vector_get(x,i);
        const double y_ = a_*pow(x_,3)+ b_*pow(x_,2)+c_*x_+d_;

        gsl_matrix_set(wyniki, i, 0, x_);
        gsl_matrix_set(wyniki, i, 1, y_);
        gsl_matrix_set(wyniki, i, 2, fun(x_));

    }
    save_matrix_to_file(wyniki,"wyniki.csv");


    return 0;
}
