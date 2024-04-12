#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>




int save_vector_to_file(const gsl_vector *b, const size_t x,const char *name)
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
int save_value_to_file(const double value,const char *name)
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
double deviation(const gsl_vector *c,const gsl_vector *b, const size_t n)
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
gsl_matrix* matrix_multiplication(const gsl_matrix *A, const gsl_matrix *B, const size_t n)
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


void fill_start_matrix(gsl_matrix* A,const size_t n)
{
    for(size_t i=0;i<n;i++)
        for(size_t j=0;j<n;j++)
        {
            gsl_matrix_set(A, i, j, pow(2+abs(i-j),-abs(i-j)/2.0f));
        }
}


gsl_matrix *tensor_product(const gsl_vector *a, const size_t n)
{
    gsl_matrix *B = gsl_matrix_calloc(n, n);
    for(size_t i=0;i<n;i++)
        for(size_t j=0;j<n;j++)
        {
            gsl_matrix_set(B,i,j,gsl_vector_get(a,i)*gsl_vector_get(a,j));
        }

    return B;
}

void calculate(gsl_matrix *A,const size_t n,const size_t M)
{
    // 1
    gsl_matrix *W0 = gsl_matrix_calloc(n, n);
    gsl_matrix_memcpy(W0, A);

    // 2
    gsl_vector *x0 = gsl_vector_calloc(n);
    gsl_matrix *B = gsl_matrix_calloc(n, n);

    gsl_matrix *VecW = gsl_matrix_calloc(n, n); //macierz wektorów własnych
    
    for(size_t k=0;k<n;k++)
    {
        for(size_t i=0;i<n;i++)
            gsl_vector_set(x0, i, 1);


        double lambda;
        gsl_vector *x_new = gsl_vector_calloc(n);

        // PKT 3 powtarzamy M razy
        for(size_t i=0;i<M;i++)
        {
            gsl_vector_memcpy(x_new, x0);
            

            // 3.1
            x_new=matrix_vector_multiplication(W0,x0,n);

            // 3.2
            lambda = vector_multiplication(x_new, x0, n) / vector_multiplication(x0, x0, n);

            // 3.3
            for (int i = 0; i < n; i++)
            {
                gsl_vector_set(x0, i, gsl_vector_get(x_new, i) / sqrt(vector_multiplication(x_new, x_new, n)));
            }
            save_value_to_file(lambda,"kolejne przyblizenia wartości wlasnych.dat");
        }

        for(size_t i=0;i<n;i++)
        {
            gsl_matrix_set(VecW, i, k, gsl_vector_get(x0, i));
        }
        

        save_value_to_file(lambda,"wartosci wlasne.dat");

        B=tensor_product(x0,n);
        
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                gsl_matrix_set(W0, i, j, gsl_matrix_get(W0, i, j)-lambda* gsl_matrix_get(B, i, j));
            }
        }
    }

    

    gsl_matrix *D = gsl_matrix_calloc(n, n);
    gsl_matrix *D1 = gsl_matrix_calloc(n, n);
    gsl_matrix_transpose(VecW);
    D1=matrix_multiplication(VecW, A, n);
    gsl_matrix_transpose(VecW);
    D=matrix_multiplication(D1, VecW, n);
   

    save_matrix_to_file(VecW,n,n,"wektory wlasne.dat");
    save_matrix_to_file(D,n,n,"macierz D.dat");
    save_matrix_to_file(W0,n,n,"macierz W0.dat");

}

int main()
{
    const size_t M=12;
    const size_t n = 7;
    gsl_matrix *A = gsl_matrix_calloc(n, n);


    fill_start_matrix(A,n);
    // save_matrix_to_file(A,n,n,"macierz A");

    calculate(A,n,M);





    return 0;
}
