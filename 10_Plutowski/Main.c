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


unsigned long U1_generate_next(unsigned long previous) {
    const int a=17;
    const unsigned long m=pow(2,13)-1;
    return (a * previous) % m;
}
double U1_normalize(unsigned long x) {
    const unsigned long m=pow(2,13)-1;
    return (double)x / (m - 1);
}

unsigned long U2_generate_next(unsigned long previous) {
    const int a=85;
    const unsigned long m=pow(2,13)-1;
    return (a * previous) % m;
}
double U2_normalize(unsigned long x) {
    const unsigned long m=pow(2,13)-1;
    return (double)x / (m - 1);
}

unsigned long U3_generate_next(unsigned long* previous) {
    const int a=1176;
    const int b=1476;
    const int c=1776;

    const unsigned long m=pow(2,32)-5;

    unsigned long returnx=(a * previous[2]+b* previous[1]+c * previous[0]) % m;

    previous[0]=previous[1];
    previous[1]=previous[2];
    previous[2]=returnx;

    return returnx;
}
double U3_normalize(unsigned long x) {
    const unsigned long m=pow(2,32)-5;
    return (double)x / (m - 1);
}

void monte()
{
    gsl_matrix *M = gsl_matrix_calloc(200,3);
    const double r = 1;
    const unsigned long seed=10;

    unsigned long *pre1 = malloc(sizeof(unsigned long)*3);
    pre1[0]=seed;
    pre1[1]=seed;
    pre1[2]=seed;
    double x=0,y=0;

    double odleglos=0;
    unsigned long suma=0;

    for(int i=1;i<20001;i++)
    {
        x = U3_normalize(U3_generate_next(pre1));
        y = U3_normalize(U3_generate_next(pre1));
    
        odleglos=sqrt(x*x+y*y);

        if(odleglos<r)
        {
            suma++;
        }

        if(i%100==0)
        {
            double pi=((double)suma/(double)i)*4.0;
            // printf("%d\n",i/100-1);

            gsl_matrix_set(M,i/100-1,0,i/100);

            gsl_matrix_set(M,i/100-1,1,log(fabs(pi-M_PI)));
            gsl_matrix_set(M,i/100-1,2,pi);


            // printf("iteracja: %d, pi: %f, suma: %ld, logarytm: %f\n",i,pi, suma, log(fabs(pi-M_PI)));
        }

    }
    save_matrix_to_file(M,"monte.csv");

}

int main()
{
    const unsigned long seed = 10;  // Można zmienić na inne ziarno

    gsl_matrix *M = gsl_matrix_calloc(2000,2);
    monte();

    unsigned long current = seed;
    for (int i = 0; i < 2000; i++) {
        current = U1_generate_next(current);
        double normalized_value = U1_normalize(current);

        gsl_matrix_set(M,i,0,i+1);
        gsl_matrix_set(M,i,1,current);
        // printf("Liczba pseudolosowa %d: %lu, Po normalizacji: %f\n", i + 1, current, normalized_value);
    }
    save_matrix_to_file(M,"U1.csv");


    current = seed;
    for (int i = 0; i < 2000; i++) {
        current = U2_generate_next(current);
        double normalized_value = U2_normalize(current);

        gsl_matrix_set(M,i,0,i+1);
        gsl_matrix_set(M,i,1,current);
        // printf("Liczba pseudolosowa %d: %lu, Po normalizacji: %f\n", i + 1, current, normalized_value);
    }
    save_matrix_to_file(M,"U2.csv");

    unsigned long *pre = malloc(sizeof(unsigned long)*3);
    pre[0]=seed;
    pre[1]=seed;
    pre[2]=seed;
    current=0;
    for (int i = 0; i < 2000; i++) {
        current = U3_generate_next(pre);
        double normalized_value = U3_normalize(current);

        gsl_matrix_set(M,i,0,i+1);
        gsl_matrix_set(M,i,1,current);
        // printf("Liczba pseudolosowa %d: %lu, Po normalizacji: %f\n", i + 1, current, normalized_value);
    }
    save_matrix_to_file(M,"U3.csv");

    return 0;
}
