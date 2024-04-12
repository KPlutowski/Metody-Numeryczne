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
    fprintf(f, "%5g\t", gsl_vector_get(b, i));
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
  {
    deviation += pow(gsl_vector_get(c, i) - gsl_vector_get(b, i), 2);
  }
  deviation = sqrt(deviation) / n;
  return deviation;
}
gsl_vector *matrix_vector_multiply(gsl_matrix *A, gsl_vector *x, const int n)
{
  gsl_vector *c = gsl_vector_calloc(n);
  double sum = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      sum += gsl_matrix_get(A, i, j) * gsl_vector_get(x, j);
    }
    gsl_vector_set(c, i, sum);
    sum = 0;
  }
  return c;
}

int main()
{
  const int n = 5;

  double tmp_matrix[5][5] = {{2, 1, 6, 9, 10}, {2, 1, 6, 9, 10}, {1, 6, 6, 8, 6}, {5, 9, 10, 7, 10}, {3, 4, 9, 7, 9}};
  double tmp_vector[5] = {10, 2, 9, 9, 3};

  gsl_matrix *A = gsl_matrix_calloc(n, n);
  gsl_vector *b = gsl_vector_calloc(n);

  // wypelnanie macierzy
  for (int i = 0; i < n; i++)
  {
    gsl_vector_set(b, i, tmp_vector[i]);
    for (int j = 0; j < n; j++)
      gsl_matrix_set(A, i, j, tmp_matrix[i][j]);
  }

  gsl_vector *x = gsl_vector_calloc(n);
  gsl_vector *c = gsl_vector_calloc(n);

  FILE *f;
  f = fopen("data.dat", "w");
  if (f == NULL)
  {
    perror("Nie udalo sie otworzyc pliku do zapisu");
    return 0;
  }

  // algo na zajeciach
  for (double q = 0.01; q < 3; q += 0.01)
  {
    gsl_matrix_set(A, 0, 0, tmp_matrix[0][0] * q);
    x = gauss_jordan(A, b, n);

    c = matrix_vector_multiply(A, x, n);

    fprintf(f, "%.3f\t%.3e\n", q, deviation(c, b, n));
  }

  // file close
  fflush(f);
  fclose(f);
  return 0;
}