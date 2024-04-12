#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>

int save_to_file(gsl_matrix *A, int n, const char *name)
{
  FILE *f;
  f = fopen(name, "w");
  if (f == NULL)
  {
    perror("Nie udalo sie otworzyc pliku do zapisu");
    return 0;
  }

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      fprintf(f, "% -g\t", gsl_matrix_get(A, i, j));
    fprintf(f, "\n");
  }
  fflush(f);
  fclose(f);
  return 1;
}

gsl_matrix *calloc_unit_matrix(int n)
{
  gsl_matrix *A = gsl_matrix_calloc(n, n);
  for (int i = 0; i < n; i++)
    gsl_matrix_set(A, i, i, 0);
  return A;
}

int main()
{
  int n = 4; // Liczba wierszy

  gsl_matrix *A = gsl_matrix_calloc(n, n);
  gsl_permutation *p = gsl_permutation_calloc(n);

  // zapisanie wartości do macierzy
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
      double value = 1.0 / (double)(i + j + 2);
      gsl_matrix_set(A, i, j, value);
    }

  gsl_matrix *Anew = gsl_matrix_calloc(n, n);
  gsl_matrix_memcpy(Anew, A);
  int signum = 1;
  gsl_linalg_LU_decomp(A, p, &signum);

  // ZAD 2 diagonale macierzy
  FILE *f;
  f = fopen("diagonale dla U i wyznacznik.dat", "w");
  if (f == NULL)
  {
    perror("Nie udalo sie otworzyc pliku do zapisu");
    return 0;
  }
  double determinant_of_the_matrix_A = 1;

  for (int i = 0; i < n; i++)
  {
    fprintf(f, "%g\t", gsl_matrix_get(A, i, i));
    fprintf(f, "\n");

    determinant_of_the_matrix_A *= gsl_matrix_get(A, i, i);
  }
  fprintf(f, "\ndet(A) = %e", determinant_of_the_matrix_A);
  fflush(f);
  fclose(f);

  // ZAD 3 macierz odwrotna
  gsl_vector *b = gsl_vector_calloc(n);
  gsl_vector *x = gsl_vector_calloc(n);
  gsl_matrix *X = gsl_matrix_calloc(n, n); // macierz odwrotna
  for (int i = 0; i < n; i++)
  {
    // wektor jednostkowy
    for (int j = 0; j < n; j++)
    {
      if (i == j)
        gsl_vector_set(b, j, 1);
      else
        gsl_vector_set(b, j, 0);
    }

    gsl_linalg_LU_solve(A, p, b, x);

    // zapis x do X
    for (int j = 0; j < n; j++)
      gsl_matrix_set(X, i, j, gsl_vector_get(x, j));
  }
  save_to_file(X, n, "Macierz odwrotna.dat");

  // ZAD 4
  double value = 0.0;
  gsl_matrix *C = gsl_matrix_calloc(n, n);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int k = 0; k < n; k++)
      {
        value += gsl_matrix_get(Anew, i, k) * gsl_matrix_get(X, k, j);
      }
      gsl_matrix_set(C, i, j, value);
      value = 0.0;
    }
  }
  save_to_file(C, n, "iloczyn macierzy A i odwrotnej do niej.dat");

  // ZAD 5 wskaźnik uwarunkowania
  double a_norm = gsl_matrix_max(Anew), a_odwr_norm = gsl_matrix_max(X);
  double wsk_uwarun = a_norm * a_odwr_norm;

  f = fopen("wskaznik uwarunkowania macierzy.dat", "w");
  if (f == NULL)
  {
    perror("Nie udalo sie otworzyc pliku do zapisu");
    return 0;
  }
  fprintf(f, "κ(A) = %f", wsk_uwarun);
  fflush(f);
  fclose(f);

  return 0;
}