#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_errno.h>


#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])

void fft_file()
{
    const int k = 11;
    const int N =  pow(2, k);
    const int colum=3;
    double y0[N];
    double dane[2 * N];
    double modul[N];

    // Wczytywanie z pliku
    FILE *data = fopen("data.csv", "r");

    char ignore[1024];
    fgets(ignore, sizeof(ignore), data);

    double value[N][4];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            fscanf(data, "%lf", &value[i][j]);
            // printf("%lf\n",value[i][j]);
        }
    }
    fclose(data);

    for (int i = 0; i < N; i++)
    {
        y0[i] = value[i][colum];
        REAL(dane, i)= y0[i];   // Część real
        IMAG(dane, i) = 0.0;    // Część imag
    }

    // Obliczanie transformaty
    gsl_fft_complex_radix2_forward(dane, 1, N);

    // Obliczenie modułów oraz progu
    double max_modul = 0.0;
    for (int i = 0; i < N; i++)
    {
        // printf("%f\n",REAL(dane, i));
        modul[i] = sqrt(REAL(dane, i) * REAL(dane, i) + IMAG(dane, i) * IMAG(dane, i));
        if (modul[i] > max_modul)
        {
            max_modul = modul[i];
        }
    }

    // Dyskryminacja współczynników
    double threshold = max_modul / 2;
    for (int i = 0; i < N; i++)
    {
        if (modul[i] < threshold)
        {
            REAL(dane, i) = 0;
            IMAG(dane, i) = 0;
        }
    }

    // Obliczenie odwrotnej transformaty
    gsl_fft_complex_radix2_backward(dane, 1, N);

    // Normalizacja danych
    for (int i = 0; i < N; i++)
    {
        REAL(dane, i) /= N;
        IMAG(dane, i) /= N;
    }


    // Zapis wyników do plików
    FILE *file_denoised = fopen("denoised_signal_d.csv", "w");
    FILE *file_original = fopen("original_signal_d.csv", "w");
    FILE *modul_file = fopen("modul_file_d.csv", "w");

    for (int i = 0; i < N; i++)
    {
        fprintf(file_original, "%f\t%f\n",value[i][0], y0[i]);// dziala

        fprintf(file_denoised, "%f\t%f\n",value[i][0], REAL(dane, i));
        fprintf(modul_file, "%f\t%f\n",value[i][0], modul[i]);
    }

    fclose(file_denoised);
    fclose(file_original);
    fclose(modul_file);




    // Częstotliwości
    double Fs = 2048; // Zakładamy, że częstotliwość próbkowania jest 1000 Hz
    double widmo[N / 2 + 1];

    // Obliczanie widma
    for (int i = 0; i <= N / 2; i++)
    {
        widmo[i] = (Fs / N) * i;
    }

    // Zapis widma do pliku
    FILE *widmo_file = fopen("widmo_file_d.csv", "w");

    for (int i = 0; i <= N / 2; i++)
    {
        fprintf(widmo_file, "%f\t%f\n", widmo[i], modul[i]);
    }

    fclose(widmo_file);



}

void fft_generated(const int k)
{
    const int N =  pow(2, k);

    const double omega = 4 * M_PI / N;
    
    double y0[N];
    double dane[2 * N];
    double modul[N];
    srand(time(NULL));
    
    // Generowanie sygnału niezaszumionego
    for (int i = 0; i < N; i++)
    {
        y0[i] = sin(omega * i) + sin(2 * omega * i) + sin(3 * omega * i);
    }

    // Tworzenie sygnału zaszumionego
    for (int i = 0; i < N; i++)
    {
        double delta = 2 * ((double)rand() / RAND_MAX) - 1; // Szum w przedziale (-1,1]
        REAL(dane, i) = y0[i] + delta; // Część rzeczywista
        IMAG(dane, i) = 0.0; // Część urojona
    }

    // Obliczanie transformaty
    gsl_fft_complex_radix2_forward(dane, 1, N);

    // Obliczenie modułów oraz progu
    double max_modul = 0.0;
    for (int i = 0; i < N; i++)
    {
        modul[i] = sqrt(REAL(dane, i) * REAL(dane, i) + IMAG(dane, i) * IMAG(dane, i));
        if (modul[i] > max_modul)
        {
            max_modul = modul[i];
        }
    }
    
    FILE *real = fopen("real.csv", "w");
    FILE *imag = fopen("imag.csv", "w");
    for (int i = 0; i < N; i++)
    {
        fprintf(real, "%d\t%f\n",i, REAL(dane, i));
        fprintf(imag, "%d\t%f\n",i, IMAG(dane, i));
    }
    fclose(real);
    fclose(imag);    
    // Dyskryminacja współczynników
    double threshold = max_modul / 2;
    for (int i = 0; i < N; i++)
    {
        if (modul[i] < threshold)
        {
            REAL(dane, i) = 0;
            IMAG(dane, i) = 0;
        }
    }


    // Obliczenie odwrotnej transformaty
    gsl_fft_complex_radix2_backward(dane, 1, N);

    // Normalizacja danych
    for (int i = 0; i < N; i++)
    {
        REAL(dane, i) /= N;
        IMAG(dane, i) /= N;
    }

    // Zapis wyników do plików
    FILE *file_noisy = fopen("noisy_signal.csv", "w");
    FILE *file_denoised = fopen("denoised_signal.csv", "w");
    FILE *file_original = fopen("original_signal.csv", "w");
    FILE *modul_file = fopen("modul_file.csv", "w");

    for (int i = 0; i < N; i++)
    {
        fprintf(file_original, "%d\t%f\n",i, y0[i]);
        fprintf(file_noisy, "%d\t%f\n", i,y0[i] + ((double)rand() / RAND_MAX) * 2 - 1);
        fprintf(file_denoised, "%d\t%f\n",i, REAL(dane, i));
        fprintf(modul_file, "%d\t%f\n",i, modul[i]);
    }

    fclose(file_noisy);
    fclose(file_denoised);
    fclose(file_original);
    fclose(modul_file);
}

int main()
{
    // fft_generated(8);
    fft_file();


    return 0;
}
