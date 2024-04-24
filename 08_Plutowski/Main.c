#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// Obliczanie gęstości powietrza na danej wysokości
double ro(double h, double T, double P0)
{
    return 0.0289644 / (8.31446 * (T - 0.0065 * h)) * P0 * pow((1 - 0.0065 * h / T), 5.255);
}
// Obliczanie oporu areodynamicznego
double F_drag(double v, double Cx, double Sm, double h, double T, double P0)
{
    return Cx * ro(h, T, P0) * Sm * pow(v, 2) / 2;
}

int wziuuum(double alfa)
{
    // Zmienne
    double x, y, v, Time;
    double pi, g, Cx, diam, Sm, M, T, P0;
    double dt, v_tmp;

    // Stałe
    pi = 3.14159265358979; // Pi, tak +-
    g = 9.807;             // Przyśpieszenie ziemskie[m / s2]

    // Parametry pocisku 155mm ERFB-BB i Armatochaubicy Krab
    //  T. Kuśnierz, METODA WYZNACZENIA WARTOŚCI WSPÓŁCZYNNIKA AERODYNAMICZNEGO
    //  POCISKÓW STABILIZOWANYCH OBROTOWO, WITU
    Cx = 0.187; // Współczynnik tarcia pocisku

    diam = 155;                          // Kaliber pocisku [mm]
    Sm = pi * pow((diam / 2 / 1000), 2); // Powierzchnia przekroju[m2]

    M = 48;    // Masa pocisku [kg]
    v = 958.0; // Początkowa prędkość pocisku [m/s] 334-968

    // Parametry środowiskowe
    T = 273;     //  Temperatura [K]
    P0 = 101300; //  Ciśnienie znormalizowane [Pa]

    // Przeliczanie kątu wystrzału [deg --> rad]
    alfa = alfa / 180 * pi; //  [deg --> rad]

    // Parametry obliczeniowe
    dt = 0.01; /*  krok czasowy [s] (Jeżeli program będzie wykonywał się zbyt długo,
               zrezygnuj z dynamicznego obliczania kroku czasowago na początku pętli */

    //  Inicjalizacja wartościami początkowymi
    x = 0.0;    // współrzędna pozioma
    y = 1.0E-5; //  wsp pionowa
    Time = 0.0; // Czas lotu pocisku

    // Obliczanie trajektorii pocisku w czasie
    while (y > 0.0)
    {
        // dt=1.0/v;    // Dynamiczny krok czasowy, dla zwiększenia precyzji
        //   Wykorzystujemy metodę Eulera do obliczenia zmiany położenia,
        //   szybkości i kąta nachylenia stycznej do toru ruchu pocisku
        x += v * cos(alfa) * dt; // Wsp. pozioma
        y += v * sin(alfa) * dt; // Wsp. pionowa (wysokość npm.)

        // Prędkość tymczasowa
        v_tmp = v + (-F_drag(v, Cx, Sm, y, T, P0) / M - g * sin(alfa)) * dt;
        alfa -= g * cos(alfa) / v * dt; // Teraz alfa jest nachyleniem toru lotu pocisku
        v = v_tmp;
        Time += dt;
    }
    // Funkcja zwraca dystans na jaki doleciał pocisk [m]
    return x;
}




double golden_section_search(int (*func)(double), const double a, const double b, const double epsilon)
{
    const double golden_ratio = (sqrt(5) - 1) / 2;
    const double lambda1 = golden_ratio * golden_ratio;
    const double lambda2 = golden_ratio;

    double xa = a;
    double xb = b;
    double x1 = xa + lambda1 * (xb - xa);
    double x2 = xa + lambda2 * (xb - xa);

    while (fabs(xb - xa) > epsilon)
    {
        if (-func(x1) > -func(x2))
            xa = x1;
        else
            xb = x2;
        x1 = xa + lambda1 * (xb - xa);
        x2 = xa + lambda2 * (xb - xa);
    }
    return (xa + xb) / 2;
}

void print_trajectory_data(double alfa)
{
    // Zmienne
    double x, y, v, Time;
    double pi, g, Cx, diam, Sm, M, T, P0;
    double dt, v_tmp;

    // Stałe
    pi = 3.14159265358979; // Pi, tak +-
    g = 9.807;             // Przyśpieszenie ziemskie[m / s2]

    // Parametry pocisku 155mm ERFB-BB i Armatochaubicy Krab
    //  T. Kuśnierz, METODA WYZNACZENIA WARTOŚCI WSPÓŁCZYNNIKA AERODYNAMICZNEGO
    //  POCISKÓW STABILIZOWANYCH OBROTOWO, WITU
    Cx = 0.187; // Współczynnik tarcia pocisku

    diam = 155;                          // Kaliber pocisku [mm]
    Sm = pi * pow((diam / 2 / 1000), 2); // Powierzchnia przekroju[m2]

    M = 48;    // Masa pocisku [kg]
    v = 958.0; // Początkowa prędkość pocisku [m/s] 334-968

    // Parametry środowiskowe
    T = 273;     //  Temperatura [K]
    P0 = 101300; //  Ciśnienie znormalizowane [Pa]

    // Przeliczanie kątu wystrzału [deg --> rad]
    alfa = alfa / 180 * pi; //  [deg --> rad]

    // Parametry obliczeniowe
    dt = 0.01; /*  krok czasowy [s] (Jeżeli program będzie wykonywał się zbyt długo,
               zrezygnuj z dynamicznego obliczania kroku czasowago na początku pętli */

    //  Inicjalizacja wartościami początkowymi
    x = 0.0;    // współrzędna pozioma
    y = 1.0E-5; //  wsp pionowa
    Time = 0.0; // Czas lotu pocisku

    FILE *file = fopen("trajectory_data.txt", "w");
    if (file == NULL)
    {
        printf("Błąd otwierania pliku.");
        return;
    }
    fprintf(file, "Czas[s]\tPozycja x[m]\tPozycja y[m]\n");

    // Obliczanie trajektorii pocisku w czasie
    while (y > 0.0)
    {
        // dt=1.0/v;    // Dynamiczny krok czasowy, dla zwiększenia precyzji
        //   Wykorzystujemy metodę Eulera do obliczenia zmiany położenia,
        //   szybkości i kąta nachylenia stycznej do toru ruchu pocisku
        x += v * cos(alfa) * dt; // Wsp. pozioma
        y += v * sin(alfa) * dt; // Wsp. pionowa (wysokość npm.)

        fprintf(file, "%.2f\t%.2f\t%.2f\n", Time, x, y); // zaips do pliku

        // Prędkość tymczasowa
        v_tmp = v + (-F_drag(v, Cx, Sm, y, T, P0) / M - g * sin(alfa)) * dt;
        alfa -= g * cos(alfa) / v * dt; // Teraz alfa jest nachyleniem toru lotu pocisku
        v = v_tmp;
        Time += dt;
    }
    fclose(file);
}

// Funkcja zwracająca logarytm modułu różnicy między rozwiązaniem dokładnym a przybliżonym
void golden_section_log_diff(int (*func)(double), const double a, const double b, const double exact_solution, const double epsilon)
{
    const double golden_ratio = (sqrt(5) - 1) / 2;
    const double lambda1 = golden_ratio * golden_ratio;
    const double lambda2 = golden_ratio;

    double xa = a;
    double xb = b;
    double x1 = xa + lambda1 * (xb - xa);
    double x2 = xa + lambda2 * (xb - xa);



    FILE *file = fopen("logarytm_modułu_różnicy.txt", "w");
    if (file == NULL)
    {
        printf("Błąd otwierania pliku.");
        return ;
    }
    
    fprintf(file, "iteracja\tlog(fabs())\n");
    int iteration = 0;
    while (fabs(xb - xa) > epsilon)
    {
        if (-func(x1) > -func(x2))
            xa = x1;
        else
            xb = x2;
        x1 = xa + lambda1 * (xb - xa);
        x2 = xa + lambda2 * (xb - xa);

        fprintf(file,"%d\t%.15f\n", iteration, log(fabs(exact_solution - (xa + xb) / 2.0)));
        // printf("%f\n", (xa + xb) / 2);

        iteration++;
    }

}


int main()
{
    const double epsilon = 1e-6;
    const double alfa_start = 20;
    const double alfa_end = 70;

    double max_range_angle = golden_section_search(wziuuum, alfa_start, alfa_end, epsilon);
    double max_range = wziuuum(max_range_angle);

    printf("%f\n", max_range_angle);

    // zasięg działa w funkcji kąta ostrzału z zaznaczonym maximum
    FILE *file = fopen("zasięg_działa_w_funkcji_kąta.txt", "w");
    if (file == NULL)
    {
        printf("Błąd otwierania pliku.");
        return 1;
    }

    fprintf(file, "kąt [stopnie]\tzasięg [m] \n");
    for (int i = alfa_start; i < alfa_end; i++)
        fprintf(file, "%.2d\t%.2d\n", i, wziuuum(i)); // zaips do pliku

    fclose(file);

    print_trajectory_data(max_range_angle); // trajektorię lotu pocisku dla maksymalnego zasięgu


    golden_section_log_diff(wziuuum, alfa_start, alfa_end,50.8057242086663, epsilon); // logarytmu modułu różnicy rozwiązania dokładnego



    // 4
    // double alfa_300m = golden_section_search(target_function, 0, M_PI / 2, 0.0001);
    // printf("Kąt ostrzału dla celu na wysokości 300 m n.p.m.: %f stopni\n",alfa_300m * 180 / M_PI);

    return 0;
}
