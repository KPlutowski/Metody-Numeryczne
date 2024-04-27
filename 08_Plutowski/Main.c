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
double wziuuum(double alfa)
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
    // dt = 0.01;
    /*  krok czasowy [s] (Jeżeli program będzie wykonywał się zbyt długo,
   zrezygnuj z dynamicznego obliczania kroku czasowago na początku pętli */

    //  Inicjalizacja wartościami początkowymi
    x = 0.0;    // współrzędna pozioma
    y = 1.0E-5; //  wsp pionowa
    Time = 0.0; // Czas lotu pocisku

    // Obliczanie trajektorii pocisku w czasie
    while (y > 0.0)
    {
        dt = 1.0 / v;            // Dynamiczny krok czasowy, dla zwiększenia precyzji
                                 //  Wykorzystujemy metodę Eulera do obliczenia zmiany położenia,
                                 //  szybkości i kąta nachylenia stycznej do toru ruchu pocisku
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

// zad 1 metodę złotego podziału do poszukiwania minimum funkcji.
double golden_section_search(double (*func)(double), const double a, const double b, const double epsilon)
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
void golden_section_log_diff(double (*func)(double), const double a, const double b, const double exact_solution, const double epsilon)
{
    const double golden_ratio = (sqrt(5.0) - 1.0) / 2.0;
    // const double lambda1 = golden_ratio * golden_ratio;
    // const double lambda2 = golden_ratio;
    const double lambda1 = 1.0 / 3.0;
    const double lambda2 = 2.0 / 3.0;

    double xa = a;
    double xb = b;


    FILE *file = fopen("logarytm_modułu_różnicy.txt", "w");
    if (file == NULL)
    {
        printf("Błąd otwierania pliku.");
        return;
    }

    fprintf(file, "iteracja\t(xa + xb) / 2\tlog(fabs())\n");
    int iteration = 0;
    while (fabs(xb - xa) > epsilon)
    {
        double x1 = xa + lambda1 * (xb - xa);
        double x2 = xa + lambda2 * (xb - xa);
        if (-wziuuum(x1) > -wziuuum(x2))
            xa = x1;
        else
            xb = x2;

        fprintf(file, "%d\t%f\t%f\n", iteration, (xa + xb) / 2, log(fabs(exact_solution - (xa + xb) / 2)));
        iteration++;
    }
    fclose(file);
}

double wziuuum2(double alfa, const double target_distance, const double target_height)
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
    // dt = 0.01;
    /*  krok czasowy [s] (Jeżeli program będzie wykonywał się zbyt długo,
   zrezygnuj z dynamicznego obliczania kroku czasowago na początku pętli */

    //  Inicjalizacja wartościami początkowymi
    x = 0.0;    // współrzędna pozioma
    y = 1.0E-5; //  wsp pionowa
    Time = 0.0; // Czas lotu pocisku

    // Obliczanie trajektorii pocisku w czasie
    while ((y > target_height && x < target_distance + 15) || (y > 0.0 && x < target_distance))
    {
        dt = 1.0 / v;            // Dynamiczny krok czasowy, dla zwiększenia precyzji
                                 //  Wykorzystujemy metodę Eulera do obliczenia zmiany położenia,
                                 //  szybkości i kąta nachylenia stycznej do toru ruchu pocisku
        x += v * cos(alfa) * dt; // Wsp. pozioma
        y += v * sin(alfa) * dt; // Wsp. pionowa (wysokość npm.)

        // Prędkość tymczasowa
        v_tmp = v + (-F_drag(v, Cx, Sm, y, T, P0) / M - g * sin(alfa)) * dt;
        alfa -= g * cos(alfa) / v * dt; // Teraz alfa jest nachyleniem toru lotu pocisku
        v = v_tmp;
        Time += dt;
        // printf("%f\t%f\n",x,target_distance);
    }
    // Funkcja zwraca dystans na jaki doleciał pocisk [m]
    return fabs(x - target_distance)+fabs(y - target_height);
}

void print_trajectory_data_2(double alfa,const double target_distance, const double target_height)
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
    while ((y > target_height && x < target_distance + 15) || (y > 0.0 && x < target_distance))
    {
        dt = 1.0 / v;            // Dynamiczny krok czasowy, dla zwiększenia precyzji
                                 //  Wykorzystujemy metodę Eulera do obliczenia zmiany położenia,
                                 //  szybkości i kąta nachylenia stycznej do toru ruchu pocisku
        x += v * cos(alfa) * dt; // Wsp. pozioma
        y += v * sin(alfa) * dt; // Wsp. pionowa (wysokość npm.)

        fprintf(file, "%.2f\t%.2f\t%.2f\n", Time, x, y); // zaips do pliku


        // Prędkość tymczasowa
        v_tmp = v + (-F_drag(v, Cx, Sm, y, T, P0) / M - g * sin(alfa)) * dt;
        alfa -= g * cos(alfa) / v * dt; // Teraz alfa jest nachyleniem toru lotu pocisku
        v = v_tmp;
        Time += dt;
        // printf("%f\t%f\n",x,target_distance);
    }


    fclose(file);
}

int main()
{
    const double epsilon = 1e-4;
    const double alfa_start = 20;
    const double alfa_end = 70;

    // zad 2
    double max_range_angle = golden_section_search(wziuuum, alfa_start, alfa_end, epsilon);
    double max_range = wziuuum(max_range_angle);
    printf("max_range_angle: %f\nmax_range: %f\n", max_range_angle, max_range);

    // zad 3.1 zasięg działa w funkcji kąta ostrzału z zaznaczonym maximum
    FILE *file = fopen("zasięg_działa_w_funkcji_kąta.txt", "w");
    if (file == NULL)
    {
        printf("Błąd otwierania pliku.");
        return 1;
    }
    fprintf(file, "kąt [stopnie]\tzasięg [m] \n");
    for (double i = alfa_start; i < alfa_end; i += 1)
        fprintf(file, "%.2f\t%.2f\n", i, wziuuum(i)); // zaips do pliku
    fclose(file);

    // zad 3.2 trajektorię lotu pocisku dla maksymalnego zasięgu
    print_trajectory_data(max_range_angle);

    // zad 3.3 logarytmu modułu różnicy rozwiązania dokładnego
    golden_section_log_diff(wziuuum, alfa_start, alfa_end, 50.8057242086663, epsilon);

    // zad 4
    // Wyznaczenie kąta ostrzału dla celu znajdującego się 30 km od stanowiska artylerii
    double distance = 30000.0; // [m]
    // Cel na wysokości 0 m n.p.m.
    double target_height_a = 300.0; // [m]

    for (double alfa = alfa_start; alfa < alfa_end; alfa += 0.1)
    {
        double wzu = wziuuum2(alfa, distance, target_height_a);
        if (wzu < 15)
        {
            printf("shooting_angle_a: %f\nwziuuum2: %f\n", alfa, wzu);
        }
    }
    print_trajectory_data_2( 26.3, distance, target_height_a);


    return 0;
}
