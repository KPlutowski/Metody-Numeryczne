#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"

#include <stdio.h>
#include <math.h>

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


void licz_Euler() {
    const double a = 1;
    const double b = 1;
    const double c = 1;
    const double d = 1;

    double K = 1;
    double Z = 0.5;

    double t0 = 0;
    double tK = 40;
    double dt = 0.1;

    double t = t0;
    for (double i = t0; i < tK; i += dt) {
        K += dt * (a * K - b * K * Z);
        Z += dt * (c * K * Z - d * Z);
        

        // printf("%lf\t%lf\t%lf\n", t, K, Z);
        save_3_value_to_file( t, K, Z,"Euler.csv");
        t = t + dt;
    }
}

void licz_Runge() {
    const double a = 1;
    const double b = 1;
    const double c = 1;
    const double d = 1;

    double K = 1;
    double Z = 0.5;

    double t0 = 0;
    double tK = 40;
    double dt = 0.1;

    double t = t0;
    for (double i = t0; i < tK; i += dt) {
        double dK = dt * (a * K - b * K * Z);
        double dZ = dt * (c * K * Z - d * Z);

        double dKa = dt * (a * (K + dK / 2.0) - b * (K + dK / 2.0) * (Z + dZ / 2.0));
        double dZet = dt * (c * (K + dK / 2.0) * (Z + dZ / 2.0) - d * (Z + dZ / 2.0));

        K += dKa;
        Z += dZet;

        // printf("%lf\t%lf\t%lf\n", t, K, Z);
        save_3_value_to_file( t, K, Z,"Runge.csv");
        t += dt;
    }
}

int main()
{
    licz_Runge();
    licz_Euler();

    return 0;
}
