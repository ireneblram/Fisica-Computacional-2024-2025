#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


// Definimos constantes
#define T 1000// Número de pasos de tiempo
#define h 0.01 // Espaciado espacial
#define PI 3.141592653589793
#define N 1000 // Número de puntos espaciales
#define I _Complex_I // Definición de la unidad imaginaria


// Función para inicializar la función de onda Phi
void inicializar_phi(double complex *Phi, double n_ciclos) {
    double k_tilde = (2.0 * PI * n_ciclos) / N; // \tilde{k}_0 = 2π n_ciclos / N
    double x0 = (N * h) / 4.0;                 // x_0 = Nh / 4
    double sigma = (N * h) / 16.0;             // σ = Nh / 16

    for (int j = 1; j < N; j++) {
        double x = j * h; // Posición espacial
        Phi[j] = cexp(I * k_tilde * j) * cexp(-pow((x - x0), 2) / (2 * sigma * sigma)); // Φ_{j,0} = e^{i\tilde{k}_0 x} e^{-(x - x_0)^2 / (2σ^2)}
    }

    Phi[0] = 0.0;
    Phi[N] = 0.0;

    // Normalización de Phi
    double norma = 0.0;
    for (int j = 0; j <= N; j++) {
        norma += pow(cabs(Phi[j]), 2) * h; // Norma = Σ |Φ[j]|^2 * h
    }
    norma = sqrt(norma);

    for (int j = 0; j <= N; j++) {
        Phi[j] /= norma; // Normalizamos Φ[j]
    }
}

// Función para inicializar el potencial V
void inicializar_potencial(double *V, double lambda, double k_tilde) {
    double V_j_h2 = lambda * k_tilde * k_tilde; // Altura del potencial
    int inicio = (2 * N) / 5;
    int fin = (3 * N) / 5;

    for (int j = 0; j <= N; j++) {
        if (j >= inicio && j <= fin) {
            V[j] = V_j_h2; // Dentro del rango del potencial
        } else {
            V[j] = 0.0; // Fuera del rango del potencial
        }
    }
}

int main() {
    // Declaramos las variables
    double complex Phi[N + 1], Phi_next[N + 1], chi[N + 1];
    double V[N + 1];
    complex double a[N + 1], b[N + 1], c[N + 1], d[N + 1];
    complex double gamma[N -1]; // Coeficientes gamma

    a[0] = 0.0; // a_0 = 0
    b[0] = 0.0; // b_0 = 0
    c[0] = 0.0; // c_0 = 0
    d[0] = 0.0; // d_0 = 0

    // Parámetros iniciales
    int n_ciclos = 1; // Número de ciclos
    double lambda = 1; // Altura del potencial
    double k_tilde = (2.0 * PI * n_ciclos) / N; // \tilde{k}_0

    // Calculamos s tilde
    double tilde_s = 1.0 / (4.0 * k_tilde * k_tilde);

    // Inicializamos la función de onda y el potencial
    inicializar_phi(Phi, n_ciclos);
    inicializar_potencial(V, lambda, k_tilde);

     // Abrimos el fichero para guardar las normas
     FILE *f_normas = fopen("normas.txt", "w");
     if (f_normas == NULL) {
         fprintf(stderr, "Error al abrir el fichero para guardar las normas.\n");
         return 1;
     }

     // Abrimos el fichero para guardar los datos de la función de onda
    FILE *f_datos = fopen("schrodinger_data.dat", "w");
    if (f_datos == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
        return 1;
    }


      for (int j = 1; j < N; j++) {
        a[j] = 1.0;
        b[j] = -2.0 + 2.0 * I / tilde_s - h * h * V[j];
        c[j] = 1.0;
      }

      // Calculamos alpha fuera del bucle temporal
      complex double alpha[N];
      alpha[N-1] = 0.0;

        for (int j = N-1; j > 0; j--) {
            gamma[j] = b[j] + c[j] * alpha[j];
            alpha[j-1] = -a[j] / gamma[j];
        }

        int j;
    // Iteración en el tiempo
    for (int n = 0; n < T; n++) {

        // Guardamos el valor de la función de onda en el fichero
        for ( j = 0; j <= N; j++) {

            fprintf(f_datos, "%d, %.10f\n", j, cabs(Phi[j]));
        }
        fprintf(f_datos, "\n"); // Separador entre instantes de tiempo

        // Construimos los coeficientes del sistema tridiagonal
        for (int j = 1; j < N; j++) {
            d[j] = 4.0 * I / tilde_s * Phi[j];
        }

        // Resolvemos el sistema tridiagonal para chi
        complex double beta[N]; // Beta aún depende del tiempo
        beta[N - 1] = 0.0 + 0.0*I;      // β_{N-1,n} = 0
        // Orden descendente
        for (int j = N - 1; j > 0; j--) {
            beta[j - 1] = (d[j] - c[j] * beta[j]) / gamma[j];
        }

        // Orden ascendente para calcular χ
        chi[0] = 0.0;
        chi[N] = 0.0;
        for (int j = 0; j < N - 1; j++) {
            chi[j + 1] = alpha[j]*chi[j] + beta[j];
        }

        // Calculamos Phi_next
        for (int j = 1; j < N; j++) {
            Phi_next[j] = chi[j] - Phi[j];
        }

        // Actualizamos Phi
        for (int j = 1; j < N; j++) {
            Phi[j] = Phi_next[j];
        }

        // Calculamos la norma de la función de onda
        double norma = 0.0;
        for (int j = 0; j <= N; j++) {
            norma += pow(cabs(Phi[j]), 2) * h; // Norma = Σ |Φ[j]|^2 * h
        }

        // Guardamos la norma en el fichero
        fprintf(f_normas, "Paso %d: Norma = %.10f\n", n, norma);

    }

    fclose(f_normas);
    fclose(f_datos);

    return 0;
}