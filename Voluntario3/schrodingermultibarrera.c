#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#define I _Complex_I
#define PI 3.141592653589793
#define h 0.01
#define N 2000
#define T 1000

typedef double complex cplx;

void phi_inicial(cplx *phi1, double ciclos) {
    double k_tilde, x0, sigma, norma;
    k_tilde = (2.0 * PI * ciclos) / N;
    x0 = (N * h) / 4.0;
    sigma = (N * h) / 16.0;

    phi1[0] = 0.0;
    phi1[N] = 0.0;

    norma = 0.0;
    for (int j = 1; j < N; j++) {
        double x = j * h;
        phi1[j] = cexp(I * k_tilde * j) * cexp(-pow((x - x0), 2) / (2 * sigma * sigma));
        norma += pow(cabs(phi1[j]), 2) * h;
    }
    norma = sqrt(norma);

    for (int j = 1; j < N; j++) {
        phi1[j] /= norma;
    }
}


void potencial_multibarrera(double *V, double lambda, double k_tilde, int n_barreras) {
    double V_tilde = lambda * k_tilde * k_tilde;
    int ancho_barrera = N / 5;
    int separacion = 2*N / 5; // Separación entre barreras
    int offset = 2*N / 5; // Desplazamiento inicial para la primera barrera

    // Inicializa todo el potencial a 0
    for (int j = 0; j <= N; j++) {
        V[j] = 0.0;
    }

    // Coloca cada barrera
    for (int k = 0; k < n_barreras; k++) {
        int inicio = offset + k * separacion;
        int fin = inicio + ancho_barrera;
        if (fin > N) fin = N;
        for (int j = inicio; j < fin; j++) {
            if (j <= N) V[j] = V_tilde;
        }
    }
}

void construir_tridiagonal(cplx *A1, cplx *A2, cplx *A3, double *V, double s_tilde) {
    for (int j = 1; j < N; j++) {
        A1[j] = 1.0;
        A2[j] = -2.0 + 2.0 * I / s_tilde - V[j];
        A3[j] = 1.0;
    }
}

void precalcular_alfa_gamma(cplx *A1, cplx *A2, cplx *A3, cplx *alfa, cplx *gamma) {
    alfa[N-1] = 0.0;
    for (int j = N-1; j > 0; j--) {
        gamma[j] = A2[j] + A3[j] * alfa[j];
        alfa[j-1] = -A1[j] / gamma[j];
    }
}

void paso_temporal(
    cplx *phi1, cplx *chi,
    cplx *A1, cplx *A2, cplx *A3, cplx *b,
    cplx *gamma, cplx *alfa, cplx *beta,
    double s_tilde,
    cplx *phi_next // Nuevo argumento
) {
    for (int j = 1; j < N; j++) {
        b[j] = 4.0 * I / s_tilde * phi1[j];
    }

    beta[N-1] = 0.0 + 0.0*I;
    for (int j = N-1; j > 0; j--) {
        beta[j-1] = (b[j] - A3[j] * beta[j]) / gamma[j];
    }

    chi[0] = 0.0;
    chi[N] = 0.0;
    for (int j = 0; j < N-1; j++) {
        chi[j+1] = alfa[j] * chi[j] + beta[j];
    }

    // Calculamos Phi_next
    for (int j = 1; j < N; j++) {
        phi_next[j] = chi[j] - phi1[j];
    }

    // Actualizamos Phi
    for (int j = 1; j < N; j++) {
        phi1[j] = phi_next[j];
    }
}

double calcular_norma(cplx *phi1) {
    double norma = 0.0;
    for (int j = 0; j <= N; j++) {
        norma += pow(cabs(phi1[j]), 2) * h;
    }
    return norma;
}

// Calcula PD(n): suma de |phi_j|^2 desde j=4N/5 hasta N
double calcular_PD(cplx *phi1) {
    double suma = 0.0;
    int ini = 4 * N / 5;
    for (int j = ini; j <= N; j++) {
        suma += pow(cabs(phi1[j]), 2) * h;
    }
    return suma;
}

int main() {
    srand((unsigned)time(NULL));
    int repeticiones = 10;
    int mT = 0;

    int ciclos = N/16;
    double lambda = 0.3;
    double k_tilde = (2.0 * PI * ciclos) / N;
    double s_tilde = 1.0 / (4.0 * k_tilde * k_tilde);

    cplx phi1[N + 1], chi[N + 1];
    cplx A1[N + 1], A2[N + 1], A3[N + 1], b[N + 1];
    cplx gamma[N - 1], alfa[N], beta[N];
    double V[N + 1];
    cplx phi_next[N + 1];

    int nD = 0;

    FILE *f_multibarrera = fopen("multibarrera.dat", "w");
    if (f_multibarrera == NULL) {
        printf("No se pudo abrir el fichero multibarrera.dat\n");
        return 1;
    }

    for (int rep = 0; rep < repeticiones; rep++) {
        A1[0] = A2[0] = A3[0] = b[0] = 0.0;

        // Número de barreras que caben en la malla con esa separación
        int n_barreras = 3;
        potencial_multibarrera(V, lambda, k_tilde, n_barreras);
        construir_tridiagonal(A1, A2, A3, V, s_tilde);
        precalcular_alfa_gamma(A1, A2, A3, alfa, gamma);
        phi_inicial(phi1, ciclos);

        if (rep == 0) {
            // Guardar estado inicial
            for (int j = 0; j <= N; j++) {
                double x = j * h;
                fprintf(f_multibarrera, "%.8f,%.10f,%.10f\n", x, cabs(phi1[j]), V[j]);
            }
            fprintf(f_multibarrera, "\n");

            // Buscar el máximo global de PD y guardar datos en cada paso
            double PD_act = calcular_PD(phi1);
            double PD_max = PD_act;
            nD = 0;
            for (int n = 1; n < T; n++) {
                paso_temporal(phi1, chi, A1, A2, A3, b, gamma, alfa, beta, s_tilde, phi_next);

                for (int j = 0; j <= N; j++) {
                    double x = j * h;
                    fprintf(f_multibarrera, "%.8f,%.10f,%.10f\n", x, cabs(phi1[j]), V[j]);
                }
                fprintf(f_multibarrera, "\n");

                PD_act = calcular_PD(phi1);
                if (PD_act > PD_max) {
                    PD_max = PD_act;
                    nD = n;
                }
            }
            fclose(f_multibarrera);
            printf("nD calculado (máximo global) = %d\n", nD);
            printf("PD(nD) = %.10f\n", PD_max);
        } else {
            // Evolucionar hasta nD
            for (int n = 0; n < nD; n++) {
                paso_temporal(phi1, chi, A1, A2, A3, b, gamma, alfa, beta, s_tilde, phi_next);
            }
            double PD = calcular_PD(phi1);
            double p = (double)rand() / (RAND_MAX + 1.0);
            printf("PD = %.10f, p = %.10f\n", PD, p);
            // Guardar los valores de phi1 en el archivo
            for (int j = 0; j <= N; j++) {
                double x = j * h;
                fprintf(f_multibarrera, "%.10f,%.10f,%.10f\n", x, creal(phi1[j]), cimag(phi1[j]));
            }
            fprintf(f_multibarrera, "\n");

            if (p < PD) {
                mT += 1;
            }
        }
    }

    printf("mT = %d\n", mT);
    double transmision = (double)mT / (repeticiones-1);
    printf("Coeficiente de transmisión estimado: %f\n", transmision);

    return 0;
}