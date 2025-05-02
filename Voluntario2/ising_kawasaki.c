#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Constantes
#define N 4 // Tamaño de la red (N x N)
#define ITERACIONES 100*N^2 // Número de iteraciones
#define T 2.5 // Temperatura inicial
#define K_BOLTZMANN 1.0 // Constante de Boltzmann (J/K)

// Función para inicializar la red con magnetización nula
void inicializarRed(int red[N][N]) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            red[i][j] = (rand() % 2) * 2 - 1; // Asigna aleatoriamente +1 o -1
        }
    }
}

// Función para imprimir la red
void imprimirRed(int red[N][N]) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%2d ", red[i][j]);
        }
        printf("\n");
    }
}

// Función para guardar la red en un archivo
void guardarRed(FILE *archivo, int red[N][N]) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            fprintf(archivo, "%d", red[i][j]);
            if (j < N - 1) {
                fprintf(archivo, ","); // Separador entre columnas
            }
        }
        fprintf(archivo, "\n"); // Nueva línea al final de cada fila
    }
    fprintf(archivo, "\n"); // Línea en blanco para separar iteraciones
}

// Función para calcular el cambio de energía al intercambiar dos espines vecinos
double calcularDeltaE(int red[N][N], int x1, int y1, int x2, int y2) {
    int suma_vecinos1 = red[(x1 + 1) % N][y1] + red[(x1 - 1 + N) % N][y1] +
                        red[x1][(y1 + 1) % N] + red[x1][(y1 - 1 + N) % N];
    int suma_vecinos2 = red[(x2 + 1) % N][y2] + red[(x2 - 1 + N) % N][y2] +
                        red[x2][(y2 + 1) % N] + red[x2][(y2 - 1 + N) % N];

    // Cambio de energía al intercambiar los espines
    return 2.0 * (red[x1][y1] - red[x2][y2]) * (suma_vecinos1 - suma_vecinos2);
}

// Algoritmo de Monte Carlo con dinámica de Kawasaki
void monteCarloKawasaki(int red[N][N], double beta, int iteraciones) {
    int i, x1, y1, x2, y2;
    double deltaE, probabilidad, r;

    FILE *archivo = fopen("kawasaki_red.txt", "w");
    if (archivo == NULL) {
        fprintf(stderr, "Error al abrir el archivo para guardar la red.\n");
        exit(1);
    }

    for (i = 0; i < iteraciones; i++) {
        // Elegir dos espines vecinos al azar
        x1 = rand() % N;
        y1 = rand() % N;

        // Elegir un vecino aleatorio
        int direccion = rand() % 4;
        if (direccion == 0) {
            x2 = (x1 + 1) % N;
            y2 = y1;
        } else if (direccion == 1) {
            x2 = (x1 - 1 + N) % N;
            y2 = y1;
        } else if (direccion == 2) {
            x2 = x1;
            y2 = (y1 + 1) % N;
        } else {
            x2 = x1;
            y2 = (y1 - 1 + N) % N;
        }

        // Calcular el cambio de energía
        deltaE = calcularDeltaE(red, x1, y1, x2, y2);

        // Calcular la probabilidad de transición
        probabilidad = exp(-beta * deltaE);
        if (probabilidad > 1.0) {
            probabilidad = 1.0;
        }

        // Generar un número aleatorio para decidir si aceptar el cambio
        r = (double)rand() / RAND_MAX;
        if (r < probabilidad) {
            // Intercambiar los espines
            int temp = red[x1][y1];
            red[x1][y1] = red[x2][y2];
            red[x2][y2] = temp;
        }
    }

    fclose(archivo);
}

int main() {
    clock_t inicio = clock(); // Medir el tiempo de inicio
    srand(time(NULL)); // Inicializar la semilla de números aleatorios

    int red[N][N];
    inicializarRed(red);

    double beta = 1.0 / (K_BOLTZMANN * T);

    printf("Configuración inicial de la red:\n");
    imprimirRed(red);

    monteCarloKawasaki(red, beta, ITERACIONES);

    printf("\nConfiguración final de la red:\n");
    imprimirRed(red);

     // Medir el tiempo de finalización
     clock_t fin = clock();
     double tiempo = (double)(fin - inicio) / CLOCKS_PER_SEC;
     printf("\nTiempo de ejecución: %.2f segundos.\n", tiempo);

    return 0;
}