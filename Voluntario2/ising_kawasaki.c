#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 32 // Dimensión de la red cuadrada
#define pasos 1000 //Número de pasos de Monte Carlo 
#define T 5.0 // Temperatura inicial
// K_BOLTZMANN 1.0 Constante de Boltzmann se ha normlizado a 1.0 tal que beta = 1/T

// Función para inicializar la red con magnetización nula
void Red_inicial(int red[N][N]) {
   int i, j;
    //Primera fila: todos -1
    for (j = 0; j < N; j++) {
        red[0][j] = -1;
    }
    // Última fila: todos +1
    for (j = 0; j < N; j++) {
        red[N-1][j] = 1;
    }
    // Filas intermedias
    for (i = 1; i < N-1; i++) {
        int mitad = N / 2;
        int temp[N];
        // Llenar mitad con +1 y mitad con -1
        for (j = 0; j < mitad; j++) temp[j] = 1;
        for (j = mitad; j < N; j++) temp[j] = -1;
        // Mezclar aleatoriamente con Fisher-Yates shuffle (algoritmo de barajado de arrays)
        for (j = N-1; j > 0; j--) {
            int k = rand() % (j+1);
            int aux = temp[j];
            temp[j] = temp[k];
            temp[k] = aux;
        }
        // Asignar a la red
        for (j = 0; j < N; j++) {
            red[i][j] = temp[j];
        }
    }
}
//void Red_inicial(int red[N][N]) {
//    int i, j;
//    for (i = 0; i < N/2; i++) {         // Mitad superior
//        for (j = 0; j < N; j++) {
//            red[i][j] = -1;
//        }
//    }
//    for (i = N/2; i < N; i++) {         // Mitad inferior
//        for (j = 0; j < N; j++) {
//            red[i][j] = +1;
//        }
//    }
//}
// Función para mostrar la red en pantalla
void Mostrar_Red(int red[N][N]) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%2d ", red[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Función para guardar la red en un archivo
void Guardar_Red(FILE *archivo, int red[N][N]) {
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
double Energia_local(int red[N][N], int x1, int y1, int x2, int y2) {
    int dx[] = {-1, 1, 0, 0}; // arriba, abajo (desplazamiento en filas)
    int dy[] = {0, 0, -1, 1}; // izquierda, derecha (desplazamiento en columnas)
    int suma1 = 0, suma2 = 0;
    int d;

    // Suma de vecinos para (x1, y1), excluyendo (x2, y2)
    for (d = 0; d < 4; d++) {
        int nx = x1 + dx[d];
        int ny = (y1 + dy[d] + N) % N; // Periódica en y

        // Fija en x: si nx sale de rango, no suma ese vecino
        if ((nx >= 0 && nx < N) && (nx != x2 || ny != y2)) {
            suma1 += red[nx][ny];
        }
    }

    // Suma de vecinos para (x2, y2), excluyendo (x1, y1)
    for (d = 0; d < 4; d++) {
        int nx = x2 + dx[d];
        int ny = (y2 + dy[d] + N) % N;

        // Fija en x: si nx sale de rango, no suma ese vecino
        if ((nx >= 0 && nx < N) && (nx != x1 || ny != y1)) {
            suma2 += red[nx][ny];
        }
    }

    double energia = -red[x1][y1] * suma1 - red[x2][y2] * suma2;
    return energia;
}

// Realiza un paso Kawasaki en la red de Ising
void paso_Kawasaki(int red[N][N], double Temperatura) {
    int x1, y1, x2, y2;
    int dx[] = {1, -1, 0, 0};
    int dy[] = {0, 0, 1, -1};
    int vecinos[4] = {0, 1, 2, 3};

    // Selecciona un sitio aleatorio
    x1 = 1 + rand() % (N - 2); // x1 ∈ [1, N-2]
    y1 = rand() % N;

    // Baraja el orden de los vecinos (Fisher-Yates)
    for (int i = 3; i > 0; i--) {
        int j = rand() % (i + 1);
        int tmp = vecinos[i];
        vecinos[i] = vecinos[j];
        vecinos[j] = tmp;
    }

    // Busca un vecino válido (espín diferente y dentro de la red)
    int encontrado = 0;
    for (int i = 0; i < 4; i++) {
        int d = vecinos[i];
        x2 = x1 + dx[d];
        y2 = (y1 + dy[d] + N) % N;
        if (x2 > 0 && x2 < N-1) { // Solo filas 1 a N-2 pueden cambiar
            if (red[x1][y1] != red[x2][y2]) {
                encontrado = 1;
                break;
            }
        }
    }
    if (!encontrado) return; // No se encontró un par válido

    // Calcular energía antes del intercambio
    double E_antes = Energia_local(red, x1, y1, x2, y2);

    // Intercambiar espines
    int temp = red[x1][y1];
    red[x1][y1] = red[x2][y2];
    red[x2][y2] = temp;

    // Calcular energía después del intercambio
    double E_despues = Energia_local(red, x1, y1, x2, y2);

    double dE = E_despues - E_antes;

    // Criterio de Metropolis
    if (dE <= 0) {
        return;
    } else {
        double r = (double)rand() / RAND_MAX;
        if (r < exp(-dE / Temperatura)) {
            return;
        } else {
            // Revertir el intercambio
            temp = red[x1][y1];
            red[x1][y1] = red[x2][y2];
            red[x2][y2] = temp;
        }
    }
}

// Magnetización promedio en la mitad superior (filas 0 a N/2 - 1)
double magnetizacion_mitad_superior(int red[N][N]) {
    int suma = 0;
    for (int i = 0; i < N/2; i++) {
        for (int j = 0; j < N; j++) {
            suma += red[i][j];
        }
    }
    return fabs((double)suma) / ((N/2) * N);
}

// Magnetización promedio en la mitad inferior (filas N/2 a N-1)
double magnetizacion_mitad_inferior(int red[N][N]) {
    int suma = 0;
    for (int i = N/2; i < N; i++) {
        for (int j = 0; j < N; j++) {
            suma += red[i][j];
        }
    }
    return fabs((double)suma) / ((N/2) * N);
}

int main() {
    clock_t inicio = clock(); // Medir el tiempo de inicio
    srand(time(NULL)); // Inicializar la semilla de números aleatorios

    int red[N][N];
    Red_inicial(red);

    printf("Configuración inicial de la red:\n");
    Mostrar_Red(red);

    // Abrir archivo para guardar la red en cada paso
    FILE *archivo = fopen("kawasaki_red.txt", "w");
    if (archivo == NULL) {
        printf("No se pudo abrir el archivo para escritura.\n");
        return 1;
    }

    // Guardar configuración inicial
    Guardar_Red(archivo, red);

    // Arreglos para almacenar magnetizaciones
    int num_medidas = pasos / 100;
    double mags_sup[num_medidas];
    double mags_inf[num_medidas];
    int medida_idx = 0;

    // Ciclo de Monte Carlo: realiza 'k' pasos Kawasaki
    int i, j;
    for (i = 0; i < pasos; i++) {
        for (j = 0; j < N*N; j++) {
            paso_Kawasaki(red, T);
        }
        Guardar_Red(archivo, red); // Guardar la red después de cada paso Monte Carlo

        // Calcular magnetización cada 100 pasos
        if ((i+1) % 100 == 0) {
            double m_sup = magnetizacion_mitad_superior(red);
            double m_inf = magnetizacion_mitad_inferior(red);
            mags_sup[medida_idx] = m_sup;
            mags_inf[medida_idx] = m_inf;
            printf("Paso %d: Magnetización superior = %f, inferior = %f\n", i+1, m_sup, m_inf);
            medida_idx++;
        }
    }

    fclose(archivo);

    // Calcular la media de las magnetizaciones
    double suma_sup = 0.0, suma_inf = 0.0;
    for (i = 0; i < medida_idx; i++) {
        suma_sup += mags_sup[i];
        suma_inf += mags_inf[i];
    }
    double media_sup = suma_sup / medida_idx;
    double media_inf = suma_inf / medida_idx;

    printf("\nMagnetización media superior (cada 100 pasos): %f\n", media_sup);
    printf("Magnetización media inferior (cada 100 pasos): %f\n", media_inf);

    printf("\nConfiguración final de la red:\n");
    Mostrar_Red(red);

    // Medir el tiempo de finalización
    clock_t fin = clock();
    double tiempo = (double)(fin - inicio) / CLOCKS_PER_SEC;
    printf("\nTiempo de ejecución: %.2f segundos.\n", tiempo);

    return 0;
}