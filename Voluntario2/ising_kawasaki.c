#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define X 4 // Dimensión de la red en la dirección x
#define Y 4 // Dimensión de la red en la dirección y
#define pasos 100 //Número de pasos de Monte Carlo 
#define T 2.5 // Temperatura inicial
#define K_BOLTZMANN 1.0 // Constante de Boltzmann

// Función para inicializar la red con magnetización nula
void Red_inicial(int red[X][Y]) {
    int i, j;
    // Primera fila: todos -1
    for (j = 0; j < Y; j++) {
        red[0][j] = -1;
    }
    // Última fila: todos +1
    for (j = 0; j < Y; j++) {
        red[X-1][j] = 1;
    }
    // Filas intermedias
    for (i = 1; i < X-1; i++) {
        int mitad = Y / 2;
        int temp[Y];
        // Llenar mitad con +1 y mitad con -1
        for (j = 0; j < mitad; j++) temp[j] = 1;
        for (j = mitad; j < Y; j++) temp[j] = -1;
        // Mezclar aleatoriamente con Fisher-Yates shuffle (algoritmo de barajado de arrays)
        for (j = Y-1; j > 0; j--) {
            int k = rand() % (j+1);
            int aux = temp[j];
            temp[j] = temp[k];
            temp[k] = aux;
        }
        // Asignar a la red
        for (j = 0; j < Y; j++) {
            red[i][j] = temp[j];
        }
    }
}

// Función para mostrar la red en pantalla
void Mostrar_Red(int red[X][Y]) {
    int i, j;
    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            printf("%2d ", red[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Función para guardar la red en un archivo
void Guardar_Red(FILE *archivo, int red[X][Y]) {
    int i, j;
    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            fprintf(archivo, "%d", red[i][j]);
            if (j < Y - 1) {
                fprintf(archivo, ","); // Separador entre columnas
            }
        }
        fprintf(archivo, "\n"); // Nueva línea al final de cada fila
    }
    fprintf(archivo, "\n"); // Línea en blanco para separar iteraciones
}


// Función para calcular el cambio de energía al intercambiar dos espines vecinos
double calcular_E(int red[X][Y], int x1, int y1, int x2, int y2) {
    int dx[] = {-1, 1, 0, 0}; // arriba, abajo (desplazamiento en filas)
    int dy[] = {0, 0, -1, 1}; // izquierda, derecha (desplazamiento en columnas)
    int suma1 = 0, suma2 = 0;
    int d;

    // Suma de vecinos para (x1, y1), excluyendo (x2, y2)
    for (d = 0; d < 4; d++) {
        int nx = x1 + dx[d];
        int ny = (y1 + dy[d] + Y) % Y; // Periódica en y

        // Fija en x: si nx sale de rango, no suma ese vecino
        if ((nx >= 0 && nx < X) && (nx != x2 || ny != y2)) {
            suma1 += red[nx][ny];
        }
    }

    // Suma de vecinos para (x2, y2), excluyendo (x1, y1)
    for (d = 0; d < 4; d++) {
        int nx = x2 + dx[d];
        int ny = (y2 + dy[d] + Y) % Y;

        // Fija en x: si nx sale de rango, no suma ese vecino
        if ((nx >= 0 && nx < X) && (nx != x1 || ny != y1)) {
            suma2 += red[nx][ny];
        }
    }

    double energia = -red[x1][y1] * suma1 - red[x2][y2] * suma2;
    return energia;
}

// Realiza un paso Kawasaki en la red de Ising
void paso_Kawasaki(int red[X][Y], double Temperatura) {
    int x1, y1, x2, y2;
    int dx[] = {1, -1, 0, 0};
    int dy[] = {0, 0, 1, -1};
    int vecinos[4] = {0, 1, 2, 3};

    // Selecciona un sitio aleatorio
    x1 = 1 + rand() % (X - 2); // x1 ∈ [1, X-2]
    y1 = rand() % Y;

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
        y2 = (y1 + dy[d] + Y) % Y;
        if (x2 > 0 && x2 < X-1) { // Solo filas 1 a X-2 pueden cambiar
            if (red[x1][y1] != red[x2][y2]) {
                encontrado = 1;
                break;
            }
        }
    }
    if (!encontrado) return; // No se encontró un par válido

    // Calcular energía antes del intercambio
    double E_antes = calcular_E(red, x1, y1, x2, y2);

    // Intercambiar espines
    int temp = red[x1][y1];
    red[x1][y1] = red[x2][y2];
    red[x2][y2] = temp;

    // Calcular energía después del intercambio
    double E_despues = calcular_E(red, x1, y1, x2, y2);

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


int main() {
    clock_t inicio = clock(); // Medir el tiempo de inicio
    srand(time(NULL)); // Inicializar la semilla de números aleatorios

    int red[X][Y];
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

    // Ciclo de Monte Carlo: realiza 'k' pasos Kawasaki
    int i, j;
    for (i = 0; i < pasos; i++) {
        for (j = 0; j < X*Y; j++) {
            pasoKawasaki(red, T);
        }
        Guardar_Red(archivo, red); // Guardar la red después de cada paso Monte Carlo
    }


    fclose(archivo);

    printf("\nConfiguración final de la red:\n");
    Mostrar_Red(red);

    // Medir el tiempo de finalización
    clock_t fin = clock();
    double tiempo = (double)(fin - inicio) / CLOCKS_PER_SEC;
    printf("\nTiempo de ejecución: %.2f segundos.\n", tiempo);

    return 0;
}