#include <iostream>
#include <cmath>
#include <array>
#include <random>
#include <limits>

// Constantes
const int N = 64; // Tamaño de la red (N x N)
const int ITERACIONES = 10000; // Número de iteraciones
const double T = 2.5; // Temperatura inicial
const double K_BOLTZMANN = 1.0; // Constante de Boltzmann (J/K)

// Implementación básica del generador Philox
class Philox {
public:
    using result_type = uint32_t;

    explicit Philox(uint64_t seed) : key(seed), counter(0) {}

    result_type operator()() {
        counter++;
        return static_cast<result_type>((key ^ counter) * 0x9E3779B97F4A7C15ULL);
    }

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

private:
    uint64_t key;
    uint64_t counter;
};

// Función para inicializar la red con magnetización nula
void inicializarRed(int red[N][N], Philox &philox) {
    std::uniform_int_distribution<int> distribucion(0, 1);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            red[i][j] = distribucion(philox) * 2 - 1; // Asigna aleatoriamente +1 o -1
        }
    }
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

// Función para imprimir la red
void imprimirRed(int red[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << (red[i][j] > 0 ? "+" : "-") << " ";
        }
        std::cout << "\n";
    }
}

// Algoritmo de Monte Carlo con dinámica de Kawasaki
void monteCarloKawasaki(int red[N][N], double beta, int iteraciones, Philox &philox) {
    std::uniform_int_distribution<int> distribucion(0, N - 1);
    std::uniform_real_distribution<double> distribucion_real(0.0, 1.0);

    for (int i = 0; i < iteraciones; i++) {
        // Elegir dos espines vecinos al azar
        int x1 = distribucion(philox);
        int y1 = distribucion(philox);

        // Elegir un vecino aleatorio
        int x2, y2;
        int direccion = distribucion(philox) % 4;
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
        double deltaE = calcularDeltaE(red, x1, y1, x2, y2);

        // Calcular la probabilidad de transición
        double probabilidad = std::exp(-beta * deltaE);
        probabilidad = std::min(1.0, probabilidad);

        // Generar un número aleatorio para decidir si aceptar el cambio
        double r = distribucion_real(philox);
        if (r < probabilidad) {
            // Intercambiar los espines
            std::swap(red[x1][y1], red[x2][y2]);
        }
    }
}

int main() {
    // Inicializar el generador Philox con una semilla
    uint64_t seed = std::random_device{}();
    Philox philox(seed);

    // Inicializar la red
    int red[N][N];
    inicializarRed(red, philox);

    double beta = 1.0 / (K_BOLTZMANN * T);

    std::cout << "Configuración inicial de la red:\n";
    imprimirRed(red);

    monteCarloKawasaki(red, beta, ITERACIONES, philox);

    std::cout << "\nConfiguración final de la red:\n";
    imprimirRed(red);

    return 0;
}