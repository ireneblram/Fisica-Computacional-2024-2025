#include <stdio.h>
#include <stdlib.h>
#include <complex.h> // Para manejar números complejos en C
#include <math.h>    // Para funciones matemáticas como pow, sqrt, M_PI
#include <time.h>    // Para la inicialización de la semilla de números aleatorios

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Define complex double para conveniencia
typedef double complex cdouble;

// Función para resolver un sistema tridiagonal (adaptada del apéndice del PDF)
// Esta versión es destructiva para los arrays 'b_diag' y 'c_diag' (la diagonal principal y superior),
// pero es la más sencilla de implementar directamente del algoritmo del PDF.
// Para un uso más general y preservación de los datos, se requeriría espacio de "scratch" adicional.
void solve_tridiagonal_destructive(cdouble x[], int N, cdouble a_diag[], cdouble b_diag[], cdouble c_diag[]) {
    int n;

    // Eliminación hacia adelante
    // c_diag[0] aquí es en realidad A_0^+ / A_0^0 según el algoritmo de recurrencia,
    // y x[0] es b_0,n / A_0^0.
    // Sin embargo, en la ecuación (18) y (19) del PDF, A_j^- = 1, A_j^+ = 1, y A_j^0 es el término central.
    // La nomenclatura de las variables 'a_diag', 'b_diag', 'c_diag' aquí se mapea a:
    // a_diag[j] = A_j^- (sub-diagonal)
    // b_diag[j] = A_j^0 (diagonal principal)
    // c_diag[j] = A_j^+ (super-diagonal)
    // x[j] = b_j,n (lado derecho del sistema)

    // El PDF presenta la solución en términos de alpha y beta.
    // Para simplificar la implementación y usar una función de resolución de tridiagonales más genérica (como la del apéndice),
    // necesitamos mapear los coeficientes de la ecuación (18) a las notaciones 'a', 'b', 'c' de un solucionador estándar.
    // La ecuación (18) es: chi_{j+1,n} + [-2 + 2i/s_tilde - V_tilde_j] * chi_{j,n} + chi_{j-1,n} = (4i/s_tilde) * Phi_j,n
    // Reordenando para el formato estándar: A_j^- * chi_{j-1,n} + A_j^0 * chi_{j,n} + A_j^+ * chi_{j+1,n} = RHS_j
    // A_j^- = 1
    // A_j^0 = -2 + 2i/s_tilde - V_tilde_j
    // A_j^+ = 1
    // RHS_j = (4i/s_tilde) * Phi_j,n

    // El algoritmo de Thomas (o de eliminación de Gauss para tridiagonales)
    // reescribe el sistema como: A_j^- * x_{j-1} + b_j * x_j + c_j * x_{j+1} = d_j
    // Donde a, b, c son los nombres genéricos para las diagonales y d el lado derecho.
    // En nuestro caso, 'a_diag' y 'c_diag' son siempre 1.0. 'b_diag' es A_j^0. 'x' contiene RHS_j.

    // Paso de eliminación hacia adelante (modifica b_diag y x)
    for (n = 1; n < N; n++) {
        cdouble m = a_diag[n] / b_diag[n - 1]; // b_diag[n-1] ha sido modificado en la iteración anterior
        b_diag[n] = b_diag[n] - m * c_diag[n - 1];
        x[n] = x[n] - m * x[n - 1];
    }

    // Paso de sustitución hacia atrás
    x[N - 1] = x[N - 1] / b_diag[N - 1];
    for (n = N - 2; n >= 0; --n) {
        x[n] = (x[n] - c_diag[n] * x[n + 1]) / b_diag[n];
    }
}


// Función para calcular la densidad de probabilidad |Phi(x)|^2
double calculate_probability_density(cdouble phi) {
    return creal(phi * conj(phi));
}

// Función para calcular P_D(n)
// P_D(n) = Suma de |Phi_j,n|^2 para j desde 4N/5 hasta N [cite: 598]
double calculate_PD(cdouble *phi, int N, double h) {
    double pd = 0.0;
    // La región del detector a la derecha va desde 4N/5 hasta N [cite: 598]
    for (int j = (int)(4.0 * N / 5.0); j <= N; j++) {
        pd += calculate_probability_density(phi[j]) * h; // Multiplicar por h para la integral numérica
    }
    return pd;
}

// Función para calcular la normalización (suma de |Phi_j,n|^2 * h)
double calculate_normalization(cdouble *phi, int N, double h) {
    double norm = 0.0;
    for (int j = 0; j <= N; j++) {
        norm += calculate_probability_density(phi[j]) * h;
    }
    return norm;
}

int main() {
    // Inicializar el generador de números aleatorios
    srand(time(NULL));

    // 1. Parámetros iniciales [cite: 593]
    int N_values[] = {500, 1000}; // N (número de puntos espaciales)
    double lambda_values[] = {0.1, 0.3, 0.5, 1.0, 5.0, 10.0}; // lambda (factor de altura del potencial)
    int n_ciclos = 5; // Número de oscilaciones completas [cite: 584, 586]
    int num_simulations = 5000; // Al menos 10^3 simulaciones [cite: 598] (aumentado para mejor estadística)

    // h = 1.0 en unidades reescaladas, como se asume implícitamente por el uso de j directamente como 'x' en algunas fórmulas reescaladas.
    // Aunque la discretización es x_j = j*h, para las constantes reescaladas del PDF, podemos tomar h=1.
    double h = 1.0;

    printf("Estudio del Coeficiente de Transmisión\n");
    printf("---------------------------------------\n");

    // Bucle sobre diferentes valores de N
    for (int n_idx = 0; n_idx < sizeof(N_values) / sizeof(N_values[0]); n_idx++) {
        int N = N_values[n_idx];
        printf("\n--- Simulacion para N = %d ---\n", N);

        // Calcular k0_tilde y s_tilde que dependen de N y n_ciclos [cite: 589, 591]
        double k0_tilde = 2.0 * M_PI * n_ciclos / N; // k_0_tilde = k_0 * h [cite: 589]
        double s_tilde = 1.0 / (4.0 * k0_tilde * k0_tilde); // s_tilde = s / h^2 [cite: 591]

        // Bucle sobre diferentes valores de lambda
        for (int lambda_idx = 0; lambda_idx < sizeof(lambda_values) / sizeof(lambda_values[0]); lambda_idx++) {
            double lambda = lambda_values[lambda_idx];
            int total_mT = 0; // Contador de detecciones a la derecha [cite: 598]
            double total_PD_at_nD = 0.0; // Suma de P_D(n_D) para promediar

            printf("\n  Lambda = %.2f:\n", lambda);

            // Reservar memoria para las variables que se reutilizarán en cada simulación
            cdouble *phi = (cdouble *)malloc((N + 1) * sizeof(cdouble));
            cdouble *chi = (cdouble *)malloc((N + 1) * sizeof(cdouble));
            double *V_tilde = (double *)malloc((N + 1) * sizeof(double));
            // Coeficientes para el solucionador tridiagonal (necesitan ser copiados si son modificados destructivamente)
            cdouble *A0_coeffs_copy = (cdouble *)malloc((N + 1) * sizeof(cdouble));
            cdouble *A_minus_coeffs_copy = (cdouble *)malloc((N + 1) * sizeof(cdouble));
            cdouble *A_plus_coeffs_copy = (cdouble *)malloc((N + 1) * sizeof(cdouble));


            // Calcular el potencial V_tilde_j [cite: 590]
            for (int j = 0; j <= N; j++) {
                if (j >= (int)(2.0 * N / 5.0) && j <= (int)(3.0 * N / 5.0)) {
                    V_tilde[j] = lambda * k0_tilde * k0_tilde;
                } else {
                    V_tilde[j] = 0.0;
                }
            }

            for (int sim = 0; sim < num_simulations; sim++) {

                // Generar la función de onda inicial Phi(x,0) [cite: 589]
                // phi_j_0 = exp(i * k0_tilde * j) * exp(-8 * (4j - N)^2 / N^2) [cite: 589]
                for (int j = 0; j <= N; j++) {
                    phi[j] = cexp(I * k0_tilde * j) * cexp(-8.0 * pow((4.0 * j - N), 2) / pow(N, 2));
                }

                // Aplicar condiciones de contorno: Phi_0,n = Phi_N,n = 0 [cite: 235, 562]
                phi[0] = 0.0;
                phi[N] = 0.0;

                // Normalizar la función de onda inicial (importante para la interpretación de probabilidad) [cite: 222]
                double norm_factor = 0.0;
                for (int j = 0; j <= N; j++) {
                    norm_factor += calculate_probability_density(phi[j]);
                }
                norm_factor = sqrt(norm_factor * h); // Multiplicar por h para la integral.

                if (norm_factor == 0.0) { // Evitar división por cero si la función de onda es nula
                    // Esto no debería ocurrir con una inicialización gaussiana adecuada, pero es una buena práctica.
                    fprintf(stderr, "Error: La normalización inicial es cero. Saltando simulación.\n");
                    continue;
                }
                for (int j = 0; j <= N; j++) {
                    phi[j] /= norm_factor;
                }

                // 3. Buscar el valor t=nD correspondiente al primer máximo local de P_D(t) [cite: 490]
                int nD = 0;
                double max_PD_val = 0.0;
                double prev_PD = 0.0;
                double current_PD = 0.0;
                // Un número máximo de pasos para evitar bucles infinitos si no hay un máximo claro.
                // Ajustar según la física del problema (ej. el tiempo que tarda la onda en cruzar la barrera).
                int max_evolution_steps = (int)(2.0 * N / k0_tilde) + 1000; // Heurística, depende de k0_tilde

                // Usar copias de los coeficientes del solucionador tridiagonal ya que `solve_tridiagonal_destructive`
                // los modifica.
                cdouble *current_phi = (cdouble *)malloc((N + 1) * sizeof(cdouble));
                for(int j=0; j<=N; j++) current_phi[j] = phi[j]; // Copia de la phi inicial

                for (int n = 0; n < max_evolution_steps; n++) {
                    // Preparar los coeficientes para el sistema tridiagonal (Ecuación 18) [cite: 572]
                    // chi_{j+1,n} + [-2 + 2i/s_tilde - V_tilde_j] * chi_{j,n} + chi_{j-1,n} = (4i/s_tilde) * Phi_j,n [cite: 572]
                    // Mapeo a A_j^- * chi_{j-1} + A_j^0 * chi_{j} + A_j^+ * chi_{j+1} = RHS_j [cite: 575]
                    // A_j^- = 1 [cite: 576]
                    // A_j^0 = -2 + 2i/s_tilde - V_tilde_j [cite: 576]
                    // A_j^+ = 1 [cite: 576]
                    // RHS_j = (4i/s_tilde) * Phi_j,n [cite: 576]

                    for (int j = 0; j <= N; j++) {
                        A_minus_coeffs_copy[j] = 1.0; // A_j^-
                        A_plus_coeffs_copy[j] = 1.0;  // A_j^+
                        A0_coeffs_copy[j] = (-2.0 + 2.0 * I / s_tilde - V_tilde[j]); // A_j^0
                        chi[j] = (4.0 * I / s_tilde) * current_phi[j]; // Lado derecho (RHS_j)
                    }

                    // Aplicar condiciones de contorno para chi [cite: 576]
                    chi[0] = 0.0;
                    chi[N] = 0.0;

                    // Resolver el sistema tridiagonal para chi.
                    // ¡Importante!: La función solve_tridiagonal_destructive modifica A0_coeffs_copy y A_plus_coeffs_copy.
                    // Por eso, se usan copias para cada paso de tiempo si la función es destructiva.
                    solve_tridiagonal_destructive(chi, N + 1, A_minus_coeffs_copy, A0_coeffs_copy, A_plus_coeffs_copy);

                    // Actualizar Phi para el siguiente paso de tiempo (Ecuación 15) [cite: 573]
                    // Phi_{j,n+1} = chi_{j,n} - Phi_{j,n} [cite: 570]
                    for (int j = 0; j <= N; j++) {
                        current_phi[j] = chi[j] - current_phi[j];
                    }

                    // Reaplicar condiciones de contorno después de la actualización [cite: 235, 562]
                    current_phi[0] = 0.0;
                    current_phi[N] = 0.0;

                    // Calcular P_D(n) para el tiempo actual [cite: 598]
                    current_PD = calculate_PD(current_phi, N, h);

                    // Lógica para encontrar el primer máximo local de P_D(t) [cite: 489, 490]
                    // Se busca un punto donde P_D(n) es mayor que el valor anterior
                    // y luego se monitorea si empieza a decrecer. Esto es una heurística
                    // para el "primer máximo local" y podría necesitar ajustes.
                    // Un enfoque más robusto podría involucrar guardar los últimos 3 valores
                    // y verificar (prev < current && current > next)
                    if (n > 0) {
                        if (current_PD > max_PD_val) {
                            max_PD_val = current_PD;
                            nD = n + 1; // nD corresponde al número de pasos
                        } else if (nD > 0 && current_PD < prev_PD * 0.95) { // Si ya encontramos un máximo y empieza a bajar significativamente
                            // Esto es una simplificación. Un pico real requiere un análisis de la "derivada" numérica.
                            // Aquí asumimos que el valor más alto encontrado hasta ahora es el "primer máximo relevante".
                            break; // Romper el bucle una vez que el máximo se considera encontrado
                        }
                    }
                    prev_PD = current_PD;
                }

                // 4. Evolucionar nD pasos (ya se hizo en el bucle anterior) [cite: 490]
                // 5. Calcular P_D(nD) (ya tenemos max_PD_val que es P_D(nD)) [cite: 490]
                //    Aunque el algoritmo del PDF dice "evolucionar nD pasos" y luego "calcular P_D(nD)",
                //    nuestra implementación ya lo hace mientras busca el máximo.
                //    Asumimos que `max_PD_val` es el P_D(nD) para esta simulación.
                total_PD_at_nD += max_PD_val; // Acumular para el promedio

                // 6. Simular el proceso de medición generando un número aleatorio p ∈ [0, 1] [cite: 491, 492]
                double p_rand = (double)rand() / RAND_MAX; // Generar un número aleatorio entre 0 y 1
                if (p_rand < max_PD_val) { // Si p < P_D(nD) habremos detectado la partícula [cite: 491, 492]
                    total_mT++; // Actualizar el valor de mT [cite: 492]
                }
                free(current_phi); // Liberar memoria para el próximo ciclo de simulación
            }

            // Liberar memoria para los buffers temporales
            free(phi);
            free(chi);
            free(V_tilde);
            free(A0_coeffs_copy);
            free(A_minus_coeffs_copy);
            free(A_plus_coeffs_copy);


            // Calcular el coeficiente de transmisión K = mT / m [cite: 598]
            double K = (double)total_mT / num_simulations;
            double avg_PD_at_nD = total_PD_at_nD / num_simulations;

            printf("    K (Coeficiente de Transmision) = %.6f\n", K);
            printf("    Promedio de P_D(nD)           = %.6f\n", avg_PD_at_nD);

            // 2. Comparar K y P_D(nD) y explicar la relación
            // El coeficiente de transmisión 'K' (obtenido experimentalmente por las simulaciones) [cite: 598]
            // debería ser aproximadamente igual a la probabilidad promedio de detección 'P_D(nD)' (la probabilidad cuántica)[cite: 598].
            // Esto se debe a que, en un gran número de repeticiones (simulaciones), la frecuencia relativa de un evento
            // (la detección de la partícula) converge a su probabilidad cuántica. [cite: 598]
            // La maximización de P_D(t) para encontrar nD asegura que la medición se realiza en un momento óptimo,
            // cuando la onda ya ha interaccionado con la barrera y la probabilidad de encontrarla al otro lado es la más alta. [cite: 489, 490]

            // 3. Estimar el error estadístico de K
            // Dado que K es la proporción de "éxitos" en una serie de "intentos" (simulaciones),
            // se puede modelar como una variable binomial.
            // La desviación estándar para una proporción K es sqrt(K * (1 - K) / num_simulations).
            double std_error_K = sqrt(K * (1.0 - K) / num_simulations);
            printf("    Error Estadistico de K         = %.6f\n", std_error_K);
        }
    }

    return 0;
}