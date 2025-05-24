#include <stdio.h>
#include <complex.h> // Esta es la cabecera para C

int main() {
    double complex z1 = 1.0 + 2.0 * I; // 'I' es la unidad imaginaria definida en complex.h
    double complex z2 = 3.0 - 4.0 * I;
    double complex sum = z1 + z2;

    printf("z1 = %.2f + %.2fi\n", creal(z1), cimag(z1));
    printf("z2 = %.2f + %.2fi\n", creal(z2), cimag(z2));
    printf("sum = %.2f + %.2fi\n", creal(sum), cimag(sum));

    return 0; 
    
}