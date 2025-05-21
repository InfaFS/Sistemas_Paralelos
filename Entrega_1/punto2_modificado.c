#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>

#define N 4096 // Tamaño de la matriz, modificar según necesidad
#define BLOCK_SIZE 256 // Tamaño del bloque para la multiplicación por bloques

// Calcula el valor máximo, mínimo y promedio de una matriz almacenada por filas
void calcular_max_min_prom(double *matriz, double *max, double *min, double *prom)
{
    *max = -DBL_MAX;
    *min = DBL_MAX;
    double suma = 0;

    for (int i = 0; i < N * N; i++)
    {
        double val = matriz[i];
        if (val > *max)
            *max = val;
        if (val < *min)
            *min = val;
        suma += val;
    }
    *prom = suma / (N * N);
}

// Transpone la matriz B antes de la multiplicación para mejorar la localidad
void transponer_matriz(double *matriz, double *transpuesta)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            transpuesta[j * N + i] = matriz[i * N + j];
        }
    }
}

// Multiplica matrices por bloques para mejorar el acceso a memoria
void multiplicar_matrices_bloques(double *A, double *B, double *resultado)
{
    for (int i = 0; i < N * N; i++)
    {
        resultado[i] = 0;
    }

    for (int bi = 0; bi < N; bi += BLOCK_SIZE)
    {
        for (int bj = 0; bj < N; bj += BLOCK_SIZE)
        {
            for (int bk = 0; bk < N; bk += BLOCK_SIZE)
            {
                for (int i = bi; i < bi + BLOCK_SIZE && i < N; i++)
                {
                    for (int j = bj; j < bj + BLOCK_SIZE && j < N; j++)
                    {
                        double sum = 0;
                        for (int k = bk; k < bk + BLOCK_SIZE && k < N; k++)
                        {
                            sum += A[i * N + k] * B[j * N + k];
                        }
                        resultado[i * N + j] += sum;
                    }
                }
            }
        }
    }
}

// Calcula la matriz R según la ecuación dada
void calcular_R(double *A, double *B, double *C, double *R, double *Bt, double *AB, double *CBt)
{
    double maxA, minA, promA, maxB, minB, promB;
    calcular_max_min_prom(A, &maxA, &minA, &promA);
    calcular_max_min_prom(B, &maxB, &minB, &promB);

    double coef = ((maxA * maxB) - (minA * minB)) / (promA * promB);

    transponer_matriz(B, Bt);
    multiplicar_matrices_bloques(A, Bt, AB);
    multiplicar_matrices_bloques(C, Bt, CBt);

    for (int i = 0; i < N * N; i++)
    {
        R[i] = coef * AB[i] + CBt[i];
    }
}

// Retorna el tiempo en segundos
double dwalltime()
{
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

int main()
{
    // Inicialización de matrices dinámicas
    double *A = (double *)malloc(N * N * sizeof(double));
    double *B = (double *)malloc(N * N * sizeof(double));
    double *C = (double *)malloc(N * N * sizeof(double));
    double *R = (double *)malloc(N * N * sizeof(double));

    // Modificamos la alocacion de memoria de las siguientes matrices
    // Pasandolo al "main" para que no influyan en el dwalltime()

    double *AB = (double *)malloc(N * N * sizeof(double));
    double *CBt = (double *)malloc(N * N * sizeof(double));
    double *Bt = (double *)malloc(N * N * sizeof(double));

    // Comprobaciones adicionales por si no se reservó memoria correctamente
    if (!A || !B || !C || !R || !AB || !CBt || !Bt)
    {
        printf("Error al reservar memoria\n");
        return 1;
    }

    // Inicialización de matrices con valores aleatorios
    // Como las transponemos despues la matriz necesaria, inicializamos de forma contigua en memoria
    for (int i = 0; i < N * N; i++)
    {
        A[i] = rand() % 10;
        B[i] = rand() % 10;
        C[i] = rand() % 10;
    }

    double tick = dwalltime();

    calcular_R(A, B, C, R, Bt, AB, CBt);
    double time = dwalltime() - tick;

    printf("Tiempo requerido: %f\n", time);

    free(A);
    free(B);
    free(C);
    free(R);
    free(AB);
    free(CBt);
    free(Bt);
    return 0;
}
