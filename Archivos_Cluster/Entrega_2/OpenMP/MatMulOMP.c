#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <sys/time.h>

#define BS 64
int N;

double dwalltime()
{
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

void initvalmat(double *mat, int n, double val, int orden)
{
    int i, j;
    if (orden == 0)
    {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[i * n + j] = val;
    }
    else
    {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[j * n + i] = val;
    }
}

// A y BT son row-major, B y C column-major
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs)
{
    int i, j, k, offsetI, offsetJ;
    double suma;

    for (i = 0; i < bs; i++)
    {
        offsetI = i * n;
        for (j = 0; j < bs; j++)
        {
            suma = 0;
            offsetJ = j * n;
            for (k = 0; k < bs; k++)
            {
                suma += ablk[offsetI + k] * bblk[offsetJ + k];
            }

            cblk[offsetI + j] += suma;
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Uso: %s <N> <T>\n", argv[0]);
        return EXIT_FAILURE;
    }

    N = atoi(argv[1]);
    int T = atoi(argv[2]);

    omp_set_num_threads(T);

    double *A = malloc(sizeof(double) * N * N);
    double *B = malloc(sizeof(double) * N * N);
    double *BT = malloc(sizeof(double) * N * N);
    double *C = malloc(sizeof(double) * N * N);
    double *R = malloc(sizeof(double) * N * N);
    double *AB = malloc(sizeof(double) * N * N);

    initvalmat(A, N, 1.0, 0);
    initvalmat(B, N, 1.0, 1);
    initvalmat(C, N, 1.0, 0);
    initvalmat(R, N, 0.0, 0);
    initvalmat(AB, N, 0.0, 0);

    if (!A || !B || !C || !R || !BT || !AB)
    {
        fprintf(stderr, "Error al alocar memoria\n");
        return EXIT_FAILURE;
    }

    double maxA = -DBL_MAX, maxB = -DBL_MAX;
    double minA = DBL_MAX, minB = DBL_MAX;
    double sumA = 0.0, sumB = 0.0;
    double escalar = 0.0;
    int i, j, k, offsetI, offsetJ;
    double posA, posB;
    double timetick = dwalltime();

// Una única sección paralela para todas las operaciones
#pragma omp parallel private(i, j, k, offsetI, offsetJ, posA, posB)
    {
// Fase 0: Transposición de matriz B
#pragma omp for schedule(static)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                BT[j * N + i] = B[i * N + j];
            }
        }

// Fase 1a: Calcular min, max, sum para matriz A - Primer for separado
#pragma omp for reduction(+ : sumA) reduction(min : minA) reduction(max : maxA) schedule(static)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                posA = A[i * N + j];
                sumA += posA;
                if (posA < minA)
                    minA = posA;
                if (posA > maxA)
                    maxA = posA;
            }
        }

// Fase 1b: Calcular min, max, sum para matriz B - Segundo for separado
#pragma omp for reduction(+ : sumB) reduction(min : minB) reduction(max : maxB) schedule(static)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                posB = B[j * N + i]; // column-major
                sumB += posB;
                if (posB < minB)
                    minB = posB;
                if (posB > maxB)
                    maxB = posB;
            }
        }

// Calcular escalar - Solo un hilo lo hace
#pragma omp single
        {
            double promA = sumA / (N * N);
            double promB = sumB / (N * N);
            escalar = (maxA * maxB - minA * minB) / (promA * promB);
        }

// Fase 2: Multiplicación A*B -> AB
#pragma omp for nowait schedule(static)
        for (i = 0; i < N; i += BS)
        {
            for (j = 0; j < N; j += BS)
            {
                for (k = 0; k < N; k += BS)
                {
                    blkmul(&A[i * N + k], &B[j * N + k], &AB[i * N + j], N, BS);
                }
            }
        }

// Fase 3: Multiplicación BT*C -> BTC
#pragma omp for schedule(static)
        for (i = 0; i < N; i += BS)
        {
            for (j = 0; j < N; j += BS)
            {
                for (k = 0; k < N; k += BS)
                {
                    blkmul(&BT[i * N + k], &C[j * N + k], &R[i * N + j], N, BS);
                }
            }
        }

// Fase 4: Suma final R = BTC + escalar*AB
#pragma omp for nowait schedule(static)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                R[i * N + j] += escalar * AB[i * N + j];
            }
        }
    }

    double workTime = dwalltime() - timetick;
    printf("Tiempo en segundos: %f\n", workTime);

    free(A);
    free(B);
    free(BT);
    free(C);
    free(R);
    free(AB);
    return 0;
}
