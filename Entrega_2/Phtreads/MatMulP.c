#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <pthread.h>

#define MAX_THREADS 8

typedef struct
{
    int id;
} thread_data_t;

void initvalmat(double *mat, int n, double val, int orden);
void *worker_thread(void *arg);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);
void transponerMatriz(double *matriz, double *resultado, int n);
double dwalltime();

// variables globales
double *A, *B, *C, *Bt, *R, *resultMatriz;
double escalar = 0.0;
pthread_mutex_t mutexCalc;
pthread_barrier_t barrier;
double minA = INT_MAX, maxA = INT_MIN, minB = INT_MAX, maxB = INT_MIN;
double promA = 0.0, promB = 0.0;
int N, T, BS;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("\nParámetros inválidos. Debe utilizar %s N BS\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    T = atoi(argv[2]);
    BS = 64;

    if ((N <= 0) || (BS <= 0) || ((N % BS) != 0) || (T > 8))
    {
        printf("\nError: N debe ser múltiplo de BS y T <=8\n");
        exit(1);
    }
    pthread_attr_t attr;
    int i;
    double timetick;
    pthread_t threads[T];
    double size;
    int ids[T];
    size = N * N;

    A = (double *)malloc(size * sizeof(double));
    B = (double *)malloc(size * sizeof(double));
    C = (double *)malloc(size * sizeof(double));
    Bt = (double *)malloc(size * sizeof(double));
    R = (double *)malloc(size * sizeof(double));
    resultMatriz = (double *)malloc(size * sizeof(double));

    if (A == NULL || B == NULL || C == NULL || Bt == NULL || R == NULL || resultMatriz == NULL)
    {
        fprintf(stderr, "Error al alocar memoria\n");
        return EXIT_FAILURE;
    }

    initvalmat(A, N, 1.0, 0); // Filas
    initvalmat(B, N, 1.0, 1); // Columnas
    initvalmat(C, N, 1.0, 0); // Filas
    initvalmat(R, N, 0.0, 0);
    initvalmat(resultMatriz, N, 0.0, 0);

    pthread_attr_init(&attr);
    pthread_mutex_init(&mutexCalc, NULL);
    pthread_barrier_init(&barrier, NULL, T);

    timetick = dwalltime();

    // Crear threads para realizar todas las operaciones
    for (i = 0; i < T; i++)
    {
        ids[i] = i;
        pthread_create(&threads[i], NULL, worker_thread, &ids[i]);
    }

    // Esperar a que todos los threads terminen
    for (i = 0; i < T; i++)
        pthread_join(threads[i], NULL);

    double workTime = dwalltime() - timetick;
    printf("Tiempo en segundos %f\n", workTime);

    free(A);
    free(B);
    free(C);
    free(Bt);
    free(R);
    free(resultMatriz);

    pthread_mutex_destroy(&mutexCalc);
    pthread_barrier_destroy(&barrier);
    pthread_attr_destroy(&attr);
    return 0;
}

// ----------------------- FUNCIONES ------------------------

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

void *worker_thread(void *arg)
{
    int id = *((int *)arg);
    double size = N * N;

    // Como nos vamos a mover por filas y columnas, usamos N para delimitar
    int start = id * (N / T);
    int end = (id == T - 1) ? N : start + (N / T);
    int i, j, k, offsetI, offsetJ;

    // Variables locales para el cálculo del escalar
    double localMinA = INT_MAX, localMaxA = INT_MIN, localMinB = INT_MAX, localMaxB = INT_MIN;
    double localSumA = 0.0, localSumB = 0.0;
    double posA, posB;

    // PARTE 0: Calcular escalar - Procesar matriz A
    for (i = start; i < end; i++)
    {
        offsetI = i * N;
        for (j = 0; j < N; j++)
        {
            posA = A[offsetI + j];

            if (posA < localMinA)
                localMinA = posA;
            if (posA > localMaxA)
                localMaxA = posA;
            localSumA += posA;
        }
    }

    // Procesar matriz B
    for (i = start; i < end; i++)
    {
        for (j = 0; j < N; j++)
        {
            offsetJ = j * N;
            posB = B[offsetJ + i];

            if (posB < localMinB)
                localMinB = posB;
            if (posB > localMaxB)
                localMaxB = posB;
            localSumB += posB;
        }
    }

    // Actualizar las variables globales con mutex
    pthread_mutex_lock(&mutexCalc);
    if (localMinA < minA)
        minA = localMinA;
    if (localMaxA > maxA)
        maxA = localMaxA;
    if (localMinB < minB)
        minB = localMinB;
    if (localMaxB > maxB)
        maxB = localMaxB;
    promA += localSumA;
    promB += localSumB;
    pthread_mutex_unlock(&mutexCalc);

    // PARTE 1: Transponer la matriz B -> Bt en paralelo antes de empezar a calcular el escalar
    // Aprovechamos tambien para usar la barrera siguiente y no tene que usar 2 (causando mayor overhead)
    for (i = start; i < end; i++)
    {
        for (j = 0; j < N; j++)
        {
            Bt[j * N + i] = B[i * N + j];
        }
    }

    // Barrera 1: Esperar a que todos los hilos terminen el cálculo
    pthread_barrier_wait(&barrier);

    // Solo el hilo 0 calcula el escalar final
    if (id == 0)
    {
        promA /= size;
        promB /= size;
        escalar = (maxA * maxB - minA * minB) / (promA * promB);
    }

    // PARTE 2: Multiplicación de matrices A*B -> resultMatriz
    for (i = start; i < end; i += BS)
    {
        offsetI = i * N;
        for (j = 0; j < N; j += BS)
        {
            offsetJ = j * N;
            for (k = 0; k < N; k += BS)
            {
                blkmul(&A[offsetI + k], &B[offsetJ + k], &resultMatriz[offsetI + j], N, BS);
            }
        }
    }

    // PARTE 3: Multiplicación de matrices C*Bt -> R
    for (i = start; i < end; i += BS)
    {
        offsetI = i * N;
        for (j = 0; j < N; j += BS)
        {
            offsetJ = j * N;
            for (k = 0; k < N; k += BS)
            {
                blkmul(&C[offsetI + k], &Bt[offsetJ + k], &R[offsetI + j], N, BS);
            }
        }
    }

    // Barrera 3: Esperar a que todas las multiplicaciones estén completas
    pthread_barrier_wait(&barrier);

    // PARTE 4: Suma final R += resultMatriz * escalar
    for (i = start; i < end; i++)
    {
        offsetI = i * N;
        for (j = 0; j < N; j++)
        {
            R[offsetI + j] += resultMatriz[offsetI + j] * escalar;
        }
    }

    pthread_exit(NULL);
}

void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs)
{
    int i, j, k, offsetI, offsetJ;
    for (i = 0; i < bs; i++)
    {
        offsetI = i * n;
        for (j = 0; j < bs; j++)
        {
            offsetJ = j * n;
            double auxiliar = 0.0;
            for (k = 0; k < bs; k++)
            {
                auxiliar += ablk[offsetI + k] * bblk[offsetJ + k];
            }
            cblk[offsetI + j] += auxiliar;
        }
    }
}

void transponerMatriz(double *matriz, double *resultado, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            resultado[j * n + i] = matriz[i * n + j];
        }
    }
}

double dwalltime()
{
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}
