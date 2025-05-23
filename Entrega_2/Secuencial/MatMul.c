#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

void initvalmat(double *mat, int n, double val, int orden);
void matmulblks(double *a, double *b, double *c, int n, int bs);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);
void transponerMatriz(double *matriz, double *resultado, int n);
double calcularEscalar(double *A, double *B, int N);

double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

int main(int argc, char *argv[])
{
    double *A, *B, *C, *Bt, *R, *resultMatriz;
    int N, BS, i, j, offsetI, offsetJ;
    double timetick;
    double minA = INT_MAX, maxA = INT_MIN, minB = INT_MAX, maxB = INT_MIN;
    double promA = 0.0, promB = 0.0;
    double escalar = 0.0;
    double posA, posB;

    N = atoi(argv[1]);
    BS = atoi(argv[2]);

    if ((argc != 3) || (N <= 0) || (BS <= 0) || ((N % BS) != 0))
    {
        printf("\nParámetros inválidos. Debe utilizar %s N BS (N debe ser múltiplo de BS)\n", argv[0]);
        exit(1);
    }

    double size = N * N;

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

    // Inicializa las matrices A, B, C y D en 1
    initvalmat(A, N, 1.0, 0); // Orden 0 para que se inicialice en orden de filas
    initvalmat(B, N, 1.0, 1); // Orden 1 para que se inicialice en orden de columnas
    initvalmat(C, N, 1.0, 0); // Orden 0 para que se inicialice en orden de filas

    initvalmat(R, N, 0.0, 0);
    initvalmat(resultMatriz, N, 0.0, 0);

    timetick = dwalltime();

    // Trasposicion de matriz incluida en dwalltime para que cuente en el tiempo
    transponerMatriz(B, Bt, N); // B traspuesta
    // Se calcula el escalar
    escalar = calcularEscalar(A, B, N);

    // Multiplica A x B y lo guarda en resultMatriz
    matmulblks(A, B, resultMatriz, N, BS);

    // Multiplica Bt x C y lo guarda en R
    matmulblks(C, Bt, R, N, BS);

    // Se suman los resultados de las multiplicaciones de matrices y al mismo tiempo se multiplica por el escalar resultMatriz
    for (i = 0; i < N; i++)
    {
        offsetI = i * N;
        for (j = 0; j < N; j++)
        {
            R[offsetI + j] += resultMatriz[offsetI + j] * escalar;
        }
    }

    double workTime = dwalltime() - timetick;

    printf("Tiempo en segundos %f\n", workTime);

    free(A);
    free(B);
    free(C);
    free(Bt);
    free(R);
    free(resultMatriz);

    return (0);
}

// Inicializa una matriz de nxn con un valor val
void initvalmat(double *mat, int n, double val, int orden)
{
    int i, j;

    if (orden == 0)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mat[i * n + j] = val;
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mat[j * n + i] = val;
            }
        }
    }
}

// Realiza la multiplicacion de matrices en bloques de tamaño bs
void matmulblks(double *a, double *b, double *c, int n, int bs)
{
    int i, j, k, offsetI, offsetJ;

    for (i = 0; i < n; i += bs)
    {
        offsetI = i * n;
        for (j = 0; j < n; j += bs)
        {
            offsetJ = j * n;
            for (k = 0; k < n; k += bs)
            {
                blkmul(&a[offsetI + k], &b[offsetJ + k], &c[offsetI + j], n, bs);
            }
        }
    }
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
            // Modificacion con variable auxiliar para mejorar performance
            cblk[offsetI + j] += auxiliar;
        }
    }
}

double calcularEscalar(double *A, double *B, int N)
{
    int i, j, offsetI, offsetJ;
    double posA, posB;
    double minA = INT_MAX, maxA = INT_MIN, minB = INT_MAX, maxB = INT_MIN;
    double promA = 0.0, promB = 0.0;
    double size = N * N;

    // Procesar matriz A
    for (i = 0; i < N; i++)
    {
        offsetI = i * N;
        for (j = 0; j < N; j++)
        {
            posA = A[offsetI + j];

            if (posA < minA)
                minA = posA;
            if (posA > maxA)
                maxA = posA;
            promA += posA;
        }
    }

    // Procesar matriz B
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            offsetJ = j * N;
            posB = B[offsetJ + i];

            if (posB < minB)
                minB = posB;
            if (posB > maxB)
                maxB = posB;
            promB += posB;
        }
    }

    promA = promA / size;
    promB = promB / size;

    return (maxA * maxB - minA * minB) / (promA * promB);
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
