#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>
#include <omp.h>

#define BS 64
#define COORDINATOR 0

void initvalmat(double *mat, int n, double val, int orden, int stripSize);

void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);

int main(int argc, char *argv[])
{
    double *A, *B, *C, *B_transpuesta, *R, *resultMatriz;
    int N, T, bs, i, j, k, offsetI, offsetJ, numProcs, rank, stripSize, provided;
    double local_min[2] = {INT_MAX, INT_MAX}, local_max[2] = {INT_MIN, INT_MIN}, min[2] = {INT_MAX, INT_MAX}, max[2] = {INT_MIN, INT_MIN};
    double local_prom[2] = {0, 0}, prom[2] = {0, 0};
    double escalar = 0.0;
    double posA, posB;
    MPI_Status status;
    double commTimes[10], maxCommTimes[10], minCommTimes[10], commTime, totalTime;

    N = atoi(argv[1]);
    T = atoi(argv[2]);

    if ((argc != 3) || (N <= 0) || (N % 64 != 0) || (T <= 0))
    {
        printf("\nParámetros inválidos. Debe utilizar %s N (N debe ser múltiplo de 64) T (cantidad de hilos por proceso)\n", argv[0]);
        exit(1);
    }

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (N % numProcs != 0)
    {
        printf("El tamaño de la matriz debe ser múltiplo del numero de procesos.\n");
        exit(1);
    }

    double size = N * N;

    stripSize = N / numProcs;

    double sizeWorker = N * stripSize;

    bs = (N / (numProcs * T) < BS ? N / (numProcs * T) : BS);

    if (rank == COORDINATOR)
    {
        A = (double *)malloc(size * sizeof(double));
        C = (double *)malloc(size * sizeof(double));
        R = (double *)malloc(size * sizeof(double));
        resultMatriz = (double *)malloc(size * sizeof(double));
    }
    else
    {
        A = (double *)malloc(sizeWorker * sizeof(double));
        C = (double *)malloc(sizeWorker * sizeof(double));
        R = (double *)malloc(sizeWorker * sizeof(double));
        resultMatriz = (double *)malloc(sizeWorker * sizeof(double));
    }

    B = (double *)malloc(size * sizeof(double));
    B_transpuesta = (double *)malloc(size * sizeof(double));

    if (rank == COORDINATOR)
    {
        // Inicializa las matrices A, B y C en 1
        initvalmat(A, N, 1.0, 0, N); // Orden 0 para que se inicialice en orden de filas
        initvalmat(B, N, 1.0, 1, N); // Orden 1 para que se inicialice en orden de columnas
        initvalmat(C, N, 1.0, 0, N); // Orden 0 para que se inicialice en orden de filas
        initvalmat(R, N, 0.0, 0, N);
        initvalmat(resultMatriz, N, 0.0, 0, N);
    }
    else
    {
        initvalmat(R, N, 0.0, 0, stripSize);
        initvalmat(resultMatriz, N, 0.0, 0, stripSize);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    commTimes[0] = MPI_Wtime();

    // Se reparten las matrices
    MPI_Scatter(A, sizeWorker, MPI_DOUBLE, A, sizeWorker, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Scatter(C, sizeWorker, MPI_DOUBLE, C, sizeWorker, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

    MPI_Bcast(B, size, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

    commTimes[1] = MPI_Wtime();

#pragma omp parallel num_threads(T) private(i, j, k, offsetI, offsetJ, posA, posB)
    {
        // Se calcula la transpuesta de B de forma distribuida y paralelizada
        // Cada proceso calcula una porción de columnas de B_transpuesta
        int start_col = rank * stripSize;
        int end_col = start_col + stripSize;

#pragma omp for schedule(static) nowait
        for (i = 0; i < N; i++)
        {
            for (j = start_col; j < end_col; j++)
            {
                B_transpuesta[j * N + i] = B[i * N + j];
            }
        }

#pragma omp master
        {
            commTimes[2] = MPI_Wtime();

            // Reunir todas las porciones de B_transpuesta
            MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                          B_transpuesta, stripSize * N, MPI_DOUBLE, MPI_COMM_WORLD);

            commTimes[3] = MPI_Wtime();
        }

#pragma omp barrier // Esperar a que el master termine el Allgather

// Se calculan los minimos, maximos y promedios de las matriz A
#pragma omp for reduction(min : local_min[0]) reduction(max : local_max[0]) reduction(+ : local_prom[0]) nowait schedule(static)
        for (i = 0; i < stripSize; i++)
        {
            offsetI = i * N;
            for (j = 0; j < N; j++)
            {
                posA = A[offsetI + j];

                if (posA < local_min[0])
                {
                    local_min[0] = posA;
                }

                if (posA > local_max[0])
                {
                    local_max[0] = posA;
                }

                local_prom[0] += posA;
            }
        }

        int inicio = stripSize * rank;
        int fin = inicio + stripSize;
// Se calculan los minimos, maximos y promedios de las matriz B
#pragma omp for reduction(min : local_min[1]) reduction(max : local_max[1]) reduction(+ : local_prom[1]) schedule(static)
        for (i = inicio; i < fin; i++)
        {
            offsetI = i * N;
            for (j = 0; j < N; j++)
            {
                posB = B[offsetI + j];

                if (posB < local_min[1])
                {
                    local_min[1] = posB;
                }

                if (posB > local_max[1])
                {
                    local_max[1] = posB;
                }

                local_prom[1] += posB;
            }
        }

// Se realiza la multiplicacion de matrices A y B
#pragma omp for nowait schedule(static)
        for (i = 0; i < stripSize; i += bs)
        {
            offsetI = i * N;
            for (j = 0; j < N; j += bs)
            {
                offsetJ = j * N;
                for (k = 0; k < N; k += bs)
                {
                    blkmul(&A[offsetI + k], &B[offsetJ + k], &resultMatriz[offsetI + j], N, bs);
                }
            }
        }

// Se realiza la multiplicacion de matrices C y B_transpuesta
#pragma omp for nowait schedule(static)
        for (i = 0; i < stripSize; i += bs)
        {
            offsetI = i * N;
            for (j = 0; j < N; j += bs)
            {
                offsetJ = j * N;
                for (k = 0; k < N; k += bs)
                {
                    blkmul(&C[offsetI + k], &B_transpuesta[offsetJ + k], &R[offsetI + j], N, bs);
                }
            }
        }

#pragma omp master // Solo el hilo master realiza la comunicacion
        {
            commTimes[4] = MPI_Wtime();

            // Se recolectan los minimos, maximos y promedios
            MPI_Reduce(&local_min, &min, 2, MPI_DOUBLE, MPI_MIN, COORDINATOR, MPI_COMM_WORLD);
            MPI_Reduce(&local_max, &max, 2, MPI_DOUBLE, MPI_MAX, COORDINATOR, MPI_COMM_WORLD);
            MPI_Reduce(&local_prom, &prom, 2, MPI_DOUBLE, MPI_SUM, COORDINATOR, MPI_COMM_WORLD);

            commTimes[5] = MPI_Wtime();

            // Se calcula el escalar
            if (rank == COORDINATOR)
            {
                prom[0] = prom[0] / (size);
                prom[1] = prom[1] / (size);

                escalar = (max[0] * max[1] - min[0] * min[1]) / (prom[0] * prom[1]);
            }

            // Se transmite el escalar
            commTimes[6] = MPI_Wtime();

            MPI_Bcast(&escalar, 1, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

            commTimes[7] = MPI_Wtime();
        }

#pragma omp barrier // Los hilos esperan a que el hilo maestro tenga el escalar

// Se suman los resultados de las multiplicaciones de matrices y al mismo tiempo se multiplica por el escalar resultMatriz
#pragma omp for nowait schedule(static)
        for (i = 0; i < stripSize; i++)
        {
            offsetI = i * N;
            for (j = 0; j < N; j++)
            {
                R[offsetI + j] += resultMatriz[offsetI + j] * escalar;
            }
        }
    }

    // Se recolectan los resultados
    commTimes[8] = MPI_Wtime();

    MPI_Gather(R, sizeWorker, MPI_DOUBLE, R, sizeWorker, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

    commTimes[9] = MPI_Wtime();

    MPI_Reduce(commTimes, minCommTimes, 10, MPI_DOUBLE, MPI_MIN, COORDINATOR, MPI_COMM_WORLD);
    MPI_Reduce(commTimes, maxCommTimes, 10, MPI_DOUBLE, MPI_MAX, COORDINATOR, MPI_COMM_WORLD);

    MPI_Finalize();

    if (rank == COORDINATOR)
    {
        totalTime = maxCommTimes[9] - minCommTimes[0];
        commTime = (maxCommTimes[1] - minCommTimes[0]) + (maxCommTimes[3] - minCommTimes[2]) + (maxCommTimes[5] - minCommTimes[4]) + (maxCommTimes[7] - minCommTimes[6]) + (maxCommTimes[9] - minCommTimes[8]);

        printf("Calculo de vector (N=%d)\tTiempo total=%lf\tTiempo comunicacion=%lf\n", N, totalTime, commTime);
    }

    free(A);
    free(B);
    free(C);
    free(B_transpuesta);
    free(R);
    free(resultMatriz);

    return (0);
}

// Inicializa una matriz de nxn con un valor val
void initvalmat(double *mat, int n, double val, int orden, int stripSize)
{
    int i, j, offsetI;

    if (orden == 0)
    {
        for (i = 0; i < stripSize; i++)
        {
            offsetI = i * n;
            for (j = 0; j < n; j++)
            {
                mat[offsetI + j] = val;
            }
        }
    }
    else
    {
        for (i = 0; i < stripSize; i++)
        {
            for (j = 0; j < n; j++)
            {
                mat[j * n + i] = val;
            }
        }
    }
}

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

