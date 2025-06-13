#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <mpi.h>

#define BS 64
int N;

double dwalltime()
{
    double t;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    t = tv.tv_sec + tv.tv_usec / 1000000.0;
    return t;
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
    int rank, size, i, j, k, offsetI, offsetJ;
    double maxA, minA, sumA, maxB, minB, sumB, escalar;
    double posA, posB, tstart, tend;

    if (argc != 2)
    {
        printf("Uso: %s <N>\n", argv[0]);
        return EXIT_FAILURE;
    }

    N = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows_per_proc = N / size;
    int extra = N % size;
    int my_rows = rows_per_proc + (rank < extra ? 1 : 0);
    int my_start = rank * rows_per_proc + (rank < extra ? rank : extra);

    double *A = malloc(N * N * sizeof(double));
    double *B = malloc(N * N * sizeof(double));
    double *BT = malloc(N * N * sizeof(double));
    double *C = malloc(N * N * sizeof(double));
    double *AB = calloc(N * my_rows, sizeof(double));
    double *CBT = calloc(N * my_rows, sizeof(double));
    double *Rlocal = calloc(N * my_rows, sizeof(double));
    double *R = NULL;

    if (rank == 0)
        R = malloc(N * N * sizeof(double));

    // Inicialización
    if (rank == 0)
    {
        for (i = 0; i < N * N; i++)
        {
            A[i] = 1.0;
            B[i] = 1.0;
            C[i] = 1.0;
        }
    }

    MPI_Bcast(A, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(C, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    maxA = -DBL_MAX;
    minA = DBL_MAX;
    sumA = 0.0;
    maxB = -DBL_MAX;
    minB = DBL_MAX;
    sumB = 0.0;

    tstart = dwalltime();

    // Transposición paralela de B en BT
    for (i = my_start; i < my_start + my_rows; i++)
        for (j = 0; j < N; j++)
            BT[j * N + i] = B[i * N + j];

    // Calcular min, max, promedio locales de A
    for (i = my_start; i < my_start + my_rows; i++)
    {
        for (j = 0; j < N; j++)
        {
            posA = A[i * N + j];

            if (posA < minA)
                minA = posA;
            if (posA > maxA)
                maxA = posA;
            sumA += posA;
        }
    }

    // Calcular min, max, promedio locales de B
    for (i = my_start; i < my_start + my_rows; i++)
    {
        for (j = 0; j < N; j++)
        {
            posB = B[j * N + i];

            if (posB < minB)
                minB = posB;
            if (posB > maxB)
                maxB = posB;
            sumB += posB;
        }
    }

    double gMaxA, gMinA, gSumA;
    double gMaxB, gMinB, gSumB;

    MPI_Reduce(&maxA, &gMaxA, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&minA, &gMinA, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sumA, &gSumA, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&maxB, &gMaxB, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&minB, &gMinB, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sumB, &gSumB, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        double promA = gSumA / (N * N);
        double promB = gSumB / (N * N);
        escalar = (gMaxA * gMaxB - gMinA * gMinB) / (promA * promB);
    }

    MPI_Bcast(&escalar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calcular AB = A × B (cada proceso sus filas)
    for (i = my_start; i < my_start + my_rows; i += BS)
        for (j = 0; j < N; j += BS)
            for (k = 0; k < N; k += BS)
                blkmul(&A[i * N + k], &B[j * N + k], &AB[(i - my_start) * N + j], N, BS);

    // Calcular CBT = C × BT (cada proceso sus filas)
    for (i = my_start; i < my_start + my_rows; i += BS)
        for (j = 0; j < N; j += BS)
            for (k = 0; k < N; k += BS)
                blkmul(&C[i * N + k], &BT[j * N + k], &CBT[(i - my_start) * N + j], N, BS);

    // Rlocal = escalar × AB + CBT
    for (i = 0; i < my_rows * N; i++)
        Rlocal[i] = escalar * AB[i] + CBT[i];

    int *recvcounts = NULL, *displs = NULL;
    if (rank == 0)
    {
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        int disp = 0;
        for (i = 0; i < size; i++)
        {
            int rows = rows_per_proc + (i < extra ? 1 : 0);
            recvcounts[i] = rows * N;
            displs[i] = disp;
            disp += rows * N;
        }
    }

    MPI_Gatherv(Rlocal, my_rows * N, MPI_DOUBLE, R, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    tend = dwalltime();

    if (rank == 0)
        printf("Tiempo en segundos: %f\n", tend - tstart);

    free(A);
    free(B);
    free(BT);
    free(C);
    free(AB);
    free(CBT);
    free(Rlocal);
    if (rank == 0)
    {
        free(R);
        free(recvcounts);
        free(displs);
    }

    MPI_Finalize();
    return 0;
}
