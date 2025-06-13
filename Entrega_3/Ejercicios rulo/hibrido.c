#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <mpi.h>
#include <sys/time.h>

#define BS 64
int N;

double dwalltime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
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
                suma += ablk[offsetI + k] * bblk[offsetJ + k];
            cblk[offsetI + j] += suma;
        }
    }
}

int main(int argc, char *argv[])
{
    int rank, size, i, j, k;

    if (argc != 3)
    {
        printf("Uso: %s <N> <T>\n", argv[0]);
        return EXIT_FAILURE;
    }

    N = atoi(argv[1]);
    int T = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    omp_set_num_threads(T);

    int rows_per_proc = N / size;
    int extra = N % size;
    int my_rows = rows_per_proc + (rank < extra ? 1 : 0);
    int my_start = rank * rows_per_proc + (rank < extra ? rank : extra);

    double *A = malloc(sizeof(double) * N * N);
    double *B = malloc(sizeof(double) * N * N);
    double *BT = malloc(sizeof(double) * N * N);
    double *C = malloc(sizeof(double) * N * N);
    double *R = malloc(sizeof(double) * N * my_rows);
    double *AB = malloc(sizeof(double) * N * my_rows);
    double *CBT = malloc(sizeof(double) * N * my_rows);

    if (rank == 0)
    {
        initvalmat(A, N, 1.0, 0);  // row-major
        initvalmat(B, N, 1.0, 1);  // column-major
        initvalmat(C, N, 1.0, 0);  // row-major
    }

    MPI_Bcast(A, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(C, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double maxA = -DBL_MAX, minA = DBL_MAX, sumA = 0.0;
    double maxB = -DBL_MAX, minB = DBL_MAX, sumB = 0.0;
    double escalar = 0.0;
    double start = dwalltime();

#pragma omp parallel private(i, j)
    {
        // Transposición de B
#pragma omp for schedule(static)
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                BT[j * N + i] = B[i * N + j];

        // Reducción local A
#pragma omp for reduction(+ : sumA) reduction(min : minA) reduction(max : maxA) schedule(static)
        for (i = my_start; i < my_start + my_rows; i++)
            for (j = 0; j < N; j++)
            {
                double val = A[i * N + j];
                sumA += val;
                if (val < minA) minA = val;
                if (val > maxA) maxA = val;
            }

        // Reducción local B
#pragma omp for reduction(+ : sumB) reduction(min : minB) reduction(max : maxB) schedule(static)
        for (i = my_start; i < my_start + my_rows; i++)
            for (j = 0; j < N; j++)
            {
                double val = B[j * N + i];  // column-major
                sumB += val;
                if (val < minB) minB = val;
                if (val > maxB) maxB = val;
            }
    }

    // Reducción global MPI
    double gSumA, gMaxA, gMinA, gSumB, gMaxB, gMinB;
    MPI_Reduce(&sumA, &gSumA, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&maxA, &gMaxA, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&minA, &gMinA, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sumB, &gSumB, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&maxB, &gMaxB, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&minB, &gMinB, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        double promA = gSumA / (N * N);
        double promB = gSumB / (N * N);
        escalar = (gMaxA * gMaxB - gMinA * gMinB) / (promA * promB);
    }

    MPI_Bcast(&escalar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#pragma omp parallel for collapse(2) schedule(static)
    for (i = my_start; i < my_start + my_rows; i += BS)
        for (j = 0; j < N; j += BS)
            for (k = 0; k < N; k += BS)
                blkmul(&A[i * N + k], &B[j * N + k], &AB[(i - my_start) * N + j], N, BS);

#pragma omp parallel for collapse(2) schedule(static)
    for (i = my_start; i < my_start + my_rows; i += BS)
        for (j = 0; j < N; j += BS)
            for (k = 0; k < N; k += BS)
                blkmul(&C[i * N + k], &BT[j * N + k], &CBT[(i - my_start) * N + j], N, BS);

#pragma omp parallel for private(j)
    for (i = 0; i < my_rows; i++)
        for (j = 0; j < N; j++)
            R[i * N + j] = escalar * AB[i * N + j] + CBT[i * N + j];

    if (rank == 0)
    {
        int *recvcounts = malloc(size * sizeof(int));
        int *displs = malloc(size * sizeof(int));
        int offset = 0;
        for (i = 0; i < size; i++)
        {
            int rows = rows_per_proc + (i < extra ? 1 : 0);
            recvcounts[i] = rows * N;
            displs[i] = offset;
            offset += rows * N;
        }

        double *Rglobal = malloc(sizeof(double) * N * N);
        MPI_Gatherv(R, my_rows * N, MPI_DOUBLE, Rglobal, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        printf("Tiempo en segundos: %f\n", dwalltime() - start);
        free(Rglobal);
        free(recvcounts);
        free(displs);
    }
    else
    {
        MPI_Gatherv(R, my_rows * N, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    free(A); free(B); free(BT); free(C); free(R); free(AB); free(CBT);

    MPI_Finalize();
    return 0;
}
