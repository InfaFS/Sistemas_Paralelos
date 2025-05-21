#include <stdio.h>
#include <omp.h>

int main()
{
    int i, j;
    int N = 4;
#pragma omp parallel for private(i, j)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("Hilo %d → i = %d, j = %d\n", omp_get_thread_num(), i, j);
        }
    }

    return 0;
}