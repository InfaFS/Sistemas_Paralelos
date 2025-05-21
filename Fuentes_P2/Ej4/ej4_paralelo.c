#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <semaphore.h>

double *V;
int N, T, X, block_size, extra;
double sum_total = 0, prom_total = 0, min_tot = 9999, max_tot = -1; // esta es ocurrencias

// estaria bien que este aca? --> duda
sem_t sem_suma, sem_max, sem_min;

void *calculo_hilo(void *ptr)
{
    double min = 9999, max = -1, suma = 0;
    int id = *(int *)ptr; // Convierte el puntero a entero y le asigna a id el valor apuntado

    int inicio = id * block_size;
    int fin = inicio + block_size;

    if (id == T - 1) // Si es el último hilo, se le asigna el resto de la división
        fin += extra;

    for (int i = inicio; i < fin; i++)
    {
        if (min > V[i])
        {
            min = V[i];
        }

        if (V[i] > max)
        {
            max = V[i];
        }

        suma += V[i];
    }

    // usamos semaforos para cada valor

    // semaforo para el promedio
    sem_wait(&sem_suma);
    sum_total += suma;
    sem_post(&sem_suma);

    // semaforo para el max
    sem_wait(&sem_max);
    if (max > max_tot)
    {
        max_tot = max;
    }
    sem_post(&sem_max);

    // semaforo para el min
    sem_wait(&sem_min);
    if (min < min_tot)
    {
        min_tot = min;
    }
    sem_post(&sem_min);

    pthread_exit(0);
}

// Para calcular tiempo
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
    // Controla los argumentos al programa
    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0))
    {
        printf("\nUsar: %s n t x \n  n: Dimension del vector(n) \n  t: Cantidad de hilos \n", argv[0]);
        exit(1);
    }

    double suma_sec = 0, min_sec = 9999, max_sec = -1, prom_sec = 0;
    int i;
    int check = 1;
    double timetick;
    int ids[T];
    block_size = N / T;
    extra = N % T;

    pthread_attr_t attr;
    pthread_t threads[T];

    // Aloca memoria para el vector
    V = (double *)malloc(sizeof(double) * N);

    srand(time(NULL)); // Semilla para el rand
    for (i = 0; i < N; i++)
    {
        V[i] = (double)rand() / RAND_MAX * 10.0; // Random doubles entre 0.0 y 10.0
    }

    // Ejecucion secuencial
    timetick = dwalltime();

    for (i = 0; i < N; i++)
    {
        if (V[i] < min_sec)
        {
            min_sec = V[i];
        }

        if (V[i] > max_sec)
        {
            max_sec = V[i];
        }

        suma_sec += V[i];
    }
    prom_sec = suma_sec / N;

    // Calculamos el tiempo
    double timesec = dwalltime() - timetick;

    // Inicializamos el sem y los attr
    sem_init(&sem_max, 0, 1);
    sem_init(&sem_min, 0, 1);
    sem_init(&sem_suma, 0, 1);
    pthread_attr_init(&attr);

    // Crea los hilos
    for (i = 0; i < T; i++)
    {
        ids[i] = i;
        pthread_create(&threads[i], &attr, calculo_hilo, &ids[i]);
    }

    // El hilo principal espera a que terminen todos los hilos para informar el tiempo
    timetick = dwalltime();

    for (i = 0; i < T; i++)
    {
        pthread_join(threads[i], NULL);
    }
    prom_total = sum_total / N;

    double time = dwalltime() - timetick;

    // printf("prom sec= %f, prom para = %f\n", prom_sec, prom_total);
    // printf("min sec= %f, min para = %f\n", min_sec, min_tot);
    // printf("max sec= %f, max para = %f\n", max_sec, max_tot);
    // chequeamos que los valores sean los correctos
    // if ((prom_total != prom_sec) || (max_tot != max_sec) || (min_tot != min_sec))
    // {
    //     printf("Los valores entre las ejecuciones no coinciden");
    //     return 1;
    // }

    printf("Tiempo para secuencial(%f) || Tiempo para paralelo(%f)", timesec, time);
    // Destruimos el semaforo y liberamos la memoria
    sem_destroy(&sem_max);
    sem_destroy(&sem_min);
    sem_destroy(&sem_suma);
    free(V);
    return (0);
}