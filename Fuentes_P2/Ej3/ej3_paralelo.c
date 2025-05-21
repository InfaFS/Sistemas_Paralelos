#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>

int *V;
int N, T, X, block_size, extra;
int total_count; // esta es ocurrencias

// estaria bien que este aca? --> duda
pthread_mutex_t counter_value_lock;

void *suma_vector(void *ptr)
{
    int inter_count = 0;
    int id = *(int *)ptr; // Convierte el puntero a entero y le asigna a id el valor apuntado

    int inicio = id * block_size;
    int fin = inicio + block_size;

    if (id == T - 1) // Si es el último hilo, se le asigna el resto de la división
        fin += extra;

    for (int i = inicio; i < fin; i++)
    {
        if (V[i] == X) // probamos que cuente los que son multiplos de 3
        {
            // aca tiene que mutearse
            inter_count++;
        }
    }

    pthread_mutex_lock(&counter_value_lock);
    total_count += inter_count;
    pthread_mutex_unlock(&counter_value_lock);

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
    // Controla los argumentos al programa, agregamos que ponga que valor X quiere buscar
    if ((argc != 4) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0) || ((X = atoi(argv[3])) <= 0))
    {
        printf("\nUsar: %s n t x \n  n: Dimension del vector(n) \n  t: Cantidad de hilos \n  x: Numero a buscar \n", argv[0]);
        exit(1);
    }

    int i, sec_count;
    int check = 1;
    double timetick;
    int ids[T];
    block_size = N / T;
    extra = N % T;
    total_count = 0;
    sec_count = 0;
    pthread_attr_t attr;
    pthread_t threads[T];

    // Aloca memoria para el vector
    V = (int *)malloc(sizeof(int) * N);

    srand(time(NULL)); // Semilla para el rand
    for (i = 0; i < N; i++)
    {
        V[i] = rand() % 10; // Random numbers entre el 0 y el 9
    }

    // Ejecucion secuencial
    timetick = dwalltime();
    for (i = 0; i < N; i++)
    {
        if (V[i] == X)
        {
            sec_count++;
        }
    }

    // Calculamos el tiempo
    double timesec = dwalltime() - timetick;

    // Inicializamos el mutex y los attr
    pthread_mutex_init(&counter_value_lock, NULL);
    pthread_attr_init(&attr);

    // Crea los hilos
    for (i = 0; i < T; i++)
    {
        ids[i] = i;
        pthread_create(&threads[i], &attr, suma_vector, &ids[i]);
    }

    // El hilo principal espera a que terminen todos los hilos para informar el tiempo
    timetick = dwalltime();

    for (i = 0; i < T; i++)
    {
        pthread_join(threads[i], NULL);
    }

    double time = dwalltime() - timetick;
    // print_m(); // Imprime el vector resultante

    // chequeamos que los valores sean los correctos
    if (total_count != sec_count)
    {
        printf("Los valores entre las ejecuciones no coinciden");
        return 1;
    }

    printf("Tiempo en segundos(secuencial) %f\nTiempo en segundos(paralelo) %f\n", timesec, time);

    // Destruimos el semaforo y liberamos la memoria
    pthread_mutex_destroy(&counter_value_lock);
    free(V);
    return (0);
}