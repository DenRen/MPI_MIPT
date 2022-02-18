MPI - Message Passing Interface
===
[source](https://habr.com/ru/post/548266/)

## Part 1

* MIMD - Multiple Sintruction Multiple Data
* SPMD - Single Program Multiple Data (по умолчанию)

```bash
sudo apt update
sudo apt install -y gcc g++ mpich
```

* *Коммуникатор* - объект через который общается определенная группа порожденных процессов.
    Его тип **MPI_Comm**.
    При старте программы все процессы работают под единым коммуникатором **MPI_COMM_WORLD**.
    Коммуникатор текущего процесса **MPI_COMM_SELF**, ни одного процесса **MPI_COMM_NULL**.

* *Сообщение* - набор данных некоторого типа, который передается при коммуникации процессов.
  * Номер процесса отправителя
  * Номер процесса получателя
  * Id сообщения
  * Коммуникатор
  * Тег (от 0 до 32767, **MPI_TAG_UB**)

В любой **MPI** программе должно быть:
```cpp
int MPI_Init (int* argc, char*** argv)
int MPI_Finalize (void);
```

Самый простой пример (*main.cpp*):
```cpp
#include <stdio.h>
#include "mpi.h"

int main (int argc, char **argv) {
  printf ("Before MPI_INIT\n");
  MPI_Init (&argc, &argv);
  printf ("Parallel sect\n");
  MPI_Finalize ();
  printf ("After MPI_FINALIZE\n");

  return 0;
}
```

```bash
mpic++ main.cpp -o main
mpiexec -n 2 ./main
```

## Part 2
```cpp
int MPI_Comm_size (MPI_Comm comm, int* size);
int MPI_Comm_rank (MPI_Comm comm, int* rank);

double MPI_Wtime (void);
double MPI_Wtick (void);

int MPI_Get_Processor_name (char* name, int* len);
```

**MPI_WTIME_IS_GLOBAL** - синхронизированы ли по времени процессоры (0 или 1)

## Part 3
```cpp
int MPI_Send (void* buf, int count, MPI_Datatype datatype, int dest, 
						  int msgtag, MPI_Comm comm);
int MPI_Recv (void* buf, int count, MPI_Datatype datatype, int source,
						  int tag, MPI_Comm comm, MPI_Status* status);
```

MPI_Datatype:
* MPI_CHAR
* MPI_SHORT
* MPI_INT
* ...
* MPI_INT16_t
* MPI_C_COMPLEX (*float _Complex*)

Пример использования *MPI_Status*:
```cpp
struct MPI_Status status;
MPI_Recv (&buffer, 1, MPI_Float, MPI_ANY_SOURCE,
          MPI_ANY_TAG, MPI_COMM_WORLD, &status);

int tag = status.MPI_SOURCE;
int source = status.MPI_SOURCE;
```
* MPI_ANY_SOURCE
* MPI_ANY_TAG

```cpp
int MPI_Get_count (MPI_Status* status, MPI_Datatype datatype, int* count);
int MPI_Get_elements (MPI_Status* status, MPI_Datatype datatype, int* count);
int MPI_Probe (int source, int tag, MPI_Comm comm, MPI_Status* status);
```

Отправить в никуда при помощи номера процесса *MPI_PROC_NULL*

Пример:
```cpp
#include <iostream>
#include "mpi.h"

using namespace std;

void show_arr(int* arr, int size)
{
	for(int i=0; i < size; i++) cout << arr[i] << " ";
	cout << endl;
}

int main(int argc, char **argv)
{
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0)
	{
		int* arr = new int[size];
		for(int i=0; i < size; i++) arr[i] = i;
		for(int i=1; i < size; i++) MPI_Send(arr, i, MPI_INT, i, 5, MPI_COMM_WORLD);
	}
	else
	{
		int count;
		MPI_Status status;

		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &count);
		int* buf = new int[count];

		MPI_Recv(buf, count, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		
		cout << "Process:" << rank << " || Count: " << count << " || Array: ";
		show_arr(buf, count);
	}
	MPI_Finalize();
	return 0;
}
```

Разные виды MPI_Send (аргументы аналогичны):
* *MPI_Bsend* - буфферезирует у получателя
* *MPI_Ssend* - блокирует, пока получатель полностью не получил сообщение
* *MPI_Rsend* - по готовности. Получатель должен быть готов получить сообщение