```cpp
// Блокирующие 
int MPI_Send (buffer, count, type, dest, tag, comm);
int MPI_Recv (buffer, count, type, source, tag, comm, status);

// Неблокирующие
int MPI_Isend (..., request)
int MPI_Irecv (..., request)

// Ожидание
int MPI_Wait (&req, &status);        // Блокирующая
int MPI_Test (&req, &flag, &status); // Неблокирующийся

int MPI_Ssend (синхронная)
int MPI_Bsend (буферезированный)

// Приём
MPI_ANY_SOURCE
MPI_ANY_TAG
MPI_STATUS_IGNORE

struct MPI_Status
double MPI_Wtime ();
double MPI_Wtick ();
```