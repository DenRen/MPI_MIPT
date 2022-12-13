#include <array>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <mpi.h>

// #define CHECK_CORRECT
// #define WRITE_RES_TO_FILE

template <std::size_t x_size>
void CheckOnCorrect(std::size_t y_size,
                    const std::vector<std::array<double, x_size>>& a)
{
    // Test correct on sequential algo
    std::vector<std::array<double, x_size>> a_seq(y_size);
    for (size_t y = 0; y < y_size; ++y)
        for (size_t x = 0; x < x_size; ++x)
            a_seq[y][x] = 10 * y + x;

    for (size_t y = 0; y < y_size; ++y)
        for (size_t x = 0; x < x_size; ++x)
            a_seq[y][x] = std::sin(0.000001 * a_seq[y][x]);

    // Find error
    for (size_t y = 0; y < y_size; ++y)
        for (size_t x = 0; x < x_size; ++x)
            if (a_seq[y][x] != a[y][x])
            {
                printf ("In parallel algorithm erorr exists\n");

                y = y_size;
                break;
            }
}

template <std::size_t x_size>
void Dump2File(const std::vector<std::array<double, x_size>>& a)
{
    // Save result to file
    FILE* file = fopen("result.txt", "w");

    const auto y_size = a.size();
    for(int i = 0; i < y_size; ++i){
        for (int j = 0; j < x_size; ++j)
            fprintf(file, "%f ", a[i][j]);

        fprintf(file, "\n");
    }

    fclose(file);
}

template <std::size_t x_size>
void CalcParallel(std::size_t Y_size, int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0, num_threads = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_threads);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    // Main worker keep full grid
    // Secondary worker keep just a part of grid
    const int y_begin = rank * Y_size / num_threads;
    const int y_end = (rank + 1) * Y_size / num_threads;

    const int y_size = y_end - y_begin;
    const int a_y_size = rank == 0 ? Y_size : y_size;    // Vertical array size
    std::vector<std::array<double, x_size>> a(a_y_size);

    const double time_begin = MPI_Wtime();
    for (size_t y = 0; y < y_size; ++y)
        for (size_t x = 0; x < x_size; ++x)
            a[y][x] = 10 * (y + y_begin) + x;


    for (int y = 0; y < y_size; ++y)
        for (int x = 0; x < x_size; ++x)
            a[y][x] = sin(0.000001 * a[y][x]);

    if (rank != 0)
    {
        MPI_Send(a[0].data(), x_size * y_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    else
    {
        for (int k = 1; k < num_threads; ++k)
        {
            const int y_begin = k * Y_size / num_threads;
            const int y_end = (k + 1) * Y_size / num_threads;
            MPI_Recv(a[y_begin].data(), x_size * (y_end - y_begin), MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
        }
    }

    const double time_end = MPI_Wtime();

    if (rank == 0)
    {
        const double dt_parallel = time_end - time_begin;
        printf("par time: %lf\n", dt_parallel);

#ifdef CHECK_CORRECT
        CheckOnCorrect<x_size>(Y_size, a);
#endif
#ifdef WRITE_RES_TO_FILE
        Dump2File<x_size>(a);
#endif
    }

    MPI_Finalize();
}

int main(int argc, char **argv)
{
    const int y_size = 6400;
    const int x_size = 10000;

    CalcParallel<x_size>(y_size, argc, argv);

    return 0;
}
