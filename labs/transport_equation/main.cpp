#include <iostream>
#include <mpi/mpi.h>
#include "lib.hpp"
/*

    unsigned M = 800;       // coord
    unsigned K = 1'000'000; // time

    double X_max = 6;       // X_min = 0
    double T_max = 6;       // T_min = 0

    O (h^2.68)

    time: 16 seconds

*/

int main (int argc, char* argv[]) {
// #define ONLY_SOLVE
//     check_single_thread_solution ();
// #undef ONLY_SOLVE

    unsigned M = 800;       // coord
    unsigned K = 1'000'000; // time

    double X_max = 6;       // X_min = 0
    double T_max = 6;       // T_min = 0

    trans_eq_task_t task {M, K, X_max, T_max};

    solve_trans_eq_parallel (task, &argc, &argv);


    /*
    
    Проблема: как передавать вычисленные блоки данных
    Возможные решения:
        1) Т.к. процесс с rank == 0 не передаёт свои блоки (все хранятся у него),
        то можно осздать второй поток, который параллельно сделает все запросы на блоки
        из других процессов. Но таких запросов м.б. очень много, поэтому появляется важная
        задача это оптимизировать.
        2) Представлять данные через MPI_Type_vector
    */
}