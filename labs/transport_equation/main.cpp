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

    unsigned M = 2000; // coord
    unsigned K = 2000; // time

    double X_max = 2;       // X_min = 0
    double T_max = 2;       // T_min = 0

    trans_eq_task_t task {M, K, X_max, T_max};

    solve_trans_eq_parallel (task, &argc, &argv);

    // auto u = solve_transport_eq (M, K, X_max, T_max);
    // { volatile double x = u[1]; }
    // print_2d_array (u.data (), K, u.size ());
}