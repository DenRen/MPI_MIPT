#include "treq_lib.hpp"

int main (int argc, char* argv[]) {
// #define ONLY_SOLVE
//     check_single_thread_solution ();
// #undef ONLY_SOLVE

    int M = 30000; // coord
    int K = 30000; // time

    double X_max = 2;       // X_min = 0
    double T_max = 2;       // T_min = 0

    treq::trans_eq_task_t task {M, K, X_max, T_max};

#if 1
    solve_trans_eq_parallel (task, &argc, &argv);
#else
    auto u = solve_transport_eq (M, K, X_max, T_max);
    // { volatile double x = u[1]; }
    print_2d_array (u.data (), K, u.size ());
#endif
}