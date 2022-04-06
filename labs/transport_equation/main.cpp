#include <iostream>

#include "lib.hpp"
/*

    unsigned M = 800;       // coord
    unsigned K = 1'000'000; // time

    double X_max = 6;       // X_min = 0
    double T_max = 6;       // T_min = 0

    O (h^2.68)

    time: 16 seconds

*/

#define ONLY_SOLVE

void
check_single_thread_solution () {
    unsigned M = 800;       // coord
    unsigned K = 1'000'000; // time

    double X_max = 6;       // X_min = 0
    double T_max = 6;       // T_min = 0


#ifdef ONLY_SOLVE
    auto u = solve_transport_eq (M, K, X_max, T_max);
    { volatile double x = u[1]; }
#else
    auto h_err = calc_h_err_vec (M, M / 10, 3, K, X_max, T_max);
    for (auto[h, err] : h_err) {
        std::cout << std::log (h) << ' ' << std::log (err) << '\n';
    }

    for (auto&[h, err] : h_err) {
        h = std::log (h);
        err = std::log (err);
    }
    std::cout << "O(h^" << linearize (h_err).first << ")\n";
#endif
}

#undef ONLY_SOLVE

int main () {
    check_single_thread_solution ();
}