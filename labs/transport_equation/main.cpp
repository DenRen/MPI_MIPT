#include <iostream>

#include "print_lib.hpp"
#include "lib.hpp"

int main () {
    std::ios_base::sync_with_stdio (false);

    unsigned M = 10000 / 10; // coord
    unsigned K = 10000 / 10; // time

    double X_max = 1; // X_min = 0
    double T_max = 2;  // T_min = 0

    auto u = solve_transport_eq (M, K, X_max, T_max);
    auto u_ref = get_solve_ref (M, K, X_max, T_max);

    std::cout << u << '\n';
    std::cout << u_ref << '\n';

    for (int i = 0; i < u.size (); ++i)
        std::cout << std::fabs (u[i] - u_ref[i]) << ' ';
}