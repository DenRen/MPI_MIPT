#include <iostream>

#include "lib.hpp"

int main () {
    std::ios_base::sync_with_stdio (false);

    unsigned M = 10; // coord
    unsigned K = 100000; // time

    double X_max = 10;  // X_min = 0
    double T_max = 1; // T_min = 0

    // auto u = solve_transport_eq (M, K, X_max, T_max);
    // auto u_ref = get_solve_ref (M, K, X_max, T_max);

    // std::cout << u << '\n';
    // std::cout << u_ref << '\n';

    auto h_err = calc_h_err_vec (M, M / 4, 20, K, X_max, T_max);
    for (auto[h, err] : h_err) {
        // std::cout << std::log (h) << ' ' << std::log (err) << '\n';
    }

    for (auto&[h, err] : h_err) {
        h = std::log (h);
        err = std::log (err);
    }
    std::cout << "O(h^" << linearize (h_err).first << ")\n";
}