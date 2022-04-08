#include <cmath>

namespace funcs {

auto u_solve = [] (double t, double x) {
    return std::sin (x + t) - std::exp (x * t);
};

auto u_t_0 = [] (double t) { return std::sin (t) - 1; };
auto u_0_x = [] (double x) { return std::sin (x) - 1; };
auto f = [] (double t, double x) {
    return 2 * cos (x + t) - (x + t) * std::exp (x * t);
};

auto u_next = [] (double dx, double dt,
                  double u_k_m, double u_k_m1, double u_k1_m, double f_kh_mh) {
    const double a = 2 * dt * dx / (dx + dt);
    const double b = (dx - dt) / (dx + dt);

    return a * f_kh_mh + u_k_m + b * (u_k_m1 - u_k1_m);
};

} // namespace funcs