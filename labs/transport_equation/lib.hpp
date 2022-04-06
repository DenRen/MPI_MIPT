#include <vector>
#include <cmath>

#include "print_lib.hpp"

#define HOST

#ifdef HOST
    #define DUMP(obj) std::cout << #obj ": " << obj << "\n"
#else
    #define DUMP(obj)
#endif

template <typename U_T_0, typename U_0_X, typename F,
          typename U_NEXT>
std::vector <double>
calc_u (int num_points_coord_x,
        int num_points_time,
        double dx,
        double dt,
        U_T_0 u_t_0,
        U_0_X u_0_x,
        F f,
        U_NEXT u_next)
{
    std::vector <double> u (num_points_time * num_points_coord_x);

    // Fill u (t = 0, x)
    for (int i = 0; i < num_points_coord_x; ++i) {
        double x = i * dx;
        u[i] = u_0_x (x);
    }

    // Calc other u
    for (int time_i = 1; time_i < num_points_time; ++time_i) {
        double t = time_i * dt;
        std::size_t cur_shift = time_i * num_points_coord_x;
        std::size_t prev_shift = cur_shift - num_points_coord_x;

        double u_left = u[cur_shift] = u_t_0 (0);
        for (int x_i = 0; x_i < num_points_coord_x - 1; ++x_i) {
            std::size_t prev_pos = prev_shift + x_i;
            std::size_t cur_pos = cur_shift + x_i;

            double u_k_m  = u[prev_pos];
            double u_k_m1 = u[prev_pos + 1];
            double u_k1_m = u_left;

            double x = x_i * dx;
            double f_kh_mh = f (t, x);

            u_left = u[cur_pos + 1] = u_next (u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }
    }

    return u;
}

std::vector <double>
solve_transport_eq (int num_points_coord_x,
                    int num_points_time,
                    double x_max,
                    double t_max)
{
    const double dx = x_max / (num_points_coord_x - 1);
    const double dt = t_max / (num_points_time - 1);

    // u (t, x) = sin (x + t) - e^(xt)
    auto u_t_0 = [] (double t) { return std::sin (t) - 1; };
    auto u_0_x = [] (double x) { return std::sin (x) - 1; };
    auto f = [] (double t, double x) {
        return 2 * cos (x + t) - (x + t) * std::exp (x * t);
    };

    auto u_next = [dx, dt] (double u_k_m, double u_k_m1, double u_k1_m, double f_kh_mh) {
        const double a = 2 * dt * dx / (dx + dt);
        const double b = (dx - dt) / (dx + dt);

        return a * f_kh_mh + u_k_m + b * (u_k_m1 - u_k1_m);
    };

    return calc_u (num_points_coord_x, num_points_time, dx, dt,
                   u_t_0, u_0_x, f, u_next);
}

template <typename U_SOLVE>
std::vector <double>
calc_u (int num_points_coord_x,
        int num_points_time,
        double dx,
        double dt,
        U_SOLVE u_solve)
{
    std::vector <double> u (num_points_coord_x * num_points_time);

    for (int t_i = 0; t_i < num_points_time; ++t_i) {
        std::size_t cur_shift = t_i * num_points_coord_x;

        double t = t_i * dt;
        for (int x_i = 0; x_i < num_points_coord_x; ++x_i) {
            double x = x_i * dx;
            u[cur_shift + x_i] = u_solve (t, x);
        }
    }

    return u;
}

std::vector <double>
get_solve_ref (int num_points_coord_x,
               int num_points_time,
               double x_max,
               double t_max)
{
    const double dx = x_max / (num_points_coord_x - 1);
    const double dt = t_max / (num_points_time - 1);

    // u (t, x) = sin (x + t) - e^(xt)
    auto u_solve = [] (double t, double x) {
        return std::sin (x + t) - std::exp (x * t);
    };

    return calc_u (num_points_coord_x, num_points_time, dx, dt,
                   u_solve);
}

template <typename C>
void
check_eq_size (const C& cont_first,
               const C& cont_second)
{
    if (cont_first.size () != cont_second.size ()) {
        throw std::runtime_error ("Sizes of containers are difference");
    }
}

template <typename T>
T
calc_err (const std::vector <T>& u_calc,
          const std::vector <T>& u_ref)
{
    check_eq_size (u_calc, u_ref);

    T max_err = 0;
    for (std::size_t i = 0; i < u_calc.size (); ++i) {
        T err = std::abs (u_calc[i] - u_ref[i]);
        max_err = std::max (max_err, err);
    }

    return max_err;
}

//                        h      err
std::vector <std::pair <double, double>>
calc_h_err_vec (int num_points_coord_x_min,
                int num_points_coord_dx,
                int num_steps,
                int num_points_time,
                double x_max,
                double t_max)
{
    std::vector <std::pair <double, double>> h_err (num_steps);

    for (int step = 0; step < num_steps; ++step) {
        int num_points_coord_x =
            num_points_coord_x_min + num_points_coord_dx * step;

        auto u = solve_transport_eq (num_points_coord_x,
                                     num_points_time, x_max, t_max);
        auto u_ref = get_solve_ref (num_points_coord_x,
                                    num_points_time, x_max, t_max);
        
        double h = x_max / (num_points_coord_x - 1);
        double err = calc_err (u, u_ref);
        h_err[step] = {h, err};
    }

    return h_err;
}

// y = k*x + b
// return k, b
template <typename T>
std::pair <T, T>
linearize (const std::vector <T>& x,
           const std::vector <T>& y)
{
    check_eq_size (x, y);

    const auto n = x.size ();

    T sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (std::size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x2 += std::pow (x[i], 2);
    }
 
    T k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - std::pow (sum_x, 2));
    T b = (sum_y - k * sum_x) / n;
    
    return {k, b};
}

template <typename T>
std::pair <T, T>
linearize (const std::vector <std::pair <T, T>>& x_y)
{
    const auto n = x_y.size ();

    T sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (std::size_t i = 0; i < n; ++i) {
        auto[x, y] = x_y[i];

        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += std::pow (x, 2);
    }
 
    T k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - std::pow (sum_x, 2));
    T b = (sum_y - k * sum_x) / n;
    
    return {k, b};
}