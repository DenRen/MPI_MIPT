#include <vector>
#include <cmath>
#include <cassert>

#include "print_lib.hpp"
#include "func_lib.hpp"

#define HOST

#ifdef HOST
    #define DUMP(obj) std::cout << #obj ": " << obj << "\n"
#else
    #define DUMP(obj)
#endif

struct trans_eq_task_t {
    int x_size;     // Number points of x coordinate
    int t_size;     // Number points of time
    double x_max;   // x coord belongs to [0, x_max]
    double t_max;   // time belongs to [0, t_max]

    void check_valid () const
    {
        if (x_size <= 1) {
            throw std::invalid_argument ("x_size must be greater then 1");
        }
        if (t_size <= 1) {
            throw std::invalid_argument ("t_size must be greater then 1");
        }
    }

    trans_eq_task_t (int num_points_coord_x,
                     int num_points_time,
                     double x_max,
                     double t_max) :
        x_size (num_points_coord_x),
        t_size (num_points_time),
        x_max (x_max),
        t_max (t_max)
    {
        check_valid ();
    }

    double
    calc_dx () const noexcept {
        return x_max / (x_size - 1);
    }

    double
    calc_dt () const noexcept {
        return t_max / (t_size - 1);
    }
};

#define CHECK_ERR(err_code) check_error (err_code, __LINE__, __FILE__)

void
check_error (int err_code, int line, const char* file_name) {
    if (err_code != MPI_SUCCESS) {
        char err_str[MPI_MAX_ERROR_STRING] = "";
        int str_len = 0;
        MPI_Error_string (err_code, err_str, &str_len);

        std::string err_msg {file_name};
        err_msg += ": " + std::to_string (line);
        err_msg += ": ";
        err_msg += err_str;

        throw std::runtime_error (err_msg);
    }
}

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

        double u_left = u[cur_shift] = u_t_0 (t);
        for (int x_i = 0; x_i < num_points_coord_x - 1; ++x_i) {
            std::size_t prev_pos = prev_shift + x_i;
            std::size_t cur_pos = cur_shift + x_i;

            double u_k_m  = u[prev_pos];
            double u_k_m1 = u[prev_pos + 1];
            double u_k1_m = u_left;

            double x = x_i * dx;
            double f_kh_mh = f (t - dt/2, x + dx/2);

            u_left = u[cur_pos + 1] = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }
    }

    return u;
}

template <typename F, typename U_NEXT>
void
calc_u (int x_size,
        int t_size,
        int x_min, int t_min,
        double dx, double dt,
        std::vector <double>& u_buf,
        F f,
        U_NEXT u_next)
{
    // Calc other u
    for (int time_i = 1; time_i < t_size; ++time_i) {
        double t = t_min + time_i * dt;
        std::size_t cur_shift = time_i * x_size;
        std::size_t prev_shift = cur_shift - x_size;

        double u_k1_m = u_buf[cur_shift];
        for (int x_i = 0; x_i < x_size - 1; ++x_i) {
            std::size_t prev_pos = prev_shift + x_i;
            std::size_t cur_pos = cur_shift + x_i;

            double u_k_m  = u_buf[prev_pos];
            double u_k_m1 = u_buf[prev_pos + 1];

            double x = x_min + x_i * dx;
            double f_kh_mh = f (t - dt/2, x + dx/2);

            u_k1_m = u_buf[cur_pos + 1]
                   = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }
    }
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
solve_transport_eq (int num_points_coord_x,
                    int num_points_time,
                    double x_max,
                    double t_max)
{
    const double dx = x_max / (num_points_coord_x - 1);
    const double dt = t_max / (num_points_time - 1);


    using funcs::u_t_0;
    using funcs::u_0_x;
    using funcs::f;
    using funcs::u_next;
    return calc_u (num_points_coord_x, num_points_time, dx, dt,
                   u_t_0, u_0_x, f, u_next);
}

std::vector <double>
get_solve_ref (int num_points_coord_x,
               int num_points_time,
               double x_max,
               double t_max)
{
    const double dx = x_max / (num_points_coord_x - 1);
    const double dt = t_max / (num_points_time - 1);

    using funcs::u_solve;
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

template <typename T,
          typename U_SOLVE>
T
calc_err_use_solve_func (const std::vector <T>& u,
                         int num_points_coord_x,
                         int num_points_time,
                         double x_max,
                         double t_max,
                         U_SOLVE u_solve)
{
    const T dx = x_max / (num_points_coord_x - 1);
    const T dt = t_max / (num_points_time - 1);

    T max_err = 0;
    for (int t_i = 0; t_i < num_points_time; ++t_i) {
        std::size_t cur_shift = t_i * num_points_coord_x;

        T t = t_i * dt;
        for (int x_i = 0; x_i < num_points_coord_x; ++x_i) {
            T x = x_i * dx;
            T u_ref = u_solve (t, x);

            T err = std::abs (u_ref - u[cur_shift + x_i]);
            max_err = std::max (max_err, err);
        }
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
    // u (t, x) = sin (x + t) - e^(xt)
    auto u_solve = [] (double t, double x) {
        return std::sin (x + t) - std::exp (x * t);
    };

    std::vector <std::pair <double, double>> h_err (num_steps);
    for (int step = 0; step < num_steps; ++step) {
        int num_points_coord_x =
            num_points_coord_x_min + num_points_coord_dx * step;

        auto u = solve_transport_eq (num_points_coord_x,
                                     num_points_time, x_max, t_max);

        double h = x_max / (num_points_coord_x - 1);
        double err = calc_err_use_solve_func (u, num_points_coord_x,
                                              num_points_time, x_max, t_max,
                                              u_solve);

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

// #define ONLY_SOLVE
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

void
solve_trans_eq_parallel_even (const trans_eq_task_t& task,
                              int x_chunk_size,
                              int t_chunk_size,
                              int num_threads,
                              int rank)
{
    assert (!(rank & 1));

    int n = num_threads / 2;
    double dx = task.calc_dx ();
    double dt = task.calc_dt ();

    using funcs::u_0_x;
    using funcs::u_t_0;
    using funcs::f;
    using funcs::u_next;

    std::vector <double> u_buf ((x_chunk_size + 1) * (t_chunk_size + 1));
    if (rank == 0) {
        for (int i = 0; i < n; ++i) {
            int x_i_min = i * x_chunk_size;
            int x_i_max = std::min (x_i_min + x_chunk_size, task.x_size);

            // Fill u (t = 0, x)
            for (int x_i = x_i_min; x_i <= x_i_max; ++x_i) {
                double x = x_i * dx;
                u_buf[x_i - x_i_min] = u_0_x (x);
            }

            // Fill u (x_min, t = 0)
            int t_i_min = i * t_chunk_size;
            int t_i_max = std::min (t_i_min + t_chunk_size, task.t_size);
            for (int t_i = t_i_min + 1; t_i <= t_i_max; ++t_i) {
                double t = t_i * dt;
                int pos = (x_chunk_size + 1) * (t_i - t_i_min);
                u_buf[pos] = u_t_0 (t);
            }

            // Calc u_buf
            double x_min = dx * x_i_min;
            double t_min = dt * t_i_min;
            calc_u (x_chunk_size + 1, t_chunk_size + 1, x_min, t_min,
                    dx, dt, u_buf, f, u_next);

            int dest = i == 0 ? rank + 1 : rank + 2;
            MPI_Send (u_buf.data () + (x_chunk_size + 1) * t_chunk_size, x_chunk_size,
                      MPI_DOUBLE, dest, MPI_ANY_TAG, MPI_COMM_WORLD);
        }
    } else {

    }

    for (int i = 0; i < n; ++i) {
        int x_i_max = std::min ((i + 1) * x_chunk_size, task.x_size);
        for (int x_i = i * x_chunk_size; x_i < x_i_max; ++x_i) {

        }
    }
}


void
solve_trans_eq_parallel_odd (const trans_eq_task_t& task,
                             int num_threads,
                             int rank)
{}

int
calc_chunk_size (int full_size,
                 int num_chunks)
{
    return  full_size / num_chunks +
           (full_size % num_chunks != 0);
}

void
solve_trans_eq_parallel (const trans_eq_task_t& task,
                         int* argc_ptr,
                         char** argv_ptr[])
{
    MPI_Init (argc_ptr, argv_ptr);

    int num_threads = 0, rank = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_threads);
    if (num_threads & 1) {
        throw std::invalid_argument ("Number threads must be only EVEN");
    }

    int x_chunk_size = calc_chunk_size (task.x_size, num_threads);
    int t_chunk_size = calc_chunk_size (task.t_size, num_threads);


    if (rank == 0) {

    }

    if (rank & 1) {
        // solve_trans_eq_parallel_odd (task, num_threads, rank);
    } else {
        solve_trans_eq_parallel_even (task, x_chunk_size, t_chunk_size,
                                      num_threads, rank);
    }

    MPI_Finalize ();
}