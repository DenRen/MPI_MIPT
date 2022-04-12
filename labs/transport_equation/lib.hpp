#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <thread>
#include <mpi/mpi.h>

#include "copy_transform.hpp"
#include "func_lib.hpp"
#include "print_lib.hpp"
// #include "solver.hpp"

#define HOST

#ifdef HOST
    #define DUMP(obj) std::cerr << #obj ": " << obj << "\n"
#else
    #define DUMP(obj)
#endif

struct trans_eq_task_t {
    int x_size;     // Number points of x coordinate
    int t_size;     // Number points of time
    double x_min = 0;
    double x_max;   // x coord belongs to [0, x_max]
    double t_min = 0;
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

    std::pair <double, double>
    calc_dx_dt () const noexcept {
        return {calc_dx (), calc_dt ()};
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
          typename U_next>
std::vector <double>
calc_u (int num_points_coord_x,
        int num_points_time,
        double dx,
        double dt,
        U_T_0 u_t_0,
        U_0_X u_0_x,
        F f,
        U_next u_next)
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

template <typename F, typename U_next>
void
calc_u_part_buf (int x_chunk_size,
                 int t_chunk_size,
                 int buf_x_size,
                 int x_arr_begin, int t_arr_begin,
                 int x_i_min, int t_i_min,
                 double dx, double dt,
                 std::vector <double>& u_buf,
                 F f,
                 U_next u_next)
{
    std::size_t buf_offset = (x_i_min - x_arr_begin) + (t_i_min - t_arr_begin) * buf_x_size;

    // Calc other u
    for (int time_i = 1; time_i < t_chunk_size; ++time_i) {
        double t = (t_i_min + time_i) * dt;
        DUMP (t);
        std::size_t cur_shift = time_i * buf_x_size + buf_offset;
        std::size_t prev_shift = cur_shift - buf_x_size;

        double u_k1_m = u_buf[cur_shift];
        for (int x_i = 0; x_i < x_chunk_size - 1; ++x_i) {
            std::size_t prev_pos = prev_shift + x_i;
            std::size_t cur_pos = cur_shift + x_i;

            double u_k_m  = u_buf[prev_pos];
            double u_k_m1 = u_buf[prev_pos + 1];

            double x = (x_i_min + x_i) * dx;
            DUMP (x);
            double f_kh_mh = f (t - dt/2, x + dx/2);

            u_k1_m = u_buf[cur_pos + 1]
                   = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }
    }
}

template <typename F, typename U_next>
void
calc_u_full_buf (int x_chunk_size,
                 int t_chunk_size,
                 int buf_x_size,
                 int x_min, int t_min,
                 double dx, double dt,
                 std::vector <double>& u_buf,
                 F f,
                 U_next u_next)
{
    calc_u_part_buf (x_chunk_size, t_chunk_size, buf_x_size,
                     0, 0, x_min, t_min, dx, dt, u_buf, f, u_next);
}

template <typename F, typename U_next>
void
calc_u (int x_size,
        int t_size,
        int x_min, int t_min,
        double dx, double dt,
        std::vector <double>& u_buf,
        F f,
        U_next u_next)
{
    return calc_u_full_buf (x_size, t_size, x_size, 0,
                            x_min, t_min, dx, dt, u_buf, f, u_next);
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

#define TAG_BORDER_COND 1
#define TAG_SAVE_ON_HOST 2

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

int
calc_chunk_size (int full_size,
                 int num_chunks)
{
    return  full_size / num_chunks +
           (full_size % num_chunks != 0);
}

template <typename T>
T
calc_end_index (T i_cur,
                T step,
                T i_end)
{
    if (i_cur + step >= i_end) {
        return i_end;
    }

    return i_cur + step;
}

template <typename T>
T
calc_max_step (T i_cur,
               T step,
               T i_end)
{
    return calc_end_index (i_cur, step, i_end) - i_cur;
}

/*
b solve_trans_eq_parallel_zero_rank
b work_non_zero_rank_grid_border
r

*/

namespace treq {

struct rect_t {
    int x_i_min, t_i_min;
    int x_chunk_size, t_chunk_size;

    int
    get_x_i_end () const noexcept {
        return x_i_min + x_chunk_size;
    }

    int
    get_t_i_end () const noexcept {
        return t_i_min + t_chunk_size;
    }

    std::size_t
    get_offset_rect (std::size_t buf_x_size) const noexcept {
        return x_i_min + t_i_min * buf_x_size;
    }

    std::size_t
    get_offset_rect (std::size_t buf_x_size,
                     const rect_t& zero_rect) const noexcept {
        return (x_i_min - zero_rect.x_i_min) + (t_i_min - zero_rect.t_i_min) * buf_x_size;
    }
}; // struct rect_t

struct area_params_t {

    rect_t zero_rect;
    int x_i_end, t_i_end;

    int
    get_zero_chunk_x_size () const noexcept {
        return std::min (x_i_end - zero_rect.x_i_min, zero_rect.x_chunk_size);
    }

    int
    get_zero_chunk_t_size () const noexcept {
        return std::min (t_i_end - zero_rect.t_i_min, zero_rect.t_chunk_size);
    }

    int
    get_zero_and_right_rect_x_size () const noexcept {
        return x_i_end - zero_rect.x_i_min;
    }

    int
    get_zero_and_right_rect_size () const noexcept {
        return get_zero_chunk_t_size () * get_zero_and_right_rect_x_size ();
    }

    int
    get_up_rect_x_size () const noexcept {
        return zero_rect.x_chunk_size;
    }

    int
    get_up_rect_size () const noexcept {
        int area_t_size = t_i_end - zero_rect.t_i_min;
        if (area_t_size <= zero_rect.t_chunk_size) {
            return 0;
        } else {
            return (area_t_size - zero_rect.t_chunk_size) * get_up_rect_x_size ();
        }
    }

    // This rectangle is assumed to exist, otherwise UB
    rect_t
    get_right_rect (int i_chunk_l) const noexcept {
        rect_t right_rect = {
            .x_i_min = zero_rect.x_i_min + i_chunk_l * zero_rect.x_chunk_size,
            .t_i_min = zero_rect.t_i_min,
            .t_chunk_size = zero_rect.t_chunk_size
        };

        if (right_rect.x_i_min + zero_rect.x_chunk_size > x_i_end) {
            right_rect.x_chunk_size = x_i_end - right_rect.x_i_min;
            if (right_rect.x_chunk_size < 0) {
                DUMP (x_i_end - right_rect.x_i_min);
            }

        } else {
            right_rect.x_chunk_size = zero_rect.x_chunk_size;
        }

        return right_rect;
    } // rect_t get_right_rect (int i_chunk)

    rect_t
    get_up_rect (int i_chunk_l) const noexcept {
        rect_t up_rect = {
            .x_i_min = zero_rect.x_i_min,
            .t_i_min = zero_rect.t_i_min + i_chunk_l * zero_rect.t_chunk_size,
            .x_chunk_size = zero_rect.x_chunk_size
        };

        if (up_rect.t_i_min + zero_rect.t_chunk_size > t_i_end) {
            up_rect.t_chunk_size = t_i_end - up_rect.t_i_min;
        } else {
            up_rect.t_chunk_size = zero_rect.t_chunk_size;
        }

        return up_rect;
    } // rect_t get_up_rect (int i_chunk)
}; // struct area_params_t

struct map_manager_t {
    int x_size;
    int t_size;

    double x_min;
    double x_max;
    double t_min;
    double t_max;

    int num_threads;
    int num_areas;

    map_manager_t (int x_size,
                   int t_size,
                   double x_min,
                   double x_max,
                   double t_min,
                   double t_max,
                   int num_threads,
                   int num_areas) :
        x_size (x_size),
        t_size (t_size),
        x_min (x_min),
        x_max (x_max),
        t_min (t_min),
        t_max (t_max),
        num_threads (num_threads),
        num_areas (num_areas)
    {
        if (x_size <= 1 || t_size <= 1 || num_threads <= 0 ||
            num_areas < num_threads) {
            throw std::invalid_argument ("Incorrect params");
        }
    }

    map_manager_t (const trans_eq_task_t& task,
                   int num_threads,
                   int num_areas) :
        map_manager_t (task.x_size, task.t_size,
                       task.x_min, task.x_max,
                       task.t_min, task.x_max,
                       num_threads, num_areas)
    {}

    area_params_t
    get_area_params (int i_area) const noexcept {
        assert (i_area >= 0);

        area_params_t params = {};
        rect_t& zero_rect = params.zero_rect;

        int x_chunk_size = calc_chunk_size (x_size, num_areas);
        int t_chunk_size = calc_chunk_size (t_size, num_areas);

        zero_rect.x_chunk_size = calc_chunk_size (x_size, num_areas);
        zero_rect.t_chunk_size = calc_chunk_size (t_size, num_areas);

        zero_rect.x_i_min = i_area * x_chunk_size;
        zero_rect.t_i_min = i_area * t_chunk_size;

        if (x_chunk_size + zero_rect.x_i_min <= x_size) {
            zero_rect.x_chunk_size = x_chunk_size;
            zero_rect.t_chunk_size = t_chunk_size;
        } else {
            zero_rect.x_chunk_size = x_size - zero_rect.x_i_min;
            zero_rect.t_chunk_size = t_size - zero_rect.t_i_min;
        }

        params.x_i_end = x_size;
        params.t_i_end = t_size;

        return params;
    }

    double
    calc_dx () const noexcept {
        return (x_max - x_min) / (x_size - 1);
    }

    double
    calc_dt () const noexcept {
        return (t_max - t_min) / (t_size - 1);
    }

    std::pair <double, double>
    calc_dx_dt () const noexcept {
        return {calc_dx (), calc_dt ()};
    }
}; // struct map_manager_t

struct trans_eq_solver {
    const map_manager_t map_mgr;

    trans_eq_solver (const trans_eq_task_t& task,
                     int num_threads,
                     int koef_mult_threads) :
        map_mgr (task, num_threads, koef_mult_threads * num_threads)
    {}

    void
    solve (int rank) {
        if (rank == 0) {
            std::vector <double> u_buf (map_mgr.x_size * map_mgr.t_size);

            if (map_mgr.num_threads == 1) {
                calc_zero_rank_func_border (u_buf);
                receive_u_bufs (u_buf);
            } else {
                auto compute_func = [&] {
                    calc_zero_rank_func_border (u_buf);
                    calc_zero_rank_grid_border (u_buf);
                };

                auto bufferize_func = [&] {
                    receive_u_bufs (u_buf);
                };

                std::thread compute_thread {compute_func};
                std::thread bufferer_thread {bufferize_func};

                compute_thread.join ();
                bufferer_thread.join ();
            }

            { volatile double tr = u_buf[3]; }

            // print_2d_array (u_buf.data (), map_mgr.x_size, u_buf.size ());
        } else {
            calc_non_zero_rank_grid_border (rank);
        }
    } // void solve (int rank)

    void
    calc_zero_rank_func_border (std::vector <double>& u_buf) {
        const area_params_t area_params = map_mgr.get_area_params (0);
        const auto[dx, dt] = map_mgr.calc_dx_dt ();
        std::size_t u_buf_x_size = map_mgr.x_size;
        std::vector <double> u_t_buf;

        // 1) Calc u, where need only functional border conditions
        // Fill u(t = 0, x)
        for (int x_i = 0; x_i < area_params.get_zero_chunk_x_size (); ++x_i) {
            double x = x_i * dx;
            u_buf[x_i] = funcs::u_0_x (x);
        }

        // Fill u(t, x = 0)
        for (int t_i = 1; t_i < area_params.get_zero_chunk_t_size (); ++t_i) {
            double t = t_i * dt;
            u_buf[t_i * u_buf_x_size] = funcs::u_t_0 (t);
        }

        // Calc u on (x_size - 1) * (t_size - 1)
        const rect_t& zero_rect = area_params.zero_rect;
        calc_u_rect (u_buf.data () + u_buf_x_size + 1, u_buf_x_size,
                     zero_rect.x_i_min + 1, zero_rect.t_i_min + 1,
                     zero_rect.x_chunk_size - 1, zero_rect.t_chunk_size - 1,
                     funcs::f, funcs::u_next);

        // Calc u rights and ups in zero area
        const int next_rank = 1;
        for (int i_chunk = 1; i_chunk < map_mgr.num_areas; ++i_chunk) {
            // Fill u(t = 0, x)
            const rect_t rect_right = area_params.get_right_rect (i_chunk);
            for (int x_i = rect_right.x_i_min; x_i < rect_right.get_x_i_end (); ++x_i) {
                double x = x_i * dx;
                u_buf[x_i] = funcs::u_0_x (x);
            }

            // Calc right rect u
            double* u_buf_right_begin = u_buf.data () + rect_right.x_i_min;
            calc_u_rect (u_buf_right_begin + u_buf_x_size, u_buf_x_size,
                         rect_right.x_i_min, rect_right.t_i_min + 1,
                         rect_right.x_chunk_size, rect_right.t_chunk_size - 1,
                         funcs::f, funcs::u_next);

            // Send up
            if (i_chunk == 1) { // right_rect.x_chunk_size + 1
                MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size - 1,
                          rect_right.x_chunk_size + 1,
                          MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            } else {
                MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size,
                          rect_right.x_chunk_size,
                          MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }

            // Fill u(t = 0, x)
            const rect_t rect_up = area_params.get_up_rect (i_chunk);
            for (int t_i = rect_up.t_i_min; t_i < rect_up.get_t_i_end (); ++t_i) {
                double t = t_i * dt;
                u_buf[t_i * u_buf_x_size] = funcs::u_t_0 (t);
            }

            // Calc up rect u
            double* u_buf_up_begin = u_buf.data () + rect_up.t_i_min * u_buf_x_size;
            calc_u_rect (u_buf_up_begin + 1, u_buf_x_size,
                         rect_up.x_i_min + 1, rect_up.t_i_min,
                         rect_up.x_chunk_size - 1, rect_up.t_chunk_size,
                         funcs::f, funcs::u_next);

            // Copy to u_t_buf
            u_t_buf.resize (rect_up.t_chunk_size);
            copy_col_2_row (u_t_buf.data (),
                            u_buf_up_begin + rect_up.x_chunk_size - 1,
                            u_buf_x_size, rect_up.t_chunk_size);
            // print_2d_array (u_buf, task.x_size);
            // Send right
            MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                      MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
        }

        // std::cerr << "Complete area: " << 0 << std::endl;
    } // void calc_zero_rank_func_border ()

    void
    calc_zero_rank_grid_border (std::vector <double>& u_buf) {
        std::size_t u_buf_x_size = map_mgr.x_size;
        std::vector <double> u_t_buf, u_x_buf;

        const int next_rank = 1, prev_rank = map_mgr.num_threads - 1;

        // 2) Calc u, where need gridden border conditions
        // area - is set of chuks: |__
        MPI_Status status = {};
        for (int i_area = map_mgr.num_threads; i_area < map_mgr.num_areas; i_area += map_mgr.num_threads) {
            const area_params_t area_params = map_mgr.get_area_params (i_area);
            const rect_t& zero_rect = area_params.zero_rect;

            u_t_buf.resize (zero_rect.t_chunk_size + 1);
            u_x_buf.resize (zero_rect.x_chunk_size + 1);

            std::size_t zero_pos = zero_rect.t_i_min * u_buf_x_size + zero_rect.x_i_min;
            double* u_buf_begin = u_buf.data () + zero_pos;

            // Get x and t array from prev rank
            MPI_Recv (u_x_buf.data (), zero_rect.x_chunk_size + 1,
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
            MPI_Recv (u_t_buf.data (), zero_rect.t_chunk_size,
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            // Calc zero area
            calc_u_rect_x_row_t_col (u_buf_begin, u_x_buf.data (), u_t_buf.data (),
                                     u_buf_x_size, zero_rect,
                                     funcs::f, funcs::u_next);

            u_t_buf[zero_rect.t_chunk_size] = u_t_buf[zero_rect.t_chunk_size - 1];
            for (int i_chunk_g = 1 + i_area; i_chunk_g < map_mgr.num_areas; ++i_chunk_g) {
                int i_chunk_l = i_chunk_g - i_area;

                const rect_t rect_right = area_params.get_right_rect (i_chunk_l);
                const rect_t rect_up = area_params.get_up_rect (i_chunk_l);

                double* u_buf_right_begin = u_buf.data () + rect_right.get_offset_rect (u_buf_x_size);
                double* u_buf_up_begin = u_buf.data () + rect_up.get_offset_rect (u_buf_x_size);

                // Calc u odd
                u_x_buf[0] = u_x_buf[zero_rect.x_chunk_size];
                MPI_Recv (u_x_buf.data () + 1, rect_right.x_chunk_size,
                          MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

                calc_u_rect_x_row (u_buf_right_begin, u_x_buf.data (), u_buf_x_size, rect_right,
                                   funcs::f, funcs::u_next);

                // Send up row
                if (i_chunk_l == 1) {
                    MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size - 1,
                              rect_right.x_chunk_size + 1,
                              MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
                } else {
                    MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size,
                              rect_right.x_chunk_size,
                              MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
                }

                // Calc u even
                u_t_buf[0] = u_t_buf[zero_rect.t_chunk_size];
                MPI_Recv (u_t_buf.data () + 1, rect_up.t_chunk_size,
                          MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

                calc_u_rect_t_col (u_buf_up_begin, u_t_buf.data (), u_buf_x_size, rect_up,
                                   funcs::f, funcs::u_next);

                // Send right col
                copy_col_2_row (u_t_buf.data (),
                                u_buf_up_begin + rect_up.x_chunk_size - 1,
                                u_buf_x_size, rect_up.t_chunk_size);
                MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                          MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }
        }
    } // void calc_zero_rank_grid_border ()

    void
    calc_non_zero_rank_grid_border (int rank) {
        std::size_t u_buf_x_size = map_mgr.x_size;

        const int next_rank = rank == map_mgr.num_threads - 1 ? 0 : rank + 1;
        const int prev_rank = rank - 1;

        std::vector <double> u_buf_right, u_buf_up, u_x_buf, u_t_buf;

        // area - is set of chuks: |__
        MPI_Status status = {};
        for (int i_area = rank; i_area < map_mgr.num_areas; i_area += map_mgr.num_threads) {
            const area_params_t area_params = map_mgr.get_area_params (i_area);
            const rect_t& zero_rect = area_params.zero_rect;

            // Get x and t array from prev rank
            u_t_buf.resize (zero_rect.t_chunk_size + 1);
            u_x_buf.resize (zero_rect.x_chunk_size + 1);

            u_buf_right.resize (area_params.get_zero_and_right_rect_size ());
            std::size_t u_buf_right_x_size = area_params.get_zero_and_right_rect_x_size ();
            std::size_t u_buf_up_x_size = area_params.get_up_rect_x_size ();

            MPI_Recv (u_x_buf.data (), zero_rect.x_chunk_size + 1,
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
            MPI_Recv (u_t_buf.data (), zero_rect.t_chunk_size,
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            // Calc zero area
            calc_u_rect_x_row_t_col (u_buf_right.data (), u_x_buf.data (), u_t_buf.data (),
                                     u_buf_right_x_size, zero_rect,
                                     funcs::f, funcs::u_next);

            if (area_params.get_up_rect_size () != 0) {
                u_buf_up.resize (area_params.get_up_rect_size () + u_buf_up_x_size);
                std::copy_n (u_buf_right.data () + (zero_rect.t_chunk_size - 1) * u_buf_right_x_size,
                             u_buf_up_x_size, u_buf_up.data ());
                u_t_buf[zero_rect.t_chunk_size] = u_t_buf[zero_rect.t_chunk_size - 1];
            }

            double* u_buf_up_begin_ = u_buf_up.data () + u_buf_up_x_size;
            for (int i_chunk_g = 1 + i_area; i_chunk_g < map_mgr.num_areas; ++i_chunk_g) {
                int i_chunk_l = i_chunk_g - i_area;

                const rect_t rect_right = area_params.get_right_rect (i_chunk_l);
                const rect_t rect_up = area_params.get_up_rect (i_chunk_l);

                double* u_buf_right_begin = u_buf_right.data () +
                                            rect_right.get_offset_rect (u_buf_right_x_size, zero_rect);
                double* u_buf_up_begin = u_buf_up_begin_ +
                                         (i_chunk_l - 1) * u_buf_up_x_size * zero_rect.t_chunk_size;

                // Calc u odd
                u_x_buf[0] = u_x_buf[zero_rect.x_chunk_size];
                MPI_Recv (u_x_buf.data () + 1, rect_right.x_chunk_size,
                          MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

                calc_u_rect_x_row (u_buf_right_begin, u_x_buf.data (),
                                   u_buf_right_x_size, rect_right,
                                   funcs::f, funcs::u_next);

                // Send up row
                if (i_chunk_l == 1) {
                    MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_right_x_size - 1,
                              rect_right.x_chunk_size + 1,
                              MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
                } else {
                    MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_right_x_size,
                              rect_right.x_chunk_size,
                              MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
                }

                // Calc u even
                u_t_buf[0] = u_t_buf[zero_rect.t_chunk_size];
                MPI_Recv (u_t_buf.data () + 1, rect_up.t_chunk_size,
                          MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

                calc_u_rect_t_col (u_buf_up_begin, u_t_buf.data (),
                                   u_buf_up_x_size, rect_up,
                                   funcs::f, funcs::u_next);

                // Send right col
                copy_col_2_row (u_t_buf.data (),
                                u_buf_up_begin + u_buf_up_x_size - 1,
                                u_buf_up_x_size, rect_up.t_chunk_size);

                MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                          MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }

            // std::cerr << "Complete area: " << i_area << std::endl;

            // Send u_buf_right and u_buf_up to process with rank 0
            // std::cout << "host, want send buf_right_size: " << area_params.get_zero_and_right_rect_size () << '\n';
            MPI_Send (u_buf_right.data (), area_params.get_zero_and_right_rect_size (),
                      MPI_DOUBLE, 0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);

            if (area_params.get_up_rect_size () != 0) {
                // std::cout << "host, want send buf_up_size: " << area_params.get_up_rect_size () << '\n';
                MPI_Send (u_buf_up_begin_, area_params.get_up_rect_size (),
                          MPI_DOUBLE, 0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);
            }
            // std::cerr << "Sended area: " << i_area << std::endl;
        }
    } // void calc_non_zero_rank_grid_border (int rank)

    void
    receive_u_bufs (std::vector <double>& u_buf) {
        const std::size_t u_buf_x_size = map_mgr.x_size;
        std::vector <double> u_buf_right, u_buf_up;

        MPI_Status status = {};
        for (int i_area = 1; i_area < map_mgr.num_areas; ++i_area) {
            int source = i_area % map_mgr.num_threads;
            if (source == 0) {
                continue;
            }
            // DUMP (source);

            const area_params_t area_params = map_mgr.get_area_params (i_area);
            const rect_t& zero_rect = area_params.zero_rect;

            const std::size_t u_buf_right_x_size = area_params.get_zero_and_right_rect_x_size ();
            const std::size_t u_buf_right_size = area_params.get_zero_and_right_rect_size ();
            const std::size_t u_buf_up_x_size = area_params.get_up_rect_x_size ();
            const std::size_t u_buf_up_size = area_params.get_up_rect_size ();

            u_buf_right.resize (u_buf_right_size);  // todo delete
            u_buf_up.resize (u_buf_up_size);        // todo delete

            // DUMP (u_buf_right_size);
            // DUMP (u_buf_right.data ());
            // sleep (1);
            // Receive right buffer
            // std::cout << "host, want get buf_right_size: " << u_buf_right_size << '\n';
            // MPI_Probe (source, TAG_SAVE_ON_HOST, MPI_COMM_WORLD, &status);
            MPI_Recv (u_buf_right.data (), u_buf_right_size,
                      MPI_DOUBLE, source, TAG_SAVE_ON_HOST, MPI_COMM_WORLD, &status);

            // int count = -1;
            // MPI_Get_count (&status, MPI_DOUBLE, &count);
            // DUMP (count);

            double* u_buf_right_begin = u_buf.data () + zero_rect.get_offset_rect (u_buf_x_size);
            copy_row_2_rect (u_buf_right_begin, u_buf_right.data (),
                             u_buf_right_x_size, u_buf_x_size, u_buf_right_size);

            if (u_buf_up_size != 0) {
                // Receive up buffer
                // std::cout << "host, want get tmp_buf_up.size (): " << u_buf_up_size << std::endl;
                MPI_Recv (u_buf_up.data (), u_buf_up_size,
                          MPI_DOUBLE, source, TAG_SAVE_ON_HOST, MPI_COMM_WORLD, &status);

                double* u_buf_up_begin = u_buf_right_begin + zero_rect.t_chunk_size * u_buf_x_size;
                copy_row_2_rect (u_buf_up_begin, u_buf_up.data (),
                                 u_buf_up_x_size, u_buf_x_size, u_buf_up_size);
            }
            // std::cerr << "Getted area: " << i_area << std::endl;
        }
    } // void receive_u_bufs ()

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row_t_col (T* u_buf,
                             const T* u_x_buf,
                             const T* u_t_buf,
                             std::size_t u_buf_x_size,
                             std::size_t x_i_min,
                             std::size_t t_i_min,
                             std::size_t x_size_calc,
                             std::size_t t_size_calc,
                             F f,
                             U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc left down point
        u_buf[0] = u_next (dx, dt, u_x_buf[0], u_x_buf[1], u_t_buf[0],
                           f (x_min - dx/2, t_min - dt/2));

        // Calc first row of u
        double u_k1_m = u_buf[0];
        for (long x_i_rel = 1; x_i_rel < x_size_calc; ++x_i_rel) {
            double u_k_m  = u_x_buf[x_i_rel];
            double u_k_m1 = u_x_buf[x_i_rel + 1];

            double x = x_min + x_i_rel * dx;
            double f_kh_mh = f (t_min - dt/2, x - dx/2);

            u_k1_m = u_buf[x_i_rel]
                   = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc first col of u
        double u_k_m1 = u_buf[0];
        for (long t_i_rel = 1; t_i_rel < t_size_calc; ++t_i_rel) {
            double u_k_m  = u_t_buf[t_i_rel - 1];
            double u_k1_m = u_t_buf[t_i_rel];

            double t = t_min + t_i_rel * dt;
            double f_kh_mh = f (t - dt/2, x_min - dx/2);

            u_k_m1 = u_buf[t_i_rel * u_buf_x_size]
                   = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc rect (x_size_calc - 1) x (t_size_calc - 1)
        calc_u_rect (u_buf + u_buf_x_size + 1, u_buf_x_size,
                     x_i_min + 1, t_i_min + 1, x_size_calc -1, t_size_calc -1,
                     funcs::f, funcs::u_next);
    } // void calc_u_rect_x_row_t_col (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row_t_col (T* u_buf,
                             const T* u_x_buf,
                             const T* u_t_buf,
                             std::size_t u_buf_x_size,
                             const rect_t& rect_calc,
                             F f,
                             U_next u_next)
    {
        calc_u_rect_x_row_t_col (u_buf, u_x_buf, u_t_buf, u_buf_x_size,
                                 rect_calc.x_i_min, rect_calc.t_i_min,
                                 rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                                 f, u_next);
    } // void calc_u_rect_x_row_t_col (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row (T* u_buf,
                       const T* u_x_buf,
                       std::size_t u_buf_x_size,
                       std::size_t x_i_min,
                       std::size_t t_i_min,
                       std::size_t x_size_calc,
                       std::size_t t_size_calc,
                       F f,
                       U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc first row of u
        for (long x_i_rel = 0; x_i_rel < x_size_calc; ++x_i_rel) {
            double u_k_m  = u_x_buf[x_i_rel];
            double u_k_m1 = u_x_buf[x_i_rel + 1];
            double u_k1_m = u_buf[x_i_rel - 1];

            double x = x_min + x_i_rel * dx;
            double f_kh_mh = f (t_min - dt/2, x - dx/2);

            u_buf[x_i_rel] = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc (rect x_size_calc) x (t_size_calc - 1) // todo
        calc_u_rect (u_buf + u_buf_x_size, u_buf_x_size, x_i_min, t_i_min + 1,
                     x_size_calc, t_size_calc - 1,
                     funcs::f, funcs::u_next);
    } // void calc_u_rect_x_row (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row (T* u_buf,
                       const T* u_x_buf,
                       std::size_t u_buf_x_size,
                       const rect_t& rect_calc,
                       F f,
                       U_next u_next)
    {
        calc_u_rect_x_row (u_buf, u_x_buf, u_buf_x_size,
                           rect_calc.x_i_min, rect_calc.t_i_min,
                           rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                           f, u_next);
    } // void calc_u_rect_x_row (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_t_col (T* u_buf,
                       const T* u_t_buf,
                       std::size_t u_buf_x_size,
                       std::size_t x_i_min,
                       std::size_t t_i_min,
                       std::size_t x_size_calc,
                       std::size_t t_size_calc,
                       F f,
                       U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc first col of u
        for (long t_i_rel = 0; t_i_rel < t_size_calc; ++t_i_rel) {
            double* u_buf_cur = u_buf + t_i_rel * u_buf_x_size;
            double* u_buf_prev = u_buf_cur - u_buf_x_size;

            double u_k_m  = u_t_buf[t_i_rel];
            double u_k1_m = u_t_buf[t_i_rel + 1];
            double u_k_m1 = *u_buf_prev;

            double t = t_min + t_i_rel * dt;
            double f_kh_mh = f (t - dt/2, x_min - dx/2);

            *u_buf_cur = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc rect of u
        calc_u_rect (u_buf + 1, u_buf_x_size, x_i_min + 1, t_i_min,
                     x_size_calc - 1, t_size_calc,
                     funcs::f, funcs::u_next);
    } // void calc_u_rect_t_col (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_t_col (T* u_buf,
                       const T* u_t_buf,
                       std::size_t u_buf_x_size,
                       const rect_t& rect_calc,
                       F f,
                       U_next u_next)
    {
        calc_u_rect_t_col (u_buf, u_t_buf, u_buf_x_size,
                           rect_calc.x_i_min, rect_calc.t_i_min,
                           rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                           f, u_next);
    } // void calc_u_rect_t_col (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect (T* u_buf,
                 std::size_t u_buf_x_size,
                 std::size_t x_i_min,
                 std::size_t t_i_min,
                 std::size_t x_size_calc,
                 std::size_t t_size_calc,
                 F f,
                 U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc rect (x_size_calc - 1) x (t_size_calc - 1)
        for (long t_i_rel = 0; t_i_rel < t_size_calc; ++t_i_rel) {
            double* u_row_cur = u_buf + t_i_rel * u_buf_x_size;
            double* u_row_prev = u_row_cur - u_buf_x_size;

            double t = t_min + t_i_rel * dt;
            for (long x_i_rel = 0; x_i_rel < x_size_calc; ++x_i_rel) {
                double u_k_m  = u_row_prev[x_i_rel - 1];
                double u_k_m1 = u_row_prev[x_i_rel];
                double u_k1_m = u_row_cur[x_i_rel - 1];

                double x = x_min + x_i_rel * dx;
                double f_kh_mh = f (t - dt/2, x - dx/2);

                u_row_cur[x_i_rel] = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
            }
        }
    } // void calc_u_rect (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect (T* u_buf,
                 std::size_t u_buf_x_size,
                 const rect_t& rect_calc,
                 F f,
                 U_next u_next)
    {
        calc_u_rect (u_buf, u_buf_x_size,
                     rect_calc.x_i_min, rect_calc.t_i_min,
                     rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                     f, u_next);
    } // void calc_u_rect (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

}; // struct trans_eq_solver

} // namespace treq

void
solve_trans_eq_parallel (const trans_eq_task_t& task,
                         int* argc_ptr,
                         char** argv_ptr[])
{
    MPI_Init (argc_ptr, argv_ptr);

    int num_threads = 0, rank = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_threads);

    int k_zone = 2;

    treq::trans_eq_solver solver {task, num_threads, k_zone};
    solver.solve (rank);

    MPI_Finalize ();
} // void solve_trans_eq_parallel

// Какая-то проблема с синхронизацией