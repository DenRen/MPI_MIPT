#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <thread>
#include <mpi/mpi.h>

#include "copy_transform.hpp"
#include "func_lib.hpp"
#include "print_lib.hpp"

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
calc_u_part_buf (int x_chunk_size,
                 int t_chunk_size,
                 int buf_x_size,
                 int x_arr_begin, int t_arr_begin,
                 int x_i_min, int t_i_min,
                 double dx, double dt,
                 std::vector <double>& u_buf,
                 F f,
                 U_NEXT u_next)
{
    std::size_t buf_offset = (x_i_min - x_arr_begin) + (t_i_min - t_arr_begin) * buf_x_size;

    // Calc other u
    for (int time_i = 1; time_i < t_chunk_size; ++time_i) {
        double t = (t_i_min + time_i) * dt;
        std::size_t cur_shift = time_i * buf_x_size + buf_offset;
        std::size_t prev_shift = cur_shift - buf_x_size;

        double u_k1_m = u_buf[cur_shift];
        for (int x_i = 0; x_i < x_chunk_size - 1; ++x_i) {
            std::size_t prev_pos = prev_shift + x_i;
            std::size_t cur_pos = cur_shift + x_i;

            double u_k_m  = u_buf[prev_pos];
            double u_k_m1 = u_buf[prev_pos + 1];

            double x = (x_i_min + x_i) * dx;
            double f_kh_mh = f (t - dt/2, x + dx/2);

            u_k1_m = u_buf[cur_pos + 1]
                   = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }
    }
}

template <typename F, typename U_NEXT>
void
calc_u_full_buf (int x_chunk_size,
                 int t_chunk_size,
                 int buf_x_size,
                 int x_min, int t_min,
                 double dx, double dt,
                 std::vector <double>& u_buf,
                 F f,
                 U_NEXT u_next)
{
    calc_u_part_buf (x_chunk_size, t_chunk_size, buf_x_size,
                     0, 0, x_min, t_min, dx, dt, u_buf, f, u_next);
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

#define TAG_BORDER_COND 15
#define TAG_SAVE_ON_HOST 16

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

void
work_zero_rank_func_border (const trans_eq_task_t& task,
                            int num_chunk_area,
                            int num_threads,
                            std::vector <double>& u_buf,
                            std::vector <double>& u_t_buf)
{
    const int x_chunk_size = calc_chunk_size (task.x_size - 1, num_chunk_area);
    const int t_chunk_size = calc_chunk_size (task.t_size - 1, num_chunk_area);

    u_t_buf.resize (t_chunk_size + 1);

    using funcs::u_0_x;
    using funcs::u_t_0;
    using funcs::f;
    using funcs::u_next;

    auto[dx, dt] = task.calc_dx_dt ();

    // 1) Calc u, where need only functional border conditions
    // Fill u(t = 0, x)
    for (int x_i = 0; x_i <= x_chunk_size; ++x_i) {
        double x = x_i * dx;
        u_buf[x_i] = u_0_x (x);
    }

    // Fill u(t, x = 0)
    for (int t_i = 1; t_i <= t_chunk_size; ++t_i) {
        double t = t_i * dt;
        u_buf[t_i * task.x_size] = u_t_0 (t);
    }

    // Calc u on (x_size - 1) * (t_size - 1)
    calc_u_full_buf (x_chunk_size + 1, t_chunk_size + 1,
                     task.x_size,
                     0, 0, dx, dt,
                     u_buf, f, u_next);

    // Calc u right and up
    const int next_rank = 1, prev_rank = num_threads - 1;
    for (int i_chunk = 1; i_chunk < num_chunk_area; ++i_chunk) {
        // Fill u(t = 0, x)
        int x_i_min = i_chunk * x_chunk_size;   // x_min - is node
        int x_i_end = calc_end_index (x_i_min, x_chunk_size + 1, task.x_size);
        for (int x_i = x_i_min + 1; x_i < x_i_end; ++x_i) {
            double x = x_i * dx;
            u_buf[x_i] = u_0_x (x);
        }

        // Calc u odd
        int x_chunk_size_corrected = x_i_end - x_i_min - 1;
        calc_u_full_buf (x_chunk_size_corrected + 1, t_chunk_size + 1,
                         task.x_size,
                         x_i_min, 0, dx, dt,
                         u_buf, f, u_next);

        // Send up
        double* buf_odd_begin = u_buf.data () + t_chunk_size * task.x_size + x_i_min;
        MPI_Send (buf_odd_begin, x_chunk_size_corrected + 1,
                  MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);

        // Fill u(t = 0, x)
        int t_i_min = i_chunk * t_chunk_size;
        int t_i_end = calc_end_index (t_i_min, t_chunk_size + 1, task.t_size);
        for (int t_i = t_i_min + 1; t_i < t_i_end; ++t_i) {
            double t = t_i * dt;
            u_buf[t_i * task.x_size] = u_t_0 (t);
        }

        // Calc u even
        int t_chunk_size_corrected = t_i_end - t_i_min - 1;
        calc_u_full_buf (x_chunk_size + 1, t_chunk_size_corrected + 1,
                         task.x_size,
                         0, t_i_min, dx, dt,
                         u_buf, f, u_next);

        // Copy to u_t_buf
        copy_col_2_row (u_t_buf.data (),
                        u_buf.data () + (t_i_min + 1) * task.x_size + x_chunk_size,
                        task.x_size, t_chunk_size_corrected);
        // Send right
        MPI_Send (u_t_buf.data (), t_chunk_size_corrected,
                  MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
    }
    // Good. Todo: overbounds
}

void
work_zero_rank_grid_border (const trans_eq_task_t& task,
                            int num_chunk_area,
                            int num_threads,
                            std::vector <double>& u_buf,
                            std::vector <double>& u_t_buf)
{
    int x_chunk_size = calc_chunk_size (task.x_size - 1, num_chunk_area);
    int t_chunk_size = calc_chunk_size (task.t_size - 1, num_chunk_area);

    u_t_buf.resize (t_chunk_size);

    const int next_rank = 1, prev_rank = num_threads - 1;
    auto[dx, dt] = task.calc_dx_dt ();

    using funcs::f;
    using funcs::u_next;

    // 2) Calc u, where need gridden border conditions
    // area - is set of chuks: |__
    MPI_Status status = {};
    for (int i_area = num_threads; i_area < num_chunk_area; i_area += num_threads) {
        int t_i_min = i_area * t_chunk_size;
        int x_i_min = i_area * x_chunk_size;
        std::size_t zero_pos = t_i_min * task.x_size + x_i_min;
        double* buf_begin = u_buf.data () + zero_pos;

        // Get x and t array from prev rank
        int x_i_end = calc_end_index (x_i_min, x_chunk_size + 1, task.x_size);
        int t_i_end = calc_end_index (t_i_min, t_chunk_size + 1, task.t_size);
        int x_chunk_size_corrected = x_i_end - x_i_min - 1;
        int t_chunk_size_corrected = t_i_end - t_i_min - 1;

        MPI_Recv (buf_begin, x_chunk_size_corrected + 1, MPI_DOUBLE,
                  prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
        MPI_Recv (u_t_buf.data (), t_chunk_size_corrected, MPI_DOUBLE,
                  prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

        copy_row_2_col (buf_begin + task.x_size, u_t_buf.data (),
                        task.x_size, t_chunk_size_corrected);

        // Calc zero area
        calc_u_full_buf (x_chunk_size_corrected, t_chunk_size_corrected,
                         task.x_size,
                         x_i_min, t_i_min, dx, dt,
                         u_buf, f, u_next);

        for (int i_chunk_g = 1 + i_area; i_chunk_g < num_chunk_area; ++i_chunk_g) {
            int i_chunk_l = i_chunk_g - i_area;
            double* buf_begin_odd  = buf_begin + i_chunk_l * x_chunk_size;
            double* buf_begin_even = buf_begin + i_chunk_l * task.x_size;

            int x_i_odd_min = x_i_min + i_chunk_l * x_chunk_size;
            int x_i_odd_end = calc_end_index (x_i_odd_min, x_chunk_size + 1, task.x_size);
            int x_chunk_size_corrected = x_i_odd_end - x_i_odd_min - 1;

            // Calc u odd
            MPI_Recv (buf_begin_odd, x_chunk_size_corrected + 1, MPI_DOUBLE,
                      prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            calc_u_full_buf (x_chunk_size_corrected + 1, t_chunk_size + 1,
                             task.x_size,
                             x_i_odd_min, t_i_min, dx, dt,
                             u_buf, f, u_next);

            // Send up
            MPI_Send (buf_begin_odd + t_chunk_size * task.x_size, x_chunk_size_corrected + 1,
                      MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);

            int t_i_odd_min = t_i_min + i_chunk_l * t_chunk_size;
            int t_i_odd_end = calc_end_index (t_i_odd_min, t_chunk_size + 1, task.t_size);
            int t_chunk_size_corrected = t_i_odd_end - t_i_odd_min - 1;

            // Calc u even
            MPI_Recv (u_t_buf.data (), t_chunk_size_corrected, MPI_DOUBLE,
                      prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
            copy_row_2_col (buf_begin_even + task.x_size, u_t_buf.data (),
                            task.x_size, t_chunk_size_corrected);

            calc_u_full_buf (x_chunk_size + 1, t_chunk_size_corrected + 1,
                             task.x_size,
                             x_i_min, t_i_min + i_chunk_l * t_chunk_size, dx, dt,
                             u_buf, f, u_next);

            // Send right
            copy_col_2_row (u_t_buf.data (),
                            buf_begin_even + task.x_size + x_chunk_size,
                            task.x_size, t_chunk_size_corrected);
            MPI_Send (u_t_buf.data (), t_chunk_size_corrected,
                      MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
        }
    }
}

// For zero rank
void
receive_u_bufs (const trans_eq_task_t& task,
                int num_chunk_area,
                int num_threads,
                std::vector <double>& u_buf)
{
    const int x_chunk_size = calc_chunk_size (task.x_size - 1, num_chunk_area);
    const int t_chunk_size = calc_chunk_size (task.t_size - 1, num_chunk_area);

    std::vector <double> tmp_buf_right, tmp_buf_up;

    MPI_Status status = {};
    for (int i_area = 1; i_area < num_chunk_area; ++i_area) {
        int source = i_area % num_threads;
        if (source == 0) {
            continue;
        }

        // Receive right buffer
        int buf_right_x_size = task.x_size - i_area * x_chunk_size;
        int buf_up_t_size    = task.t_size - (i_area + 1) * t_chunk_size;
        int x_i_min = i_area * x_chunk_size;
        int t_i_min = i_area * t_chunk_size;

        double* buf_right_begin = u_buf.data () + x_i_min + t_i_min * task.x_size;
        std::size_t buf_right_size = buf_right_x_size * (t_chunk_size + 1);
        tmp_buf_right.resize (buf_right_size);
        MPI_Recv (tmp_buf_right.data (), buf_right_size, MPI_DOUBLE,
                  source, TAG_SAVE_ON_HOST, MPI_COMM_WORLD, &status);
        copy_row_2_rect (buf_right_begin, tmp_buf_right.data (),
                         buf_right_x_size, task.x_size, buf_right_size);

        // Receive up buffer
        std::size_t buf_up_size = buf_up_t_size * (t_chunk_size + 1);
        tmp_buf_up.resize (buf_up_size);
        MPI_Recv (tmp_buf_up.data (), buf_up_size, MPI_DOUBLE,
                  source, TAG_SAVE_ON_HOST, MPI_COMM_WORLD, &status);

        double* buf_up_begin = buf_right_begin + t_chunk_size * task.x_size;
        copy_row_2_rect (buf_up_begin, tmp_buf_up.data (),
                         t_chunk_size + 1, task.x_size, tmp_buf_up.size ());
    }
}

void
solve_trans_eq_parallel_zero_rank (const trans_eq_task_t& task,
                                   int num_chunk_area,
                                   int num_threads,
                                   std::vector <double>& u_buf)
{
    std::vector <double> u_t_buf;

    auto compute_func = [&] {
        work_zero_rank_func_border (task, num_chunk_area, num_threads, u_buf, u_t_buf);
        work_zero_rank_grid_border (task, num_chunk_area, num_threads, u_buf, u_t_buf);
    };

    auto bufferize_func = [&] {
        receive_u_bufs (task, num_chunk_area,num_threads, u_buf);
    };

    std::thread compute_thread {compute_func};
    std::thread bufferer_thread {bufferize_func};

    compute_thread.join ();
    bufferer_thread.join ();
}

void
work_non_zero_rank_grid_border (const trans_eq_task_t& task,
                                int rank,
                                int num_chunk_area,
                                int num_threads)
{
    int x_chunk_size = calc_chunk_size (task.x_size - 1, num_chunk_area);
    int t_chunk_size = calc_chunk_size (task.t_size - 1, num_chunk_area);

    int buf_right_x_size = task.x_size - rank * x_chunk_size;
    int buf_up_t_size    = task.t_size - (rank + 1) * t_chunk_size;
    std::vector <double> u_buf_right (buf_right_x_size * (t_chunk_size + 1)),
                         u_buf_up (buf_up_t_size * (x_chunk_size + 1)),
                         u_t_buf (t_chunk_size + 1);

    const int next_rank = rank == num_threads - 1 ? 0 : rank - 1;
    const int prev_rank = rank - 1;
    auto[dx, dt] = task.calc_dx_dt ();

    using funcs::f;
    using funcs::u_next;

    // area - is set of chuks: |__
    MPI_Status status = {};
    for (int i_area = rank; i_area < num_chunk_area; i_area += num_threads) {
        int t_i_min = i_area * t_chunk_size;
        int x_i_min = i_area * x_chunk_size;

        // Get x and t array from prev rank
        int x_i_end = calc_end_index (x_i_min, x_chunk_size + 1, task.x_size);
        int t_i_end = calc_end_index (t_i_min, t_chunk_size + 1, task.t_size);
        int x_chunk_size_corrected = x_i_end - x_i_min - 1;
        int t_chunk_size_corrected = t_i_end - t_i_min - 1;

        MPI_Recv (u_buf_right.data (), x_chunk_size_corrected + 1, MPI_DOUBLE,
                  prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
        MPI_Recv (u_t_buf.data (), t_chunk_size_corrected, MPI_DOUBLE,
                  prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

        copy_row_2_col (u_buf_right.data () + buf_right_x_size, u_t_buf.data (),
                        buf_right_x_size, t_chunk_size_corrected);

        // Calc zero area
        calc_u_part_buf (x_chunk_size_corrected + 1, t_chunk_size_corrected + 1,
                         buf_right_x_size,
                         x_i_min, t_i_min,
                         x_i_min, t_i_min, dx, dt,
                         u_buf_right, f, u_next);
        std::copy_n (u_buf_right.cbegin () + t_chunk_size * buf_right_x_size,
                     x_chunk_size_corrected + 1, u_buf_up.begin ());

        buf_right_x_size = task.x_size;
        buf_up_t_size = task.t_size;

        for (int i_chunk_g = 1 + i_area; i_chunk_g < num_chunk_area; ++i_chunk_g) {
            int i_chunk_l = i_chunk_g - i_area;
            double* buf_begin_odd  = u_buf_right.data () + i_chunk_l * x_chunk_size;
            double* buf_begin_even = u_buf_up.data ()
                                     + (i_chunk_l - 1) * (t_chunk_size + 1) * (x_chunk_size + 1);

            int x_i_odd_min = x_i_min + i_chunk_l * x_chunk_size;
            int x_i_odd_end = calc_end_index (x_i_odd_min, x_chunk_size + 1, task.x_size);
            int x_chunk_size_corrected = x_i_odd_end - x_i_odd_min - 1;

            // Resize buffers
            buf_right_x_size -= num_threads * x_chunk_size;
            buf_up_t_size -= num_threads * t_chunk_size;

            u_buf_right.resize (buf_right_x_size * (t_chunk_size + 1));
            u_buf_up.resize (buf_up_t_size * (x_chunk_size + 1));

            // Calc u odd
            MPI_Recv (buf_begin_odd, x_chunk_size_corrected + 1, MPI_DOUBLE,
                      prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            calc_u_part_buf (x_chunk_size_corrected + 1, t_chunk_size + 1,
                             buf_right_x_size,
                             x_i_min, t_i_min,
                             x_i_odd_min, t_i_min, dx, dt,
                             u_buf_right, f, u_next);

            // Send up
            MPI_Send (buf_begin_odd + t_chunk_size * buf_right_x_size, x_chunk_size_corrected + 1,
                      MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);

            int t_i_odd_min = t_i_min + i_chunk_l * t_chunk_size;
            int t_i_odd_end = calc_end_index (t_i_odd_min, t_chunk_size, task.t_size - 1);
            int t_chunk_size_corrected = t_i_odd_end - t_i_odd_min;

            // Calc u even
            MPI_Recv (u_t_buf.data (), t_chunk_size_corrected, MPI_DOUBLE,
                      prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
            copy_row_2_col (buf_begin_even + x_chunk_size + 1, u_t_buf.data (),
                            x_chunk_size + 1, t_chunk_size_corrected);

            calc_u_part_buf (x_chunk_size + 1, t_chunk_size_corrected + 1,
                             x_chunk_size,
                             x_i_min, t_i_min,
                             x_i_min, t_i_min + i_chunk_l * t_chunk_size, dx, dt,
                             u_buf_up, f, u_next);

            // Send right
            copy_col_2_row (u_t_buf.data (),
                            buf_begin_even + (x_chunk_size + 1) + x_chunk_size,
                            task.x_size, t_chunk_size_corrected);
            MPI_Send (u_t_buf.data (), t_chunk_size_corrected,
                      MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
        }

        // print_2d_array (u_buf_right, x_chunk_size_corrected + 1);
        // Send u_buf_right and u_buf_up to process with rank 0
        MPI_Send (u_buf_right.data (), u_buf_right.size (), MPI_DOUBLE,
                  0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);
        MPI_Send (u_buf_up.data (), u_buf_up.size (), MPI_DOUBLE,
                  0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);
    }
}

void
solve_trans_eq_parallel_non_zero_rank (const trans_eq_task_t& task,
                                       int num_chunk_area,
                                       int num_threads,
                                       int rank)
{
    assert (rank > 0);
    work_non_zero_rank_grid_border (task, rank, num_chunk_area, num_threads);
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

    int k_zone = 4/4;
    int num_chunk_area = num_threads * k_zone;

    if (rank == 0) {
        std::vector <double> u (task.t_size * task.x_size);
        
        solve_trans_eq_parallel_zero_rank (task, num_chunk_area, num_threads, u);
        print_2d_array (u.data (), task.x_size, u.size ());
    } else {
        work_non_zero_rank_grid_border (task, rank, num_chunk_area, num_threads);
    }

    MPI_Finalize ();
}

/*
b solve_trans_eq_parallel_zero_rank
b work_non_zero_rank_grid_border
r

*/