#pragma once

#include "copy_transform.hpp"
#include "func_lib.hpp"
#include "print_lib.hpp"

namespace treq_seq {

#define CHECK_ERR(err_code) check_error (err_code, __LINE__, __FILE__)

void
check_error (int err_code, int line, const char* file_name);

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
                    double t_max);

std::vector <double>
get_solve_ref (int num_points_coord_x,
               int num_points_time,
               double x_max,
               double t_max);

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
                double t_max);

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

void
check_single_thread_solution ();

/*
b solve_trans_eq_parallel_zero_rank
b work_non_zero_rank_grid_border
r

*/
} // namespace treq_seq

namespace treq {

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

    static std::size_t
    calc_chunk_size (std::size_t full_size,
                    std::size_t num_chunks) noexcept
    {
        return  full_size / num_chunks +
               (full_size % num_chunks != 0);
    }

    area_params_t
    get_area_params (int i_area) const noexcept;

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
    const int TAG_BORDER_COND = 1;
    const int TAG_SAVE_ON_HOST = 2;

    const map_manager_t map_mgr;

    trans_eq_solver (const trans_eq_task_t& task,
                     int num_threads,
                     int koef_mult_threads) :
        map_mgr (task, num_threads, koef_mult_threads * num_threads)
    {}

    void solve (int rank);

    void calc_zero_rank_func_border (std::vector <double>& u_buf);
    void calc_zero_rank_grid_border (std::vector <double>& u_buf);
    void calc_non_zero_rank_grid_border (int rank);
    void receive_u_bufs (std::vector <double>& u_buf);

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
solve_trans_eq_parallel (const treq::trans_eq_task_t& task,
                         int* argc_ptr,
                         char** argv_ptr[]);