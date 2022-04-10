#include "lib.hpp"

namespace treq {

struct rect_t {
    int x_chunk_size, t_chunk_size;
    int x_i_min, t_i_min;
}; // struct rect_t

struct area_params_t {

    rect_t zero_rect;
    int x_i_end, t_i_end;

    int
    get_x_zero_chunk_size () const noexcept {
        return std::min (x_i_end - zero_rect.x_i_min, zero_rect.x_chunk_size);
    }

    int
    get_t_zero_chunk_size () const noexcept {
        return std::min (t_i_end - zero_rect.t_i_min, zero_rect.t_chunk_size);
    }

    int
    get_zero_right_rect_size () const noexcept {
        return get_t_zero_chunk_size () * (x_i_end - zero_rect.x_i_min);
    }

    int
    get_up_rect_size () const noexcept {
        int area_t_size = t_i_end - zero_rect.t_i_min;
        if (area_t_size <= zero_rect.t_chunk_size) {
            return 0;
        } else {
            return (area_t_size - zero_rect.t_chunk_size) * zero_rect.x_chunk_size;
        }
    }

    area_params_t
    get_right_rect (int i_chunk) const noexcept {
        rect_t right_rect = {};
        right_rect.x_i_min = zero_rect.x_i_min + i_chunk * zero_rect.x_chunk_size;
        right_rect.t_i_min = zero_rect.t_i_min + i_chunk * zero_rect.t_chunk_size;

        if (right_rect.x_i_min + zero_rect.x_chunk_size > x_i_end) {
            right_rect.x_chunk_size = x_i_end - right_rect.x_i_min;
            right_rect.t_chunk_size = t_i_end - right_rect.t_i_min;
        } else {
            right_rect.x_chunk_size = zero_rect.x_chunk_size;
            right_rect.t_chunk_size = zero_rect.t_chunk_size;
        }
    }
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

        params.zero_rect.x_chunk_size = calc_chunk_size (x_size, num_areas);
        params.zero_rect.t_chunk_size = calc_chunk_size (t_size, num_areas);

        params.zero_rect.x_i_min = i_area * params.zero_rect.x_chunk_size;
        params.zero_rect.t_i_min = i_area * params.zero_rect.t_chunk_size;

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

    std::vector <double> u_buf;

    void
    solve (int rank) {
        if (rank == 0) {
            calc_zero_rank_func_border ();
            calc_zero_rank_grid_border ();
        } else {
            calc_non_zero_rank_grid_border (rank);
        }
    }

    void
    calc_zero_rank_func_border () {
        // fill U (x_min, [t_min, t_min + dt * t_chuck_size])
        // fill U (t_min, [x_min, x_min + dx * x_chuck_size])
        // calc U ([x_min, x_min + dx * x_chuck_size],
        //         [t_min, t_min + dt * t_chuck_size])

        /*
        double x_right_min, t_right_min;
        double x_up_min, t_up_min;

        for (int i_area = 1; i_area < map_mgr.num_area; ++i_area) {
            fill U (t, [x])
            calc U_right
            send up of U_right

            fill U (x, [t])
            calc U_up
            send right of U_up
        }
        */
    }

    void
    calc_zero_rank_grid_border () {

    }

    void
    calc_non_zero_rank_grid_border (int rank) {

    }

    const map_manager_t map_mgr;

    trans_eq_solver (const trans_eq_task_t& task,
                     int num_threads,
                     int koef_mult_threads) :
        map_mgr (task, num_threads, koef_mult_threads * num_threads)
    {}

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect (T* u_buf,
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
        // todo
    }

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect (T* u_buf,
                 const T* u_x_buf,
                 const T* u_t_buf,
                 std::size_t u_buf_x_size,
                 const rect_t& rect_calc,
                 F f,
                 U_next u_next)
    {
        calc_u_rect (u_buf, u_x_buf, u_t_buf, u_buf_x_size,
                     rect_calc.x_i_min, rect_calc.t_i_min,
                     rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                     f, u_next);
    }

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
        // todo
    }

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
    }

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
        // todo
    }

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
    }

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
        // todo
    }

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
    }

}; // struct trans_eq_solver

} // namespace treq