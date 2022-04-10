#include "lib.hpp"

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
    get_zero_and_right_rect_size () const noexcept {
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

    rect_t
    get_right_rect (int i_chunk) const noexcept {
        rect_t right_rect = {
            .x_i_min = zero_rect.x_i_min + i_chunk * zero_rect.x_chunk_size,
            .t_i_min = zero_rect.t_i_min
        };

        if (right_rect.x_i_min + zero_rect.x_chunk_size > x_i_end) {
            right_rect.x_chunk_size = x_i_end - right_rect.x_i_min;
            right_rect.t_chunk_size = t_i_end - right_rect.t_i_min;
        } else {
            right_rect.x_chunk_size = zero_rect.x_chunk_size;
            right_rect.t_chunk_size = zero_rect.t_chunk_size;
        }

        return right_rect;
    } // rect_t get_right_rect (int i_chunk)

    rect_t
    get_up_rect (int i_chunk) const noexcept {
        rect_t up_rect = {
            .x_i_min = zero_rect.x_i_min,
            .t_i_min = zero_rect.t_i_min + i_chunk * zero_rect.t_chunk_size,
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
    const map_manager_t map_mgr;
    std::vector <double> u_buf;

    trans_eq_solver (const trans_eq_task_t& task,
                     int num_threads,
                     int koef_mult_threads) :
        map_mgr (task, num_threads, koef_mult_threads * num_threads)
    {}

    void
    solve (int rank) {
        if (rank == 0) {
            auto compute_func = [&] {
                calc_zero_rank_func_border ();
                calc_zero_rank_grid_border ();
            };

            auto bufferize_func = [&] {
                receive_u_bufs ();
            };

            std::thread compute_thread {compute_func};
            std::thread bufferer_thread {bufferize_func};

            compute_thread.join ();
            bufferer_thread.join ();
        } else {
            calc_non_zero_rank_grid_border (rank);
        }
    } // void solve

    void
    calc_zero_rank_func_border () {
        const area_params_t area_params = map_mgr.get_area_params (0);
        const auto[dx, dt] = map_mgr.calc_dx_dt ();
        std::size_t u_buf_x_size = map_mgr.x_size;
        std::vector <double> u_t_buf;

        // 1) Calc u, where need only functional border conditions
        // Fill u(t = 0, x)
        for (int x_i = 0; x_i < area_params.get_x_zero_chunk_size (); ++x_i) {
            double x = x_i * dx;
            u_buf[x_i] = funcs::u_0_x (x);
        }

        // Fill u(t, x = 0)
        for (int t_i = 1; t_i < area_params.get_t_zero_chunk_size (); ++t_i) {
            double t = t_i * dt;
            u_buf[t_i * u_buf_x_size] = funcs::u_t_0 (t);
        }

        // Calc u on (x_size - 1) * (t_size - 1)
        calc_u_rect (u_buf.data (), u_buf_x_size,
                     area_params.zero_rect,
                     funcs::f, funcs::u_next);

        // Calc u rights and ups in zero area
        const int next_rank = 1, prev_rank = map_mgr.num_threads - 1;
        for (int i_chunk = 1; i_chunk < map_mgr.num_areas; ++i_chunk) {
            // Fill u(t = 0, x)
            const rect_t rect_right = area_params.get_right_rect (i_chunk);
            for (int x_i = rect_right.x_i_min; x_i < rect_right.get_x_i_end (); ++x_i) {
                double x = x_i * dx;
                u_buf[x_i] = funcs::u_0_x (x);
            }

            // Calc right rect u
            double* u_buf_right_begin = u_buf.data () + rect_right.x_i_min;
            calc_u_rect (u_buf_right_begin, u_buf_x_size,
                         rect_right,
                         funcs::f, funcs::u_next);

            // Send up
            MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size - 1,
                      rect_right.x_chunk_size + 1,
                      MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);

            // Fill u(t = 0, x)
            const rect_t rect_up = area_params.get_up_rect (i_chunk);
            for (int t_i = rect_up.t_i_min; t_i < rect_up.get_t_i_end (); ++t_i) {
                double t = t_i * dt;
                u_buf[t_i * u_buf_x_size] = funcs::u_t_0 (t);
            }

            // Calc up rect u
            double* u_buf_up_begin = u_buf.data () + rect_up.t_i_min * u_buf_x_size;
            calc_u_rect (u_buf_up_begin, u_buf_x_size,
                         rect_up,
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
    } // void calc_zero_rank_func_border ()

    void
    calc_zero_rank_grid_border () {
        const auto[dx, dt] = map_mgr.calc_dx_dt ();
        std::size_t u_buf_x_size = map_mgr.x_size;
        std::vector <double> u_t_buf, u_x_buf;

        const int next_rank = 1, prev_rank = map_mgr.num_threads - 1;

        // 2) Calc u, where need gridden border conditions
        // area - is set of chuks: |__
        MPI_Status status = {};
        for (int i_area = map_mgr.num_threads;
             i_area < map_mgr.num_areas;
             i_area += map_mgr.num_threads)
        {
            const area_params_t area_params = map_mgr.get_area_params (i_area);
            const rect_t& zero_rect = area_params.zero_rect;

            u_t_buf.resize (area_params.zero_rect.t_chunk_size);
            u_x_buf.resize (area_params.zero_rect.x_chunk_size + 1);

            std::size_t zero_pos = zero_rect.t_i_min * u_buf_x_size + zero_rect.x_i_min;
            double* u_buf_begin = u_buf.data () + zero_pos;

            // Get x and t array from prev rank
            MPI_Recv (u_x_buf.data (), u_x_buf.size (),
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
            MPI_Recv (u_t_buf.data (), u_t_buf.size (),
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            // Calc zero area
            calc_u_rect_x_row_t_col (u_buf_begin, u_x_buf.data (), u_t_buf.data (),
                                     u_buf_x_size, area_params.zero_rect,
                                     funcs::f, funcs::u_next);

            for (int i_chunk_g = 1 + i_area; i_chunk_g < num_chunk_area; ++i_chunk_g) {
                int i_chunk_l = i_chunk_g - i_area;
                double* buf_begin_odd  = u_buf_begin + i_chunk_l * x_chunk_size;
                double* buf_begin_even = u_buf_begin + i_chunk_l * task.x_size;

                int x_i_odd_min = x_i_min + i_chunk_l * x_chunk_size;
                int x_i_odd_end = calc_end_index (x_i_odd_min, x_chunk_size, task.x_size);
                int x_chunk_size_corrected = x_i_odd_end - x_i_odd_min;

                // Calc u odd
                MPI_Recv (buf_begin_odd - task.x_size - 1, x_chunk_size_corrected + 1,
                        MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

                calc_u_full_buf (x_chunk_size_corrected + 1, t_chunk_size + 1,
                                task.x_size,
                                x_i_odd_min - 1 , t_i_min - 1, dx, dt,
                                u_buf, f, u_next);

                // Send up
                MPI_Send (buf_begin_odd + (t_chunk_size - 1) * task.x_size - 1, x_chunk_size_corrected + 1,
                        MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);

                int t_i_odd_min = t_i_min + i_chunk_l * t_chunk_size;
                int t_i_odd_end = calc_end_index (t_i_odd_min, t_chunk_size, task.t_size);
                int t_chunk_size_corrected = t_i_odd_end - t_i_odd_min;

                // Calc u even
                MPI_Recv (u_t_buf.data (), t_chunk_size_corrected, MPI_DOUBLE,
                        prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
                copy_row_2_col (buf_begin_even - 1, u_t_buf.data (),
                                task.x_size, t_chunk_size_corrected);

                calc_u_full_buf (x_chunk_size + 1, t_chunk_size_corrected + 1,
                                task.x_size,
                                x_i_min - 1, t_i_odd_min - 1, dx, dt,
                                u_buf, f, u_next);

                // Send right
                copy_col_2_row (u_t_buf.data (),
                                buf_begin_even + x_chunk_size - 1,
                                task.x_size, t_chunk_size_corrected);
                MPI_Send (u_t_buf.data (), t_chunk_size_corrected,
                        MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }
        }
    } // void calc_zero_rank_grid_border ()

    void
    calc_non_zero_rank_grid_border (int rank) {

    } // void calc_non_zero_rank_grid_border (int rank)

    void
    receive_u_bufs () {

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
        // todo
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
        // todo
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
        // todo
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
        // todo
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